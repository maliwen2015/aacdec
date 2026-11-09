/*
 * conceal_check - verify error concealment and recovery.
 *
 *   conceal_check [stream.aac] [reference.s16] [config.asc]
 *
 * Three scenarios are checked:
 *   1. A corrupt ADTS frame is replaced by a faded copy of the previous frame
 *      and the following frames decode normally again.
 *   2. A caller-signalled lost access unit (raw transport) can be concealed
 *      with aac_dec_conceal_frame().
 *   3. AAC_DEC_FLAG_DISABLE_CONCEALMENT makes the decoder report an error
 *      instead of producing a substitute frame.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "aacdec.h"

#define MAX_FRAMES 64

typedef struct
{
    size_t offset;
    size_t length;
    size_t header;
} adts_frame;

static unsigned char* read_file(const char* path, size_t* size)
{
    FILE* f = fopen(path, "rb");
    unsigned char* data;
    long n;

    if (f == NULL) {
        perror(path);
        return NULL;
    }
    fseek(f, 0, SEEK_END);
    n = ftell(f);
    fseek(f, 0, SEEK_SET);
    if (n <= 0) {
        fclose(f);
        return NULL;
    }
    data = (unsigned char*)malloc((size_t)n);
    if (data == NULL || fread(data, 1, (size_t)n, f) != (size_t)n) {
        free(data);
        fclose(f);
        return NULL;
    }
    fclose(f);
    *size = (size_t)n;
    return data;
}

static int parse_adts(const unsigned char* data, size_t size, adts_frame* out,
                      int max_frames)
{
    size_t pos = 0;
    int count = 0;

    while (pos + 7 <= size && count < max_frames) {
        size_t len;
        size_t hdr;

        if (data[pos] != 0xFF || (data[pos + 1] & 0xF6) != 0xF0)
            return -1;
        len = ((size_t)(data[pos + 3] & 3) << 11) |
              ((size_t)data[pos + 4] << 3) |
              ((size_t)data[pos + 5] >> 5);
        hdr = (data[pos + 1] & 1) ? 7 : 9;
        if (len < hdr || pos + len > size)
            return -1;
        out[count].offset = pos;
        out[count].length = len;
        out[count].header = hdr;
        count++;
        pos += len;
    }
    return count;
}

static int collect(void* dec, short* pcm, size_t* samples, int* errors,
                   int* frames, int* failed)
{
    aac_stream_info info;
    aac_frame_t frame;
    int ret;

    for (;;) {
        int i;
        size_t count;

        memset(&info, 0, sizeof(info));
        memset(&frame, 0, sizeof(frame));
        ret = aac_dec_recv_frame(dec, &info, &frame);
        if (ret == AAC_DEC_GOT_HEADER)
            continue;
        if (ret == AAC_DEC_GOT_FRAME) {
            count = (size_t)frame.channels * (size_t)frame.frame_size;
            if (*frames >= MAX_FRAMES) {
                aac_dec_frame_unref(dec, &frame);
                return -1;
            }
            errors[*frames] = frame.error;
            for (i = 0; i < (int)count; i++)
                pcm[*samples + (size_t)i] = ((const short*)frame.planes[0])[i];
            *samples += count;
            (*frames)++;
            aac_dec_frame_unref(dec, &frame);
            continue;
        }
        if (ret == AAC_DEC_DEFAULT)
            return 0;
        *failed = ret;
        return -1;
    }
}

static int flush_frames(void* dec, short* pcm, size_t* samples, int* errors,
                        int* frames)
{
    aac_stream_info info;
    aac_frame_t frame;
    int ret;

    for (;;) {
        int i;
        size_t count;

        memset(&info, 0, sizeof(info));
        memset(&frame, 0, sizeof(frame));
        ret = aac_dec_flush(dec, &info, &frame);
        if (ret == AAC_DEC_GOT_HEADER)
            continue;
        if (ret == AAC_DEC_GOT_FRAME) {
            count = (size_t)frame.channels * (size_t)frame.frame_size;
            if (*frames >= MAX_FRAMES) {
                aac_dec_frame_unref(dec, &frame);
                return -1;
            }
            errors[*frames] = frame.error;
            for (i = 0; i < (int)count; i++)
                pcm[*samples + (size_t)i] = ((const short*)frame.planes[0])[i];
            *samples += count;
            (*frames)++;
            aac_dec_frame_unref(dec, &frame);
            continue;
        }
        if (ret == AAC_DEC_END || ret == AAC_DEC_DEFAULT)
            return 0;
        return -1;
    }
}

static int max_diff(const short* a, const short* b, size_t count)
{
    size_t i;
    int worst = 0;
    for (i = 0; i < count; i++) {
        int d = (int)a[i] - (int)b[i];
        if (d < 0)
            d = -d;
        if (d > worst)
            worst = d;
    }
    return worst;
}

static int compare_stream(const short* got, const short* ref, int frames,
                          int first_frame, int tolerance)
{
    const size_t frame_samples = 1024 * 2;
    int f;
    for (f = first_frame; f < frames; f++) {
        int d = max_diff(got + (size_t)f * frame_samples,
                         ref + (size_t)f * frame_samples, frame_samples);
        if (d > tolerance) {
            fprintf(stderr, "frame %d mismatch: max error %d\n", f, d);
            return -1;
        }
    }
    return 0;
}

static int test_corrupt_adts(const char* path, const char* ref_path)
{
    size_t size;
    size_t ref_size;
    unsigned char* data = read_file(path, &size);
    unsigned char* ref_data = read_file(ref_path, &ref_size);
    adts_frame frames[MAX_FRAMES];
    int frame_count;
    short* ref = NULL;
    short* got = NULL;
    int errors[MAX_FRAMES];
    size_t got_samples = 0;
    size_t ref_samples = ref_size / sizeof(short);
    int got_frames = 0;
    int failed = 0;
    aac_param_t param;
    aac_packet_t pkt;
    void* dec;
    int i;

    if (data == NULL || ref_data == NULL)
        goto fail;
    ref = (short*)malloc(ref_size);
    got = (short*)malloc(ref_size + 8192 * sizeof(short));
    if (ref == NULL || got == NULL)
        goto fail;
    memcpy(ref, ref_data, ref_size);
    frame_count = parse_adts(data, size, frames, MAX_FRAMES);
    if (frame_count < 8) {
        fprintf(stderr, "not enough ADTS frames\n");
        goto fail;
    }

    memset(&param, 0, sizeof(param));
    param.transport = AAC_TRANSPORT_AUTO;
    param.output_format = AAC_OUTPUT_S16;
    dec = aac_dec_open(&param);
    if (dec == NULL) {
        fprintf(stderr, "aac_dec_open failed\n");
        goto fail;
    }

    /* Corrupt the fifth frame: the payload turns into an empty
     * raw_data_block (ID_END), which the core must reject. */
    for (i = 0; i < (int)frames[4].length - (int)frames[4].header; i++)
        data[frames[4].offset + frames[4].header + (size_t)i] = 0xFF;

    memset(&pkt, 0, sizeof(pkt));
    pkt.data = data;
    pkt.len = (int)size;
    if (aac_dec_send_packet(dec, &pkt) < 0) {
        fprintf(stderr, "send_packet failed\n");
        aac_dec_close(dec);
        goto fail;
    }
    if (collect(dec, got, &got_samples, errors, &got_frames, &failed) < 0 ||
        flush_frames(dec, got, &got_samples, errors, &got_frames) < 0) {
        fprintf(stderr, "decode failed: %s\n", aac_dec_strerror(failed));
        aac_dec_close(dec);
        goto fail;
    }
    aac_dec_close(dec);

    if (got_frames != frame_count || got_samples != ref_samples) {
        fprintf(stderr, "frame count mismatch: got %d/%zu, expected %d/%zu\n",
                got_frames, got_samples, frame_count, ref_samples);
        goto fail;
    }
    {
        int error_count = 0;
        for (i = 0; i < got_frames; i++) {
            if (errors[i])
                error_count++;
        }
        if (error_count != 1 || !errors[4]) {
            fprintf(stderr, "unexpected concealment flags\n");
            goto fail;
        }
    }
    if (max_diff(got + 4 * 2048, got + 3 * 2048, 2048) != 0) {
        fprintf(stderr, "concealed frame is not a copy of the previous frame\n");
        goto fail;
    }
    /* Frame 5 still uses the overlap from frame 3; frame 6 and later are
     * identical to the pristine decode again. */
    if (compare_stream(got, ref, got_frames, 6, 8) < 0)
        goto fail;

    free(ref);
    free(got);
    free(data);
    free(ref_data);
    return 0;

fail:
    free(ref);
    free(got);
    free(data);
    free(ref_data);
    return -1;
}

static int test_signalled_loss(const char* path, const char* ref_path,
                               const char* asc_path)
{
    size_t size;
    size_t ref_size;
    size_t asc_size;
    unsigned char* data = read_file(path, &size);
    unsigned char* ref_data = read_file(ref_path, &ref_size);
    unsigned char* asc = read_file(asc_path, &asc_size);
    adts_frame frames[MAX_FRAMES];
    int frame_count;
    short* ref = NULL;
    short* got = NULL;
    int errors[MAX_FRAMES];
    size_t got_samples = 0;
    int got_frames = 0;
    int failed = 0;
    aac_param_t param;
    aac_packet_t pkt;
    void* dec;
    int i;

    if (data == NULL || ref_data == NULL || asc == NULL)
        goto fail;
    ref = (short*)malloc(ref_size);
    got = (short*)malloc(ref_size + 8192 * sizeof(short));
    if (ref == NULL || got == NULL)
        goto fail;
    memcpy(ref, ref_data, ref_size);
    frame_count = parse_adts(data, size, frames, MAX_FRAMES);
    if (frame_count < 8) {
        fprintf(stderr, "not enough ADTS frames\n");
        goto fail;
    }

    memset(&param, 0, sizeof(param));
    param.transport = AAC_TRANSPORT_RAW;
    param.config = asc;
    param.config_len = (int)asc_size;
    param.output_format = AAC_OUTPUT_S16;
    dec = aac_dec_open(&param);
    if (dec == NULL) {
        fprintf(stderr, "aac_dec_open failed\n");
        goto fail;
    }

    for (i = 0; i < frame_count; i++) {
        if (i == 4) {
            aac_stream_info info;
            aac_frame_t frame;
            int ret;
            memset(&info, 0, sizeof(info));
            memset(&frame, 0, sizeof(frame));
            ret = aac_dec_conceal_frame(dec, &info, &frame);
            if (ret != AAC_DEC_GOT_FRAME || !frame.error) {
                fprintf(stderr, "conceal_frame failed: %d\n", ret);
                aac_dec_close(dec);
                goto fail;
            }
            {
                size_t count = (size_t)frame.channels * frame.frame_size;
                int k;
                errors[got_frames] = 1;
                for (k = 0; k < (int)count; k++)
                    got[got_samples + (size_t)k] =
                        ((const short*)frame.planes[0])[k];
                got_samples += count;
                got_frames++;
            }
            aac_dec_frame_unref(dec, &frame);
        } else {
            memset(&pkt, 0, sizeof(pkt));
            pkt.data = data + frames[i].offset + frames[i].header;
            pkt.len = (int)(frames[i].length - frames[i].header);
            if (aac_dec_send_packet(dec, &pkt) < 0) {
                fprintf(stderr, "send_packet failed\n");
                aac_dec_close(dec);
                goto fail;
            }
            if (collect(dec, got, &got_samples, errors, &got_frames,
                        &failed) < 0) {
                fprintf(stderr, "decode failed: %s\n",
                        aac_dec_strerror(failed));
                aac_dec_close(dec);
                goto fail;
            }
        }
    }
    if (flush_frames(dec, got, &got_samples, errors, &got_frames) < 0) {
        aac_dec_close(dec);
        goto fail;
    }
    aac_dec_close(dec);

    if (got_frames != frame_count || got_samples != ref_size / 2) {
        fprintf(stderr, "frame count mismatch: got %d/%zu, expected %d\n",
                got_frames, got_samples, frame_count);
        goto fail;
    }
    if (!errors[4]) {
        fprintf(stderr, "signalled loss was not concealed\n");
        goto fail;
    }
    if (max_diff(got + 4 * 2048, got + 3 * 2048, 2048) != 0) {
        fprintf(stderr, "concealed frame is not a copy of the previous frame\n");
        goto fail;
    }
    if (compare_stream(got, ref, got_frames, 6, 8) < 0)
        goto fail;

    free(ref);
    free(got);
    free(data);
    free(ref_data);
    free(asc);
    return 0;

fail:
    free(ref);
    free(got);
    free(data);
    free(ref_data);
    free(asc);
    return -1;
}

static int test_disabled(const char* path)
{
    size_t size;
    unsigned char* data = read_file(path, &size);
    adts_frame frames[MAX_FRAMES];
    int frame_count;
    aac_param_t param;
    aac_packet_t pkt;
    void* dec;
    int ret;
    int i;

    if (data == NULL)
        return -1;
    frame_count = parse_adts(data, size, frames, MAX_FRAMES);
    if (frame_count < 6) {
        free(data);
        return -1;
    }
    memset(&param, 0, sizeof(param));
    param.transport = AAC_TRANSPORT_AUTO;
    param.output_format = AAC_OUTPUT_S16;
    param.flags = AAC_DEC_FLAG_DISABLE_CONCEALMENT;
    dec = aac_dec_open(&param);
    if (dec == NULL) {
        free(data);
        return -1;
    }
    for (i = 0; i < (int)frames[4].length - (int)frames[4].header; i++)
        data[frames[4].offset + frames[4].header + (size_t)i] = 0xFF;
    memset(&pkt, 0, sizeof(pkt));
    pkt.data = data;
    pkt.len = (int)size;
    aac_dec_send_packet(dec, &pkt);

    for (;;) {
        aac_stream_info info;
        aac_frame_t frame;
        memset(&info, 0, sizeof(info));
        memset(&frame, 0, sizeof(frame));
        ret = aac_dec_recv_frame(dec, &info, &frame);
        if (ret == AAC_DEC_GOT_HEADER)
            continue;
        if (ret == AAC_DEC_GOT_FRAME) {
            aac_dec_frame_unref(dec, &frame);
            continue;
        }
        break;
    }
    aac_dec_close(dec);
    free(data);
    if (ret != AAC_DEC_ERROR) {
        fprintf(stderr, "disabled concealment did not report an error\n");
        return -1;
    }
    return 0;
}

int main(int argc, char** argv)
{
    const char* path = (argc > 1) ? argv[1] : "tests/data/tiny_stereo.aac";
    const char* ref_path =
        (argc > 2) ? argv[2] : "tests/data/tiny_stereo.s16";
    const char* asc_path =
        (argc > 3) ? argv[3] : "tests/data/tiny_stereo.asc";

    if (test_corrupt_adts(path, ref_path) < 0)
        return 1;
    if (test_signalled_loss(path, ref_path, asc_path) < 0)
        return 1;
    if (test_disabled(path) < 0)
        return 1;

    printf("ok: error concealment, signalled loss and disable flag\n");
    return 0;
}
