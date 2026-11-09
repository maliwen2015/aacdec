/*
 * crc_check - verify ADTS header CRC detection and recovery.
 *
 *   crc_check [stream.aac] [reference.s16]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "aacdec.h"

#define MAX_FRAMES 64

typedef struct
{
    size_t offset;
    size_t length;
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

        if (data[pos] != 0xFF || (data[pos + 1] & 0xF6) != 0xF0)
            return -1;
        len = ((size_t)(data[pos + 3] & 3) << 11) |
              ((size_t)data[pos + 4] << 3) |
              ((size_t)data[pos + 5] >> 5);
        if (len < 7 || pos + len > size)
            return -1;
        out[count].offset = pos;
        out[count].length = len;
        count++;
        pos += len;
    }
    return count;
}

static int collect(void* dec, short* pcm, size_t* samples, int* errors,
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
        return -1;
    }
}

int main(int argc, char** argv)
{
    const char* path = (argc > 1) ? argv[1] : "tests/data/tiny_adts_crc4.aac";
    const char* ref_path =
        (argc > 2) ? argv[2] : "tests/data/tiny_adts_crc_ref.s16";
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
    int got_frames = 0;
    aac_param_t param;
    aac_packet_t pkt;
    void* dec;
    int i;
    int result = 1;

    if (data == NULL || ref_data == NULL)
        goto cleanup;
    ref = (short*)malloc(ref_size);
    got = (short*)malloc(ref_size + 8192 * sizeof(short));
    if (ref == NULL || got == NULL)
        goto cleanup;
    memcpy(ref, ref_data, ref_size);
    frame_count = parse_adts(data, size, frames, MAX_FRAMES);
    if (frame_count < 3) {
        fprintf(stderr, "not enough ADTS frames\n");
        goto cleanup;
    }

    /* Corrupt the header CRC of the last transport frame. */
    {
        size_t off = frames[frame_count - 1].offset;
        int blocks = (data[off + 6] & 3) + 1;
        size_t crc_pos = off + 7 + 2 * (size_t)(blocks - 1);
        int good_frames = (frame_count - 1) * blocks;
        data[crc_pos + 1] ^= 0x01;

        memset(&param, 0, sizeof(param));
        param.transport = AAC_TRANSPORT_AUTO;
        param.output_format = AAC_OUTPUT_S16;
        dec = aac_dec_open(&param);
        if (dec == NULL) {
            fprintf(stderr, "aac_dec_open failed\n");
            goto cleanup;
        }
        memset(&pkt, 0, sizeof(pkt));
        pkt.data = data;
        pkt.len = (int)size;
        if (aac_dec_send_packet(dec, &pkt) < 0) {
            fprintf(stderr, "send_packet failed\n");
            aac_dec_close(dec);
            goto cleanup;
        }
        if (collect(dec, got, &got_samples, errors, &got_frames) < 0) {
            fprintf(stderr, "decode failed\n");
            aac_dec_close(dec);
            goto cleanup;
        }
        aac_dec_close(dec);

        if (got_frames != good_frames + 1 ||
            got_samples != (size_t)(good_frames + 1) * 2048) {
            fprintf(stderr, "frame count mismatch: got %d/%zu, expected %d\n",
                    got_frames, got_samples, good_frames + 1);
            goto cleanup;
        }
        if (!errors[good_frames]) {
            fprintf(stderr, "header CRC error was not detected\n");
            goto cleanup;
        }
        for (i = 0; i < good_frames; i++) {
            size_t k;
            for (k = 0; k < 2048; k++) {
                int d = (int)got[(size_t)i * 2048 + k] -
                        (int)ref[(size_t)i * 2048 + k];
                if (d < 0)
                    d = -d;
                if (d > 8) {
                    fprintf(stderr, "frame %d mismatch after CRC loss: %d\n",
                            i, d);
                    goto cleanup;
                }
            }
        }
    }

    printf("ok: ADTS header CRC detection and recovery\n");
    result = 0;

cleanup:
    free(ref);
    free(got);
    free(data);
    free(ref_data);
    return result;
}
