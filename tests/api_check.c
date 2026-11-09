/*
 * api_check - exercise the packet API with different packet sizes and output
 * formats.
 *
 *   api_check input.aac
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "aacdec.h"

typedef struct
{
    uint8_t* data;
    size_t len;
    size_t cap;
} byte_buffer;

static int buffer_append(byte_buffer* b, const void* data, size_t len)
{
    if (b->len + len > b->cap) {
        size_t cap = b->cap ? b->cap : 65536;
        uint8_t* p;
        while (cap < b->len + len)
            cap *= 2;
        p = (uint8_t*)realloc(b->data, cap);
        if (p == NULL)
            return -1;
        b->data = p;
        b->cap = cap;
    }
    memcpy(b->data + b->len, data, len);
    b->len += len;
    return 0;
}

static int decode_file(const char* path, int chunk, int format,
                       byte_buffer* pcm, aac_stream_info* first_info,
                       int* frames)
{
    FILE* f;
    aac_param_t param;
    aac_packet_t pkt;
    aac_stream_info info;
    aac_frame_t frame;
    void* dec;
    uint8_t* buf;
    int64_t pts = 0;
    int eof = 0;
    int headers = 0;
    int ret;

    f = fopen(path, "rb");
    if (f == NULL)
        return -1;
    memset(&param, 0, sizeof(param));
    param.transport = AAC_TRANSPORT_AUTO;
    param.output_format = format;
    dec = aac_dec_open(&param);
    if (dec == NULL) {
        fclose(f);
        return -1;
    }
    buf = (uint8_t*)malloc((size_t)chunk);
    if (buf == NULL) {
        aac_dec_close(dec);
        fclose(f);
        return -1;
    }

    while (!eof) {
        size_t n = fread(buf, 1, (size_t)chunk, f);
        if (n == 0)
            eof = 1;
        else {
            memset(&pkt, 0, sizeof(pkt));
            pkt.data = buf;
            pkt.len = (int)n;
            pkt.pts = pts;
            pkt.dts = pts;
            pts += 1024;
            if (aac_dec_send_packet(dec, &pkt) < 0) {
                free(buf);
                aac_dec_close(dec);
                fclose(f);
                return -1;
            }
        }
        for (;;) {
            memset(&info, 0, sizeof(info));
            memset(&frame, 0, sizeof(frame));
            ret = aac_dec_recv_frame(dec, &info, &frame);
            if (ret == AAC_DEC_GOT_HEADER) {
                if (headers == 0 && first_info != NULL)
                    *first_info = info;
                headers++;
                continue;
            }
            if (ret == AAC_DEC_GOT_FRAME) {
                size_t bytes = (size_t)frame.frame_size *
                               (size_t)frame.channels *
                               (size_t)frame.bytes_per_sample;
                if (buffer_append(pcm, frame.planes[0], bytes) < 0) {
                    aac_dec_frame_unref(dec, &frame);
                    free(buf);
                    aac_dec_close(dec);
                    fclose(f);
                    return -1;
                }
                (*frames)++;
                aac_dec_frame_unref(dec, &frame);
                continue;
            }
            if (ret == AAC_DEC_DEFAULT || ret == AAC_DEC_END)
                break;
            free(buf);
            aac_dec_close(dec);
            fclose(f);
            return -1;
        }
    }

    for (;;) {
        memset(&info, 0, sizeof(info));
        memset(&frame, 0, sizeof(frame));
        ret = aac_dec_flush(dec, &info, &frame);
        if (ret == AAC_DEC_GOT_HEADER) {
            if (headers == 0 && first_info != NULL)
                *first_info = info;
            headers++;
            continue;
        }
        if (ret == AAC_DEC_GOT_FRAME) {
            size_t bytes = (size_t)frame.frame_size *
                           (size_t)frame.channels *
                           (size_t)frame.bytes_per_sample;
            if (buffer_append(pcm, frame.planes[0], bytes) < 0) {
                aac_dec_frame_unref(dec, &frame);
                free(buf);
                aac_dec_close(dec);
                fclose(f);
                return -1;
            }
            (*frames)++;
            aac_dec_frame_unref(dec, &frame);
            continue;
        }
        if (ret == AAC_DEC_END || ret == AAC_DEC_DEFAULT)
            break;
        free(buf);
        aac_dec_close(dec);
        fclose(f);
        return -1;
    }

    free(buf);
    aac_dec_close(dec);
    fclose(f);
    return 0;
}

static int run(const char* path, int chunk, int format, byte_buffer* pcm,
               aac_stream_info* info, int* frames)
{
    pcm->data = NULL;
    pcm->len = 0;
    pcm->cap = 0;
    *frames = 0;
    return decode_file(path, chunk, format, pcm, info, frames);
}

int main(int argc, char** argv)
{
    const char* path = (argc > 1) ? argv[1] : "tests/data/tiny_stereo.aac";
    byte_buffer ref, cur;
    aac_stream_info info_ref, info_cur;
    int frames_ref = 0, frames_cur = 0;
    int chunks[] = { 1, 7, 4096 };
    int i;

    if (aac_dec_version() == NULL) {
        fprintf(stderr, "version string is NULL\n");
        return 1;
    }

    if (run(path, 4096, AAC_OUTPUT_S16, &ref, &info_ref, &frames_ref) < 0) {
        fprintf(stderr, "failed to decode %s\n", path);
        return 1;
    }
    if (frames_ref <= 0 || ref.len == 0) {
        fprintf(stderr, "no frames decoded\n");
        free(ref.data);
        return 1;
    }
    if (info_ref.sample_rate <= 0 || info_ref.channels <= 0 ||
        info_ref.frame_size <= 0) {
        fprintf(stderr, "invalid stream info\n");
        free(ref.data);
        return 1;
    }
    if (info_ref.aot != AAC_AOT_AAC_LC) {
        fprintf(stderr, "unexpected AOT %d\n", info_ref.aot);
        free(ref.data);
        return 1;
    }

    for (i = 0; i < (int)(sizeof(chunks) / sizeof(chunks[0])); i++) {
        if (run(path, chunks[i], AAC_OUTPUT_S16, &cur, &info_cur,
                &frames_cur) < 0) {
            fprintf(stderr, "chunk %d failed\n", chunks[i]);
            free(ref.data);
            return 1;
        }
        if (frames_cur != frames_ref || cur.len != ref.len ||
            memcmp(cur.data, ref.data, ref.len) != 0) {
            fprintf(stderr,
                    "chunk %d mismatch: frames %d/%d bytes %lu/%lu\n",
                    chunks[i], frames_cur, frames_ref,
                    (unsigned long)cur.len, (unsigned long)ref.len);
            free(ref.data);
            free(cur.data);
            return 1;
        }
        if (info_cur.sample_rate != info_ref.sample_rate ||
            info_cur.channels != info_ref.channels ||
            info_cur.frame_size != info_ref.frame_size) {
            fprintf(stderr, "chunk %d stream info mismatch\n", chunks[i]);
            free(ref.data);
            free(cur.data);
            return 1;
        }
        free(cur.data);
    }

    /* Float and planar outputs should contain the same number of samples. */
    {
        aac_stream_info info_f;
        int frames_f = 0;
        byte_buffer f32;
        if (run(path, 4096, AAC_OUTPUT_F32, &f32, &info_f, &frames_f) < 0) {
            fprintf(stderr, "float output failed\n");
            free(ref.data);
            return 1;
        }
        if (frames_f != frames_ref ||
            f32.len != (size_t)frames_ref * info_ref.frame_size *
                           info_ref.channels * 4) {
            fprintf(stderr, "float output size mismatch\n");
            free(ref.data);
            free(f32.data);
            return 1;
        }
        free(f32.data);

        if (run(path, 4096, AAC_OUTPUT_S16_PLANAR, &cur, &info_cur,
                &frames_cur) < 0) {
            fprintf(stderr, "planar output failed\n");
            free(ref.data);
            return 1;
        }
        if (frames_cur != frames_ref || cur.len != ref.len) {
            fprintf(stderr, "planar output size mismatch\n");
            free(ref.data);
            free(cur.data);
            return 1;
        }
        free(cur.data);
    }

    printf("ok: %s, %d frames, %d Hz, %d ch\n", path, frames_ref,
           info_ref.sample_rate, info_ref.channels);
    free(ref.data);
    return 0;
}
