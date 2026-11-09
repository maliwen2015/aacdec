/*
 * compare_check - decode an AAC file and compare the S16 interleaved PCM
 * against a reference dump.
 *
 *   compare_check input.aac reference.s16 [tolerance]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "aacdec.h"

static int parse_transport(const char* s)
{
    if (s == NULL || strcmp(s, "auto") == 0)
        return AAC_TRANSPORT_AUTO;
    if (strcmp(s, "adts") == 0)
        return AAC_TRANSPORT_ADTS;
    if (strcmp(s, "adif") == 0)
        return AAC_TRANSPORT_ADIF;
    if (strcmp(s, "latm") == 0 || strcmp(s, "loas") == 0)
        return AAC_TRANSPORT_LATM;
    if (strcmp(s, "raw") == 0)
        return AAC_TRANSPORT_RAW;
    if (strcmp(s, "mcp0") == 0)
        return AAC_TRANSPORT_LATM_MCP0;
    if (strcmp(s, "mcp1") == 0)
        return AAC_TRANSPORT_LATM_MCP1;
    return -1;
}

static int decode(const char* path, const char* transport_name,
                  const char* config_path, int16_t** out, size_t* samples)
{
    FILE* f;
    aac_param_t param;
    aac_packet_t pkt;
    aac_frame_t frame;
    aac_stream_info info;
    void* dec;
    uint8_t buf[4096];
    uint8_t* config = NULL;
    int config_len = 0;
    int16_t* pcm = NULL;
    size_t count = 0, cap = 0;
    int64_t pts = 0;
    int eof = 0;
    int ret;

    if (config_path != NULL) {
        FILE* cf = fopen(config_path, "rb");
        long n;
        if (cf == NULL)
            return -1;
        fseek(cf, 0, SEEK_END);
        n = ftell(cf);
        rewind(cf);
        if (n > 0 && n < 4096) {
            config = (uint8_t*)malloc((size_t)n);
            if (config == NULL || fread(config, 1, (size_t)n, cf) != (size_t)n) {
                free(config);
                fclose(cf);
                return -1;
            }
            config_len = (int)n;
        }
        fclose(cf);
    }

    f = fopen(path, "rb");
    if (f == NULL) {
        free(config);
        return -1;
    }
    memset(&param, 0, sizeof(param));
    param.transport = parse_transport(transport_name);
    if (param.transport < 0) {
        free(config);
        fclose(f);
        return -1;
    }
    param.config = config;
    param.config_len = config_len;
    dec = aac_dec_open(&param);
    free(config);
    if (dec == NULL) {
        fclose(f);
        return -1;
    }

    while (!eof) {
        size_t n = fread(buf, 1, sizeof(buf), f);
        if (n == 0)
            eof = 1;
        else {
            memset(&pkt, 0, sizeof(pkt));
            pkt.data = buf;
            pkt.len = (int)n;
            pkt.pts = pts;
            pkt.dts = pts;
            pts += 1024;
            if (aac_dec_send_packet(dec, &pkt) < 0)
                goto fail;
        }
        for (;;) {
            memset(&frame, 0, sizeof(frame));
            ret = aac_dec_recv_frame(dec, &info, &frame);
            if (ret == AAC_DEC_GOT_HEADER)
                continue;
            if (ret == AAC_DEC_GOT_FRAME) {
                size_t n2 = (size_t)frame.frame_size *
                            (size_t)frame.channels;
                if (count + n2 > cap) {
                    size_t newcap = cap ? cap * 2 : 65536;
                    int16_t* p;
                    while (newcap < count + n2)
                        newcap *= 2;
                    p = (int16_t*)realloc(pcm, newcap * 2);
                    if (p == NULL) {
                        aac_dec_frame_unref(dec, &frame);
                        goto fail;
                    }
                    pcm = p;
                    cap = newcap;
                }
                memcpy(pcm + count, frame.planes[0], n2 * 2);
                count += n2;
                aac_dec_frame_unref(dec, &frame);
                continue;
            }
            if (ret == AAC_DEC_DEFAULT || ret == AAC_DEC_END)
                break;
            goto fail;
        }
    }
    for (;;) {
        memset(&frame, 0, sizeof(frame));
        ret = aac_dec_flush(dec, &info, &frame);
        if (ret == AAC_DEC_GOT_HEADER)
            continue;
        if (ret == AAC_DEC_GOT_FRAME) {
            size_t n2 = (size_t)frame.frame_size * (size_t)frame.channels;
            if (count + n2 > cap) {
                size_t newcap = cap ? cap * 2 : 65536;
                int16_t* p;
                while (newcap < count + n2)
                    newcap *= 2;
                p = (int16_t*)realloc(pcm, newcap * 2);
                if (p == NULL) {
                    aac_dec_frame_unref(dec, &frame);
                    goto fail;
                }
                pcm = p;
                cap = newcap;
            }
            memcpy(pcm + count, frame.planes[0], n2 * 2);
            count += n2;
            aac_dec_frame_unref(dec, &frame);
            continue;
        }
        if (ret == AAC_DEC_END || ret == AAC_DEC_DEFAULT)
            break;
        goto fail;
    }
    aac_dec_close(dec);
    fclose(f);
    *out = pcm;
    *samples = count;
    return 0;

fail:
    free(pcm);
    aac_dec_close(dec);
    fclose(f);
    return -1;
}

int main(int argc, char** argv)
{
    int16_t* ours = NULL;
    int16_t* ref = NULL;
    size_t ours_n = 0, ref_n = 0;
    long max_diff = 0;
    double sum = 0.0;
    size_t i;
    int tolerance = 32;
    const char* transport = NULL;
    const char* config_path = NULL;
    FILE* f;
    long ref_bytes;

    if (argc < 3) {
        fprintf(stderr,
                "usage: %s input reference.s16 [tolerance] [transport] [asc]\n",
                argv[0]);
        return 2;
    }
    if (argc > 3)
        tolerance = atoi(argv[3]);
    if (argc > 4)
        transport = argv[4];
    if (argc > 5)
        config_path = argv[5];
    if (decode(argv[1], transport, config_path, &ours, &ours_n) < 0) {
        fprintf(stderr, "decode failed\n");
        return 1;
    }
    f = fopen(argv[2], "rb");
    if (f == NULL) {
        perror(argv[2]);
        free(ours);
        return 1;
    }
    fseek(f, 0, SEEK_END);
    ref_bytes = ftell(f);
    rewind(f);
    ref_n = (size_t)ref_bytes / 2;
    ref = (int16_t*)malloc(ref_n * 2);
    if (ref == NULL || fread(ref, 2, ref_n, f) != ref_n) {
        fprintf(stderr, "cannot read reference\n");
        fclose(f);
        free(ours);
        free(ref);
        return 1;
    }
    fclose(f);

    if (ours_n != ref_n) {
        fprintf(stderr, "sample count mismatch: %lu/%lu\n",
                (unsigned long)ours_n, (unsigned long)ref_n);
        free(ours);
        free(ref);
        return 1;
    }
    for (i = 0; i < ours_n; i++) {
        long d = (long)ours[i] - (long)ref[i];
        if (d < 0)
            d = -d;
        if (d > max_diff)
            max_diff = d;
        sum += (double)d;
    }
    free(ours);
    free(ref);
    if (max_diff > tolerance) {
        fprintf(stderr, "mismatch: max=%ld avg=%g tolerance=%d\n", max_diff,
                sum / (double)ours_n, tolerance);
        return 1;
    }
    printf("ok: %s max=%ld avg=%g\n", argv[1], max_diff,
           sum / (double)ours_n);
    return 0;
}
