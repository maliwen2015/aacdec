/*
 * eld_check - decode an AAC-ELD LOAS stream and compare against FDK reference.
 *
 *   eld_check input.latm reference.pcm
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "aacdec.h"

int main(int argc, char** argv)
{
    const char* path = (argc > 1) ? argv[1] : "tests/data/tiny_eld.latm";
    const char* ref_path =
        (argc > 2) ? argv[2] : "tests/data/tiny_eld_ref.pcm";
    FILE* f;
    FILE* rf;
    aac_param_t param;
    aac_packet_t pkt;
    aac_stream_info info;
    aac_frame_t frame;
    void* dec;
    unsigned char buf[4096];
    short* ref;
    float* got = NULL;
    size_t ref_samples = 0;
    size_t got_samples = 0;
    size_t got_cap = 0;
    int frames = 0;
    int saw_header = 0;
    int rate = 0;
    int frame_size = 0;
    int channels = 0;
    int aot = 0;
    int eof = 0;
    int ret;
    double sig = 0.0;
    double err = 0.0;
    double max_err = 0.0;

    f = fopen(path, "rb");
    if (f == NULL) {
        perror(path);
        return 1;
    }
    rf = fopen(ref_path, "rb");
    if (rf == NULL) {
        perror(ref_path);
        fclose(f);
        return 1;
    }
    fseek(rf, 0, SEEK_END);
    ref_samples = (size_t)ftell(rf) / sizeof(short);
    fseek(rf, 0, SEEK_SET);
    ref = (short*)malloc(ref_samples * sizeof(short));
    if (ref == NULL || fread(ref, sizeof(short), ref_samples, rf) !=
                           ref_samples) {
        fprintf(stderr, "cannot read reference PCM\n");
        fclose(rf);
        fclose(f);
        free(ref);
        return 1;
    }
    fclose(rf);

    memset(&param, 0, sizeof(param));
    param.transport = AAC_TRANSPORT_AUTO;
    param.output_format = AAC_OUTPUT_F32;
    dec = aac_dec_open(&param);
    if (dec == NULL) {
        fclose(f);
        free(ref);
        return 1;
    }

    while (!eof) {
        size_t n = fread(buf, 1, sizeof(buf), f);
        if (n == 0)
            eof = 1;
        memset(&pkt, 0, sizeof(pkt));
        pkt.data = buf;
        pkt.len = (int)n;
        if (n > 0 && aac_dec_send_packet(dec, &pkt) < 0) {
            fprintf(stderr, "send_packet failed\n");
            goto fail;
        }
        for (;;) {
            memset(&info, 0, sizeof(info));
            memset(&frame, 0, sizeof(frame));
            ret = aac_dec_recv_frame(dec, &info, &frame);
            if (ret == AAC_DEC_GOT_HEADER) {
                if (!saw_header) {
                    saw_header = 1;
                    rate = info.sample_rate;
                    frame_size = info.frame_size;
                    channels = info.channels;
                    aot = info.aot;
                }
                continue;
            }
            if (ret == AAC_DEC_GOT_FRAME) {
                const float* pcm = (const float*)frame.planes[0];
                int count = frame.channels * frame.frame_size;
                int i;
                if (got_samples + (size_t)count > got_cap) {
                    size_t cap = got_cap ? got_cap * 2 : 65536;
                    float* p = (float*)realloc(got, cap * sizeof(float));
                    if (p == NULL)
                        goto fail;
                    got = p;
                    got_cap = cap;
                }
                for (i = 0; i < count; i++)
                    got[got_samples++] = pcm[i] * 32768.0f;
                frames++;
                aac_dec_frame_unref(dec, &frame);
                continue;
            }
            if (ret == AAC_DEC_DEFAULT || ret == AAC_DEC_END)
                break;
            fprintf(stderr, "decode failed: %s\n", aac_dec_strerror(ret));
            goto fail;
        }
    }
    for (;;) {
        memset(&info, 0, sizeof(info));
        memset(&frame, 0, sizeof(frame));
        ret = aac_dec_flush(dec, &info, &frame);
        if (ret == AAC_DEC_GOT_HEADER)
            continue;
        if (ret == AAC_DEC_GOT_FRAME) {
            const float* pcm = (const float*)frame.planes[0];
            int count = frame.channels * frame.frame_size;
            int i;
            if (got_samples + (size_t)count > got_cap) {
                size_t cap = got_cap ? got_cap * 2 : 65536;
                float* p = (float*)realloc(got, cap * sizeof(float));
                if (p == NULL)
                    goto fail;
                got = p;
                got_cap = cap;
            }
            for (i = 0; i < count; i++)
                got[got_samples++] = pcm[i] * 32768.0f;
            frames++;
            aac_dec_frame_unref(dec, &frame);
            continue;
        }
        break;
    }

    aac_dec_close(dec);
    dec = NULL;
    fclose(f);

    if (!saw_header || frames < 2) {
        fprintf(stderr, "no AAC-ELD frames decoded\n");
        goto fail;
    }
    if (rate != 48000 || frame_size != 512 || channels != 2 || aot != 39) {
        fprintf(stderr,
                "unexpected AAC-ELD format: %d Hz, %d samples, %d ch, AOT %d\n",
                rate, frame_size, channels, aot);
        goto fail;
    }
    if (got_samples != ref_samples) {
        fprintf(stderr, "sample count mismatch: got %zu, expected %zu\n",
                got_samples, ref_samples);
        goto fail;
    }

    {
        /* The reference decoder fades in its first frames; compare the
         * steady-state waveform. */
        size_t skip = (size_t)frame_size * (size_t)channels * 10;
        size_t i;
        if (skip >= ref_samples)
            skip = 0;
        for (i = skip; i < ref_samples; i++) {
            double r = (double)ref[i];
            double d = (double)got[i] - r;
            sig += r * r;
            err += d * d;
            if (fabs(d) > max_err)
                max_err = fabs(d);
        }
    }
    if (!(sig > 1.0) || 10.0 * log10(sig / (err + 1e-30)) < 60.0 ||
        max_err > 4.0) {
        fprintf(stderr,
                "AAC-ELD reference mismatch: SNR %.2f dB, max error %.3f\n",
                10.0 * log10(sig / (err + 1e-30)), max_err);
        goto fail;
    }

    printf("ok: %s, %d frames, %d Hz, %d samples, %d ch, SNR %.2f dB\n",
           path, frames, rate, frame_size, channels,
           10.0 * log10(sig / (err + 1e-30)));
    free(got);
    free(ref);
    return 0;

fail:
    if (dec != NULL)
        aac_dec_close(dec);
    fclose(f);
    free(got);
    free(ref);
    return 1;
}
