/*
 * he_check - decode an HE-AAC ADTS stream and verify the SBR output format.
 *
 *   he_check input.aac
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "aacdec.h"

int main(int argc, char** argv)
{
    const char* path = (argc > 1) ? argv[1] : "tests/data/he_stereo.aac";
    FILE* f;
    aac_param_t param;
    aac_packet_t pkt;
    aac_stream_info info;
    aac_frame_t frame;
    void* dec;
    unsigned char buf[4096];
    int frames = 0;
    int saw_header = 0;
    int rate = 0;
    int frame_size = 0;
    int channels = 0;
    int eof = 0;
    double energy = 0.0;
    int ret;

    f = fopen(path, "rb");
    if (f == NULL) {
        perror(path);
        return 1;
    }
    memset(&param, 0, sizeof(param));
    param.transport = AAC_TRANSPORT_AUTO;
    param.output_format = AAC_OUTPUT_F32;
    dec = aac_dec_open(&param);
    if (dec == NULL) {
        fclose(f);
        return 1;
    }

    while (!eof) {
        size_t n = fread(buf, 1, sizeof(buf), f);
        int i;
        if (n == 0)
            eof = 1;
        memset(&pkt, 0, sizeof(pkt));
        pkt.data = buf;
        pkt.len = (int)n;
        if (n > 0 && aac_dec_send_packet(dec, &pkt) < 0) {
            fprintf(stderr, "send_packet failed\n");
            aac_dec_close(dec);
            fclose(f);
            return 1;
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
                }
                continue;
            }
            if (ret == AAC_DEC_GOT_FRAME) {
                const float* pcm = (const float*)frame.planes[0];
                int count = frame.channels * frame.frame_size;
                int i;
                for (i = 0; i < count; i++) {
                    if (!isfinite(pcm[i])) {
                        fprintf(stderr, "non-finite PCM sample\n");
                        aac_dec_frame_unref(dec, &frame);
                        aac_dec_close(dec);
                        fclose(f);
                        return 1;
                    }
                    energy += (double)pcm[i] * pcm[i];
                }
                frames++;
                aac_dec_frame_unref(dec, &frame);
                continue;
            }
            if (ret == AAC_DEC_DEFAULT || ret == AAC_DEC_END)
                break;
            fprintf(stderr, "decode failed: %s\n", aac_dec_strerror(ret));
            aac_dec_close(dec);
            fclose(f);
            return 1;
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
            for (i = 0; i < count; i++)
                energy += (double)pcm[i] * pcm[i];
            frames++;
            aac_dec_frame_unref(dec, &frame);
            continue;
        }
        break;
    }

    aac_dec_close(dec);
    fclose(f);

    if (frames < 2 || !saw_header) {
        fprintf(stderr, "no HE-AAC frames decoded\n");
        return 1;
    }
    if (rate != 44100 || frame_size != 2048 ||
        (channels != 1 && channels != 2)) {
        fprintf(stderr, "unexpected SBR format: %d Hz, %d samples, %d ch\n",
                rate, frame_size, channels);
        return 1;
    }
    if (!(energy > 1.0)) {
        fprintf(stderr, "decoded audio is silent\n");
        return 1;
    }
    printf("ok: %s, %d frames, %d Hz, %d samples, %d ch\n", path, frames,
           rate, frame_size, channels);
    return 0;
}
