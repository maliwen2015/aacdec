/*
 * ps_check - verify that HE-AACv2 Parametric Stereo produces a stereo image.
 *
 *   ps_check input.aac
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "aacdec.h"

int main(int argc, char** argv)
{
    const char* path = (argc > 1) ? argv[1] : "tests/data/he_v2.aac";
    FILE* f;
    aac_param_t param;
    aac_packet_t pkt;
    aac_stream_info info;
    aac_frame_t frame;
    void* dec;
    unsigned char buf[4096];
    double sl = 0.0, sr = 0.0, slr = 0.0, sd = 0.0;
    long count = 0;
    int frames = 0;
    int channels = 0;
    int rate = 0;
    int frame_size = 0;
    int eof = 0;
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
        if (n == 0)
            eof = 1;
        memset(&pkt, 0, sizeof(pkt));
        pkt.data = buf;
        pkt.len = (int)n;
        if (n > 0 && aac_dec_send_packet(dec, &pkt) < 0)
            break;
        for (;;) {
            memset(&info, 0, sizeof(info));
            memset(&frame, 0, sizeof(frame));
            ret = aac_dec_recv_frame(dec, &info, &frame);
            if (ret == AAC_DEC_GOT_HEADER) {
                if (frame_size == 0) {
                    channels = info.channels;
                    rate = info.sample_rate;
                    frame_size = info.frame_size;
                }
                continue;
            }
            if (ret == AAC_DEC_GOT_FRAME) {
                const float* pcm = (const float*)frame.planes[0];
                int i;
                for (i = 0; i < frame.frame_size; i++) {
                    double l = pcm[(size_t)i * frame.channels + 0];
                    double r = pcm[(size_t)i * frame.channels + 1];
                    sl += l * l;
                    sr += r * r;
                    slr += l * r;
                    sd += (l - r) * (l - r);
                    count++;
                }
                frames++;
                aac_dec_frame_unref(dec, &frame);
                continue;
            }
            if (ret == AAC_DEC_DEFAULT || ret == AAC_DEC_END)
                break;
            aac_dec_close(dec);
            fclose(f);
            return 1;
        }
    }

    aac_dec_close(dec);
    fclose(f);

    if (frames < 2 || count == 0 || channels != 2 || rate != 44100 ||
        frame_size != 2048) {
        fprintf(stderr, "unexpected output: %d frames, %d Hz, %d ch, %d\n",
                frames, rate, channels, frame_size);
        return 1;
    }
    {
        double corr = slr / sqrt(sl * sr + 1e-30);
        double side = sqrt(sd / (double)count) * 0.5;
        if (corr > 0.5 || side < 1e-3) {
            fprintf(stderr, "PS did not produce a stereo image: corr=%g side=%g\n",
                    corr, side);
            return 1;
        }
        printf("ok: %s, %d frames, %d Hz, LR corr=%g, side=%g\n", path,
               frames, rate, corr, side);
    }
    return 0;
}
