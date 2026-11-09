/*
 * usac_check - decode an xHE-AAC (USAC, AOT 42) LOAS stream and verify the
 * output format and that every access unit decodes without errors.
 *
 *   usac_check input.loas
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "aacdec.h"

int main(int argc, char** argv)
{
    const char* path = (argc > 1) ? argv[1] : "tests/data/usac_stereo.loas";
    FILE* f;
    aac_param_t param;
    aac_packet_t pkt;
    aac_stream_info info;
    aac_frame_t frame;
    void* dec;
    unsigned char buf[4096];
    int frames = 0;
    int bad = 0;
    int saw_header = 0;
    int aot = 0;
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
    param.transport = AAC_TRANSPORT_LOAS;
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
                    aot = info.aot;
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
                if (frame.error)
                    bad++;
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

    aac_dec_close(dec);
    fclose(f);

    if (frames < 2 || !saw_header) {
        fprintf(stderr, "no USAC frames decoded\n");
        return 1;
    }
    if (aot != 42 || rate != 48000 || frame_size != 2048 || channels != 2) {
        fprintf(stderr, "unexpected USAC format: AOT %d, %d Hz, %d samples, "
                        "%d ch\n",
                aot, rate, frame_size, channels);
        return 1;
    }
    if (bad != 0) {
        fprintf(stderr, "%d of %d USAC frames reported errors\n", bad, frames);
        return 1;
    }
    if (!(energy > 1.0)) {
        fprintf(stderr, "decoded audio is silent\n");
        return 1;
    }
    printf("ok: %s, %d frames, AOT %d, %d Hz, %d samples, %d ch\n", path,
           frames, aot, rate, frame_size, channels);
    return 0;
}
