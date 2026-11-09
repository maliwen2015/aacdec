/*
 * er_check - decode a minimal ER AAC-LC (AOT 17) access unit.
 *
 *   er_check
 *
 * The stream is built in memory: an AudioSpecificConfig with AOT 17 and one
 * silent SCE.  It exercises the error-resilient raw_data_block path.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "aacdec.h"

typedef struct
{
    unsigned char data[32];
    size_t bits;
} bit_writer;

static void bw_put(bit_writer* w, unsigned value, int n)
{
    int i;
    for (i = n - 1; i >= 0; i--) {
        size_t byte = w->bits >> 3;
        int shift = 7 - (int)(w->bits & 7u);
        if ((value >> i) & 1u)
            w->data[byte] |= (unsigned char)(1u << shift);
        w->bits++;
    }
}

static size_t bw_bytes(const bit_writer* w)
{
    return (w->bits + 7u) / 8u;
}

static size_t build_asc(unsigned char* out)
{
    bit_writer w;

    memset(&w, 0, sizeof(w));
    bw_put(&w, 17, 5);      /* AOT_ER_AAC_LC */
    bw_put(&w, 4, 4);       /* 44100 Hz */
    bw_put(&w, 1, 4);       /* mono */
    bw_put(&w, 0, 1);       /* frameLengthFlag */
    bw_put(&w, 0, 1);       /* dependsOnCoreCoder */
    bw_put(&w, 1, 1);       /* extensionFlag */
    bw_put(&w, 0, 3);       /* resilience flags */
    bw_put(&w, 0, 1);       /* extensionFlag3 */
    bw_put(&w, 0, 2);       /* epConfig */
    while (w.bits & 7u)
        bw_put(&w, 0, 1);
    memcpy(out, w.data, bw_bytes(&w));
    return bw_bytes(&w);
}

static size_t build_au(unsigned char* out)
{
    bit_writer w;

    memset(&w, 0, sizeof(w));
    bw_put(&w, 0, 4);       /* element_instance_tag */
    bw_put(&w, 0, 8);       /* global_gain */
    bw_put(&w, 0, 1);       /* reserved */
    bw_put(&w, 0, 2);       /* ONLY_LONG_SEQUENCE */
    bw_put(&w, 0, 1);       /* window_shape */
    bw_put(&w, 0, 6);       /* max_sfb */
    bw_put(&w, 0, 1);       /* predictor_data_present */
    bw_put(&w, 0, 1);       /* pulse_data_present */
    bw_put(&w, 0, 1);       /* tns_data_present */
    bw_put(&w, 0, 1);       /* gain_control_data_present */
    while (w.bits & 7u)
        bw_put(&w, 0, 1);
    memcpy(out, w.data, bw_bytes(&w));
    return bw_bytes(&w);
}

int main(void)
{
    unsigned char asc[8];
    unsigned char au[16];
    size_t asc_len = build_asc(asc);
    size_t au_len = build_au(au);
    aac_param_t param;
    aac_packet_t pkt;
    aac_stream_info info;
    aac_frame_t frame;
    void* dec;
    int ret;
    int saw_header = 0;
    int saw_frame = 0;
    int rate = 0;
    int frame_size = 0;
    int channels = 0;
    int aot = 0;

    memset(&param, 0, sizeof(param));
    param.transport = AAC_TRANSPORT_RAW;
    param.config = asc;
    param.config_len = (int)asc_len;
    param.output_format = AAC_OUTPUT_S16;
    dec = aac_dec_open(&param);
    if (dec == NULL) {
        fprintf(stderr, "aac_dec_open failed\n");
        return 1;
    }

    memset(&pkt, 0, sizeof(pkt));
    pkt.data = au;
    pkt.len = (int)au_len;
    if (aac_dec_send_packet(dec, &pkt) < 0) {
        fprintf(stderr, "send_packet failed\n");
        aac_dec_close(dec);
        return 1;
    }

    for (;;) {
        memset(&info, 0, sizeof(info));
        memset(&frame, 0, sizeof(frame));
        ret = aac_dec_recv_frame(dec, &info, &frame);
        if (ret == AAC_DEC_GOT_HEADER) {
            saw_header = 1;
            rate = info.sample_rate;
            frame_size = info.frame_size;
            channels = info.channels;
            aot = info.aot;
            continue;
        }
        if (ret == AAC_DEC_GOT_FRAME) {
            saw_frame = 1;
            aac_dec_frame_unref(dec, &frame);
            continue;
        }
        if (ret == AAC_DEC_DEFAULT)
            break;
        fprintf(stderr, "decode failed: %s\n", aac_dec_strerror(ret));
        aac_dec_close(dec);
        return 1;
    }
    aac_dec_close(dec);

    if (!saw_header || !saw_frame || rate != 44100 || frame_size != 1024 ||
        channels != 1 || aot != 17) {
        fprintf(stderr,
                "unexpected ER format: header=%d frame=%d %d Hz %d samples "
                "%d ch AOT %d\n",
                saw_header, saw_frame, rate, frame_size, channels, aot);
        return 1;
    }

    printf("ok: ER AAC-LC raw_data_block, %d Hz, %d samples, %d ch, AOT %d\n",
           rate, frame_size, channels, aot);
    return 0;
}
