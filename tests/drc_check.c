/*
 * drc_check - verify MPEG-4 dynamic range control parsing.
 *
 *   drc_check
 *
 * The first part feeds a synthetic dynamic_range_info() payload to the core
 * parser.  The second part builds a minimal raw_data_block containing an SCE
 * and a FIL/DRC element and checks that the public stream info exports the
 * program reference level.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "aacdec.h"
#include "aac_bits.h"
#include "aac_core.h"

typedef struct
{
    unsigned char data[128];
    size_t bits;
} bit_writer;

static void bw_put(bit_writer* w, unsigned value, int n)
{
    int i;
    for (i = n - 1; i >= 0; i--) {
        size_t byte = w->bits >> 3;
        int shift = 7 - (int)(w->bits & 7u);
        if (byte >= sizeof(w->data))
            return;
        if ((value >> i) & 1u)
            w->data[byte] |= (unsigned char)(1u << shift);
        w->bits++;
    }
}

static size_t bw_bytes(const bit_writer* w)
{
    return (w->bits + 7u) / 8u;
}

static void put_drc_payload(bit_writer* w)
{
    bw_put(w, 1, 1);            /* pce_tag_present */
    bw_put(w, 3, 4);            /* pce_instance_tag */
    bw_put(w, 0, 4);            /* tag_reserved_bits */
    bw_put(w, 0, 1);            /* excluded_chns_present */
    bw_put(w, 1, 1);            /* drc_bands_present */
    bw_put(w, 1, 4);            /* band_incr => 2 bands */
    bw_put(w, 2, 4);            /* interpolation_scheme */
    bw_put(w, 10, 8);           /* band_top[0] */
    bw_put(w, 20, 8);           /* band_top[1] */
    bw_put(w, 1, 1);            /* prog_ref_level_present */
    bw_put(w, 80, 7);           /* prog_ref_level */
    bw_put(w, 0, 1);            /* prog_ref_level_reserved_bits */
    bw_put(w, 0, 1);            /* dyn_rng_sgn[0] */
    bw_put(w, 40, 7);           /* dyn_rng_ctl[0] */
    bw_put(w, 1, 1);            /* dyn_rng_sgn[1] */
    bw_put(w, 20, 7);           /* dyn_rng_ctl[1] */
}

static size_t build_au(unsigned char* out, int with_drc)
{
    bit_writer w;
    size_t bits;

    memset(&w, 0, sizeof(w));
    bw_put(&w, 0, 3);           /* ID_SCE */
    bw_put(&w, 0, 4);           /* element_instance_tag */
    bw_put(&w, 0, 8);           /* global_gain */
    bw_put(&w, 0, 1);           /* reserved */
    bw_put(&w, 0, 2);           /* ONLY_LONG_SEQUENCE */
    bw_put(&w, 0, 1);           /* window_shape */
    bw_put(&w, 0, 6);           /* max_sfb */
    bw_put(&w, 0, 1);           /* predictor_data_present */
    bw_put(&w, 0, 1);           /* pulse_data_present */
    bw_put(&w, 0, 1);           /* tns_data_present */
    bw_put(&w, 0, 1);           /* gain_control_data_present */

    if (with_drc) {
        bw_put(&w, 6, 3);       /* ID_FIL */
        bw_put(&w, 8, 4);       /* count: 4-bit extension type + 60 bits */
        bw_put(&w, 11, 4);      /* EXT_DYNAMIC_RANGE */
        put_drc_payload(&w);
    }

    bw_put(&w, 7, 3);           /* ID_END */
    bits = w.bits;
    while (w.bits & 7u)
        bw_put(&w, 0, 1);
    memcpy(out, w.data, bw_bytes(&w));
    (void)bits;
    return bw_bytes(&w);
}

static int test_parser(void)
{
    aac_core core;
    aac_bits b;
    bit_writer w;

    memset(&w, 0, sizeof(w));
    put_drc_payload(&w);

    if (aac_core_init(&core) < 0)
        return -1;
    aac_bits_init(&b, w.data, bw_bytes(&w));
    if (aac_core_parse_dynamic_range(&core, &b) != AAC_CORE_OK) {
        fprintf(stderr, "dynamic_range_info parse failed\n");
        aac_core_close(&core);
        return -1;
    }
    if (!core.cfg.drc_present || core.cfg.drc_num_bands != 2 ||
        core.cfg.drc_prog_ref_level != 80 ||
        core.cfg.drc_pce_instance_tag != 3 ||
        core.cfg.drc_dyn_rng_sgn[0] != 0 || core.cfg.drc_dyn_rng_ctl[0] != 40 ||
        core.cfg.drc_dyn_rng_sgn[1] != 1 || core.cfg.drc_dyn_rng_ctl[1] != 20) {
        fprintf(stderr, "unexpected DRC values\n");
        aac_core_close(&core);
        return -1;
    }
    aac_core_close(&core);
    return 0;
}

static int test_stream_info(void)
{
    static const unsigned char asc[] = { 0x12, 0x08 };
    unsigned char au[64];
    aac_param_t param;
    aac_packet_t pkt;
    aac_stream_info info;
    aac_stream_info saved;
    aac_frame_t frame;
    void* dec;
    int ret;
    int saw_header = 0;
    int saw_frame = 0;
    size_t au_len = build_au(au, 1);

    memset(&saved, 0, sizeof(saved));

    memset(&param, 0, sizeof(param));
    param.transport = AAC_TRANSPORT_RAW;
    param.config = asc;
    param.config_len = sizeof(asc);
    param.output_format = AAC_OUTPUT_S16;
    dec = aac_dec_open(&param);
    if (dec == NULL) {
        fprintf(stderr, "aac_dec_open failed\n");
        return -1;
    }

    memset(&pkt, 0, sizeof(pkt));
    pkt.data = au;
    pkt.len = (int)au_len;
    if (aac_dec_send_packet(dec, &pkt) < 0) {
        fprintf(stderr, "send_packet failed\n");
        aac_dec_close(dec);
        return -1;
    }

    for (;;) {
        memset(&info, 0, sizeof(info));
        memset(&frame, 0, sizeof(frame));
        ret = aac_dec_recv_frame(dec, &info, &frame);
        if (ret == AAC_DEC_GOT_HEADER) {
            saw_header = 1;
            saved = info;
            continue;
        }
        if (ret == AAC_DEC_GOT_FRAME) {
            saw_frame = 1;
            if (!saw_header)
                saved = info;
            aac_dec_frame_unref(dec, &frame);
            continue;
        }
        if (ret == AAC_DEC_DEFAULT)
            break;
        fprintf(stderr, "decode failed: %s\n", aac_dec_strerror(ret));
        aac_dec_close(dec);
        return -1;
    }

    if (!saw_header || !saw_frame ||
        saved.channels != 1 || saved.frame_size != 1024 ||
        saved.sample_rate != 44100 || saved.drc_program_ref_level != 80) {
        fprintf(stderr, "unexpected stream info: %d ch, %d samples, %d Hz, DRC %d\n",
                saved.channels, saved.frame_size, saved.sample_rate,
                saved.drc_program_ref_level);
        aac_dec_close(dec);
        return -1;
    }

    aac_dec_close(dec);
    return 0;
}

int main(void)
{
    if (test_parser() < 0 || test_stream_info() < 0)
        return 1;
    printf("ok: MPEG-4 DRC parsing and stream info export\n");
    return 0;
}
