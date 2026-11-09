/*
 * aacdec - packet based public API and transport demultiplexer.
 */
#include "aacdec.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "aac_bits.h"
#include "aac_core.h"
#include "aac_tables.h"


#define AAC_DETECT_PROBE 8192
#define AAC_TS_MAX_ENTRIES (1u << 18)
#define AAC_LATM_MAX_SUBFRAMES 64

#ifdef AAC_USAC_DEBUG
#include <stdio.h>
#define USDBG(...) do { fprintf(stderr, __VA_ARGS__); } while (0)
#else
#define USDBG(...)
#endif

/* ------------------------------------------------------------------------- */
/* Growable byte buffer                                                      */
/* ------------------------------------------------------------------------- */

typedef struct
{
    uint8_t* data;
    size_t len;
    size_t cap;
    size_t pos;                 /* first unconsumed byte */
} aac_buffer;

typedef struct
{
    int64_t pts;
    int64_t dts;
    void* opaque;
} aac_ts_entry;

typedef struct
{
    aac_ts_entry* v;
    size_t head;
    size_t count;
    size_t cap;
} aac_ts_fifo;

typedef struct
{
    size_t start;
    size_t len;
    int64_t pts;
    int64_t dts;
    void* opaque;
} aac_send_record;

typedef struct
{
    uint8_t* buf;
    size_t len;
    size_t cap;
    int in_frame;
    size_t frame_pos;
    size_t frame_len;
    int frame_rdb;
    int64_t frame_pts;
    int64_t frame_dts;
    void* frame_opaque;
} aac_scanner;

typedef struct aac_decoder
{
    aac_core core;
    aac_param_t param;
    uint8_t* asc;
    int asc_len;
    int asc_set;

    int transport;
    int opened;
    int fatal;
    int started;
    int flushed;

    aac_buffer in;
    aac_ts_fifo ts;
    aac_scanner scan;
    aac_send_record* records;
    int record_count;

    int adif_header_done;
    int adif_bitstream_type;
    int adif_bitrate;

    int latm_config_valid;
    int latm_num_subframes;
    int latm_frame_length_type;
    int latm_frame_length;
    int latm_audio_mux_version_a;
    int latm_other_data_bits;
    int latm_crc_present;

    uint8_t* au_buf;
    size_t au_len;
    size_t au_pos;
    int au_blocks_left;
    int au_block_crc;           /* per-raw_data_block CRC bytes (ADTS) */

    /* Decoded AudioMuxElement subframes (LATM/LOAS). */
    uint8_t* latm_au_buf;
    size_t latm_au_cap;
    size_t latm_au_pos;
    size_t latm_sub_off[AAC_LATM_MAX_SUBFRAMES];
    size_t latm_sub_len[AAC_LATM_MAX_SUBFRAMES];
    int latm_sub_count;
    int latm_sub_index;

    float* pcm_scratch;
    float* conceal_pcm;         /* last good frame used for concealment */
    int conceal_valid;
    int conceal_channels;
    int conceal_frame_size;
    int conceal_sample_rate;
    int conceal_run;

    int info_valid;
    aac_stream_info info;
    aac_frame_t pending_frame;
    int has_pending;

    int64_t last_pts;
    int64_t last_dts;
    void* last_opaque;

    int64_t frames_out;
    int64_t bad_frames;
    int64_t total_bytes;
} aac_decoder;

/* ------------------------------------------------------------------------- */
/* Small helpers                                                             */
/* ------------------------------------------------------------------------- */

static int buffer_append(aac_decoder* d, const uint8_t* data, size_t len)
{
    aac_buffer* b = &d->in;
    size_t need;

    if (len == 0)
        return 0;
    need = b->len + len;
    if (need > b->cap) {
        size_t cap = b->cap ? b->cap : 4096;
        uint8_t* p;
        while (cap < need)
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

static size_t buffer_avail(const aac_buffer* b)
{
    return b->len - b->pos;
}

static uint8_t* buffer_ptr(aac_buffer* b)
{
    return b->data + b->pos;
}

static void buffer_consume(aac_buffer* b, size_t n)
{
    if (n > b->len - b->pos)
        n = b->len - b->pos;
    b->pos += n;
    if (b->pos == b->len) {
        b->pos = 0;
        b->len = 0;
    }
}

static void buffer_compact(aac_buffer* b)
{
    if (b->pos > 0 && b->pos < b->len) {
        memmove(b->data, b->data + b->pos, b->len - b->pos);
        b->len -= b->pos;
        b->pos = 0;
    } else if (b->pos >= b->len) {
        b->pos = 0;
        b->len = 0;
    }
}

static int ts_reserve(aac_ts_fifo* f, size_t extra)
{
    size_t need = f->count + extra;
    size_t cap;
    aac_ts_entry* v;
    size_t i;

    if (need <= f->cap)
        return 0;
    cap = f->cap ? f->cap : 4096;
    while (cap < need)
        cap *= 2;
    if (cap > AAC_TS_MAX_ENTRIES)
        cap = AAC_TS_MAX_ENTRIES;
    v = (aac_ts_entry*)malloc(cap * sizeof(*v));
    if (v == NULL)
        return -1;
    for (i = 0; i < f->count; i++)
        v[i] = f->v[(f->head + i) % f->cap];
    free(f->v);
    f->v = v;
    f->cap = cap;
    f->head = 0;
    return 0;
}

static int ts_push(aac_decoder* d, int64_t pts, int64_t dts, void* opaque)
{
    aac_ts_fifo* f = &d->ts;
    if (f->count == AAC_TS_MAX_ENTRIES) {
        f->head = (f->head + 1) % f->cap;
        f->count--;
    }
    if (ts_reserve(f, 1) < 0)
        return -1;
    {
        size_t idx = (f->head + f->count) % f->cap;
        f->v[idx].pts = pts;
        f->v[idx].dts = dts;
        f->v[idx].opaque = opaque;
        f->count++;
    }
    return 0;
}

static int ts_pop(aac_ts_fifo* f, aac_ts_entry* out)
{
    if (f->count == 0)
        return 0;
    *out = f->v[f->head];
    f->head = (f->head + 1) % f->cap;
    f->count--;
    if (f->count == 0)
        f->head = 0;
    return 1;
}

static void ts_free(aac_ts_fifo* f)
{
    free(f->v);
    memset(f, 0, sizeof(*f));
}

static void frame_release(aac_frame_t* frame)
{
    if (frame == NULL)
        return;
    free(frame->dec_frame);
    memset(frame, 0, sizeof(*frame));
}

/* ------------------------------------------------------------------------- */
/* Timestamp scanner for ADTS / LOAS                                         */
/* ------------------------------------------------------------------------- */

static int scan_reserve(aac_scanner* s, size_t extra)
{
    size_t need = s->len + extra;
    size_t cap;
    uint8_t* p;

    if (need <= s->cap)
        return 0;
    cap = s->cap ? s->cap : 4096;
    while (cap < need)
        cap *= 2;
    p = (uint8_t*)realloc(s->buf, cap);
    if (p == NULL)
        return -1;
    s->buf = p;
    s->cap = cap;
    return 0;
}

static int scan_push(aac_decoder* d, const uint8_t* data, size_t len,
                     int64_t pts, int64_t dts, void* opaque)
{
    aac_scanner* s = &d->scan;
    const int is_adts = (d->transport == AAC_TRANSPORT_ADTS);
    size_t pos = 0;

    if (len > 0) {
        if (scan_reserve(s, len) < 0)
            return -1;
        memcpy(s->buf + s->len, data, len);
        s->len += len;
    }

    while (pos < s->len) {
        if (!s->in_frame) {
            size_t i = pos;
            int found = 0;
            for (; i + 1 < s->len; i++) {
                if (is_adts) {
                    if (s->buf[i] == 0xFF && (s->buf[i + 1] & 0xF6) == 0xF0) {
                        found = 1;
                        break;
                    }
                } else {
                    if (s->buf[i] == 0x56 && (s->buf[i + 1] & 0xE0) == 0xE0) {
                        found = 1;
                        break;
                    }
                }
            }
            if (!found) {
                size_t keep = 0;
                if (s->len > pos) {
                    uint8_t last = s->buf[s->len - 1];
                    if ((is_adts && last == 0xFF) ||
                        (!is_adts && last == 0x56))
                        keep = 1;
                }
                pos = s->len - keep;
                break;
            }
            pos = i;

            if (is_adts) {
                size_t hdr;
                unsigned layer, flen, sfidx;
                if (s->len - pos < 7)
                    break;
                hdr = (s->buf[pos + 1] & 1) ? 7 : 9;
                if (s->len - pos < hdr)
                    break;
                layer = (unsigned)((s->buf[pos + 1] >> 1) & 3);
                flen = ((unsigned)(s->buf[pos + 3] & 3) << 11) |
                       ((unsigned)s->buf[pos + 4] << 3) |
                       ((unsigned)s->buf[pos + 5] >> 5);
                sfidx = (unsigned)((s->buf[pos + 2] >> 2) & 0xF);
                if (layer != 0 || flen < hdr || flen > 8191 || sfidx == 0xF) {
                    pos++;
                    continue;
                }
                s->frame_len = flen;
                s->frame_rdb = (s->buf[pos + 6] & 3) + 1;
            } else {
                unsigned aml, flen;
                if (s->len - pos < 3)
                    break;
                aml = ((unsigned)(s->buf[pos + 1] & 0x0F) << 8) |
                      (unsigned)s->buf[pos + 2];
                flen = 3 + aml;
                if (flen < 3 || flen > 8194) {
                    pos++;
                    continue;
                }
                s->frame_len = flen;
                s->frame_rdb = 1;
            }
            s->in_frame = 1;
            s->frame_pos = pos;
            s->frame_pts = pts;
            s->frame_dts = dts;
            s->frame_opaque = opaque;
        }

        if (s->len - s->frame_pos < s->frame_len)
            break;

        {
            int k;
            for (k = 0; k < s->frame_rdb; k++) {
                if (ts_push(d, s->frame_pts, s->frame_dts,
                            s->frame_opaque) < 0)
                    return -1;
            }
        }
        pos = s->frame_pos + s->frame_len;
        s->in_frame = 0;
        s->frame_pos = 0;
        s->frame_len = 0;
    }

    if (pos > s->len)
        pos = s->len;
    if (pos > 0) {
        if (s->in_frame) {
            if (s->frame_pos >= pos)
                s->frame_pos -= pos;
            else
                s->frame_pos = 0;
        }
        memmove(s->buf, s->buf + pos, s->len - pos);
        s->len -= pos;
    }
    return 0;
}

/* ------------------------------------------------------------------------- */
/* Transport probing                                                         */
/* ------------------------------------------------------------------------- */

static int probe_transport(const uint8_t* buf, size_t len)
{
    size_t i;

    if (len >= 4 && memcmp(buf, "ADIF", 4) == 0)
        return AAC_TRANSPORT_ADIF;

    for (i = 0; i + 1 < len; i++) {
        if (buf[i] == 0xFF && (buf[i + 1] & 0xF6) == 0xF0) {
            unsigned flen, sfidx;
            if (i + 7 > len)
                return 0;
            flen = ((unsigned)(buf[i + 3] & 3) << 11) |
                   ((unsigned)buf[i + 4] << 3) |
                   ((unsigned)buf[i + 5] >> 5);
            sfidx = (unsigned)((buf[i + 2] >> 2) & 0xF);
            if (((buf[i + 1] >> 1) & 3) == 0 && flen >= 7 && flen <= 8191 &&
                sfidx != 0xF)
                return AAC_TRANSPORT_ADTS;
            continue;
        }
        if (buf[i] == 0x56 && (buf[i + 1] & 0xE0) == 0xE0) {
            if (i + 3 > len)
                return 0;
            return AAC_TRANSPORT_LATM;
        }
    }

    if (len < AAC_DETECT_PROBE)
        return 0;
    return -1;
}

/* ------------------------------------------------------------------------- */
/* LATM / LOAS                                                               */
/* ------------------------------------------------------------------------- */

static int latm_read_value(aac_bits* b, uint32_t* value)
{
    int length = (int)aac_bits_read(b, 2);
    *value = aac_bits_read(b, (length + 1) * 8);
    return b->error ? -1 : 0;
}

static int latm_read_payload_length(aac_decoder* d, aac_bits* b, int* length)
{
    if (d->latm_frame_length_type == 0) {
        int total = 0;
        int v;
        do {
            v = (int)aac_bits_read(b, 8);
            if (b->error)
                return AAC_CORE_NEED_MORE;
            total += v;
        } while (v == 255);
        *length = total;
        return AAC_CORE_OK;
    }
    if (d->latm_frame_length_type == 1) {
        *length = d->latm_frame_length;
        return AAC_CORE_OK;
    }
    return AAC_CORE_ERROR;
}

static int latm_parse_stream_mux_config(aac_decoder* d, aac_bits* b)
{
    int audio_mux_version;
    int audio_mux_version_a = 0;
    int all_same;
    int num_subframes;
    int num_program;
    int num_layer;
    int frame_length_type;
    int frame_length = 0;
    int ret;

    audio_mux_version = (int)aac_bits_read(b, 1);
    if (audio_mux_version)
        audio_mux_version_a = (int)aac_bits_read(b, 1);
    if (audio_mux_version_a)
        return AAC_CORE_ERROR;
    if (audio_mux_version) {
        uint32_t v;
        if (latm_read_value(b, &v) < 0)
            return AAC_CORE_NEED_MORE;
    }
    all_same = (int)aac_bits_read(b, 1);
    num_subframes = (int)aac_bits_read(b, 6);
    num_program = (int)aac_bits_read(b, 4);
    num_layer = (int)aac_bits_read(b, 3);

    if (all_same != 1 || num_program != 0 || num_layer != 0)
        return AAC_CORE_ERROR;

    if (!audio_mux_version) {
        ret = aac_core_parse_asc_bits(&d->core, b);
        if (ret < 0)
            return b->error ? AAC_CORE_NEED_MORE : AAC_CORE_ERROR;
    } else {
        uint32_t asc_len;
        size_t start;
        size_t end;
        aac_bits sub;

        if (latm_read_value(b, &asc_len) < 0)
            return AAC_CORE_NEED_MORE;
        start = b->bitpos;
        end = start + (size_t)asc_len;
        if (end > b->size * 8u)
            return AAC_CORE_NEED_MORE;

        sub.data = b->data;
        sub.size = (end + 7u) / 8u;
        sub.bitpos = start;
        sub.error = 0;
        ret = aac_core_parse_asc_bits(&d->core, &sub);
        b->bitpos = end;
        if (ret < 0)
            return b->error ? AAC_CORE_NEED_MORE : AAC_CORE_ERROR;
    }

    frame_length_type = (int)aac_bits_read(b, 3);
    USDBG("smc: after asc bitpos=%zu flt=%d\n", b->bitpos, frame_length_type);
    if (frame_length_type == 0) {
        (void)aac_bits_read(b, 8); /* latmBufferFullness */
    } else if (frame_length_type == 1) {
        frame_length = (int)aac_bits_read(b, 9);
    } else {
        return AAC_CORE_ERROR;
    }

    d->latm_other_data_bits = 0;
    d->latm_crc_present = 0;

    if (aac_bits_read(b, 1)) { /* otherDataPresent */
        if (audio_mux_version) {
            uint32_t v;
            if (latm_read_value(b, &v) < 0)
                return AAC_CORE_NEED_MORE;
            d->latm_other_data_bits = (int)v;
        } else {
            int esc;
            int bits = 0;
            do {
                esc = (int)aac_bits_read(b, 1);
                bits = (bits << 8) | (int)aac_bits_read(b, 8);
                if (b->error)
                    return AAC_CORE_NEED_MORE;
            } while (esc);
            d->latm_other_data_bits = bits;
        }
    }
    d->latm_crc_present = (int)aac_bits_read(b, 1);
    if (d->latm_crc_present)
        (void)aac_bits_read(b, 8);

    if (b->error)
        return AAC_CORE_NEED_MORE;

    USDBG("smc: flt=%d len=%d subs=%d other=%d crc=%d\n", frame_length_type, frame_length, num_subframes, d->latm_other_data_bits, d->latm_crc_present);
    d->latm_config_valid = 1;
    d->latm_num_subframes = num_subframes;
    d->latm_frame_length_type = frame_length_type;
    d->latm_frame_length = frame_length;
    d->latm_audio_mux_version_a = audio_mux_version_a;
    return AAC_CORE_OK;
}

/* Parse one AudioMuxElement and queue all of its subframe payloads. */
static void copy_bits_to_bytes(aac_bits* b, uint8_t* dst, size_t nbits)
{
    size_t i;
    memset(dst, 0, (nbits + 7u) >> 3);
    for (i = 0; i < nbits; i++) {
        if (aac_bits_read(b, 1))
            dst[i >> 3] |= (uint8_t)(1u << (7u - (i & 7u)));
    }
}

static int latm_parse_element(aac_decoder* d, aac_bits* b, int mux_config_present)
{
    int count;
    int i;
    int ret;
    size_t total = 0;

    if (mux_config_present) {
        int use_same = (int)aac_bits_read(b, 1);
        if (!use_same) {
            ret = latm_parse_stream_mux_config(d, b);
            if (ret < 0)
                return ret;
        } else if (!d->latm_config_valid) {
            return AAC_CORE_ERROR;
        }
    } else if (!d->latm_config_valid) {
        return AAC_CORE_ERROR;
    }

    if (d->latm_audio_mux_version_a != 0)
        return AAC_CORE_ERROR;
    count = d->latm_num_subframes + 1;
    if (count <= 0 || count > AAC_LATM_MAX_SUBFRAMES)
        return AAC_CORE_ERROR;

    for (i = 0; i < count; i++) {
        int payload_len;
        size_t need;

        ret = latm_read_payload_length(d, b, &payload_len);
        USDBG("latm: sub %d payload_len=%d left=%zu\n", i, payload_len, aac_bits_left(b));
        if (ret < 0)
            return ret;
        if (payload_len < 0 || aac_bits_left(b) < (size_t)payload_len * 8u)
            return AAC_CORE_NEED_MORE;

        need = total + (size_t)payload_len;
        if (need > d->latm_au_cap) {
            size_t cap = d->latm_au_cap ? d->latm_au_cap : 4096;
            uint8_t* p;
            while (cap < need)
                cap *= 2;
            p = (uint8_t*)realloc(d->latm_au_buf, cap);
            if (p == NULL)
                return AAC_CORE_ERROR;
            d->latm_au_buf = p;
            d->latm_au_cap = cap;
        }
        d->latm_sub_off[i] = total;
        d->latm_sub_len[i] = (size_t)payload_len;
        copy_bits_to_bytes(b, d->latm_au_buf + total,
                           (size_t)payload_len * 8u);
        total += (size_t)payload_len;
    }

    if (d->latm_other_data_bits > 0)
        aac_bits_skip(b, (size_t)d->latm_other_data_bits);

    d->latm_au_pos = total;
    d->latm_sub_count = count;
    d->latm_sub_index = 0;
    d->au_blocks_left = 0;
    d->au_block_crc = 0;
    return b->error ? AAC_CORE_NEED_MORE : AAC_CORE_OK;
}

/* ------------------------------------------------------------------------- */
/* Stream information and output frames                                      */
/* ------------------------------------------------------------------------- */

static void fill_info(aac_decoder* d, aac_stream_info* out, int channels,
                      int sample_rate)
{
    aac_core* core = &d->core;
    aac_config* c = &core->cfg;
    int i;

    memset(out, 0, sizeof(*out));
    out->sample_rate = sample_rate;
    out->frame_size = aac_core_output_frame_size(core);
    out->channels = channels;
    out->aot = c->aot;
    out->ext_aot = c->ext_aot;
    out->ext_sample_rate = c->ext_sample_rate;
    out->channel_config = c->channel_config;
    out->profile = c->profile;
    out->bitrate = d->adif_bitrate;
    out->aac_sample_rate = c->sample_rate;
    out->aac_frame_size = c->frame_size;
    out->aac_channels = c->channels;
    out->ep_config = c->ep_config ? c->ep_config : -1;
    out->output_delay = (unsigned)c->frame_size;
    out->num_total_bytes = d->total_bytes;
    out->num_bad_bytes = 0;
    out->num_total_access_units = d->frames_out;
    out->num_bad_access_units = d->bad_frames;
    out->drc_program_ref_level =
        (c->drc_present && c->drc_prog_ref_present) ? c->drc_prog_ref_level : -1;
    out->drc_presentation_mode = -1;
    out->output_loudness = -1;
    for (i = 0; i < AAC_DEC_MAX_CHANNELS; i++) {
        out->channel_types[i] = c->channel_types[i];
        out->channel_indices[i] = c->channel_indices[i];
    }
}

static int info_equal(const aac_stream_info* a, const aac_stream_info* b)
{
    if (a->sample_rate != b->sample_rate ||
        a->frame_size != b->frame_size ||
        a->channels != b->channels ||
        a->aot != b->aot ||
        a->ext_aot != b->ext_aot ||
        a->ext_sample_rate != b->ext_sample_rate ||
        a->channel_config != b->channel_config ||
        a->profile != b->profile ||
        a->aac_sample_rate != b->aac_sample_rate ||
        a->aac_frame_size != b->aac_frame_size ||
        a->aac_channels != b->aac_channels)
        return 0;
    if (memcmp(a->channel_types, b->channel_types,
               sizeof(a->channel_types)) != 0)
        return 0;
    if (memcmp(a->channel_indices, b->channel_indices,
               sizeof(a->channel_indices)) != 0)
        return 0;
    return 1;
}

static int build_frame(aac_decoder* d, const float* pcm, int channels,
                       int frame_size, int sample_rate, int error,
                       aac_frame_t* out)
{
    const int format = d->param.output_format;
    const int planar = (format == AAC_OUTPUT_S16_PLANAR ||
                        format == AAC_OUTPUT_F32_PLANAR);
    const int is_float = (format == AAC_OUTPUT_F32 ||
                          format == AAC_OUTPUT_F32_PLANAR);
    const int bps = is_float ? 4 : 2;
    const size_t total = (size_t)channels * (size_t)frame_size;
    uint8_t* buf;
    aac_ts_entry ts;
    int have_ts;
    int c, i;

    if (channels <= 0 || channels > AAC_DEC_MAX_CHANNELS ||
        frame_size <= 0 || frame_size > AAC_DEC_MAX_FRAME_SIZE)
        return -1;

    buf = (uint8_t*)malloc(total * (size_t)bps);
    if (buf == NULL)
        return -1;

    if (!is_float) {
        int16_t* dst = (int16_t*)buf;
        for (i = 0; i < (int)total; i++) {
            int v = (int)lrintf(pcm[i]);
            if (v > 32767)
                v = 32767;
            if (v < -32768)
                v = -32768;
            dst[i] = (int16_t)v;
        }
    } else {
        float* dst = (float*)buf;
        for (i = 0; i < (int)total; i++)
            dst[i] = pcm[i] * (1.0f / 32768.0f);
    }

    if (planar) {
        if (!is_float) {
            int16_t* src = (int16_t*)buf;
            int16_t* tmp = (int16_t*)malloc(total * sizeof(int16_t));
            if (tmp == NULL) {
                free(buf);
                return -1;
            }
            for (c = 0; c < channels; c++)
                for (i = 0; i < frame_size; i++)
                    tmp[(size_t)c * frame_size + i] =
                        src[(size_t)i * channels + c];
            memcpy(buf, tmp, total * sizeof(int16_t));
            free(tmp);
        } else {
            float* src = (float*)buf;
            float* tmp = (float*)malloc(total * sizeof(float));
            if (tmp == NULL) {
                free(buf);
                return -1;
            }
            for (c = 0; c < channels; c++)
                for (i = 0; i < frame_size; i++)
                    tmp[(size_t)c * frame_size + i] =
                        src[(size_t)i * channels + c];
            memcpy(buf, tmp, total * sizeof(float));
            free(tmp);
        }
    }

    memset(out, 0, sizeof(*out));
    out->dec_frame = buf;
    out->channels = channels;
    out->frame_size = frame_size;
    out->sample_rate = sample_rate;
    out->format = format;
    out->bytes_per_sample = bps;
    out->is_planar = planar;
    out->error = error;

    if (planar) {
        out->num_planes = channels;
        for (c = 0; c < channels; c++) {
            out->planes[c] = buf + (size_t)c * frame_size * bps;
            out->strides[c] = bps;
        }
    } else {
        out->num_planes = 1;
        out->planes[0] = buf;
        out->strides[0] = channels * bps;
    }

    have_ts = ts_pop(&d->ts, &ts);
    if (have_ts) {
        out->pts = ts.pts;
        out->dts = ts.dts;
        out->opaque = ts.opaque;
        d->last_pts = ts.pts;
        d->last_dts = ts.dts;
        d->last_opaque = ts.opaque;
    } else {
        out->pts = d->last_pts;
        out->dts = d->last_dts;
        out->opaque = d->last_opaque;
    }
    return 0;
}

/* ------------------------------------------------------------------------- */
/* Decoding one access unit                                                  */
/* ------------------------------------------------------------------------- */

static int finish_frame(aac_decoder* d, const float* pcm, int channels,
                        int sample_rate, int error, aac_stream_info* info,
                        aac_frame_t* frame)
{
    aac_stream_info new_info;
    int changed;

    fill_info(d, &new_info, channels, sample_rate);
    changed = !d->info_valid || !info_equal(&d->info, &new_info);
    d->info = new_info;
    d->info_valid = 1;

    if (build_frame(d, pcm, channels, aac_core_output_frame_size(&d->core),
                    sample_rate, error, frame) < 0) {
        d->fatal = 1;
        return AAC_DEC_ERROR;
    }
    d->frames_out++;

    if (changed) {
        d->pending_frame = *frame;
        d->has_pending = 1;
        memset(frame, 0, sizeof(*frame));
        if (info != NULL)
            *info = d->info;
        return AAC_DEC_GOT_HEADER;
    }

    if (info != NULL)
        *info = d->info;
    return AAC_DEC_GOT_FRAME;
}

/* ------------------------------------------------------------------------- */
/* Error concealment                                                         */
/* ------------------------------------------------------------------------- */

static void remember_pcm(aac_decoder* d, const float* pcm, int channels,
                         int sample_rate)
{
    const int frame_size = aac_core_output_frame_size(&d->core);
    const size_t total = (size_t)channels * (size_t)frame_size;

    if (d->conceal_pcm == NULL || channels <= 0 || frame_size <= 0 ||
        total > (size_t)AAC_DEC_MAX_CHANNELS * AAC_DEC_MAX_FRAME_SIZE)
        return;
    memcpy(d->conceal_pcm, pcm, total * sizeof(float));
    d->conceal_valid = 1;
    d->conceal_channels = channels;
    d->conceal_frame_size = frame_size;
    d->conceal_sample_rate = sample_rate;
    d->conceal_run = 0;
}

static int emit_concealed(aac_decoder* d, aac_stream_info* info,
                          aac_frame_t* frame)
{
    float gain;
    int total;
    int i;

    if (!d->conceal_valid || !d->info_valid || d->conceal_pcm == NULL ||
        d->conceal_frame_size != aac_core_output_frame_size(&d->core))
        return AAC_DEC_ERROR;

    /* Repeat the last good frame with an exponential fade; after four
     * consecutive losses emit silence to avoid a stuck tonal artefact. */
    total = d->conceal_channels * d->conceal_frame_size;
    gain = (d->conceal_run >= 4) ? 0.0f
                                 : 1.0f / (float)(1 << d->conceal_run);
    for (i = 0; i < total; i++)
        d->pcm_scratch[i] = d->conceal_pcm[i] * gain;
    d->conceal_run++;

    return finish_frame(d, d->pcm_scratch, d->conceal_channels,
                        d->conceal_sample_rate, 1, info, frame);
}

static int conceal_or_error(aac_decoder* d, aac_stream_info* info,
                            aac_frame_t* frame)
{
    int ret = AAC_DEC_ERROR;

    if (!(d->param.flags & AAC_DEC_FLAG_DISABLE_CONCEALMENT))
        ret = emit_concealed(d, info, frame);
    d->bad_frames++;
    return ret;
}

static void latm_next_subframe(aac_decoder* d)
{
    d->latm_sub_index++;
    if (d->latm_sub_index >= d->latm_sub_count) {
        d->latm_sub_count = 0;
        d->latm_sub_index = 0;
    }
}

static int emit_decoded(aac_decoder* d, const uint8_t* data, size_t len,
                        size_t* consumed, aac_stream_info* info,
                        aac_frame_t* frame)
{
    float* pcm = d->pcm_scratch;
    int channels = 0;
    int sample_rate = 0;
    int latm_au = 0;
    int ret;

    if (d->au_blocks_left > 0) {
        data = d->au_buf + d->au_pos;
        len = d->au_len - d->au_pos;
        *consumed = 0;
    } else if (d->latm_sub_count > 0) {
        if (d->latm_sub_index < 0 ||
            d->latm_sub_index >= d->latm_sub_count) {
            d->latm_sub_count = 0;
            d->latm_sub_index = 0;
        } else {
            data = d->latm_au_buf + d->latm_sub_off[d->latm_sub_index];
            len = d->latm_sub_len[d->latm_sub_index];
            *consumed = 0;
            latm_au = 1;
        }
    }

    ret = aac_core_decode(&d->core, data, len, consumed, pcm, &channels,
                          &sample_rate);
    if (ret == AAC_CORE_NEED_MORE &&
        (latm_au || d->au_blocks_left > 0 ||
         d->transport == AAC_TRANSPORT_RAW)) {
        /* A complete access unit was truncated: drop it. */
        ret = AAC_CORE_ERROR;
    }
    if (ret == AAC_CORE_NEED_MORE)
        return AAC_DEC_DEFAULT;
    if (ret == AAC_CORE_ERROR) {
        if (latm_au)
            latm_next_subframe(d);
        else if (d->au_blocks_left > 0) {
            d->au_blocks_left = 0;
            d->au_pos = 0;
        } else {
            /* Packetized transports deliver one access unit per packet. */
            *consumed = len;
        }
        return conceal_or_error(d, info, frame);
    }

    if (latm_au) {
        latm_next_subframe(d);
        *consumed = 0;
    } else if (d->au_blocks_left > 0) {
        d->au_pos += *consumed;
        if (d->au_block_crc > 0 &&
            d->au_pos + (size_t)d->au_block_crc <= d->au_len)
            d->au_pos += (size_t)d->au_block_crc;
        d->au_blocks_left--;
        *consumed = 0;
    } else {
        d->total_bytes += (int64_t)*consumed;
    }

    remember_pcm(d, pcm, channels, sample_rate);
    return finish_frame(d, pcm, channels, sample_rate, 0, info, frame);
}

/* ------------------------------------------------------------------------- */
/* ADTS                                                                      */
/* ------------------------------------------------------------------------- */

static int find_adts_sync(aac_buffer* b, size_t* offset)
{
    size_t n = buffer_avail(b);
    uint8_t* p = buffer_ptr(b);
    size_t i;
    for (i = 0; i + 6 < n; i++) {
        if (p[i] == 0xFF && (p[i + 1] & 0xF6) == 0xF0) {
            unsigned flen = ((unsigned)(p[i + 3] & 3) << 11) |
                            ((unsigned)p[i + 4] << 3) |
                            ((unsigned)p[i + 5] >> 5);
            unsigned sfidx = (unsigned)((p[i + 2] >> 2) & 0xF);
            if (((p[i + 1] >> 1) & 3) == 0 && flen >= 7 && flen <= 8191 &&
                sfidx != 0xF) {
                *offset = i;
                return 1;
            }
        }
    }
    return 0;
}

/* CRC-16 used by the ADTS header (ISO/IEC 13818-7, polynomial 0x8005). */
static uint16_t crc16_8005(const uint8_t* data, size_t bits)
{
    uint16_t crc = 0xFFFF;
    size_t i;

    for (i = 0; i < bits; i++) {
        uint16_t bit = (uint16_t)((data[i >> 3] >> (7u - (i & 7u))) & 1u);
        uint16_t msb = (uint16_t)((crc >> 15) & 1u);
        crc = (uint16_t)(crc << 1);
        if (msb ^ bit)
            crc ^= 0x8005;
    }
    return crc;
}

static int decode_adts(aac_decoder* d, aac_stream_info* info,
                       aac_frame_t* frame)
{
    size_t avail;
    uint8_t* p;
    size_t hdr, frame_len;
    int profile, sfidx, channel_config, mpeg2;
    int blocks;
    int protection;
    size_t payload_len;
    int ret;

    if (d->au_blocks_left > 0)
        return emit_decoded(d, NULL, 0, &payload_len, info, frame);

    buffer_compact(&d->in);
    avail = buffer_avail(&d->in);
    if (avail < 7)
        return AAC_DEC_DEFAULT;
    p = buffer_ptr(&d->in);

    if (!(p[0] == 0xFF && (p[1] & 0xF6) == 0xF0)) {
        size_t off;
        if (find_adts_sync(&d->in, &off)) {
            buffer_consume(&d->in, off);
            return AAC_DEC_DEFAULT;
        }
        if (avail > 1024)
            buffer_consume(&d->in, avail - 1);
        return AAC_DEC_DEFAULT;
    }

    blocks = (p[6] & 3) + 1;
    protection = !(p[1] & 1);
    hdr = protection ? 7 + 2 * (size_t)blocks : 7;
    if (avail < hdr)
        return AAC_DEC_DEFAULT;
    frame_len = ((size_t)(p[3] & 3) << 11) | ((size_t)p[4] << 3) |
                ((size_t)p[5] >> 5);
    if (frame_len < hdr) {
        buffer_consume(&d->in, 1);
        return AAC_DEC_DEFAULT;
    }
    if (frame_len > avail)
        return AAC_DEC_DEFAULT;

    profile = (p[2] >> 6) & 3;
    sfidx = (p[2] >> 2) & 0xF;
    channel_config = ((p[2] & 1) << 2) | (p[3] >> 6);
    mpeg2 = (p[1] >> 3) & 1;

    if (protection && blocks > 1) {
        size_t crc_bits = 56 + 16 * (size_t)(blocks - 1);
        size_t crc_pos = 7 + 2 * (size_t)(blocks - 1);
        uint16_t stored = (uint16_t)((p[crc_pos] << 8) | p[crc_pos + 1]);
        if (crc16_8005(p, crc_bits) != stored) {
            buffer_consume(&d->in, frame_len);
            return conceal_or_error(d, info, frame);
        }
    }
    d->au_block_crc = (protection && blocks > 1) ? 2 : 0;

    if (aac_core_set_adts(&d->core, profile, sfidx, channel_config,
                          mpeg2) < 0) {
        buffer_consume(&d->in, frame_len);
        return conceal_or_error(d, info, frame);
    }

    payload_len = frame_len - hdr;
    if (payload_len > d->au_len) {
        uint8_t* q = (uint8_t*)realloc(d->au_buf, payload_len);
        if (q == NULL) {
            d->fatal = 1;
            return AAC_DEC_ERROR;
        }
        d->au_buf = q;
        d->au_len = payload_len;
    }
    memcpy(d->au_buf, p + hdr, payload_len);
    d->au_pos = 0;
    d->au_blocks_left = blocks;
    buffer_consume(&d->in, frame_len);

    ret = emit_decoded(d, NULL, 0, &payload_len, info, frame);
    return ret;
}

/* ------------------------------------------------------------------------- */
/* ADIF                                                                      */
/* ------------------------------------------------------------------------- */

static int decode_adif(aac_decoder* d, aac_stream_info* info,
                       aac_frame_t* frame)
{
    size_t avail;
    uint8_t* p;
    size_t consumed;

    if (!d->adif_header_done) {
        aac_bits b;
        int copyright, bitstream_type, num_pce, i, ret;
        size_t header_len;

        buffer_compact(&d->in);
        avail = buffer_avail(&d->in);
        if (avail < 4)
            return AAC_DEC_DEFAULT;
        p = buffer_ptr(&d->in);
        if (memcmp(p, "ADIF", 4) != 0) {
            d->fatal = 1;
            return AAC_DEC_ERROR;
        }

        aac_bits_init(&b, p, avail);
        aac_bits_skip(&b, 32);
        copyright = (int)aac_bits_read(&b, 1);
        if (copyright)
            aac_bits_skip(&b, 72);
        (void)aac_bits_read(&b, 1); /* original_copy */
        (void)aac_bits_read(&b, 1); /* home */
        bitstream_type = (int)aac_bits_read(&b, 1);
        d->adif_bitrate = (int)aac_bits_read(&b, 23);
        num_pce = (int)aac_bits_read(&b, 4) + 1;
        if (bitstream_type == 0)
            aac_bits_skip(&b, 20);
        if (b.error)
            return AAC_DEC_DEFAULT;

        d->core.cfg.mpeg2 = 1;
        d->core.cfg.profile = -1;
        d->core.cfg.aot = 0;
        d->core.cfg.core_aot = 0;
        d->core.cfg.sample_rate_index = -1;
        d->core.cfg.sample_rate = 0;
        d->core.cfg.frame_size = 1024;
        d->core.cfg.short_len = 128;
        d->core.cfg.pce_valid = 0;

        for (i = 0; i < num_pce; i++) {
            ret = aac_core_parse_pce_bits(&d->core, &b);
            if (ret < 0) {
                if (b.error)
                    return AAC_DEC_DEFAULT;
                d->fatal = 1;
                return AAC_DEC_ERROR;
            }
        }
        aac_bits_align(&b);
        if (b.error)
            return AAC_DEC_DEFAULT;
        if (aac_core_apply_config(&d->core, 1) < 0) {
            d->fatal = 1;
            return AAC_DEC_ERROR;
        }
        header_len = aac_bits_bytes_consumed(&b);
        buffer_consume(&d->in, header_len);
        d->adif_header_done = 1;
        d->adif_bitstream_type = bitstream_type;
    }

    buffer_compact(&d->in);
    avail = buffer_avail(&d->in);
    if (avail == 0)
        return AAC_DEC_DEFAULT;
    p = buffer_ptr(&d->in);

    {
        int channels = 0;
        int sample_rate = 0;
        int core_ret = aac_core_decode(&d->core, p, avail, &consumed,
                                       d->pcm_scratch, &channels,
                                       &sample_rate);
        if (core_ret == AAC_CORE_NEED_MORE)
            return AAC_DEC_DEFAULT;
        if (core_ret == AAC_CORE_ERROR) {
            d->bad_frames++;
            d->fatal = 1;
            return AAC_DEC_ERROR;
        }
        buffer_consume(&d->in, consumed);
        d->total_bytes += (int64_t)consumed;
        return finish_frame(d, d->pcm_scratch, channels, sample_rate, 0, info,
                            frame);
    }
}

/* ------------------------------------------------------------------------- */
/* LOAS / LATM                                                               */
/* ------------------------------------------------------------------------- */

static int find_loas_sync(aac_buffer* b, size_t* offset)
{
    size_t n = buffer_avail(b);
    uint8_t* p = buffer_ptr(b);
    size_t i;
    for (i = 0; i + 2 < n; i++) {
        if (p[i] == 0x56 && (p[i + 1] & 0xE0) == 0xE0) {
            *offset = i;
            return 1;
        }
    }
    return 0;
}

static int decode_loas(aac_decoder* d, aac_stream_info* info,
                       aac_frame_t* frame)
{
    size_t avail;
    uint8_t* p;
    size_t frame_len;
    aac_bits b;
    int ret;
    size_t payload_len = 0;
    size_t consumed_bits;

    if (d->au_blocks_left > 0 || d->latm_sub_count > 0)
        return emit_decoded(d, NULL, 0, &payload_len, info, frame);

    buffer_compact(&d->in);
    avail = buffer_avail(&d->in);
    if (avail < 3)
        return AAC_DEC_DEFAULT;
    p = buffer_ptr(&d->in);

    if (!(p[0] == 0x56 && (p[1] & 0xE0) == 0xE0)) {
        size_t off;
        if (find_loas_sync(&d->in, &off)) {
            buffer_consume(&d->in, off);
            return AAC_DEC_DEFAULT;
        }
        if (avail > 1024)
            buffer_consume(&d->in, avail - 1);
        return AAC_DEC_DEFAULT;
    }

    frame_len = 3 + (((size_t)(p[1] & 0x1F) << 8) | (size_t)p[2]);
    if (frame_len < 3 || frame_len > 8194) {
        buffer_consume(&d->in, 1);
        return AAC_DEC_DEFAULT;
    }
    if (avail < frame_len)
        return AAC_DEC_DEFAULT;

    aac_bits_init(&b, p + 3, frame_len - 3);
    ret = latm_parse_element(d, &b, 1);
    USDBG("loas: frame_len=%zu parse ret=%d subs=%d\n", frame_len, ret,
          d->latm_sub_count);
    if (ret == AAC_CORE_NEED_MORE) {
        buffer_consume(&d->in, frame_len);
        return conceal_or_error(d, info, frame);
    }
    if (ret < 0) {
        buffer_consume(&d->in, frame_len);
        return conceal_or_error(d, info, frame);
    }
    consumed_bits = b.bitpos;
    (void)consumed_bits;
    buffer_consume(&d->in, frame_len);

    return emit_decoded(d, NULL, 0, &payload_len, info, frame);
}

/* ------------------------------------------------------------------------- */
/* Raw / packetized LATM                                                     */
/* ------------------------------------------------------------------------- */

static int decode_packetized(aac_decoder* d, aac_stream_info* info,
                             aac_frame_t* frame)
{
    size_t avail;
    uint8_t* p;
    size_t consumed = 0;
    int ret;
    size_t payload_len = 0;

    if (d->au_blocks_left > 0 || d->latm_sub_count > 0)
        return emit_decoded(d, NULL, 0, &payload_len, info, frame);

    buffer_compact(&d->in);
    avail = buffer_avail(&d->in);
    if (avail == 0)
        return AAC_DEC_DEFAULT;
    p = buffer_ptr(&d->in);

    if (d->transport == AAC_TRANSPORT_RAW) {
        int f = emit_decoded(d, p, avail, &consumed, info, frame);
        if (f == AAC_DEC_GOT_FRAME || f == AAC_DEC_GOT_HEADER) {
            buffer_consume(&d->in, consumed);
        } else if (f == AAC_DEC_ERROR) {
            buffer_consume(&d->in, consumed > 0 ? consumed : 1);
        }
        return f;
    }

    {
        aac_bits b;
        int mux_config_present = (d->transport == AAC_TRANSPORT_LATM_MCP1);
        aac_bits_init(&b, p, avail);
        ret = latm_parse_element(d, &b, mux_config_present);
        if (ret == AAC_CORE_NEED_MORE) {
            if (avail > 65536) {
                d->bad_frames++;
                return AAC_DEC_ERROR;
            }
            return AAC_DEC_DEFAULT;
        }
        if (ret < 0) {
            d->bad_frames++;
            d->fatal = 1;
            return AAC_DEC_ERROR;
        }
        consumed = aac_bits_bytes_consumed(&b);
        if (consumed == 0)
            consumed = 1;
        buffer_consume(&d->in, consumed);
    }

    return emit_decoded(d, NULL, 0, &payload_len, info, frame);
}

/* ------------------------------------------------------------------------- */
/* Transport resolution                                                      */
/* ------------------------------------------------------------------------- */

static int process_records(aac_decoder* d)
{
    int i;
    for (i = 0; i < d->record_count; i++) {
        aac_send_record* r = &d->records[i];
        if (d->transport == AAC_TRANSPORT_ADTS ||
            d->transport == AAC_TRANSPORT_LATM) {
            if (scan_push(d, d->in.data + r->start, r->len, r->pts, r->dts,
                          r->opaque) < 0)
                return -1;
        } else {
            if (ts_push(d, r->pts, r->dts, r->opaque) < 0)
                return -1;
        }
    }
    free(d->records);
    d->records = NULL;
    d->record_count = 0;
    return 0;
}

static int resolve_transport(aac_decoder* d)
{
    int t;

    if (d->opened)
        return 1;

    t = d->param.transport;
    if (t == AAC_TRANSPORT_AUTO) {
        if (d->asc_set) {
            t = AAC_TRANSPORT_RAW;
        } else {
            size_t avail = buffer_avail(&d->in);
            int detected = probe_transport(buffer_ptr(&d->in), avail);
            if (detected == 0)
                return 0;
            if (detected < 0) {
                d->fatal = 1;
                return -1;
            }
            t = detected;
        }
    }

    d->transport = t;
    d->opened = 1;

    if (t == AAC_TRANSPORT_RAW) {
        if (d->asc_set &&
            aac_core_set_asc(&d->core, d->asc, (size_t)d->asc_len) < 0) {
            d->fatal = 1;
            return -1;
        }
    } else if (t == AAC_TRANSPORT_LATM_MCP0) {
        if (d->asc_set) {
            aac_bits b;
            int ret;
            aac_bits_init(&b, d->asc, (size_t)d->asc_len);
            ret = latm_parse_stream_mux_config(d, &b);
            if (ret < 0) {
                d->fatal = 1;
                return -1;
            }
        }
    } else if (t == AAC_TRANSPORT_LATM || t == AAC_TRANSPORT_LATM_MCP1) {
        if (d->asc_set) {
            aac_bits b;
            aac_bits_init(&b, d->asc, (size_t)d->asc_len);
            (void)latm_parse_stream_mux_config(d, &b);
        }
    }

    if (d->record_count > 0 && process_records(d) < 0) {
        d->fatal = 1;
        return -1;
    }
    return 1;
}

/* ------------------------------------------------------------------------- */
/* Public API                                                                */
/* ------------------------------------------------------------------------- */

static int valid_transport(int t)
{
    return t == AAC_TRANSPORT_AUTO || t == AAC_TRANSPORT_RAW ||
           t == AAC_TRANSPORT_ADTS || t == AAC_TRANSPORT_ADIF ||
           t == AAC_TRANSPORT_LATM || t == AAC_TRANSPORT_LATM_MCP0 ||
           t == AAC_TRANSPORT_LATM_MCP1;
}

void* aac_dec_open(const aac_param_t* param)
{
    aac_decoder* d;
    aac_param_t defaults;

    memset(&defaults, 0, sizeof(defaults));
    d = (aac_decoder*)calloc(1, sizeof(*d));
    if (d == NULL)
        return NULL;

    if (param != NULL) {
        d->param = *param;
        if (param->config_len < 0 ||
            (param->config_len > 0 && param->config == NULL) ||
            !valid_transport(param->transport) ||
            param->output_format < AAC_OUTPUT_S16 ||
            param->output_format > AAC_OUTPUT_F32_PLANAR ||
            param->channel_order < AAC_CHANNEL_ORDER_WAV ||
            param->channel_order > AAC_CHANNEL_ORDER_MPEG ||
            (param->flags & ~AAC_DEC_FLAG_DISABLE_CONCEALMENT) != 0) {
            free(d);
            return NULL;
        }
    } else {
        d->param = defaults;
    }
    d->param.config = NULL;
    d->param.config_len = 0;

    if (aac_core_init(&d->core) < 0) {
        free(d);
        return NULL;
    }
    d->core.mpeg_order =
        (d->param.channel_order == AAC_CHANNEL_ORDER_MPEG) ? 1 : 0;
    d->pcm_scratch = (float*)malloc((size_t)AAC_DEC_MAX_CHANNELS *
                                    AAC_DEC_MAX_FRAME_SIZE * sizeof(float));
    d->conceal_pcm = (float*)malloc((size_t)AAC_DEC_MAX_CHANNELS *
                                    AAC_DEC_MAX_FRAME_SIZE * sizeof(float));
    if (d->pcm_scratch == NULL || d->conceal_pcm == NULL) {
        free(d->pcm_scratch);
        free(d->conceal_pcm);
        aac_core_close(&d->core);
        free(d);
        return NULL;
    }

    if (param != NULL && param->config_len > 0) {
        d->asc = (uint8_t*)malloc((size_t)param->config_len);
        if (d->asc == NULL) {
            free(d->pcm_scratch);
            free(d->conceal_pcm);
            aac_core_close(&d->core);
            free(d);
            return NULL;
        }
        memcpy(d->asc, param->config, (size_t)param->config_len);
        d->asc_len = param->config_len;
        d->asc_set = 1;
    }

    if (d->param.transport != AAC_TRANSPORT_AUTO || d->asc_set) {
        if (resolve_transport(d) < 0) {
            aac_dec_close(d);
            return NULL;
        }
    }
    return d;
}

int aac_dec_set_config(void* decoder, const void* config, int len)
{
    aac_decoder* d = (aac_decoder*)decoder;

    if (d == NULL || d->fatal)
        return -1;
    if (len < 0 || (len > 0 && config == NULL))
        return -1;
    if (d->started || d->frames_out > 0 || d->has_pending)
        return -1;

    free(d->asc);
    d->asc = NULL;
    d->asc_len = 0;
    d->asc_set = 0;
    if (len > 0) {
        d->asc = (uint8_t*)malloc((size_t)len);
        if (d->asc == NULL)
            return -1;
        memcpy(d->asc, config, (size_t)len);
        d->asc_len = len;
        d->asc_set = 1;
    }

    if (!d->opened) {
        if (d->param.transport == AAC_TRANSPORT_AUTO && !d->asc_set)
            return -1;
        if (resolve_transport(d) < 0)
            return -1;
        return 0;
    }

    if (d->transport == AAC_TRANSPORT_RAW) {
        if (aac_core_set_asc(&d->core, d->asc, (size_t)d->asc_len) < 0)
            return -1;
        return 0;
    }
    if (d->transport == AAC_TRANSPORT_LATM_MCP0 ||
        d->transport == AAC_TRANSPORT_LATM ||
        d->transport == AAC_TRANSPORT_LATM_MCP1) {
        aac_bits b;
        aac_bits_init(&b, d->asc, (size_t)d->asc_len);
        if (latm_parse_stream_mux_config(d, &b) < 0)
            return -1;
        return 0;
    }
    return -1;
}

int aac_dec_send_packet(void* decoder, const aac_packet_t* packet)
{
    aac_decoder* d = (aac_decoder*)decoder;
    size_t start;

    if (d == NULL || d->fatal)
        return AAC_DEC_ERROR;
    if (packet == NULL)
        return AAC_DEC_ERROR;
    if (packet->len < 0 || (packet->len > 0 && packet->data == NULL))
        return AAC_DEC_ERROR;

    d->flushed = 0;
    if (packet->len == 0)
        return 0;

    start = d->in.len;
    if (buffer_append(d, packet->data, (size_t)packet->len) < 0) {
        d->fatal = 1;
        return AAC_DEC_ERROR;
    }

    if (!d->opened) {
        aac_send_record* rec;
        aac_send_record* p = (aac_send_record*)realloc(
            d->records, (size_t)(d->record_count + 1) * sizeof(*p));
        int r;
        if (p == NULL) {
            d->fatal = 1;
            return AAC_DEC_ERROR;
        }
        d->records = p;
        rec = &d->records[d->record_count++];
        rec->start = start;
        rec->len = (size_t)packet->len;
        rec->pts = packet->pts;
        rec->dts = packet->dts;
        rec->opaque = packet->opaque;

        r = resolve_transport(d);
        if (r < 0)
            return AAC_DEC_ERROR;
        if (r == 0)
            return 0;
        d->started = 1;
        return 0;
    }

    d->started = 1;
    if (d->transport == AAC_TRANSPORT_ADTS ||
        d->transport == AAC_TRANSPORT_LATM) {
        if (scan_push(d, packet->data, (size_t)packet->len, packet->pts,
                      packet->dts, packet->opaque) < 0) {
            d->fatal = 1;
            return AAC_DEC_ERROR;
        }
    } else {
        if (ts_push(d, packet->pts, packet->dts, packet->opaque) < 0) {
            d->fatal = 1;
            return AAC_DEC_ERROR;
        }
    }
    return 0;
}

static int decode_next(aac_decoder* d, aac_stream_info* info,
                       aac_frame_t* frame)
{
    int ret;

    if (d->has_pending) {
        *frame = d->pending_frame;
        memset(&d->pending_frame, 0, sizeof(d->pending_frame));
        d->has_pending = 0;
        if (info != NULL)
            *info = d->info;
        return AAC_DEC_GOT_FRAME;
    }

    if (d->flushed)
        return AAC_DEC_END;

    if (!d->opened) {
        ret = resolve_transport(d);
        if (ret < 0)
            return AAC_DEC_ERROR;
        if (ret == 0)
            return AAC_DEC_DEFAULT;
    }

    switch (d->transport) {
    case AAC_TRANSPORT_ADTS:
        return decode_adts(d, info, frame);
    case AAC_TRANSPORT_ADIF:
        return decode_adif(d, info, frame);
    case AAC_TRANSPORT_LATM:
        return decode_loas(d, info, frame);
    case AAC_TRANSPORT_RAW:
    case AAC_TRANSPORT_LATM_MCP0:
    case AAC_TRANSPORT_LATM_MCP1:
        return decode_packetized(d, info, frame);
    default:
        return AAC_DEC_ERROR;
    }
}

int aac_dec_recv_frame(void* decoder, aac_stream_info* info,
                       aac_frame_t* frame)
{
    aac_decoder* d = (aac_decoder*)decoder;

    if (d == NULL || frame == NULL)
        return AAC_DEC_ERROR;
    if (d->fatal)
        return AAC_DEC_ERROR;
    memset(frame, 0, sizeof(*frame));
    return decode_next(d, info, frame);
}

int aac_dec_flush(void* decoder, aac_stream_info* info, aac_frame_t* frame)
{
    aac_decoder* d = (aac_decoder*)decoder;
    int ret;

    if (d == NULL || frame == NULL)
        return AAC_DEC_ERROR;
    if (d->fatal)
        return AAC_DEC_ERROR;
    memset(frame, 0, sizeof(*frame));

    if (d->has_pending) {
        *frame = d->pending_frame;
        memset(&d->pending_frame, 0, sizeof(d->pending_frame));
        d->has_pending = 0;
        if (info != NULL)
            *info = d->info;
        return AAC_DEC_GOT_FRAME;
    }

    if (!d->opened) {
        ret = resolve_transport(d);
        if (ret <= 0) {
            d->flushed = 1;
            return AAC_DEC_END;
        }
    }

    ret = decode_next(d, info, frame);
    if (ret == AAC_DEC_GOT_FRAME || ret == AAC_DEC_GOT_HEADER)
        return ret;
    d->flushed = 1;
    if (ret == AAC_DEC_ERROR)
        return AAC_DEC_ERROR;
    return AAC_DEC_END;
}

int aac_dec_conceal_frame(void* decoder, aac_stream_info* info,
                          aac_frame_t* frame)
{
    aac_decoder* d = (aac_decoder*)decoder;
    int ret;

    if (d == NULL || frame == NULL || d->fatal)
        return AAC_DEC_ERROR;
    memset(frame, 0, sizeof(*frame));
    if (d->param.flags & AAC_DEC_FLAG_DISABLE_CONCEALMENT)
        return AAC_DEC_ERROR;
    if (!d->info_valid || !d->conceal_valid)
        return AAC_DEC_ERROR;

    ret = emit_concealed(d, info, frame);
    d->bad_frames++;
    return ret;
}

void aac_dec_frame_unref(void* decoder, aac_frame_t* frame)
{
    (void)decoder;
    frame_release(frame);
}

int aac_dec_get_stream_info(void* decoder, aac_stream_info* info)
{
    aac_decoder* d = (aac_decoder*)decoder;
    if (d == NULL || info == NULL || !d->info_valid)
        return -1;
    *info = d->info;
    return 0;
}

void aac_dec_close(void* decoder)
{
    aac_decoder* d = (aac_decoder*)decoder;
    if (d == NULL)
        return;
    aac_core_close(&d->core);
    frame_release(&d->pending_frame);
    free(d->in.data);
    free(d->asc);
    free(d->scan.buf);
    ts_free(&d->ts);
    free(d->records);
    free(d->au_buf);
    free(d->latm_au_buf);
    free(d->pcm_scratch);
    free(d->conceal_pcm);
    free(d);
}

const char* aac_dec_version(void)
{
    return "aacdec 1.7.0 (clean-room AAC-LC / HE-AAC / HE-AACv2 / AAC-LD / "
           "AAC-ELD / USAC decoder)";
}

const char* aac_dec_strerror(int code)
{
    switch (code) {
    case AAC_DEC_ERROR:
        return "decoding error";
    case AAC_DEC_DEFAULT:
        return "no output available";
    case AAC_DEC_GOT_FRAME:
        return "frame available";
    case AAC_DEC_GOT_HEADER:
        return "configuration available";
    case AAC_DEC_END:
        return "end of stream";
    default:
        return "unknown error";
    }
}
