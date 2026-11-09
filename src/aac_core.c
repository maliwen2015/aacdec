/*
 * aacdec - clean-room MPEG-2 / MPEG-4 AAC-LC / ER AAC-LC / AAC-LD / AAC-ELD
 * / USAC decoder core.
 *
 * The decoder implements the raw_data_block syntax, noiseless decoding,
 * inverse quantization, M/S and intensity stereo, PNS, TNS, the IMDCT
 * filterbank with sine/KBD windows and the standard channel configurations,
 * plus the low-delay filterbank used by AAC-LD and AAC-ELD and the USAC
 * (xHE-AAC) arithmetic-coded spectral path.  LD-SBR, ELDv2, MPS212, LPD and
 * time warping are not part of this core; the transport layer reports
 * unsupported extensions as errors.
 */
#include "aac_core.h"
#include "aac_eld_tables.h"
#include "aac_eld_ld_tables.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>


#include "aac_bits.h"
#include "aac_huffman.h"
#include "aac_tables.h"
#include "aac_tables_data.h"

#define AAC_PI 3.14159265358979323846

#ifdef AAC_USAC_DEBUG
#include <stdio.h>
#define USDBG(...) do { fprintf(stderr, __VA_ARGS__); } while (0)
#else
#define USDBG(...)
#endif

enum
{
    ONLY_LONG_SEQUENCE = 0,
    LONG_START_SEQUENCE = 1,
    EIGHT_SHORT_SEQUENCE = 2,
    LONG_STOP_SEQUENCE = 3
};

enum
{
    ID_SCE = 0,
    ID_CPE = 1,
    ID_CCE = 2,
    ID_LFE = 3,
    ID_DSE = 4,
    ID_PCE = 5,
    ID_FIL = 6,
    ID_END = 7
};

enum
{
    ZERO_HCB = 0,
    NOISE_HCB = 13,
    INTENSITY_HCB2 = 14,
    INTENSITY_HCB = 15
};

enum
{
    ACT_NONE = 0x00,
    ACT_FRONT = 0x01,
    ACT_SIDE = 0x02,
    ACT_BACK = 0x03,
    ACT_LFE = 0x04,
    ACT_FRONT_TOP = 0x11,
    ACT_SIDE_TOP = 0x12,
    ACT_BACK_TOP = 0x13,
    ACT_FRONT_BOTTOM = 0x21,
    ACT_SIDE_BOTTOM = 0x22,
    ACT_BACK_BOTTOM = 0x23
};

#define ELEMENT_SCE 0
#define ELEMENT_CPE 1
#define ELEMENT_LFE 3

/* ------------------------------------------------------------------------- */
/* IMDCT / filterbank                                                        */
/* ------------------------------------------------------------------------- */

static float cos_long_1024[8192];
static float cos_long_960[7680];
static float cos_long_512[4096];
static float cos_long_480[3840];
static float cos_short_128[1024];
static float cos_short_120[960];
static float ld_low_half_512[512];
static float ld_low_half_480[480];
static int synth_ready;

static void init_cos_table(float* table, int M)
{
    int i;
    int n = 8 * M;
    for (i = 0; i < n; i++)
        table[i] = (float)cos(AAC_PI * (double)i / (double)(4 * M));
}

/* Low-overlap AAC-LD window (128-sample overlap).  The returned table is the
 * rising half of the full 2*M window; the falling half is its mirror. */
static void init_ld_low_window(float* table, int M)
{
    int i;
    const int overlap = 128;
    for (i = 0; i < M; i++) {
        if (i < overlap)
            table[i] = (float)sin(3.14159265358979323846 *
                                  ((double)i + 0.5) / (double)(2 * overlap));
        else
            table[i] = 1.0f;
    }
}

static void synth_init(void)
{
    if (synth_ready)
        return;
    init_cos_table(cos_long_1024, 1024);
    init_cos_table(cos_long_960, 960);
    init_cos_table(cos_long_512, 512);
    init_cos_table(cos_long_480, 480);
    init_cos_table(cos_short_128, 128);
    init_cos_table(cos_short_120, 120);
    init_ld_low_window(ld_low_half_512, 512);
    init_ld_low_window(ld_low_half_480, 480);
    synth_ready = 1;
}

static const float* cos_table_for(int M)
{
    if (M == 1024)
        return cos_long_1024;
    if (M == 960)
        return cos_long_960;
    if (M == 512)
        return cos_long_512;
    if (M == 480)
        return cos_long_480;
    if (M == 128)
        return cos_short_128;
    return cos_short_120;
}

/* Direct DCT-IV: y[m] = sum_k in[k] * cos(pi/M*(m+0.5)*(k+0.5)). */
static void imdct_dct4(const float* in, int M, float* y)
{
    const float* tab = cos_table_for(M);
    const int mod = 8 * M;
    const float scale = 1.0f / (float)M;
    int m;

    for (m = 0; m < M; m++) {
        int idx = 2 * m + 1;
        int step = 4 * m + 2;
        float sum = 0.0f;
        int k;

        for (k = 0; k < M; k++) {
            sum += in[k] * tab[idx];
            idx += step;
            if (idx >= mod)
                idx -= mod;
        }
        y[m] = sum * scale;
    }
}

/* Full IMDCT of length 2M from the DCT-IV extension. */
static void imdct_full(const float* in, int M, float* out)
{
    float y[1024];
    int half = M / 2;
    int n;

    imdct_dct4(in, M, y);
    for (n = 0; n < 2 * M; n++) {
        int m = n + half;
        if (m < M)
            out[n] = y[m];
        else if (m < 2 * M)
            out[n] = -y[2 * M - 1 - m];
        else
            out[n] = -y[m - 2 * M];
    }
}

static const float* ld_window_half(int frame_size, int shape){
    if (shape == 0)
        return aac_window_half(frame_size, 0, 0);
    return (frame_size == 480) ? ld_low_half_480 : ld_low_half_512;
}

static void synth_channel(aac_channel_state* st, const float* spec, int seq,
                          int shape, int frame_size, float* out, int ld)
{
    const int sh = frame_size / 8;      /* short half length (128/120) */
    const int short_off = 7 * sh / 2;   /* 448/420 */
    const int prev_shape = st->initialized ? st->prev_window_shape : shape;
    const float* long_prev = ld ? ld_window_half(frame_size, prev_shape)
                                : aac_window_half(frame_size, 0, prev_shape);
    const float* long_cur = ld ? ld_window_half(frame_size, shape)
                               : aac_window_half(frame_size, 0, shape);
    const float* short_prev = aac_window_half(frame_size, 1, prev_shape);
    const float* short_cur = aac_window_half(frame_size, 1, shape);
    int i;

    if (!ld && seq == EIGHT_SHORT_SEQUENCE) {
        float y[2048];
        int w;

        memset(y, 0, sizeof(y));
        for (w = 0; w < 8; w++) {
            float t[256];
            imdct_full(spec + w * sh, sh, t);
            for (i = 0; i < sh; i++) {
                float wl = (w == 0) ? short_prev[i] : short_cur[i];
                y[short_off + w * sh + i] += t[i] * wl;
                y[short_off + w * sh + sh + i] +=
                    t[sh + i] * short_cur[sh - 1 - i];
            }
        }
        for (i = 0; i < frame_size; i++)
            out[i] = st->overlap[i] + y[i];
        memcpy(st->overlap, y + frame_size, (size_t)frame_size * sizeof(float));
    } else {
        float t[2048];
        float win[2048];

        imdct_full(spec, frame_size, t);

        if (seq == ONLY_LONG_SEQUENCE) {
            for (i = 0; i < frame_size; i++) {
                win[i] = long_prev[i];
                win[frame_size + i] = long_cur[frame_size - 1 - i];
            }
        } else if (seq == LONG_START_SEQUENCE) {
            for (i = 0; i < frame_size; i++)
                win[i] = long_prev[i];
            for (i = 0; i < short_off; i++)
                win[frame_size + i] = 1.0f;
            for (i = 0; i < sh; i++)
                win[frame_size + short_off + i] =
                    short_cur[sh - 1 - i];
            for (i = 0; i < short_off; i++)
                win[frame_size + short_off + sh + i] = 0.0f;
        } else { /* LONG_STOP_SEQUENCE */
            for (i = 0; i < short_off; i++)
                win[i] = 0.0f;
            for (i = 0; i < sh; i++)
                win[short_off + i] = short_prev[i];
            for (i = 0; i < short_off; i++)
                win[short_off + sh + i] = 1.0f;
            for (i = 0; i < frame_size; i++)
                win[frame_size + i] = long_cur[frame_size - 1 - i];
        }

        for (i = 0; i < 2 * frame_size; i++)
            t[i] *= win[i];
        for (i = 0; i < frame_size; i++)
            out[i] = st->overlap[i] + t[i];
        memcpy(st->overlap, t + frame_size,
               (size_t)frame_size * sizeof(float));
    }
    st->prev_window_shape = shape;
    st->initialized = 1;
}

/* AAC-ELD low-delay synthesis filterbank.  The transform is mapped to a
 * conventional IMDCT; the windowing and overlap-add follow the low-delay
 * structure defined by ISO/IEC 14496-3. */
static void synth_channel_eld(aac_channel_state* st, const float* spec,
                              int frame_size, float* out)
{
    const int n = frame_size;
    const int n2 = n >> 1;
    const int n4 = n >> 2;
    const float* fb =
        (n == 480) ? aac_eld_ld_window_480 : aac_eld_ld_window_512;
    float x[512];
    float* z = st->eld_z[0];
    int i;

    imdct_dct4(spec, n, x);

    for (i = 0; i < n4; i++) {
        float z2 = x[n2 + i];
        float z0 = z2 + 0.25f * z[n2 + i] * fb[2 * n + i];
        z[n2 + i] = x[n2 - 1 - i] +
                    0.25f * z[n + i] * fb[2 * n + n2 + i];
        out[3 * n4 - 1 - i] = z[n2 + i] * fb[n + n2 - 1 - i] +
                              z[i] * fb[n + n2 + i];
        z[i] = z0;
        z[n + i] = z2;
    }
    for (i = n4; i < n2; i++) {
        float z2 = x[n2 + i];
        float z0 = z2 + 0.25f * z[n2 + i] * fb[2 * n + i];
        z[n2 + i] = x[n2 - 1 - i] +
                    0.25f * z[n + i] * fb[2 * n + n2 + i];
        out[i - n4] = 2.0f * z[n2 + i] * fb[n2 - 1 - i] +
                      2.0f * z[i] * fb[n2 + i];
        out[3 * n4 - 1 - i] = z[n2 + i] * fb[n + n2 - 1 - i] +
                              z[i] * fb[n + n2 + i];
        z[i] = z0;
        z[n + i] = z2;
    }
    for (i = 0; i < n4; i++)
        out[3 * n4 + i] = 2.0f * z[i] * fb[n2 + i];
}

/* ------------------------------------------------------------------------- */
/* Channel layout                                                            */
/* ------------------------------------------------------------------------- */

static void set_element_channels(aac_config* c, int ch, int type, int index)
{
    if (ch < AAC_CORE_MAX_CHANNELS) {
        c->channel_types[ch] = type;
        c->channel_indices[ch] = index;
    }
}

static void set_standard_layout(aac_config* c, int cc)
{
    int i;

    c->pce_valid = 1;
    c->channels = 0;
    for (i = 0; i < AAC_CORE_MAX_CHANNELS; i++) {
        c->channel_types[i] = ACT_NONE;
        c->channel_indices[i] = 0;
        c->output_order[i] = i;
    }

    switch (cc) {
    case 1: /* mono */
        c->channels = 1;
        set_element_channels(c, 0, ACT_FRONT, 0);
        c->output_order[0] = 0;
        break;
    case 2: /* stereo */
        c->channels = 2;
        set_element_channels(c, 0, ACT_FRONT, 0);
        set_element_channels(c, 1, ACT_FRONT, 1);
        c->output_order[0] = 0;
        c->output_order[1] = 1;
        break;
    case 3: /* C, L, R */
        c->channels = 3;
        set_element_channels(c, 0, ACT_FRONT, 0);
        set_element_channels(c, 1, ACT_FRONT, 1);
        set_element_channels(c, 2, ACT_FRONT, 2);
        c->output_order[0] = 1;
        c->output_order[1] = 2;
        c->output_order[2] = 0;
        break;
    case 4: /* C, L, R, Cs */
        c->channels = 4;
        set_element_channels(c, 0, ACT_FRONT, 0);
        set_element_channels(c, 1, ACT_FRONT, 1);
        set_element_channels(c, 2, ACT_FRONT, 2);
        set_element_channels(c, 3, ACT_BACK, 0);
        c->output_order[0] = 1;
        c->output_order[1] = 2;
        c->output_order[2] = 0;
        c->output_order[3] = 3;
        break;
    case 5: /* C, L, R, Ls, Rs */
        c->channels = 5;
        set_element_channels(c, 0, ACT_FRONT, 0);
        set_element_channels(c, 1, ACT_FRONT, 1);
        set_element_channels(c, 2, ACT_FRONT, 2);
        set_element_channels(c, 3, ACT_SIDE, 0);
        set_element_channels(c, 4, ACT_SIDE, 1);
        c->output_order[0] = 1;
        c->output_order[1] = 2;
        c->output_order[2] = 0;
        c->output_order[3] = 3;
        c->output_order[4] = 4;
        break;
    case 6: /* C, L, R, Ls, Rs, LFE */
        c->channels = 6;
        set_element_channels(c, 0, ACT_FRONT, 0);
        set_element_channels(c, 1, ACT_FRONT, 1);
        set_element_channels(c, 2, ACT_FRONT, 2);
        set_element_channels(c, 3, ACT_SIDE, 0);
        set_element_channels(c, 4, ACT_SIDE, 1);
        set_element_channels(c, 5, ACT_LFE, 0);
        c->output_order[0] = 1;
        c->output_order[1] = 2;
        c->output_order[2] = 0;
        c->output_order[3] = 5;
        c->output_order[4] = 3;
        c->output_order[5] = 4;
        break;
    case 7: /* C, L, R, Ls, Rs, Lrs, Rrs, LFE */
        c->channels = 8;
        set_element_channels(c, 0, ACT_FRONT, 0);
        set_element_channels(c, 1, ACT_FRONT, 1);
        set_element_channels(c, 2, ACT_FRONT, 2);
        set_element_channels(c, 3, ACT_SIDE, 0);
        set_element_channels(c, 4, ACT_SIDE, 1);
        set_element_channels(c, 5, ACT_BACK, 0);
        set_element_channels(c, 6, ACT_BACK, 1);
        set_element_channels(c, 7, ACT_LFE, 0);
        c->output_order[0] = 1;
        c->output_order[1] = 2;
        c->output_order[2] = 0;
        c->output_order[3] = 7;
        c->output_order[4] = 3;
        c->output_order[5] = 4;
        c->output_order[6] = 5;
        c->output_order[7] = 6;
        break;
    default:
        c->channels = 0;
        c->pce_valid = 0;
        break;
    }
}

static int wav_rank(int type, int index)
{
    switch (type) {
    case ACT_FRONT:
        if (index == 1)
            return 10;
        if (index == 2)
            return 11;
        if (index == 0)
            return 12;
        return 13 + index;
    case ACT_LFE:
        return 20 + index;
    case ACT_SIDE:
        return 30 + index;
    case ACT_BACK:
        return 40 + index;
    case ACT_FRONT_TOP:
        return 50 + index;
    case ACT_SIDE_TOP:
        return 55 + index;
    case ACT_BACK_TOP:
        return 60 + index;
    case ACT_FRONT_BOTTOM:
        return 65 + index;
    case ACT_SIDE_BOTTOM:
        return 70 + index;
    case ACT_BACK_BOTTOM:
        return 75 + index;
    default:
        return 100 + index;
    }
}

static void build_wav_order(aac_config* c)
{
    int i, j;

    for (i = 0; i < c->channels; i++)
        c->output_order[i] = i;

    for (i = 1; i < c->channels; i++) {
        int key = c->output_order[i];
        int key_rank = wav_rank(c->channel_types[key], c->channel_indices[key]);
        j = i - 1;
        while (j >= 0) {
            int prev = c->output_order[j];
            int prev_rank =
                wav_rank(c->channel_types[prev], c->channel_indices[prev]);
            if (prev_rank <= key_rank)
                break;
            c->output_order[j + 1] = prev;
            j--;
        }
        c->output_order[j + 1] = key;
    }
}

/* ------------------------------------------------------------------------- */
/* PCE                                                                       */
/* ------------------------------------------------------------------------- */

static int read_pce(aac_core* core, aac_bits* b, int from_asc)
{
    aac_config* c = &core->cfg;
    int tag, object_type, sf_index;
    int num_front, num_side, num_back, num_lfe;
    int num_assoc, num_cc;
    int i;
    int front_sce_count = 0;
    int ch = 0;

    tag = (int)aac_bits_read(b, 4);
    object_type = (int)aac_bits_read(b, 2);
    sf_index = (int)aac_bits_read(b, 4);
    num_front = (int)aac_bits_read(b, 4);
    num_side = (int)aac_bits_read(b, 4);
    num_back = (int)aac_bits_read(b, 4);
    num_lfe = (int)aac_bits_read(b, 2);
    num_assoc = (int)aac_bits_read(b, 3);
    num_cc = (int)aac_bits_read(b, 4);

    if (aac_bits_read(b, 1))          /* mono mixdown */
        aac_bits_read(b, 4);
    if (aac_bits_read(b, 1))          /* stereo mixdown */
        aac_bits_read(b, 4);
    if (aac_bits_read(b, 1)) {        /* matrix mixdown */
        aac_bits_read(b, 2);
        aac_bits_read(b, 1);
    }

    if (b->error)
        return AAC_CORE_NEED_MORE;

    if (!from_asc || c->sample_rate_index < 0) {
        if (sf_index < AAC_SR_INDEX_MAX) {
            c->sample_rate_index = sf_index;
            c->sample_rate = aac_sample_rates[sf_index];
        }
    }

    if (c->core_aot == 0) {
        c->core_aot = object_type + 1;
        c->aot = c->core_aot;
        if (object_type == 1) {
            c->mpeg2 = 1;
            c->profile = 1;
        }
    }

    if (c->core_aot != 2) {
        /* Only AAC-LC is implemented by this core. */
        return AAC_CORE_ERROR;
    }

    /* Read the element lists. */
    c->pce_element_count = 0;
    for (i = 0; i < num_front; i++) {
        int is_cpe = (int)aac_bits_read(b, 1);
        int etag = (int)aac_bits_read(b, 4);
        int idx = c->pce_element_count++;
        c->pce_elements[idx].type = is_cpe ? ELEMENT_CPE : ELEMENT_SCE;
        c->pce_elements[idx].is_cpe = is_cpe;
        c->pce_elements[idx].tag = etag;
        if (!is_cpe)
            front_sce_count++;
    }
    for (i = 0; i < num_side; i++) {
        int is_cpe = (int)aac_bits_read(b, 1);
        int etag = (int)aac_bits_read(b, 4);
        int idx = c->pce_element_count++;
        c->pce_elements[idx].type = is_cpe ? ELEMENT_CPE : ELEMENT_SCE;
        c->pce_elements[idx].is_cpe = is_cpe;
        c->pce_elements[idx].tag = etag;
    }
    for (i = 0; i < num_back; i++) {
        int is_cpe = (int)aac_bits_read(b, 1);
        int etag = (int)aac_bits_read(b, 4);
        int idx = c->pce_element_count++;
        c->pce_elements[idx].type = is_cpe ? ELEMENT_CPE : ELEMENT_SCE;
        c->pce_elements[idx].is_cpe = is_cpe;
        c->pce_elements[idx].tag = etag;
    }
    for (i = 0; i < num_lfe; i++) {
        int etag = (int)aac_bits_read(b, 4);
        int idx = c->pce_element_count++;
        c->pce_elements[idx].type = ELEMENT_LFE;
        c->pce_elements[idx].is_cpe = 0;
        c->pce_elements[idx].tag = etag;
    }
    for (i = 0; i < num_assoc; i++)
        aac_bits_read(b, 4);
    for (i = 0; i < num_cc; i++) {
        aac_bits_read(b, 1);
        aac_bits_read(b, 4);
    }

    aac_bits_align(b);
    if (b->error)
        return AAC_CORE_NEED_MORE;
    {
        int comment_bytes = (int)aac_bits_read(b, 8);
        if (b->error)
            return AAC_CORE_NEED_MORE;
        aac_bits_skip(b, (size_t)comment_bytes * 8u);
    }
    if (b->error)
        return AAC_CORE_NEED_MORE;

    if (c->pce_element_count > AAC_CORE_MAX_CHANNELS)
        return AAC_CORE_ERROR;

    /* Assign speaker positions in element order. */
    {
        int front_sce_seen = 0;
        int front_cpe_idx = (front_sce_count > 0) ? 1 : 0;
        int side_sce_idx = 0, side_cpe_idx = 0;
        int back_sce_idx = 0, back_cpe_idx = 0;
        int lfe_idx = 0;

        ch = 0;
        for (i = 0; i < c->pce_element_count; i++) {
            int is_front = (i < num_front);
            int is_side = (i >= num_front && i < num_front + num_side);

            if (c->pce_elements[i].type == ELEMENT_CPE) {
                int t, idx;
                if (is_front) {
                    t = ACT_FRONT;
                    idx = front_cpe_idx;
                    front_cpe_idx += 2;
                } else if (is_side) {
                    t = ACT_SIDE;
                    idx = side_cpe_idx;
                    side_cpe_idx += 2;
                } else {
                    t = ACT_BACK;
                    idx = back_cpe_idx;
                    back_cpe_idx += 2;
                }
                c->pce_elements[i].channel_type[0] = t;
                c->pce_elements[i].channel_index[0] = idx;
                c->pce_elements[i].channel_type[1] = t;
                c->pce_elements[i].channel_index[1] = idx + 1;
                set_element_channels(c, ch, t, idx);
                set_element_channels(c, ch + 1, t, idx + 1);
                ch += 2;
            } else if (c->pce_elements[i].type == ELEMENT_SCE) {
                int t, idx;
                if (is_front) {
                    t = ACT_FRONT;
                    if (front_sce_seen == 0 && front_sce_count > 0) {
                        idx = 0;
                        front_sce_seen = 1;
                    } else {
                        idx = 3 + front_sce_seen;
                        front_sce_seen++;
                    }
                } else if (is_side) {
                    t = ACT_SIDE;
                    idx = side_sce_idx++;
                } else {
                    t = ACT_BACK;
                    idx = back_sce_idx++;
                }
                c->pce_elements[i].channel_type[0] = t;
                c->pce_elements[i].channel_index[0] = idx;
                set_element_channels(c, ch, t, idx);
                ch += 1;
            } else { /* LFE */
                c->pce_elements[i].channel_type[0] = ACT_LFE;
                c->pce_elements[i].channel_index[0] = lfe_idx++;
                set_element_channels(c, ch, ACT_LFE, lfe_idx - 1);
                ch += 1;
            }
        }
    }

    if (b->error)
        return AAC_CORE_NEED_MORE;

    c->channels = ch;
    c->channel_config = 0;
    c->pce_valid = 1;
    build_wav_order(c);
    (void)tag;
    return AAC_CORE_OK;
}

/* ------------------------------------------------------------------------- */
/* Small syntax helpers                                                      */
/* ------------------------------------------------------------------------- */

static int codebook_dim(int cb)
{
    return (cb <= 4) ? 4 : 2;
}

static int codebook_is_unsigned(int cb)
{
    return (cb >= 1 && cb <= 11) && !(cb <= 2) && !(cb == 5 || cb == 6);
}

static void codebook_vector(int cb, int idx, int* v)
{
    int i;

    if (cb <= 2) {
        for (i = 3; i >= 0; i--) {
            v[i] = (idx % 3) - 1;
            idx /= 3;
        }
    } else if (cb <= 4) {
        for (i = 3; i >= 0; i--) {
            v[i] = idx % 3;
            idx /= 3;
        }
    } else if (cb <= 6) {
        v[1] = (idx % 9) - 4;
        v[0] = (idx / 9) - 4;
    } else if (cb <= 8) {
        v[1] = idx % 8;
        v[0] = idx / 8;
    } else if (cb <= 10) {
        v[1] = idx % 13;
        v[0] = idx / 13;
    } else {
        v[1] = idx % 17;
        v[0] = idx / 17;
    }
}

static int decode_ics_info(aac_core* core, aac_bits* b, aac_ics* ics,
                           int common_window)
{
    aac_config* c = &core->cfg;
    int short_window;
    int eld = (c->aot == 39);
    int i;

    if (eld) {
        ics->window_sequence = ONLY_LONG_SEQUENCE;
        ics->window_shape = 0;
    } else if (!common_window) {
        int reserved = (int)aac_bits_read(b, 1);
        (void)reserved;
        ics->window_sequence = (int)aac_bits_read(b, 2);
        ics->window_shape = (int)aac_bits_read(b, 1);
    }

    if (ics->window_sequence < ONLY_LONG_SEQUENCE ||
        ics->window_sequence > LONG_STOP_SEQUENCE)
        return AAC_CORE_ERROR;
    if ((c->aot == 23 || eld) &&
        ics->window_sequence != ONLY_LONG_SEQUENCE)
        return AAC_CORE_ERROR;

    short_window = (ics->window_sequence == EIGHT_SHORT_SEQUENCE);
    ics->num_windows = short_window ? 8 : 1;
    ics->num_groups = 1;
    ics->group_len[0] = 1;

    if (short_window) {
        ics->max_sfb = (int)aac_bits_read(b, 4);
        for (i = 0; i < 7; i++) {
            if (aac_bits_read(b, 1))
                ics->group_len[ics->num_groups - 1]++;
            else {
                ics->num_groups++;
                ics->group_len[ics->num_groups - 1] = 1;
            }
        }
        ics->swb_offset = aac_swb_offsets(c->frame_size, 1,
                                          c->sample_rate_index);
        ics->num_swb = aac_num_swb(c->frame_size, 1, c->sample_rate_index);
        ics->tns_max_bands =
            aac_tns_max_bands(c->frame_size, 1, c->sample_rate_index);
    } else {
        ics->max_sfb = (int)aac_bits_read(b, 6);
        ics->swb_offset = aac_swb_offsets(c->frame_size, 0,
                                          c->sample_rate_index);
        ics->num_swb = aac_num_swb(c->frame_size, 0, c->sample_rate_index);
        ics->tns_max_bands =
            aac_tns_max_bands(c->frame_size, 0, c->sample_rate_index);
        if (!eld) {
            int predictor = (int)aac_bits_read(b, 1);
            if (predictor)
                return AAC_CORE_ERROR; /* Main prediction / LTP unsupported */
        }
    }
    if (b->error)
        return AAC_CORE_NEED_MORE;
    if (ics->max_sfb > ics->num_swb)
        return AAC_CORE_ERROR;
    return AAC_CORE_OK;
}

static int decode_sections(aac_bits* b, aac_ics* ics,
                           aac_channel_work* ch)
{
    int bits = (ics->window_sequence == EIGHT_SHORT_SEQUENCE) ? 3 : 5;
    int g, sfb;

    for (g = 0; g < ics->num_groups; g++) {
        int k = 0;
        for (sfb = 0; sfb < ics->max_sfb; sfb++)
            ch->ics.band_type[g][sfb] = 0;

        while (k < ics->max_sfb) {
            int cb = (int)aac_bits_read(b, 4);
            int len = 0;
            int inc;

            if (cb == 12)
                return AAC_CORE_ERROR;
            do {
                inc = (int)aac_bits_read(b, bits);
                len += inc;
                if (len > ics->max_sfb - k)
                    return AAC_CORE_ERROR;
            } while (inc == ((1 << bits) - 1));

            if (len <= 0)
                return AAC_CORE_ERROR;
            for (sfb = 0; sfb < len; sfb++)
                ch->ics.band_type[g][k + sfb] = cb;
            k += len;
        }
        if (b->error)
            return AAC_CORE_NEED_MORE;
    }
    return AAC_CORE_OK;
}

static int decode_scalefactors(aac_bits* b, aac_channel_work* ch)
{
    aac_ics* ics = &ch->ics;
    int factor = ics->global_gain;
    int position = 0;
    int noise_active = 0;
    int noise_energy = 0;
    int g, sfb;

    for (g = 0; g < ics->num_groups; g++) {
        for (sfb = 0; sfb < ics->max_sfb; sfb++) {
            int cb = ch->ics.band_type[g][sfb];
            int delta;

            if (cb == ZERO_HCB) {
                ch->ics.sf[g][sfb] = 0;
                continue;
            }
            if (cb == NOISE_HCB) {
                if (!noise_active) {
                    int start = (int)aac_bits_read(b, 9);
                    noise_active = 1;
                    noise_energy = ics->global_gain - 90 + start - 256;
                } else {
                    delta = aac_huff_decode_sf(b);
                    if (delta < 0)
                        return AAC_CORE_ERROR;
                    noise_energy += delta - 60;
                }
                ch->noise_used[g][sfb] = 1;
                ch->ics.noise_energy[g][sfb] = noise_energy;
            } else if (cb == INTENSITY_HCB || cb == INTENSITY_HCB2) {
                delta = aac_huff_decode_sf(b);
                if (delta < 0)
                    return AAC_CORE_ERROR;
                position += delta - 60;
                if (position < -155)
                    position = -155;
                if (position > 100)
                    position = 100;
                ch->intensity_used[g][sfb] = 1;
                ch->ics.is_position[g][sfb] = position;
            } else {
                delta = aac_huff_decode_sf(b);
                if (delta < 0)
                    return AAC_CORE_ERROR;
                factor += delta - 60;
                if (factor < 0)
                    factor = 0;
                if (factor > 255)
                    factor = 255;
                ch->ics.sf[g][sfb] = factor - 100;
            }
        }
    }
    if (b->error)
        return AAC_CORE_NEED_MORE;
    return AAC_CORE_OK;
}

static int decode_pulse(aac_bits* b, aac_ics* ics)
{
    int i;
    int pos;

    if (ics->window_sequence == EIGHT_SHORT_SEQUENCE)
        return AAC_CORE_ERROR;

    ics->num_pulse = (int)aac_bits_read(b, 2) + 1;
    pos = (int)aac_bits_read(b, 6);
    if (pos >= ics->num_swb)
        return AAC_CORE_ERROR;
    pos = ics->swb_offset[pos] + (int)aac_bits_read(b, 5);
    for (i = 0; i < ics->num_pulse; i++) {
        if (i > 0)
            pos += (int)aac_bits_read(b, 5);
        if (pos < 0 || pos >= ics->swb_offset[ics->num_swb])
            return AAC_CORE_ERROR;
        ics->pulse_pos[i] = pos;
        ics->pulse_amp[i] = (int)aac_bits_read(b, 4);
    }
    return b->error ? AAC_CORE_NEED_MORE : AAC_CORE_OK;
}

static int decode_tns(aac_bits* b, aac_channel_work* ch)
{
    aac_ics* ics = &ch->ics;
    int short_window = (ics->window_sequence == EIGHT_SHORT_SEQUENCE);
    int max_order = short_window ? 7 : 12;
    int w;

    for (w = 0; w < ics->num_windows; w++) {
        int n_filt = (int)aac_bits_read(b, short_window ? 1 : 2);
        ics->tns_nfilt[w] = n_filt;
        if (n_filt) {
            int coef_res = (int)aac_bits_read(b, 1);
            int f;
            for (f = 0; f < n_filt; f++) {
                int length = (int)aac_bits_read(b, short_window ? 4 : 6);
                int order = (int)aac_bits_read(b, short_window ? 3 : 5);
                int dir = 0, comp = 0;
                int coef_len, i;
                float* coef = ics->tns_coef[w][f];

                if (order > max_order)
                    return AAC_CORE_ERROR;
                ics->tns_length[w][f] = length;
                ics->tns_order[w][f] = order;
                if (order == 0)
                    continue;
                dir = (int)aac_bits_read(b, 1);
                comp = (int)aac_bits_read(b, 1);
                ics->tns_direction[w][f] = dir;
                coef_len = 3 + coef_res - comp;
                for (i = 0; i < order; i++) {
                    int code = (int)aac_bits_read(b, coef_len);
                    if (coef_res == 0 && comp == 0)
                        coef[i] = aac_tns_coef_3[code & 7];
                    else if (coef_res == 1 && comp == 0)
                        coef[i] = aac_tns_coef_4[code & 15];
                    else if (coef_res == 0 && comp == 1) {
                        static const float tab[4] = {
                            0.0f, -0.43388373f, 0.64278758f, 0.34202015f
                        };
                        coef[i] = tab[code & 3];
                    } else {
                        static const float tab[8] = {
                            0.0f, -0.20791170f, -0.40673664f, -0.58778524f,
                            0.67369562f, 0.52643216f, 0.36124167f, 0.18374951f
                        };
                        coef[i] = tab[code & 7];
                    }
                }
            }
        }
        if (b->error)
            return AAC_CORE_NEED_MORE;
    }
    return AAC_CORE_OK;
}

static int decode_spectral_data(aac_bits* b, aac_channel_work* ch,
                                int frame_size)
{
    aac_ics* ics = &ch->ics;
    int short_window = (ics->window_sequence == EIGHT_SHORT_SEQUENCE);
    int stride = short_window ? frame_size / 8 : frame_size;
    int window_base = 0;
    int g, sfb;

    memset(ch->q, 0, sizeof(ch->q));

    for (g = 0; g < ics->num_groups; g++) {
        for (sfb = 0; sfb < ics->max_sfb; sfb++) {
            int cb = ics->band_type[g][sfb];
            int width = ics->swb_offset[sfb + 1] - ics->swb_offset[sfb];
            int dim, vecs, w;

            if (cb == ZERO_HCB || cb == NOISE_HCB ||
                cb == INTENSITY_HCB || cb == INTENSITY_HCB2)
                continue;
            if (cb < 1 || cb > 11)
                return AAC_CORE_ERROR;
            dim = codebook_dim(cb);
            if (width % dim != 0)
                return AAC_CORE_ERROR;
            vecs = width / dim;

            for (w = 0; w < ics->group_len[g]; w++) {
                int base = (window_base + w) * stride +
                           ics->swb_offset[sfb];
                int vi;
                for (vi = 0; vi < vecs; vi++) {
                    int sym = aac_huff_decode(b, cb);
                    int vals[4];
                    int signs[4];
                    int nsign = 0;
                    int sign_index = 0;
                    int j;
                    if (sym < 0)
                        return AAC_CORE_ERROR;
                    codebook_vector(cb, sym, vals);

                    if (codebook_is_unsigned(cb)) {
                        for (j = 0; j < dim; j++) {
                            if (vals[j] != 0) {
                                signs[nsign++] = (int)aac_bits_read(b, 1);
                            }
                        }
                    }

                    for (j = 0; j < dim; j++) {
                        int val = vals[j];
                        int q;

                        if (cb == 11 && val == 16) {
                            int nbits = 0;
                            while (aac_bits_read(b, 1)) {
                                nbits++;
                                if (nbits > 8)
                                    return AAC_CORE_ERROR;
                            }
                            q = (1 << (nbits + 4)) +
                                (int)aac_bits_read(b, nbits + 4);
                        } else {
                            q = val;
                        }

                        if (codebook_is_unsigned(cb) && q != 0) {
                            if (sign_index >= nsign)
                                return AAC_CORE_ERROR;
                            if (signs[sign_index++])
                                q = -q;
                        }
                        ch->q[base + vi * dim + j] = q;
                    }
                }
            }
            if (b->error)
                return AAC_CORE_NEED_MORE;
        }
        window_base += ics->group_len[g];
    }
    return AAC_CORE_OK;
}

static void dequantize_channel(aac_channel_work* ch, int frame_size)
{
    aac_ics* ics = &ch->ics;
    int short_window = (ics->window_sequence == EIGHT_SHORT_SEQUENCE);
    int stride = short_window ? frame_size / 8 : frame_size;
    int window_base = 0;
    int g, sfb;

    memset(ch->spec, 0, sizeof(ch->spec));

    for (g = 0; g < ics->num_groups; g++) {
        for (sfb = 0; sfb < ics->max_sfb; sfb++) {
            int cb = ics->band_type[g][sfb];
            int start = ics->swb_offset[sfb];
            int width = ics->swb_offset[sfb + 1] - start;
            float scale;
            int w, i;

            if (cb == ZERO_HCB || cb == NOISE_HCB ||
                cb == INTENSITY_HCB || cb == INTENSITY_HCB2)
                continue;
            scale = aac_pow2_quarter(ics->sf[g][sfb]);

            for (w = 0; w < ics->group_len[g]; w++) {
                int base = (window_base + w) * stride + start;
                for (i = 0; i < width; i++) {
                    int q = ch->q[base + i];
                    if (q > 0)
                        ch->spec[base + i] = aac_pow43(q) * scale;
                    else if (q < 0)
                        ch->spec[base + i] = -aac_pow43(-q) * scale;
                }
            }
        }
        window_base += ics->group_len[g];
    }
}

static void apply_pulses(aac_channel_work* ch)
{
    aac_ics* ics = &ch->ics;
    int i;

    for (i = 0; i < ics->num_pulse; i++) {
        int pos = ics->pulse_pos[i];
        if (pos >= 0 && pos < 1024) {
            int q = ch->q[pos];
            if (q != 0) {
                int mag = q < 0 ? -q : q;
                mag -= ics->pulse_amp[i];
                if (mag < 0)
                    mag = 0;
                ch->q[pos] = (q < 0) ? -mag : mag;
            }
        }
    }
}

/* ------------------------------------------------------------------------- */
/* PNS, M/S, intensity stereo, TNS                                           */
/* ------------------------------------------------------------------------- */

static void noise_vector(aac_core* core, float* dst, int n)
{
    float sum = 0.0f;
    int i;

    for (i = 0; i < n; i++) {
        uint32_t s;
        core->pns_seed = core->pns_seed * 1664525u + 1013904223u;
        s = core->pns_seed;
        dst[i] = (float)(int32_t)s * (1.0f / 2147483648.0f);
        sum += dst[i] * dst[i];
    }
    if (sum < 1e-20f) {
        dst[0] = 1.0f;
        return;
    }
    sum = 1.0f / (float)sqrt((double)sum);
    for (i = 0; i < n; i++)
        dst[i] *= sum;
}

static void apply_pns_sce(aac_core* core, aac_channel_work* ch,
                          int frame_size)
{
    aac_ics* ics = &ch->ics;
    int short_window = (ics->window_sequence == EIGHT_SHORT_SEQUENCE);
    int stride = short_window ? frame_size / 8 : frame_size;
    int window_base = 0;
    int g, sfb;

    for (g = 0; g < ics->num_groups; g++) {
        for (sfb = 0; sfb < ics->max_sfb; sfb++) {
            if (ch->noise_used[g][sfb]) {
                int start = ics->swb_offset[sfb];
                int width = ics->swb_offset[sfb + 1] - start;
                float scale = aac_pow2_quarter(ics->noise_energy[g][sfb]);
                int w;
                for (w = 0; w < ics->group_len[g]; w++) {
                    float n[128];
                    int base = (window_base + w) * stride + start;
                    int i;
                    noise_vector(core, n, width);
                    for (i = 0; i < width; i++)
                        ch->spec[base + i] = n[i] * scale;
                }
            }
        }
        window_base += ics->group_len[g];
    }
}

static void apply_pns_cpe(aac_core* core, aac_channel_work* ch0,
                          aac_channel_work* ch1, int frame_size)
{
    aac_ics* ics = &ch0->ics;
    int short_window = (ics->window_sequence == EIGHT_SHORT_SEQUENCE);
    int stride = short_window ? frame_size / 8 : frame_size;
    int window_base = 0;
    int g, sfb;

    for (g = 0; g < ics->num_groups; g++) {
        for (sfb = 0; sfb < ics->max_sfb; sfb++) {
            int n0 = ch0->noise_used[g][sfb];
            int n1 = ch1->noise_used[g][sfb];
            int corr = n0 && n1 && ics->ms_used[g][sfb];
            int start = ics->swb_offset[sfb];
            int width = ics->swb_offset[sfb + 1] - start;
            float s0 = aac_pow2_quarter(ics->noise_energy[g][sfb]);
            float s1 = aac_pow2_quarter(ch1->ics.noise_energy[g][sfb]);
            int w;

            if (!n0 && !n1)
                continue;

            for (w = 0; w < ics->group_len[g]; w++) {
                int base = (window_base + w) * stride + start;
                int i;
                if (corr) {
                    float n[128];
                    noise_vector(core, n, width);
                    for (i = 0; i < width; i++) {
                        ch0->spec[base + i] = n[i] * s0;
                        ch1->spec[base + i] = n[i] * s1;
                    }
                } else {
                    if (n0) {
                        float n[128];
                        noise_vector(core, n, width);
                        for (i = 0; i < width; i++)
                            ch0->spec[base + i] = n[i] * s0;
                    }
                    if (n1) {
                        float n[128];
                        noise_vector(core, n, width);
                        for (i = 0; i < width; i++)
                            ch1->spec[base + i] = n[i] * s1;
                    }
                }
            }
        }
        window_base += ics->group_len[g];
    }
}

static void apply_ms(aac_channel_work* ch0, aac_channel_work* ch1,
                     int frame_size)
{
    aac_ics* ics = &ch0->ics;
    int short_window = (ics->window_sequence == EIGHT_SHORT_SEQUENCE);
    int stride = short_window ? frame_size / 8 : frame_size;
    int window_base = 0;
    int g, sfb;

    for (g = 0; g < ics->num_groups; g++) {
        for (sfb = 0; sfb < ics->max_sfb; sfb++) {
            int cb0 = ch0->ics.band_type[g][sfb];
            int cb1 = ch1->ics.band_type[g][sfb];
            if (ics->ms_used[g][sfb] && cb0 < NOISE_HCB &&
                cb1 < NOISE_HCB) {
                int start = ics->swb_offset[sfb];
                int width = ics->swb_offset[sfb + 1] - start;
                int w, i;
                for (w = 0; w < ics->group_len[g]; w++) {
                    int base = (window_base + w) * stride + start;
                    for (i = 0; i < width; i++) {
                        float l = ch0->spec[base + i];
                        float r = ch1->spec[base + i];
                        ch0->spec[base + i] = l + r;
                        ch1->spec[base + i] = l - r;
                    }
                }
            }
        }
        window_base += ics->group_len[g];
    }
}

static void apply_is(aac_channel_work* ch0, aac_channel_work* ch1,
                     int frame_size)
{
    aac_ics* ics = &ch1->ics;
    int short_window = (ics->window_sequence == EIGHT_SHORT_SEQUENCE);
    int stride = short_window ? frame_size / 8 : frame_size;
    int window_base = 0;
    int g, sfb;

    for (g = 0; g < ics->num_groups; g++) {
        for (sfb = 0; sfb < ics->max_sfb; sfb++) {
            int cb = ch1->ics.band_type[g][sfb];
            if (cb == INTENSITY_HCB || cb == INTENSITY_HCB2) {
                float scale =
                    aac_pow2_quarter(-ch1->ics.is_position[g][sfb]);
                int start = ics->swb_offset[sfb];
                int width = ics->swb_offset[sfb + 1] - start;
                int w, i;

                if (cb == INTENSITY_HCB2)
                    scale = -scale;
                if (ch0->ics.ms_used[g][sfb])
                    scale = -scale;

                for (w = 0; w < ics->group_len[g]; w++) {
                    int base = (window_base + w) * stride + start;
                    for (i = 0; i < width; i++)
                        ch1->spec[base + i] = ch0->spec[base + i] * scale;
                }
            }
        }
        window_base += ics->group_len[g];
    }
}

static void reflection_to_lpc(const float* refl, int order, float* lpc)
{
    float prev[20];
    float cur[20];
    int i, j;

    memset(prev, 0, sizeof(prev));
    for (i = 0; i < order; i++) {
        float r = -refl[i];
        cur[i] = r;
        for (j = 0; j < (i + 1) / 2; j++) {
            float f = prev[j];
            float b = prev[i - 1 - j];
            cur[j] = f + r * b;
            cur[i - 1 - j] = b + r * f;
        }
        for (j = 0; j <= i; j++)
            prev[j] = cur[j];
    }
    for (i = 0; i < order; i++)
        lpc[i] = prev[i];
}

static void apply_tns(aac_channel_work* ch, int frame_size)
{
    aac_ics* ics = &ch->ics;
    int short_window = (ics->window_sequence == EIGHT_SHORT_SEQUENCE);
    int stride = short_window ? frame_size / 8 : frame_size;
    int mmm = ics->tns_max_bands;
    int w;

    if (mmm > ics->max_sfb)
        mmm = ics->max_sfb;
    if (mmm <= 0)
        return;

    for (w = 0; w < ics->num_windows; w++) {
        int bottom = ics->num_swb;
        int f;
        for (f = 0; f < ics->tns_nfilt[w]; f++) {
            int top = bottom;
            int order = ics->tns_order[w][f];
            int start, end, size, inc, m, i;
            float lpc[20];

            bottom = top - ics->tns_length[w][f];
            if (bottom < 0)
                bottom = 0;
            if (order <= 0)
                continue;

            reflection_to_lpc(ics->tns_coef[w][f], order, lpc);
            start = ics->swb_offset[bottom < mmm ? bottom : mmm];
            end = ics->swb_offset[top < mmm ? top : mmm];
            size = end - start;
            if (size <= 0)
                continue;

            if (ics->tns_direction[w][f]) {
                inc = -1;
                start = end - 1;
            } else {
                inc = 1;
            }
            start += w * stride;

            for (m = 0; m < size; m++, start += inc) {
                int limit = (m < order) ? m : order;
                float v = ch->spec[start];
                for (i = 1; i <= limit; i++)
                    v -= ch->spec[start - i * inc] * lpc[i - 1];
                ch->spec[start] = v;
            }
        }
    }
}

/* ------------------------------------------------------------------------- */
/* Dynamic range control                                                     */
/* ------------------------------------------------------------------------- */

/* dynamic_range_info() from ISO/IEC 14496-3 (table 4.52).  The payload is
 * carried in a FIL element with extension type 11.  Values are kept in the
 * stream configuration so the transport layer can export them. */
int aac_core_parse_dynamic_range(aac_core* core, void* bits_ctx)
{
    aac_bits* b = (aac_bits*)bits_ctx;
    aac_config* c = &core->cfg;
    int num_bands = 1;
    int i;

    if (core == NULL || b == NULL)
        return -1;

    memset(c->drc_dyn_rng_sgn, 0, sizeof(c->drc_dyn_rng_sgn));
    memset(c->drc_dyn_rng_ctl, 0, sizeof(c->drc_dyn_rng_ctl));
    c->drc_present = 1;
    c->drc_pce_instance_tag = -1;
    c->drc_prog_ref_present = 0;
    c->drc_prog_ref_level = -1;
    c->drc_band_incr = 0;
    c->drc_interpolation_scheme = 0;
    c->drc_num_bands = 1;

    if (aac_bits_read(b, 1)) {          /* pce_tag_present */
        c->drc_pce_instance_tag = (int)aac_bits_read(b, 4);
        aac_bits_skip(b, 4);            /* tag_reserved_bits */
    }

    if (aac_bits_read(b, 1)) {          /* excluded_chns_present */
        int num_excl_chan = 0;
        do {
            for (i = 0; i < 7; i++)
                (void)aac_bits_read(b, 1);
            num_excl_chan += 7;
        } while (num_excl_chan < 57 && aac_bits_read(b, 1));
    }

    if (aac_bits_read(b, 1)) {          /* drc_bands_present */
        c->drc_band_incr = (int)aac_bits_read(b, 4);
        c->drc_interpolation_scheme = (int)aac_bits_read(b, 4);
        num_bands = 1 + c->drc_band_incr;
        for (i = 0; i < num_bands; i++)
            (void)aac_bits_read(b, 8);  /* band_top */
    }

    if (aac_bits_read(b, 1)) {          /* prog_ref_level_present */
        c->drc_prog_ref_level = (int)aac_bits_read(b, 7);
        aac_bits_skip(b, 1);            /* prog_ref_level_reserved_bits */
        c->drc_prog_ref_present = 1;
    }

    c->drc_num_bands = num_bands;
    for (i = 0; i < num_bands; i++) {
        int sgn = (int)aac_bits_read(b, 1);
        int ctl = (int)aac_bits_read(b, 7);
        if (i < 17) {
            c->drc_dyn_rng_sgn[i] = sgn;
            c->drc_dyn_rng_ctl[i] = ctl;
        }
    }

    return b->error ? AAC_CORE_NEED_MORE : AAC_CORE_OK;
}

/* ------------------------------------------------------------------------- */
/* Individual channel streams                                                */
/* ------------------------------------------------------------------------- */

static int decode_ics_channel(aac_core* core, aac_bits* b,
                              aac_channel_work* ch, int common_window,
                              int lfe)
{
    aac_ics* ics = &ch->ics;
    int pulse_present = 0, tns_present = 0, gain_control_present = 0;
    int er_syntax = (core->cfg.aot == 17 || core->cfg.aot == 23);
    int eld = (core->cfg.aot == 39);
    int ret;

    (void)lfe;
    ics->global_gain = (int)aac_bits_read(b, 8);

    if (!common_window) {
        ret = decode_ics_info(core, b, ics, 0);
        if (ret < 0)
            return ret;
    }

    ret = decode_sections(b, ics, ch);
    if (ret < 0)
        return ret;
    ret = decode_scalefactors(b, ch);
    if (ret < 0)
        return ret;

    if (eld) {
        tns_present = (int)aac_bits_read(b, 1);
        ics->tns_present = tns_present;
        memset(ics->tns_nfilt, 0, sizeof(ics->tns_nfilt));
        if (tns_present) {
            ret = decode_tns(b, ch);
            if (ret < 0)
                return ret;
        }
    } else {
        pulse_present = (int)aac_bits_read(b, 1);
        ics->pulse_present = pulse_present;
        if (pulse_present) {
            ret = decode_pulse(b, ics);
            if (ret < 0)
                return ret;
        }

        tns_present = (int)aac_bits_read(b, 1);
        ics->tns_present = tns_present;
        memset(ics->tns_nfilt, 0, sizeof(ics->tns_nfilt));
        if (tns_present && !er_syntax) {
            ret = decode_tns(b, ch);
            if (ret < 0)
                return ret;
        }

        gain_control_present = (int)aac_bits_read(b, 1);
        if (gain_control_present)
            return AAC_CORE_ERROR; /* SSR gain control is not implemented */

        /* In the ER syntax the TNS payload follows the gain control flag. */
        if (tns_present && er_syntax) {
            ret = decode_tns(b, ch);
            if (ret < 0)
                return ret;
        }
    }

    ret = decode_spectral_data(b, ch, core->cfg.frame_size);
    if (ret < 0)
        return ret;

    if (pulse_present)
        apply_pulses(ch);
    dequantize_channel(ch, core->cfg.frame_size);
    return AAC_CORE_OK;
}

static int frame_add_element(aac_frame_work* frame, int type, int tag,
                             int channels)
{
    int e;
    if (frame->element_count >= 16 ||
        frame->channels + channels > AAC_CORE_MAX_CHANNELS)
        return -1;
    e = frame->element_count++;
    frame->element_type[e] = type;
    frame->element_tag[e] = tag;
    frame->element_first_channel[e] = frame->channels;
    frame->element_channels[e] = channels;
    frame->channels += channels;
    return e;
}

static int decode_sce(aac_core* core, aac_bits* b, aac_frame_work* frame,
                      int lfe)
{
    int eld = (core->cfg.aot == 39);
    int tag = eld ? 0 : (int)aac_bits_read(b, 4);
    int chi = frame->channels;
    aac_channel_work* ch;
    int ret;

    if (chi >= AAC_CORE_MAX_CHANNELS)
        return AAC_CORE_ERROR;
    ch = &frame->ch[chi];
    memset(ch, 0, sizeof(*ch));
    ch->element_type = lfe ? ELEMENT_LFE : ELEMENT_SCE;
    ch->element_tag = tag;
    ch->ics.window_sequence = ONLY_LONG_SEQUENCE;
    ch->ics.window_shape = 0;
    ch->ics.ms_mask_present = 0;

    ret = decode_ics_channel(core, b, ch, 0, lfe);
    if (ret < 0)
        return ret;
    if (frame_add_element(frame, lfe ? ELEMENT_LFE : ELEMENT_SCE, tag, 1) < 0)
        return AAC_CORE_ERROR;
    return AAC_CORE_OK;
}

static void copy_common_ics(aac_ics* dst, const aac_ics* src)
{
    dst->window_sequence = src->window_sequence;
    dst->window_shape = src->window_shape;
    dst->max_sfb = src->max_sfb;
    dst->num_windows = src->num_windows;
    dst->num_groups = src->num_groups;
    memcpy(dst->group_len, src->group_len, sizeof(dst->group_len));
    dst->num_swb = src->num_swb;
    dst->tns_max_bands = src->tns_max_bands;
    dst->swb_offset = src->swb_offset;
    dst->ms_mask_present = src->ms_mask_present;
    memcpy(dst->ms_used, src->ms_used, sizeof(dst->ms_used));
}

static int decode_cpe(aac_core* core, aac_bits* b, aac_frame_work* frame)
{
    int eld = (core->cfg.aot == 39);
    int tag = eld ? 0 : (int)aac_bits_read(b, 4);
    int common = eld ? 1 : (int)aac_bits_read(b, 1);
    int chi = frame->channels;
    aac_channel_work* ch0;
    aac_channel_work* ch1;
    int ret;

    if (chi + 1 >= AAC_CORE_MAX_CHANNELS)
        return AAC_CORE_ERROR;
    ch0 = &frame->ch[chi];
    ch1 = &frame->ch[chi + 1];
    memset(ch0, 0, sizeof(*ch0));
    memset(ch1, 0, sizeof(*ch1));
    ch0->element_type = ELEMENT_CPE;
    ch0->element_tag = tag;
    ch1->element_type = ELEMENT_CPE;
    ch1->element_tag = tag;
    ch0->ics.window_sequence = ONLY_LONG_SEQUENCE;
    ch1->ics.window_sequence = ONLY_LONG_SEQUENCE;

    if (common) {
        int ms_present;
        ret = decode_ics_info(core, b, &ch0->ics, 0);
        if (ret < 0)
            return ret;
        copy_common_ics(&ch1->ics, &ch0->ics);

        ms_present = (int)aac_bits_read(b, 2);
        ch0->ics.ms_mask_present = ms_present;
        ch1->ics.ms_mask_present = ms_present;
        if (ms_present == 3)
            return AAC_CORE_ERROR;
        if (ms_present == 1) {
            int g, sfb;
            for (g = 0; g < ch0->ics.num_groups; g++) {
                for (sfb = 0; sfb < ch0->ics.max_sfb; sfb++) {
                    int used = (int)aac_bits_read(b, 1);
                    ch0->ics.ms_used[g][sfb] = used;
                    ch1->ics.ms_used[g][sfb] = used;
                }
            }
        } else if (ms_present == 2) {
            int g, sfb;
            for (g = 0; g < ch0->ics.num_groups; g++)
                for (sfb = 0; sfb < ch0->ics.max_sfb; sfb++)
                    ch0->ics.ms_used[g][sfb] = ch1->ics.ms_used[g][sfb] = 1;
        }
        if (b->error)
            return AAC_CORE_NEED_MORE;

        ret = decode_ics_channel(core, b, ch0, 1, 0);
        if (ret < 0)
            return ret;
        ret = decode_ics_channel(core, b, ch1, 1, 0);
        if (ret < 0)
            return ret;
    } else {
        ret = decode_ics_channel(core, b, ch0, 0, 0);
        if (ret < 0)
            return ret;
        ret = decode_ics_channel(core, b, ch1, 0, 0);
        if (ret < 0)
            return ret;
        /* Intensity stereo requires a common window. */
        {
            int g, sfb;
            for (g = 0; g < ch1->ics.num_groups; g++)
                for (sfb = 0; sfb < ch1->ics.max_sfb; sfb++)
                    if (ch1->ics.band_type[g][sfb] == INTENSITY_HCB ||
                        ch1->ics.band_type[g][sfb] == INTENSITY_HCB2)
                        return AAC_CORE_ERROR;
        }
    }

    if (frame_add_element(frame, ELEMENT_CPE, tag, 2) < 0)
        return AAC_CORE_ERROR;
    return AAC_CORE_OK;
}

static int skip_dse(aac_bits* b)
{
    int tag = (int)aac_bits_read(b, 4);
    int align = (int)aac_bits_read(b, 1);
    int count = (int)aac_bits_read(b, 8);
    (void)tag;
    if (count == 255) {
        int esc;
        do {
            esc = (int)aac_bits_read(b, 8);
            count += esc;
        } while (esc == 255);
    }
    if (align)
        aac_bits_align(b);
    aac_bits_skip(b, (size_t)count * 8u);
    return b->error ? AAC_CORE_NEED_MORE : AAC_CORE_OK;
}

static void free_sbr_states(aac_core* core)
{
    int i;
    for (i = 0; i < AAC_CORE_MAX_CHANNELS; i++) {
        if (core->sbr[i] != NULL) {
            aac_sbr_close(core->sbr[i]);
            free(core->sbr[i]);
            core->sbr[i] = NULL;
        }
    }
}

static aac_sbr_decoder* ensure_sbr(aac_core* core, int first, int channels)
{
    aac_sbr_decoder* sbr;
    int out_rate;

    if (first < 0 || first >= AAC_CORE_MAX_CHANNELS || channels <= 0 ||
        channels > 2)
        return NULL;
    out_rate = 2 * core->cfg.sample_rate;
    sbr = core->sbr[first];
    if (sbr != NULL) {
        if (sbr->channels == channels && sbr->core_rate == core->cfg.sample_rate &&
            sbr->frame_size == core->cfg.frame_size)
            return sbr;
        aac_sbr_close(sbr);
        free(sbr);
        core->sbr[first] = NULL;
    }
    sbr = (aac_sbr_decoder*)malloc(sizeof(*sbr));
    if (sbr == NULL)
        return NULL;
    if (aac_sbr_open(sbr, channels, core->cfg.sample_rate, out_rate,
                     core->cfg.frame_size) < 0) {
        free(sbr);
        return NULL;
    }
    core->sbr[first] = sbr;
    return sbr;
}

/* Decode one fill element.  SBR extension payloads (types 13/14) are parsed
 * and associated with the most recent audio element. */
static int decode_fil(aac_core* core, aac_bits* b, aac_frame_work* frame,
                      int last_element)
{
    int count = (int)aac_bits_read(b, 4);
    size_t end_bit;
    size_t start_bit;

    if (count == 15) {
        int esc = (int)aac_bits_read(b, 8);
        count += esc - 1;
    }
    if (b->error)
        return AAC_CORE_NEED_MORE;

    start_bit = b->bitpos;
    end_bit = start_bit + (size_t)count * 8u;
    if (end_bit > b->size * 8u) {
        b->error = 1;
        return AAC_CORE_NEED_MORE;
    }

    if (count > 0) {
        int ext_type = (int)aac_bits_read(b, 4);
        if ((ext_type == 13 || ext_type == 14) && last_element >= 0 &&
            last_element < frame->element_count) {
            int first = frame->element_first_channel[last_element];
            int ech = frame->element_channels[last_element];
            int etype = (ech == 2) ? 1 : 0;
            aac_sbr_decoder* sbr = ensure_sbr(core, first, ech);
            if (sbr != NULL) {
                aac_bits sub;
                int len_bits = count * 8 - 4;
                int ret;
                sub.data = b->data;
                /* The SBR syntax is byte-aligned only at the fill element
                 * level; allow the parser to read the last few padding bits of
                 * the access unit in the same way the reference decoders do. */
                sub.size = b->size;
                sub.bitpos = b->bitpos;
                sub.error = 0;
                ret = aac_sbr_parse(sbr, &sub, len_bits, etype,
                                    ext_type == 14);
                if (ret == 0) {
                    frame->element_sbr[last_element] = 1;
                } else {
                    /* A malformed SBR payload must not discard the audio
                     * frame; fall back to QMF upsampling for this frame. */
                    sbr->active = 0;
                }
                if (!core->cfg.has_sbr) {
                    core->cfg.has_sbr = 1;
                    core->cfg.ext_aot = 5;
                    if (core->cfg.ext_sample_rate <= 0)
                        core->cfg.ext_sample_rate =
                            2 * core->cfg.sample_rate;
                }
                if (ret == 0 && sbr->ps_present) {
                    core->cfg.has_ps = 1;
                    core->cfg.ext_aot = 29;
                }
            }
        } else if (ext_type == 11) {
            aac_bits sub;
            sub.data = b->data;
            sub.size = (end_bit + 7u) / 8u;
            sub.bitpos = b->bitpos;
            sub.error = 0;
            (void)aac_core_parse_dynamic_range(core, &sub);
        }
    }

    b->bitpos = end_bit;
    return b->error ? AAC_CORE_NEED_MORE : AAC_CORE_OK;
}

static int decode_er_raw_data_block(aac_core* core, aac_bits* b,
                                    aac_frame_work* frame)
{
    static const int layout[8][5] = {
        {0, 0, 0, 0, 0},
        {ELEMENT_SCE, 0, 0, 0, 0},
        {ELEMENT_CPE, 0, 0, 0, 0},
        {ELEMENT_SCE, ELEMENT_CPE, 0, 0, 0},
        {ELEMENT_SCE, ELEMENT_CPE, ELEMENT_SCE, 0, 0},
        {ELEMENT_SCE, ELEMENT_CPE, ELEMENT_CPE, 0, 0},
        {ELEMENT_SCE, ELEMENT_CPE, ELEMENT_CPE, ELEMENT_LFE, 0},
        {ELEMENT_SCE, ELEMENT_CPE, ELEMENT_CPE, ELEMENT_CPE, ELEMENT_LFE}
    };
    int cc = core->cfg.channel_config;
    int i, count;

    if (cc < 1 || cc > 7)
        return AAC_CORE_ERROR;
    count = (cc == 1) ? 1 : (cc == 2) ? 1 : (cc == 3) ? 2 : (cc == 4) ? 3 :
            (cc == 5) ? 3 : (cc == 6) ? 4 : 5;

    for (i = 0; i < count; i++) {
        int type = layout[cc][i];
        int r;
        if (type == ELEMENT_CPE)
            r = decode_cpe(core, b, frame);
        else
            r = decode_sce(core, b, frame, type == ELEMENT_LFE);
        if (r < 0)
            return r;
    }

    aac_bits_align(b);
    if (b->error)
        return AAC_CORE_NEED_MORE;
    return AAC_CORE_OK;
}

static int decode_raw_data_block(aac_core* core, aac_bits* b,
                                 aac_frame_work* frame)
{
    int id;
    int last_element = -1;

    memset(frame, 0, sizeof(*frame));
    frame->frame_size = core->cfg.frame_size;
    frame->sample_rate = core->cfg.sample_rate;

    /* AAC-LD uses the error resilient syntax where the element types are
     * implied by the channel configuration instead of an id_syn_ele field. */
    if (core->cfg.aot == 17 || core->cfg.aot == 23 ||
        core->cfg.aot == 39)
        return decode_er_raw_data_block(core, b, frame);

    for (;;) {
        id = (int)aac_bits_read(b, 3);
        if (b->error)
            return AAC_CORE_NEED_MORE;

        switch (id) {
        case ID_SCE: {
            int before = frame->element_count;
            int r = decode_sce(core, b, frame, 0);
            if (r < 0)
                return r;
            last_element = before;
            break;
        }
        case ID_CPE: {
            int before = frame->element_count;
            int r = decode_cpe(core, b, frame);
            if (r < 0)
                return r;
            last_element = before;
            break;
        }
        case ID_CCE:
            /* Coupling channel elements are not implemented. */
            return AAC_CORE_ERROR;
        case ID_LFE: {
            int before = frame->element_count;
            int r = decode_sce(core, b, frame, 1);
            if (r < 0)
                return r;
            last_element = before;
            break;
        }
        case ID_DSE: {
            int r = skip_dse(b);
            if (r < 0)
                return r;
            last_element = -1;
            break;
        }
        case ID_PCE: {
            int r = read_pce(core, b, 0);
            if (r < 0)
                return r;
            last_element = -1;
            break;
        }
        case ID_FIL: {
            int r = decode_fil(core, b, frame, last_element);
            if (r < 0)
                return r;
            break;
        }
        case ID_END:
            aac_bits_align(b);
            if (b->error)
                return AAC_CORE_NEED_MORE;
            return AAC_CORE_OK;
        default:
            return AAC_CORE_ERROR;
        }
    }
}

/* ------------------------------------------------------------------------- */
/* Configuration                                                             */
/* ------------------------------------------------------------------------- */

static void config_clear(aac_config* c)
{
    memset(c, 0, sizeof(*c));
    c->profile = -1;
    c->sample_rate_index = -1;
    c->core_aot = 0;
    c->aot = 0;
    c->frame_size = 1024;
    c->short_len = 128;
    c->drc_pce_instance_tag = -1;
    c->drc_prog_ref_level = -1;
    c->drc_num_bands = 1;
}

static int config_apply(aac_core* core, int reset)
{
    aac_config* c = &core->cfg;
    int old_channels = core->cfg.channels;
    int old_rate = core->cfg.sample_rate;
    int old_frame = core->cfg.frame_size;
    (void)old_channels;
    (void)old_rate;
    (void)old_frame;

    if (c->frame_size <= 0)
        c->frame_size = 1024;
    c->short_len = c->frame_size / 8;

    if (c->channel_config >= 1 && c->channel_config <= 7)
        set_standard_layout(c, c->channel_config);
    else if (c->channel_config == 0 && c->pce_valid) {
        /* The PCE parser has already filled channels and the output order. */
    } else {
        c->channels = 0;
    }

    if (c->sample_rate <= 0)
        return AAC_CORE_ERROR;

    if (reset) {
        memset(core->state, 0, sizeof(core->state));
        core->pns_seed = 0x1f2e3d4cu;
        core->frames = 0;
        c->drc_present = 0;
        c->drc_pce_instance_tag = -1;
        c->drc_prog_ref_present = 0;
        c->drc_prog_ref_level = -1;
        c->drc_num_bands = 1;
        memset(c->drc_dyn_rng_sgn, 0, sizeof(c->drc_dyn_rng_sgn));
        memset(c->drc_dyn_rng_ctl, 0, sizeof(c->drc_dyn_rng_ctl));
        free_sbr_states(core);
    }
    core->configured = 1;
    return AAC_CORE_OK;
}

static int find_sr_index(int rate)
{
    int i;
    for (i = 0; i < AAC_SR_INDEX_MAX; i++)
        if (aac_sample_rates[i] == rate)
            return i;
    return 4; /* 44100 fallback */
}

/* ------------------------------------------------------------------------- */
/* USAC (xHE-AAC, ISO/IEC 23003-3) frame decoding (audioObjectType 42)       */
/* ------------------------------------------------------------------------- */

static const uint8_t usac_tns_max_bands_1024[13] = {
    31, 31, 34, 40, 42, 51, 47, 47, 43, 43, 43, 40, 40
};
static const uint8_t usac_tns_max_bands_128[13] = {
    9, 9, 10, 14, 14, 14, 15, 15, 15, 15, 15, 15, 15
};
static const uint8_t usac_noise_fill_start[2][2] = {
    { 160, 20 },
    { 120, 15 }
};

/* Per-element state that is shared between the stereo information parser and
 * the per-channel decoding loops of one frame. */
typedef struct
{
    int common_window;
    int ms_mask_mode;
    int tns_on_lr;
    int max_sfb_ste;
    int core_mode[2];
    int tns_data_present[2];

    /* Complex (inter-channel) stereo prediction, ISO/IEC 23003-3 7.3.4.2 */
    int cplx_pred;
    int pred_dir;
    int complex_coef;
    int use_prev_frame;
    int delta_code_time;
    uint8_t pred_used[AAC_CORE_MAX_GROUPS * AAC_CORE_MAX_SFB];
    int alpha_q_re[AAC_CORE_MAX_GROUPS * AAC_CORE_MAX_SFB];
    int alpha_q_im[AAC_CORE_MAX_GROUPS * AAC_CORE_MAX_SFB];
} usac_cpe_state;

/* usac_stereo_cplx(): complex stereo prediction parameters.  The alpha/beta
 * prediction coefficients are differentially coded with the same Huffman
 * table as the scale factors.  The DPCM direction is either within the frame
 * (frequency) or across frames/groups (time), signalled by delta_code_time. */
static int usac_decode_stereo_cplx(aac_core* core, aac_bits* b,
                                   usac_cpe_state* st, int indep,
                                   int num_groups, int window_sequence)
{
    const int stride = st->max_sfb_ste;
    int g, sfb;

    memset(st->pred_used, 0, sizeof(st->pred_used));
    memset(st->alpha_q_re, 0, sizeof(st->alpha_q_re));
    memset(st->alpha_q_im, 0, sizeof(st->alpha_q_im));

    if (!aac_bits_read(b, 1)) { /* cplx_pred_all */
        for (g = 0; g < num_groups; g++) {
            for (sfb = 0; sfb < st->max_sfb_ste; sfb += 2) {
                int val = (int)aac_bits_read(b, 1);
                st->pred_used[g * stride + sfb] = (uint8_t)val;
                if (sfb + 1 < st->max_sfb_ste)
                    st->pred_used[g * stride + sfb + 1] = (uint8_t)val;
            }
        }
    } else {
        for (g = 0; g < num_groups; g++)
            for (sfb = 0; sfb < st->max_sfb_ste; sfb++)
                st->pred_used[g * stride + sfb] = 1;
    }

    st->pred_dir = (int)aac_bits_read(b, 1);
    st->complex_coef = (int)aac_bits_read(b, 1);

    st->use_prev_frame = 0;
    if (st->complex_coef && !indep)
        st->use_prev_frame = (int)aac_bits_read(b, 1);

    st->delta_code_time = 0;
    if (!indep)
        st->delta_code_time = (int)aac_bits_read(b, 1);

    /* A long <-> short window transition invalidates the temporal predictors. */
    if (st->delta_code_time &&
        ((window_sequence == EIGHT_SHORT_SEQUENCE) !=
         (core->usac_cplx_win_seq_prev == EIGHT_SHORT_SEQUENCE))) {
        memset(core->usac_cplx_alpha_re_prev, 0,
               sizeof(core->usac_cplx_alpha_re_prev));
        memset(core->usac_cplx_alpha_im_prev, 0,
               sizeof(core->usac_cplx_alpha_im_prev));
    }

    for (g = 0; g < num_groups; g++) {
        for (sfb = 0; sfb < st->max_sfb_ste; sfb += 2) {
            int idx = g * stride + sfb;
            int last_re = 0, last_im = 0;

            if (st->delta_code_time) {
                if (g > 0) {
                    last_re = st->alpha_q_re[(g - 1) * stride + sfb];
                    last_im = st->alpha_q_im[(g - 1) * stride + sfb];
                } else if (window_sequence == EIGHT_SHORT_SEQUENCE &&
                           core->usac_cplx_win_seq_prev ==
                               EIGHT_SHORT_SEQUENCE) {
                    int pg = core->usac_cplx_win_groups_prev - 1;
                    if (pg < 0)
                        pg = 0;
                    last_re = core->usac_cplx_alpha_re_prev[pg * stride + sfb];
                    last_im = core->usac_cplx_alpha_im_prev[pg * stride + sfb];
                } else {
                    last_re = core->usac_cplx_alpha_re_prev[sfb];
                    last_im = core->usac_cplx_alpha_im_prev[sfb];
                }
            } else if (sfb > 0) {
                last_re = st->alpha_q_re[idx - 1];
                last_im = st->alpha_q_im[idx - 1];
            }

            if (st->pred_used[idx]) {
                int v = aac_huff_decode_sf(b);
                if (v < 0)
                    return -1;
                st->alpha_q_re[idx] = (-v + 60) + last_re;
                if (st->complex_coef) {
                    v = aac_huff_decode_sf(b);
                    if (v < 0)
                        return -1;
                    st->alpha_q_im[idx] = (-v + 60) + last_im;
                } else {
                    st->alpha_q_im[idx] = 0;
                }
            } else {
                st->alpha_q_re[idx] = 0;
                st->alpha_q_im[idx] = 0;
            }

            if (sfb + 1 < st->max_sfb_ste) {
                st->alpha_q_re[idx + 1] = st->alpha_q_re[idx];
                st->alpha_q_im[idx + 1] = st->alpha_q_im[idx];
            }
            core->usac_cplx_alpha_re_prev[idx] = st->alpha_q_re[idx];
            core->usac_cplx_alpha_im_prev[idx] = st->alpha_q_im[idx];
        }
    }

    core->usac_cplx_valid = 1;

    return b->error ? -1 : 0;
}

static void usac_apply_grouping(aac_ics* ics, int sfg)
{
    int j;

    ics->num_groups = 0;
    for (j = 0; j < 7; j++) {
        ics->group_len[j] = 1;
        if (sfg & (1 << (6 - j)))
            ics->group_len[ics->num_groups] += 1;
        else
            ics->num_groups++;
    }
    ics->group_len[7] = 1;
    ics->num_groups++;
    ics->num_windows = 8;
}

static int usac_ics_setup(aac_channel_work* ch, int sr_index, int frame_size,
                          int sfg)
{
    aac_ics* ics = &ch->ics;
    int short_window = (ics->window_sequence == EIGHT_SHORT_SEQUENCE);
    int g, sfb;

    if (short_window) {
        ics->swb_offset = aac_swb_offsets(frame_size, 1, sr_index);
        ics->num_swb = aac_num_swb(frame_size, 1, sr_index);
        ics->tns_max_bands = usac_tns_max_bands_128[sr_index];
        usac_apply_grouping(ics, sfg);
    } else {
        ics->swb_offset = aac_swb_offsets(frame_size, 0, sr_index);
        ics->num_swb = aac_num_swb(frame_size, 0, sr_index);
        ics->tns_max_bands = usac_tns_max_bands_1024[sr_index];
        ics->group_len[0] = 1;
        ics->num_groups = 1;
        ics->num_windows = 1;
    }
    if (ics->max_sfb > ics->num_swb)
        return -1;
    for (g = 0; g < AAC_CORE_MAX_GROUPS; g++)
        for (sfb = 0; sfb < AAC_CORE_MAX_SFB; sfb++)
            ics->band_type[g][sfb] = 1; /* uniformly quantized bands */
    return 0;
}

/* USAC tns_data(): identical to AAC except that the filter order field is
 * one bit shorter for long windows and there is no maximum order check. */
static int usac_decode_tns(aac_bits* b, aac_channel_work* ch)
{
    static const float coef_tab4[4] = {
        0.0f, -0.43388373f, 0.64278758f, 0.34202015f
    };
    static const float coef_tab8[8] = {
        0.0f, -0.20791170f, -0.40673664f, -0.58778524f,
        0.67369562f, 0.52643216f, 0.36124167f, 0.18374951f
    };
    aac_ics* ics = &ch->ics;
    int short_window = (ics->window_sequence == EIGHT_SHORT_SEQUENCE);
    int w, f, i;

    for (w = 0; w < ics->num_windows; w++) {
        int n_filt = (int)aac_bits_read(b, short_window ? 1 : 2);
        ics->tns_nfilt[w] = n_filt;
        if (n_filt) {
            int coef_res = (int)aac_bits_read(b, 1);
            for (f = 0; f < n_filt; f++) {
                int length = (int)aac_bits_read(b, short_window ? 4 : 6);
                int order = (int)aac_bits_read(b, short_window ? 3 : 4);
                float* coef = ics->tns_coef[w][f];

                ics->tns_length[w][f] = length;
                ics->tns_order[w][f] = order;
                if (order == 0)
                    continue;
                if (order > AAC_CORE_MAX_TNS_ORDER)
                    return -1;
                ics->tns_direction[w][f] = (int)aac_bits_read(b, 1);
                {
                    int comp = (int)aac_bits_read(b, 1);
                    int coef_len = coef_res + 3 - comp;
                    for (i = 0; i < order; i++) {
                        int code = (int)aac_bits_read(b, coef_len);
                        if (coef_res == 0 && comp == 0)
                            coef[i] = aac_tns_coef_3[code & 7];
                        else if (coef_res == 1 && comp == 0)
                            coef[i] = aac_tns_coef_4[code & 15];
                        else if (coef_res == 0 && comp == 1)
                            coef[i] = coef_tab4[code & 3];
                        else
                            coef[i] = coef_tab8[code & 7];
                    }
                }
            }
        }
        if (b->error)
            return -1;
    }
    return 0;
}

static void usac_copy_tns(aac_ics* dst, const aac_ics* src)
{
    int w, f, i;

    dst->tns_present = src->tns_present;
    for (w = 0; w < AAC_CORE_MAX_WINDOWS; w++) {
        dst->tns_nfilt[w] = src->tns_nfilt[w];
        for (f = 0; f < AAC_CORE_MAX_TNS_FILTERS; f++) {
            dst->tns_length[w][f] = src->tns_length[w][f];
            dst->tns_order[w][f] = src->tns_order[w][f];
            dst->tns_direction[w][f] = src->tns_direction[w][f];
            for (i = 0; i < AAC_CORE_MAX_TNS_ORDER; i++)
                dst->tns_coef[w][f][i] = src->tns_coef[w][f][i];
        }
    }
}

static int usac_decode_scalefactors(aac_bits* b, aac_channel_work* ch,
                                    int global_gain)
{
    aac_ics* ics = &ch->ics;
    int offset = global_gain;
    int g, sfb;

    for (g = 0; g < ics->num_groups; g++) {
        for (sfb = 0; sfb < ics->max_sfb; sfb++) {
            if (g || sfb) {
                int delta = aac_huff_decode_sf(b);
                if (delta < 0)
                    return -1;
                offset += delta - 60;
                if (offset < 0 || offset > 255)
                    return -1;
            }
            ics->sf[g][sfb] = offset - 100;
        }
    }
    return b->error ? -1 : 0;
}

/* Random signs for noise filling (ISO/IEC 23003-3 7.2.4). */
static float usac_noise_sign(uint32_t* seed)
{
    uint32_t s = *seed = (*seed) * 69069u + 5u;
    return (s & 0x10000u) ? -1.0f : 1.0f;
}

/* facData() of an FD channel stream (gain + codebook indices).  The decoded
 * values only matter for LPD transitions; parsing keeps the bitstream in
 * sync. */
static int usac_fac_codebook_idx(aac_bits* b)
{
    int qn = 0, nk, n, i;

    /* unary value, stop bit = 0 */
    for (;;) {
        int bit = (int)aac_bits_read(b, 1);
        if (b->error)
            return -1;
        if (!bit)
            break;
        qn++;
    }
    if (qn)
        qn++;

    if (qn > 4) {
        nk = (qn - 3) / 2;
        n = qn - nk * 2;
    } else {
        nk = 0;
        n = qn;
    }

    aac_bits_skip(b, (size_t)(4 * n));
    if (nk > 0)
        for (i = 0; i < 8; i++)
            (void)aac_bits_read(b, nk);
    return b->error ? -1 : 0;
}

static int usac_parse_fac_data(aac_bits* b, int fac_len)
{
    int i;

    (void)aac_bits_read(b, 7); /* facGain */
    if (fac_len / 8 > 32)
        return -1;
    for (i = 0; i < fac_len / 8; i++)
        if (usac_fac_codebook_idx(b) < 0)
            return -1;
    return b->error ? -1 : 0;
}

/* Noise filling (noiseFillingData) and band scaling.  USAC spectra are
 * stored window interleaved with a stride of 128 coefficients. */
static void usac_spectrum_scale(aac_core* core, aac_channel_work* ch,
                                int frame_size)
{
    aac_ics* ics = &ch->ics;
    int short_window = (ics->window_sequence == EIGHT_SHORT_SEQUENCE);
    int g, sfb, w, i;
    int wbase;

    if (ch->usac_noise_level) {
        float noise_val = powf(2.0f, ((float)ch->usac_noise_level - 14.0f) /
                                          3.0f);
        int noise_offset = ch->usac_noise_offset - 16;
        int band_off =
            usac_noise_fill_start[frame_size == 768][short_window];

        wbase = 0;
        for (g = 0; g < ics->num_groups; g++) {
            for (sfb = 0; sfb < ics->max_sfb; sfb++) {
                int start = ics->swb_offset[sfb];
                int width = ics->swb_offset[sfb + 1] - start;
                int band_quantized_to_zero = 1;

                if (start < band_off)
                    continue;
                for (w = 0; w < ics->group_len[g]; w++) {
                    int base = (wbase + w) * 128 + start;
                    for (i = 0; i < width; i++) {
                        if (ch->spec[base + i] == 0.0f)
                            ch->spec[base + i] =
                                usac_noise_sign(&core->usac_noise_seed
                                                    [ch - core->work.ch]) *
                                noise_val;
                        else
                            band_quantized_to_zero = 0;
                    }
                }
                if (band_quantized_to_zero) {
                    int sf = ics->sf[g][sfb] + noise_offset;
                    if (sf < -200)
                        sf = -200;
                    ics->sf[g][sfb] = sf;
                }
            }
            wbase += ics->group_len[g];
        }
    }

    wbase = 0;
    for (g = 0; g < ics->num_groups; g++) {
        for (sfb = 0; sfb < ics->max_sfb; sfb++) {
            int start = ics->swb_offset[sfb];
            int width = ics->swb_offset[sfb + 1] - start;
            float sf = aac_pow2_quarter(ics->sf[g][sfb]);
            for (w = 0; w < ics->group_len[g]; w++) {
                int base = (wbase + w) * 128 + start;
                for (i = 0; i < width; i++)
                    ch->spec[base + i] *= sf;
            }
        }
        wbase += ics->group_len[g];
    }
}

/* usac_stereo_data(): common window, M/S masks and common TNS. */
static int usac_stereo_info(aac_core* core, aac_bits* b, int ch0, int indep,
                            usac_cpe_state* st)
{
    aac_frame_work* frame = &core->work;
    aac_channel_work* w0 = &frame->ch[ch0];
    aac_channel_work* w1 = &frame->ch[ch0 + 1];
    aac_ics* ics1 = &w0->ics;
    aac_ics* ics2 = &w1->ics;
    int tns_active;
    int sr_index = core->cfg.sample_rate_index;
    int frame_size = core->cfg.frame_size;
    int g, sfb;

    st->common_window = 0;
    memset(ics1->ms_used, 0, sizeof(ics1->ms_used));
    memset(ics2->ms_used, 0, sizeof(ics2->ms_used));

    /* Rotate the complex-prediction history: the values of the previous
     * frame stay available through *_prev while this frame is decoded. */
    if (core->usac_cplx_valid) {
        core->usac_cplx_win_seq_prev = core->usac_cplx_win_seq_cur;
        core->usac_cplx_win_shape_prev = core->usac_cplx_win_shape_cur;
        core->usac_cplx_win_groups_prev = core->usac_cplx_win_groups_cur;
    }

    if (st->core_mode[0] || st->core_mode[1])
        return AAC_CORE_OK;

    tns_active = (int)aac_bits_read(b, 1);
    st->common_window = (int)aac_bits_read(b, 1);

    if (st->common_window) {
        int sfg = 0;
        int seq = (int)aac_bits_read(b, 2);
        int kb = (int)aac_bits_read(b, 1);

        ics1->window_sequence = ics2->window_sequence = seq;
        ics1->window_shape = ics2->window_shape = kb;

        if (seq == EIGHT_SHORT_SEQUENCE) {
            ics1->max_sfb = ics2->max_sfb = (int)aac_bits_read(b, 4);
            sfg = (int)aac_bits_read(b, 7);
        } else {
            ics1->max_sfb = ics2->max_sfb = (int)aac_bits_read(b, 6);
        }

        if (!(int)aac_bits_read(b, 1) /* common_max_sfb */) {
            ics2->max_sfb = (seq == EIGHT_SHORT_SEQUENCE)
                                ? (int)aac_bits_read(b, 4)
                                : (int)aac_bits_read(b, 6);
        }

        if (usac_ics_setup(w0, sr_index, frame_size, sfg) < 0 ||
            usac_ics_setup(w1, sr_index, frame_size, sfg) < 0)
            return AAC_CORE_ERROR;

        st->max_sfb_ste = ics1->max_sfb > ics2->max_sfb ? ics1->max_sfb
                                                        : ics2->max_sfb;
        core->usac_cplx_win_groups_cur = ics1->num_groups;
        core->usac_cplx_win_seq_cur = ics1->window_sequence;
        core->usac_cplx_win_shape_cur = ics1->window_shape;

        st->ms_mask_mode = (int)aac_bits_read(b, 2);
        if (st->ms_mask_mode == 1) {
            for (g = 0; g < ics1->num_groups; g++)
                for (sfb = 0; sfb < st->max_sfb_ste; sfb++)
                    ics1->ms_used[g][sfb] = (int)aac_bits_read(b, 1);
        } else if (st->ms_mask_mode == 2) {
            memset(ics1->ms_used, 1, sizeof(ics1->ms_used));
        } else if (st->ms_mask_mode == 3) {
            if (usac_decode_stereo_cplx(core, b, st, indep,
                                        ics1->num_groups,
                                        ics1->window_sequence) < 0) {
                USDBG("usac: complex stereo prediction parse failed\n");
                return AAC_CORE_ERROR;
            }
            st->cplx_pred = 1;
#ifdef AAC_USAC_DEBUG
            {
                int cnt = 0, k, asum = 0;
                for (k = 0; k < ics1->num_groups * st->max_sfb_ste; k++) {
                    cnt += st->pred_used[k];
                    asum += st->alpha_q_re[k] < 0 ? -st->alpha_q_re[k]
                                                  : st->alpha_q_re[k];
                }
                fprintf(stderr,
                        "usac: cplx pred dir=%d coef=%d prev=%d dt=%d used=%d "
                        "asum=%d\n",
                        st->pred_dir, st->complex_coef, st->use_prev_frame,
                        st->delta_code_time, cnt, asum);
            }
#endif
        }
        memcpy(ics2->ms_used, ics1->ms_used, sizeof(ics1->ms_used));
        ics1->ms_mask_present = ics2->ms_mask_present = st->ms_mask_mode;
        if (st->ms_mask_mode != 3) {
            /* No complex prediction in this frame: drop the temporal state. */
            memset(core->usac_cplx_alpha_re_prev, 0,
                   sizeof(core->usac_cplx_alpha_re_prev));
            memset(core->usac_cplx_alpha_im_prev, 0,
                   sizeof(core->usac_cplx_alpha_im_prev));
            core->usac_cplx_prev_len = 0;
        }
    }

    st->tns_on_lr = 0;
    st->tns_data_present[0] = st->tns_data_present[1] = 0;
    ics1->tns_present = ics2->tns_present = 0;
    if (tns_active) {
        int common_tns = 0;
        if (st->common_window)
            common_tns = (int)aac_bits_read(b, 1);

        st->tns_on_lr = (int)aac_bits_read(b, 1);
        if (common_tns) {
            if (usac_decode_tns(b, w0) < 0)
                return AAC_CORE_ERROR;
            ics1->tns_present = ics2->tns_present = 1;
            usac_copy_tns(ics2, ics1);
        } else if ((int)aac_bits_read(b, 1)) {
            st->tns_data_present[0] = st->tns_data_present[1] = 1;
        } else {
            st->tns_data_present[1] = (int)aac_bits_read(b, 1);
            st->tns_data_present[0] = !st->tns_data_present[1];
        }
    }
    (void)indep;
    return b->error ? AAC_CORE_NEED_MORE : AAC_CORE_OK;
}

/* MDST interpolation kernels used by the complex stereo prediction,
 * ISO/IEC 23003-3 tables 124 (current) and 125 (previous window). */
static const float usac_mdst_cur[4][4][7] = {
    {{0.000000f, 0.000000f, 0.500000f, 0.000000f, -0.500000f, 0.000000f,
      0.000000f},
     {0.091497f, 0.000000f, 0.581427f, 0.000000f, -0.581427f, 0.000000f,
      -0.091497f},
     {0.045748f, 0.057238f, 0.540714f, 0.000000f, -0.540714f, -0.057238f,
      -0.045748f},
     {0.045748f, -0.057238f, 0.540714f, 0.000000f, -0.540714f, 0.057238f,
      -0.045748f}},
    {{0.102658f, 0.103791f, 0.567149f, 0.000000f, -0.567149f, -0.103791f,
      -0.102658f},
     {0.150512f, 0.047969f, 0.608574f, 0.000000f, -0.608574f, -0.047969f,
      -0.150512f},
     {0.104763f, 0.105207f, 0.567861f, 0.000000f, -0.567861f, -0.105207f,
      -0.104763f},
     {0.148406f, 0.046553f, 0.607863f, 0.000000f, -0.607863f, -0.046553f,
      -0.148406f}},
    {{0.102658f, -0.103791f, 0.567149f, 0.000000f, -0.567149f, 0.103791f,
      -0.102658f},
     {0.150512f, -0.047969f, 0.608574f, 0.000000f, -0.608574f, 0.047969f,
      -0.150512f},
     {0.148406f, -0.046553f, 0.607863f, 0.000000f, -0.607863f, 0.046553f,
      -0.148406f},
     {0.104763f, -0.105207f, 0.567861f, 0.000000f, -0.567861f, 0.105207f,
      -0.104763f}},
    {{0.205316f, 0.000000f, 0.634298f, 0.000000f, -0.634298f, 0.000000f,
      -0.205316f},
     {0.209526f, 0.000000f, 0.635722f, 0.000000f, -0.635722f, 0.000000f,
      -0.209526f},
     {0.207421f, 0.001416f, 0.635010f, 0.000000f, -0.635010f, -0.001416f,
      -0.207421f},
     {0.207421f, -0.001416f, 0.635010f, 0.000000f, -0.635010f, 0.001416f,
      -0.207421f}}};

static const float usac_mdst_prev[2][2][7] = {
    {{0.000000f, 0.106103f, 0.250000f, 0.318310f, 0.250000f, 0.106103f,
      0.000000f},
     {0.059509f, 0.123714f, 0.186579f, 0.213077f, 0.186579f, 0.123714f,
      0.059509f}},
    {{0.038498f, 0.039212f, 0.039645f, 0.039790f, 0.039645f, 0.039212f,
      0.038498f},
     {0.026142f, 0.026413f, 0.026577f, 0.026631f, 0.026577f, 0.026413f,
      0.026142f}}};

/* Interpolate the MDST (imaginary part) from the MDCT (real part) using the
 * 7-tap kernel f[].  len is the number of coefficients to process. */
static void usac_mdst_filter(float* im, const float* re, const float f[7],
                             int len, int factor_even, int factor_odd)
{
    int i;
    float s;

    if (len < 7)
        return;

    s = f[6] * re[2] + f[5] * re[1] + f[4] * re[0] + f[3] * re[0] +
        f[2] * re[1] + f[1] * re[2] + f[0] * re[3];
    im[0] += s * (float)factor_even;

    s = f[6] * re[1] + f[5] * re[0] + f[4] * re[0] + f[3] * re[1] +
        f[2] * re[2] + f[1] * re[3] + f[0] * re[4];
    im[1] += s * (float)factor_odd;

    s = f[6] * re[0] + f[5] * re[0] + f[4] * re[1] + f[3] * re[2] +
        f[2] * re[3] + f[1] * re[4] + f[0] * re[5];
    im[2] += s * (float)factor_even;

    for (i = 3; i < len - 4; i += 2) {
        s = f[6] * re[i - 3] + f[5] * re[i - 2] + f[4] * re[i - 1] +
            f[3] * re[i] + f[2] * re[i + 1] + f[1] * re[i + 2] +
            f[0] * re[i + 3];
        im[i] += s * (float)factor_odd;

        s = f[6] * re[i - 2] + f[5] * re[i - 1] + f[4] * re[i] +
            f[3] * re[i + 1] + f[2] * re[i + 2] + f[1] * re[i + 3] +
            f[0] * re[i + 4];
        im[i + 1] += s * (float)factor_even;
    }

    i = len - 3;
    s = f[6] * re[i - 3] + f[5] * re[i - 2] + f[4] * re[i - 1] +
        f[3] * re[i] + f[2] * re[i + 1] + f[1] * re[i + 2] +
        f[0] * re[i + 2];
    im[i] += s * (float)factor_odd;

    i = len - 2;
    s = f[6] * re[i - 3] + f[5] * re[i - 2] + f[4] * re[i - 1] +
        f[3] * re[i] + f[2] * re[i + 1] + f[1] * re[i + 1] +
        f[0] * re[i];
    im[i] += s * (float)factor_even;

    i = len - 1;
    s = f[6] * re[i - 3] + f[5] * re[i - 2] + f[4] * re[i - 1] +
        f[3] * re[i] + f[2] * re[i] + f[1] * re[i - 1] +
        f[0] * re[i - 2];
    im[i] += s * (float)factor_odd;
}

/* Complex stereo prediction upmix (ISO/IEC 23003-3 7.7.2.3.3). */
static void usac_apply_cplx(aac_core* core, aac_channel_work* w0,
                            aac_channel_work* w1, int frame_size,
                            usac_cpe_state* st)
{
    aac_ics* ics = &w0->ics;
    int seq = ics->window_sequence;
    int shape = ics->window_shape;
    int prev_seq = core->usac_cplx_win_seq_prev;
    int prev_shape = core->usac_cplx_win_shape_prev;
    int short_win = (seq == EIGHT_SHORT_SEQUENCE);
    int window_len = short_win ? frame_size / 8 : frame_size;
    int num_windows = short_win ? 8 : 1;
    int max_sfb = st->max_sfb_ste;
    const uint16_t* swb = ics->swb_offset;
    int stride = max_sfb;
    int sign = st->pred_dir ? -1 : 1;
    float dmx_re[AAC_CORE_MAX_FRAME_SIZE];
    float dmx_im[AAC_CORE_MAX_FRAME_SIZE];
    float dmx_prev[AAC_CORE_MAX_FRAME_SIZE];
    const float* filt_cur;
    const float* filt_prev;
    int shape_off;
    int g, gw, band, i, w = 0;

    if (prev_shape == 0 && shape == 0)
        shape_off = 0;
    else if (prev_shape == 1 && shape == 1)
        shape_off = 1;
    else if (prev_shape == 0 && shape == 1)
        shape_off = 2;
    else
        shape_off = 3;

    if (seq == LONG_START_SEQUENCE)
        filt_cur = usac_mdst_cur[3][shape_off]; /* stop-start table */
    else if (seq == LONG_STOP_SEQUENCE)
        filt_cur = usac_mdst_cur[2][shape_off];
    else
        filt_cur = usac_mdst_cur[0][shape_off];

    if (seq == LONG_START_SEQUENCE &&
        prev_seq != EIGHT_SHORT_SEQUENCE && prev_seq != LONG_START_SEQUENCE)
        filt_cur = usac_mdst_cur[1][shape_off];

    filt_prev = (seq == LONG_STOP_SEQUENCE ||
                 (seq == LONG_START_SEQUENCE &&
                  (prev_seq == EIGHT_SHORT_SEQUENCE ||
                   prev_seq == LONG_START_SEQUENCE)))
                    ? usac_mdst_prev[1][prev_shape]
                    : usac_mdst_prev[0][prev_shape];

    memset(dmx_re, 0, sizeof(dmx_re));
    memset(dmx_im, 0, sizeof(dmx_im));
    memset(dmx_prev, 0, sizeof(dmx_prev));

    /* Build the previous frame's downmix used by the temporal MDST term.  For
     * a short current frame the last window of the previous frame is used. */
    if (st->complex_coef && st->use_prev_frame) {
        int off = short_win ? window_len * 7 : 0;
        int n = core->usac_cplx_prev_len;
        for (i = 0; i < window_len && off + i < n; i++)
            dmx_prev[i] = 0.5f * (core->usac_cplx_prev_l[off + i] +
                                  (float)sign * core->usac_cplx_prev_r[off + i]);
    }

    /* Remember this frame's spectra for the next frame's MDST prediction.
     * This must happen before the prediction is applied. */
    if (st->complex_coef) {
        memcpy(core->usac_cplx_prev_l, w0->spec,
               sizeof(float) * (size_t)window_len * (size_t)num_windows);
        memcpy(core->usac_cplx_prev_r, w1->spec,
               sizeof(float) * (size_t)window_len * (size_t)num_windows);
        core->usac_cplx_prev_len = window_len * num_windows;
    }

    for (g = 0; g < ics->num_groups; g++) {
        for (gw = 0; gw < ics->group_len[g]; gw++, w++) {
            int base = w * window_len;
            int max_off = swb[max_sfb];
            int ext;

            for (band = 0; band < max_sfb; band++) {
                for (i = swb[band]; i < swb[band + 1]; i++) {
                    if (st->pred_used[g * stride + band])
                        dmx_re[base + i] = w0->spec[base + i];
                    else
                        dmx_re[base + i] =
                            0.5f * (w0->spec[base + i] +
                                    (float)sign * w1->spec[base + i]);
                }
            }
            for (i = max_off; i < window_len; i++)
                dmx_re[base + i] = 0.0f;

            if (st->complex_coef) {
                ext = max_off + 7;
                if (ext > window_len)
                    ext = window_len;
                if (short_win && w != 0) {
                    usac_mdst_filter(&dmx_im[base], &dmx_re[base], filt_cur,
                                     ext, 1, 1);
                    usac_mdst_filter(&dmx_im[base], &dmx_re[base - window_len],
                                     filt_prev, ext, -1, 1);
                } else {
                    usac_mdst_filter(&dmx_im[base], dmx_re, filt_cur, ext, 1,
                                     1);
                    if (st->use_prev_frame)
                        usac_mdst_filter(&dmx_im[base], dmx_prev, filt_prev,
                                         ext, -1, 1);
                }
            }

            for (band = 0; band < max_sfb; band++) {
                float are, aim;

                if (!st->pred_used[g * stride + band])
                    continue;
                are = (float)st->alpha_q_re[g * stride + band] * 0.1f;
                aim = (float)st->alpha_q_im[g * stride + band] * 0.1f;
                for (i = swb[band]; i < swb[band + 1]; i++) {
                    int k = base + i;
                    float L = w0->spec[k];
                    float R = w1->spec[k];
                    float side = R - are * dmx_re[k] - aim * dmx_im[k];

                    w0->spec[k] = L + side;
                    w1->spec[k] = (L - side) * (float)sign;
                }
            }
        }
    }
}

/* Decode one UsacSingleChannelElement / UsacChannelPairElement payload
 * (UsacCoreCoderData).  SBR data parsing, noise filling, TNS, M/S and the
 * filterbank are applied here as well, mirroring the USAC syntax order. */
static int usac_core_coder(aac_core* core, aac_bits* b,
                           aac_usac_elem_config* e, int ch0, int nb_channels,
                           int indep, int elem_index)
{
    aac_usac_config* uc = &core->usac.config;
    aac_frame_work* frame = &core->work;
    int frame_size = uc->core_frame_len;
    int sr_index = core->cfg.sample_rate_index;
    usac_cpe_state cpe;
    int core_nb = nb_channels;
    int ch, ret;

    memset(&cpe, 0, sizeof(cpe));

    for (ch = 0; ch < nb_channels; ch++) {
        aac_channel_work* w = &frame->ch[ch0 + ch];
        memset(w, 0, sizeof(*w));
        w->element_type = (nb_channels == 2) ? ELEMENT_CPE : ELEMENT_SCE;
        cpe.core_mode[ch] = (int)aac_bits_read(b, 1);
        if (b->error)
            return AAC_CORE_NEED_MORE;
    }

    if (nb_channels > 1 && e->stereo_config_index == 1)
        core_nb = 1; /* MPS212 carries a single core channel */

    if (core_nb == 2) {
        ret = usac_stereo_info(core, b, ch0, indep, &cpe);
        if (ret < 0) {
            USDBG("usac: stereo info failed (ret=%d)\n", ret);
            return ret;
        }
    }

    for (ch = 0; ch < core_nb; ch++) {
        aac_channel_work* w = &frame->ch[ch0 + ch];
        aac_ics* ics = &w->ics;
        int tns_data_present;
        int global_gain;
        int arith_reset;
        int win;
        int sfg = -1;

        if (cpe.core_mode[ch]) {
            USDBG("usac: ch%d core_mode LPD\n", ch);
            return AAC_CORE_ERROR; /* LPD is not supported */
        }

        if (core_nb == 1 || cpe.core_mode[0] != cpe.core_mode[1])
            tns_data_present = (int)aac_bits_read(b, 1);
        else
            tns_data_present = cpe.tns_data_present[ch];

        global_gain = (int)aac_bits_read(b, 8);
        USDBG("usac: gain ch=%d gg=%d bit=%zu\n", ch, global_gain, b->bitpos);
        if (e->noise_fill) {
            w->usac_noise_level = (int)aac_bits_read(b, 3);
            w->usac_noise_offset = (int)aac_bits_read(b, 5);
        }

        if (!(nb_channels == 2 && cpe.common_window)) {
            /* fd_channel_ics_info() */
            ics->window_sequence = (int)aac_bits_read(b, 2);
            ics->window_shape = (int)aac_bits_read(b, 1);
            if (ics->window_sequence == EIGHT_SHORT_SEQUENCE) {
                ics->max_sfb = (int)aac_bits_read(b, 4);
                sfg = (int)aac_bits_read(b, 7);
            } else {
                ics->max_sfb = (int)aac_bits_read(b, 6);
            }
            if (usac_ics_setup(w, sr_index, frame_size, sfg) < 0) {
                USDBG("usac: ch%d ics setup failed (max_sfb=%d)\n", ch, ics->max_sfb);
                return AAC_CORE_ERROR;
            }
        }

        if (e->tw_mdct && (int)aac_bits_read(b, 1))
            return AAC_CORE_ERROR; /* time warping is not supported */

        if (usac_decode_scalefactors(b, w, global_gain) < 0) {
            USDBG("usac: ch%d scalefactor decode failed\n", ch);
            return AAC_CORE_ERROR;
        }
        USDBG("usac: ch%d post-sf bit=%zu\n", ch, b->bitpos);

        if (tns_data_present) {
            if (usac_decode_tns(b, w) < 0) {
                USDBG("usac: ch%d tns decode failed\n", ch);
                return AAC_CORE_ERROR;
            }
            ics->tns_present = 1;
        }

        /* ac_spectral_data */
        arith_reset = indep;
        if (!arith_reset)
            arith_reset = (int)aac_bits_read(b, 1);
        USDBG("usac: ch%d pre-arith bit=%zu seq=%d sfb=%d ng=%d\n", ch,
              b->bitpos, ics->window_sequence, ics->max_sfb, ics->num_groups);

        memset(w->spec, 0, sizeof(w->spec));
        for (win = 0; win < ics->num_windows; win++) {
            int lg = ics->swb_offset[ics->max_sfb];
            int N = (ics->window_sequence == EIGHT_SHORT_SEQUENCE)
                        ? frame_size / 8
                        : frame_size;
#ifdef AAC_USAC_DEBUG
            size_t wstart = b->bitpos;
#endif
            if (aac_usac_ac_decode(&core->usac_ac[ch0 + ch],
                                   b, w->spec + win * 128, lg, N,
                                   arith_reset && win == 0) < 0) {
                USDBG("usac: ch%d win%d arith decode failed (lg=%d N=%d)\n", ch, win, lg, N);
                return AAC_CORE_ERROR;
            }
#ifdef AAC_USAC_DEBUG
            fprintf(stderr, "usac: ch%d win%d arith bits=%zu lg=%d\n", ch,
                    win, b->bitpos - wstart, lg);
#endif
        }

        if ((int)aac_bits_read(b, 1) /* fac_data_present */) {
            int fac_len = (ics->window_sequence == EIGHT_SHORT_SEQUENCE)
                              ? frame_size / 16
                              : frame_size / 8;
            if (usac_parse_fac_data(b, fac_len) < 0) {
                USDBG("usac: ch%d fac data parse failed\n", ch);
                return AAC_CORE_ERROR;
            }
        }
        USDBG("usac: ch%d post-fac bit=%zu\n", ch, b->bitpos);
        if (b->error)
            return AAC_CORE_NEED_MORE;
    }

    /* sbrData() */
    if (e->sbr_ratio > 0) {
        int sbr_ch = nb_channels;
        aac_sbr_decoder* sbr;
        aac_sbr_usac_defaults dflt;

        if (nb_channels == 2 && !(e->stereo_config_index == 0 ||
                                  e->stereo_config_index == 3))
            sbr_ch = 1;

        dflt.start_freq = e->sbr_start_freq;
        dflt.stop_freq = e->sbr_stop_freq;
        dflt.freq_scale = e->sbr_freq_scale;
        dflt.alter_scale = e->sbr_alter_scale;
        dflt.noise_bands = e->sbr_noise_bands;
        dflt.limiter_bands = e->sbr_limiter_bands;
        dflt.limiter_gains = e->sbr_limiter_gains;
        dflt.interpol_freq = e->sbr_interpol_freq;
        dflt.smoothing_mode = e->sbr_smoothing_mode;

        sbr = ensure_sbr(core, ch0, nb_channels);
        if (sbr == NULL ||
            aac_sbr_parse_usac(sbr, b, sbr_ch, indep, &dflt) < 0) {
            USDBG("usac: sbr data parse failed (sbr_ch=%d indep=%d)\n", sbr_ch, indep);
            return AAC_CORE_ERROR;
        }
        frame->element_sbr[elem_index] = 1;
    }

    /* Noise filling, scaling, TNS and M/S, then the filterbank. */
    for (ch = 0; ch < core_nb; ch++)
        usac_spectrum_scale(core, &frame->ch[ch0 + ch], frame_size);

    if (nb_channels == 2 && cpe.common_window) {
        aac_channel_work* w0 = &frame->ch[ch0];
        aac_channel_work* w1 = &frame->ch[ch0 + 1];
        int save = w0->ics.max_sfb;

        if (w0->ics.tns_present && !cpe.tns_on_lr) {
            apply_tns(w0, frame_size);
            apply_tns(w1, frame_size);
        }
        if (cpe.ms_mask_mode == 3 && cpe.cplx_pred) {
            usac_apply_cplx(core, w0, w1, frame_size, &cpe);
        } else if (cpe.ms_mask_mode > 0) {
            w0->ics.max_sfb = cpe.max_sfb_ste;
            apply_ms(w0, w1, frame_size);
            w0->ics.max_sfb = save;
        }
        if (w0->ics.tns_present && cpe.tns_on_lr) {
            apply_tns(w0, frame_size);
            apply_tns(w1, frame_size);
        }
    } else {
        for (ch = 0; ch < core_nb; ch++) {
            aac_channel_work* w = &frame->ch[ch0 + ch];
            if (w->ics.tns_present)
                apply_tns(w, frame_size);
        }
    }

    for (ch = 0; ch < nb_channels; ch++) {
        aac_channel_work* w = &frame->ch[ch0 + ch];
        synth_channel(&core->state[ch0 + ch], w->spec,
                      w->ics.window_sequence, w->ics.window_shape,
                      frame_size, w->pcm, 0);
    }

    return b->error ? AAC_CORE_NEED_MORE : AAC_CORE_OK;
}

/* UsacExtElement(): filler and UniDrc payloads are skipped, audio preroll
 * frames are not supported. */
static int usac_ext_element(aac_core* core, aac_bits* b,
                            aac_usac_elem_config* e)
{
    int len;
    int start = 1;
    int stop = 1;

    (void)core;
    if (!(int)aac_bits_read(b, 1) /* usacExtElementPresent */)
        return b->error ? AAC_CORE_NEED_MORE : AAC_CORE_OK;
    USDBG("usac: ext present type=%d bit=%zu\n", e->ext_type, b->bitpos);

    if ((int)aac_bits_read(b, 1) /* useDefaultLength */) {
        len = e->ext_default_len;
    } else {
        len = (int)aac_bits_read(b, 8);
        if (len == 255)
            len += (int)aac_bits_read(b, 16) - 2;
    }
    if (b->error)
        return AAC_CORE_NEED_MORE;
    if (len <= 0)
        return AAC_CORE_OK;
    USDBG("usac: ext len=%d bit=%zu\n", len, b->bitpos);

    if (e->ext_payload_frag) {
        start = (int)aac_bits_read(b, 1);
        stop = (int)aac_bits_read(b, 1);
        if (b->error)
            return AAC_CORE_NEED_MORE;
    }
    if (!(start && stop))
        return AAC_CORE_ERROR; /* fragmented payloads are not supported */

    /* FILL, UniDrc, elevation, audio preroll and unknown payloads carry no
     * PCM data for the current frame; the payload is skipped.  (Audio preroll
     * frames only prime the decoder state for a clean start, so skipping them
     * is a valid decoder choice.) */
    aac_bits_skip(b, (size_t)len * 8u);
    return b->error ? AAC_CORE_NEED_MORE : AAC_CORE_OK;
}

/* Decode one UsacFrame(). */
static int usac_decode_frame(aac_core* core, const uint8_t* data, size_t len,
                             size_t* consumed, float* pcm, int* channels,
                             int* sample_rate)
{
    aac_usac_config* uc = &core->usac.config;
    aac_frame_work* frame = &core->work;
    aac_bits b;
    int indep;
    int ch_slot = 0;
    int i, e;
    int core_F, F, out_channels;
    int covered[AAC_CORE_MAX_CHANNELS];
    int o, src;

    aac_bits_init(&b, data, len);
    memset(frame, 0, sizeof(*frame));
    frame->frame_size = core->cfg.frame_size;
    frame->sample_rate = core->cfg.sample_rate;

    indep = (int)aac_bits_read(&b, 1);
    if (b.error)
        return AAC_CORE_NEED_MORE;
    USDBG("usac: frame start indep=%d elems=%d len=%zu bytes=%02X %02X %02X "
          "%02X\n",
          indep, uc->nb_elements, len, data[0], data[1], data[2], data[3]);

    for (i = 0; i < uc->nb_elements; i++) {
        aac_usac_elem_config* elem = &uc->elems[i];
        int ret;

        switch (elem->type) {
        case AAC_USAC_ID_SCE:
        case AAC_USAC_ID_LFE: {
            int idx = frame_add_element(frame, ELEMENT_SCE, 0, 1);
            if (idx < 0)
                return AAC_CORE_ERROR;
            ret = usac_core_coder(core, &b, elem, ch_slot, 1, indep, idx);
            if (ret < 0)
                return ret;
            ch_slot += 1;
            break;
        }
        case AAC_USAC_ID_CPE: {
            int idx = frame_add_element(frame, ELEMENT_CPE, 0, 2);
            if (idx < 0)
                return AAC_CORE_ERROR;
            ret = usac_core_coder(core, &b, elem, ch_slot, 2, indep, idx);
            if (ret < 0)
                return ret;
            ch_slot += 2;
            break;
        }
        case AAC_USAC_ID_EXT:
            ret = usac_ext_element(core, &b, elem);
            if (ret < 0)
                return ret;
            break;
        default:
            return AAC_CORE_ERROR;
        }
    }

    USDBG("usac: frame done, consumed=%zu bits err=%d\n", (size_t)b.bitpos,
          b.error);
    *consumed = aac_bits_bytes_consumed(&b);
    if (b.error)
        return AAC_CORE_NEED_MORE;

    /* Assemble the output, applying SBR upsampling where configured. */
    core_F = core->cfg.frame_size;
    F = core->cfg.has_sbr ? core_F * 2 : core_F;
    out_channels = core->cfg.channels;
    if (out_channels > frame->channels || out_channels <= 0)
        out_channels = frame->channels;
    memset(covered, 0, sizeof(covered));

    if (core->cfg.has_sbr) {
        for (e = 0; e < frame->element_count; e++) {
            int first = frame->element_first_channel[e];
            int count = frame->element_channels[e];
            aac_sbr_decoder* sbr;
            const float* srcp[2];
            float* dst[2];
            int c;

            if (first < 0 || first + count > frame->channels ||
                first >= AAC_CORE_MAX_CHANNELS)
                continue;
            sbr = core->sbr[first];
            srcp[0] = frame->ch[first].pcm;
            srcp[1] = count > 1 ? frame->ch[first + 1].pcm : NULL;
            dst[0] = core->pcm_hi[first];
            dst[1] = count > 1 ? core->pcm_hi[first + 1] : NULL;
            if (sbr != NULL &&
                aac_sbr_apply(sbr, srcp, core_F, dst, F) == 0) {
                for (c = 0; c < count; c++)
                    covered[first + c] = 1;
            }
        }
        /* Elements without SBR state are upsampled by linear interpolation. */
        for (i = 0; i < frame->channels; i++) {
            int c;
            if (covered[i])
                continue;
            for (c = 0; c < core_F; c++) {
                float v = frame->ch[i].pcm[c];
                float n = (c + 1 < core_F) ? frame->ch[i].pcm[c + 1] : v;
                core->pcm_hi[i][2 * c] = v;
                core->pcm_hi[i][2 * c + 1] = 0.5f * (v + n);
            }
        }
    } else {
        for (i = 0; i < frame->channels; i++)
            memcpy(core->pcm_hi[i], frame->ch[i].pcm,
                   (size_t)core_F * sizeof(float));
    }

    memset(pcm, 0, (size_t)out_channels * (size_t)F * sizeof(float));
    for (o = 0; o < out_channels; o++) {
        src = core->mpeg_order ? o : core->cfg.output_order[o];
        if (src < 0 || src >= frame->channels)
            src = o;
        for (i = 0; i < F; i++)
            pcm[(size_t)i * out_channels + o] = core->pcm_hi[src][i];
    }
    *channels = out_channels;
    *sample_rate = aac_core_output_sample_rate(core);

    core->frames++;
    return AAC_CORE_OK;
}

int aac_core_init(aac_core* core)
{
    if (core == NULL)
        return -1;
    memset(core, 0, sizeof(*core));
    config_clear(&core->cfg);
    core->pns_seed = 0x1f2e3d4cu;
    if (aac_tables_init() < 0)
        return -1;
    if (aac_huffman_init() < 0)
        return -1;
    synth_init();
    return 0;
}

void aac_core_close(aac_core* core)
{
    if (core == NULL)
        return;
    free_sbr_states(core);
}

void aac_core_reset(aac_core* core)
{
    if (core == NULL)
        return;
    memset(core->state, 0, sizeof(core->state));
    memset(core->usac_ac, 0, sizeof(core->usac_ac));
    memset(core->usac_cplx_alpha_re_prev, 0,
           sizeof(core->usac_cplx_alpha_re_prev));
    memset(core->usac_cplx_alpha_im_prev, 0,
           sizeof(core->usac_cplx_alpha_im_prev));
    core->usac_cplx_valid = 0;
    core->usac_cplx_prev_len = 0;
    core->pns_seed = 0x1f2e3d4cu;
    core->frames = 0;
    free_sbr_states(core);
}

int aac_core_frame_size(const aac_core* core)
{
    return core ? core->cfg.frame_size : 0;
}

int aac_core_output_frame_size(const aac_core* core)
{
    if (core == NULL)
        return 0;
    return core->cfg.has_sbr ? core->cfg.frame_size * 2 : core->cfg.frame_size;
}

int aac_core_output_sample_rate(const aac_core* core)
{
    if (core == NULL)
        return 0;
    if (core->cfg.has_sbr)
        return core->cfg.ext_sample_rate > 0 ? core->cfg.ext_sample_rate :
                                               2 * core->cfg.sample_rate;
    return core->cfg.sample_rate;
}

int aac_core_set_adts(aac_core* core, int profile, int sf_index,
                      int channel_config, int mpeg2)
{
    aac_config* c;
    int old_rate, old_frame, old_cc;

    if (core == NULL)
        return -1;
    c = &core->cfg;
    old_rate = c->sample_rate;
    old_frame = c->frame_size;
    old_cc = c->channel_config;

    if (profile != 1) /* only LC */
        return -1;
    if (sf_index < 0 || sf_index >= AAC_SR_INDEX_MAX)
        return -1;

    c->aot = 2;
    c->core_aot = 2;
    c->profile = mpeg2 ? profile : -1;
    c->mpeg2 = mpeg2 ? 1 : 0;
    c->sample_rate_index = sf_index;
    c->sample_rate = aac_sample_rates[sf_index];
    c->channel_config = channel_config;
    c->frame_length_flag = 0;
    c->frame_size = 1024;
    c->short_len = 128;
    if (channel_config != 0 || old_cc != 0)
        c->pce_valid = 0;

    return config_apply(core,
                        old_rate != c->sample_rate ||
                        old_frame != c->frame_size ||
                        old_cc != c->channel_config);
}

int aac_core_set_adif_pce(aac_core* core, int profile, int sf_index,
                          const uint8_t* pce, size_t pce_len)
{
    aac_bits b;
    int ret;
    aac_config* c;

    if (core == NULL || pce == NULL)
        return -1;
    c = &core->cfg;
    c->aot = profile + 1;
    c->core_aot = profile + 1;
    c->profile = profile;
    c->mpeg2 = 1;
    c->sample_rate_index = sf_index;
    c->sample_rate = aac_sample_rates[sf_index];
    c->frame_length_flag = 0;
    c->frame_size = 1024;
    c->short_len = 128;

    aac_bits_init(&b, pce, pce_len);
    ret = read_pce(core, &b, 1);
    if (ret < 0)
        return -1;
    return config_apply(core, 1);
}

static int parse_asc_bits(aac_core* core, aac_bits* b, int configure)
{
    aac_config* c;
    int aot, sf_index, channel_config;
    int original_aot;
    int ext_index = -1;
    int ext_rate = 0;
    int frame_flag;
    int depends, extension;
    int ret;
    int old_rate, old_frame, old_cc;

    if (core == NULL || b == NULL)
        return -1;

    c = &core->cfg;
    old_rate = c->sample_rate;
    old_frame = c->frame_size;
    old_cc = c->channel_config;

    aot = (int)aac_bits_read(b, 5);
    if (aot == 31)
        aot = 32 + (int)aac_bits_read(b, 6);
    original_aot = aot;
    sf_index = (int)aac_bits_read(b, 4);
    if (sf_index == 15) {
        int rate = (int)aac_bits_read(b, 24);
        sf_index = find_sr_index(rate);
        c->sample_rate = rate;
    } else if (sf_index >= 0 && sf_index < AAC_SR_INDEX_MAX) {
        c->sample_rate = aac_sample_rates[sf_index];
    } else {
        return -1;
    }
    channel_config = (int)aac_bits_read(b, 4);

    c->ext_aot = 0;
    c->has_sbr = 0;
    c->has_ps = 0;
    c->ext_sample_rate = 0;
    c->drc_present = 0;
    c->drc_pce_instance_tag = -1;
    c->drc_prog_ref_present = 0;
    c->drc_prog_ref_level = -1;
    c->drc_band_incr = 0;
    c->drc_interpolation_scheme = 0;
    c->drc_num_bands = 1;
    memset(c->drc_dyn_rng_sgn, 0, sizeof(c->drc_dyn_rng_sgn));
    memset(c->drc_dyn_rng_ctl, 0, sizeof(c->drc_dyn_rng_ctl));

    if (aot == 5 || aot == 29) {
        ext_index = (int)aac_bits_read(b, 4);
        if (ext_index == 15)
            ext_rate = (int)aac_bits_read(b, 24);
        else if (ext_index < AAC_SR_INDEX_MAX)
            ext_rate = aac_sample_rates[ext_index];
        else
            return -1;
        aot = (int)aac_bits_read(b, 5);
        if (aot == 31)
            aot = 32 + (int)aac_bits_read(b, 6);
        c->ext_aot = original_aot;
        c->has_sbr = 1;
        c->has_ps = (original_aot == 29);
        c->ext_sample_rate = ext_rate;
    }

    if (aot != 2 && aot != 17 && aot != 23 && aot != 39 && aot != 42) {
        /* Only AAC-LC, ER AAC-LC, AAC-LD, AAC-ELD and USAC are implemented. */
        return -1;
    }

    if (aot == 42) {
        /* USAC: AudioSpecificConfig carries UsacConfig() after the common
         * 4 bit samplingFrequencyIndex and channelConfiguration fields. */
        aac_usac_config* uc;
        int i;

        core->usac_valid = 0;
        ret = aac_usac_parse_asc_bits_state(b, &core->usac);
        if (ret < 0)
            return -1;
        uc = &core->usac.config;
        if (uc->core_frame_len != 1024)
            return -1; /* 768 frame length is not supported yet */
        if (uc->sbr_ratio != 0 && uc->sbr_ratio != 3)
            return -1; /* only the 2x SBR ratio is supported */

        c->aot = 42;
        c->core_aot = 42;
        c->profile = -1;
        c->mpeg2 = 0;
        c->sample_rate_index = sf_index;
        c->channel_config = channel_config;
        c->pce_valid = 0;
        c->ep_config = 0;
        c->frame_length_flag = 0;
        c->frame_size = uc->core_frame_len;
        c->short_len = c->frame_size / 8;
        if (uc->sbr_ratio > 0) {
            c->has_sbr = 1;
            c->ext_sample_rate = uc->sample_rate;
            c->sample_rate = uc->sample_rate / 2;
        } else {
            c->has_sbr = 0;
            c->ext_sample_rate = 0;
            c->sample_rate = uc->sample_rate;
        }
        c->sample_rate_index = find_sr_index(c->sample_rate);

        if (uc->channel_config_index >= 1 && uc->channel_config_index <= 7) {
            c->channel_config = uc->channel_config_index;
        } else {
            /* Custom (index 0) and 22.2 layouts are not mapped yet. */
            return -1;
        }

        /* Reset the USAC persistent state and initialize the noise seeds
         * per element channel. */
        memset(core->usac_ac, 0, sizeof(core->usac_ac));
        memset(core->usac_noise_seed, 0, sizeof(core->usac_noise_seed));
        {
            int slot = 0;
            for (i = 0; i < uc->nb_elements && slot < AAC_CORE_MAX_CHANNELS;
                 i++) {
                if (uc->elems[i].type == AAC_USAC_ID_EXT)
                    continue;
                core->usac_noise_seed[slot] = 0x3039u;
                if (uc->elems[i].type == AAC_USAC_ID_CPE &&
                    slot + 1 < AAC_CORE_MAX_CHANNELS)
                    core->usac_noise_seed[slot + 1] = 0x10932u;
                slot += (uc->elems[i].type == AAC_USAC_ID_CPE) ? 2 : 1;
            }
        }
        core->usac_valid = 1;

        if (!configure)
            return AAC_CORE_OK;

        return config_apply(core,
                            old_rate != c->sample_rate ||
                                old_frame != c->frame_size ||
                                old_cc != c->channel_config);
    }

    c->aot = aot;
    c->core_aot = aot;
    c->profile = -1;
    c->mpeg2 = 0;
    c->sample_rate_index = sf_index;
    c->channel_config = channel_config;
    c->pce_valid = 0;
    c->ep_config = 0;

    if (aot == 39) {
        int ldsbr;
        frame_flag = (int)aac_bits_read(b, 1);
        c->frame_length_flag = frame_flag;
        c->frame_size = frame_flag ? 480 : 512;
        c->short_len = c->frame_size / 8;

        {
            int res_flags = (int)aac_bits_read(b, 3);
            ldsbr = (int)aac_bits_read(b, 1);
            if (res_flags)
                return -1;
        }
        if (ldsbr) {
            /* LD-SBR is not implemented yet. */
            return -1;
        }
        for (;;) {
            int ext_type = (int)aac_bits_read(b, 4);
            int len;
            if (ext_type == 0)
                break;
            len = (int)aac_bits_read(b, 4);
            if (len == 15)
                len += (int)aac_bits_read(b, 8);
            if (len == 15 + 255)
                len += (int)aac_bits_read(b, 16);
            aac_bits_skip(b, (size_t)len * 8u);
            if (b->error)
                return -1;
        }
        c->ep_config = (int)aac_bits_read(b, 2);
        if (c->ep_config != 0 || b->error)
            return -1;
    } else {
        frame_flag = (int)aac_bits_read(b, 1);
        depends = (int)aac_bits_read(b, 1);
        if (depends)
            aac_bits_read(b, 14);
        extension = (int)aac_bits_read(b, 1);

        c->frame_length_flag = frame_flag;
        if (aot == 23) {
            c->frame_size = frame_flag ? 480 : 512;
        } else {
            c->frame_size = frame_flag ? 960 : 1024;
        }
        c->short_len = c->frame_size / 8;

        if (b->error)
            return -1;

        if (channel_config == 0) {
            ret = read_pce(core, b, 1);
            if (ret < 0)
                return -1;
        }

        if (aot == 17 || aot == 23) {
            if (extension) {
                int res_flags = (int)aac_bits_read(b, 3);
                (void)aac_bits_read(b, 1); /* extensionFlag3 */
                if (res_flags)
                    return -1;
            }
            c->ep_config = (int)aac_bits_read(b, 2);
            if (c->ep_config != 0)
                return -1;
            if (b->error)
                return -1;
        }
    }

    if (!configure)
        return AAC_CORE_OK;

    return config_apply(core,
                        old_rate != c->sample_rate ||
                        old_frame != c->frame_size ||
                        old_cc != c->channel_config);
}

int aac_core_set_asc(aac_core* core, const uint8_t* asc, size_t len)
{
    aac_bits b;
    if (core == NULL || asc == NULL || len == 0)
        return -1;
    aac_bits_init(&b, asc, len);
    return parse_asc_bits(core, &b, 1);
}

int aac_core_parse_asc_bits(aac_core* core, void* bits_ctx)
{
    return parse_asc_bits(core, (aac_bits*)bits_ctx, 1);
}

int aac_core_parse_pce_bits(aac_core* core, void* bits_ctx)
{
    if (core == NULL || bits_ctx == NULL)
        return -1;
    return read_pce(core, (aac_bits*)bits_ctx, 0);
}

int aac_core_apply_config(aac_core* core, int reset)
{
    if (core == NULL)
        return -1;
    return config_apply(core, reset);
}

/* ------------------------------------------------------------------------- */
/* Frame decoding                                                            */
/* ------------------------------------------------------------------------- */

int aac_core_decode(aac_core* core, const uint8_t* data, size_t len,
                    size_t* consumed, float* pcm, int* channels,
                    int* sample_rate)
{
    aac_bits b;
    aac_frame_work* frame;
    int ret;
    int e;

    if (core == NULL || data == NULL || consumed == NULL || pcm == NULL ||
        channels == NULL || sample_rate == NULL)
        return AAC_CORE_ERROR;
    if (!core->configured || core->cfg.sample_rate <= 0)
        return AAC_CORE_ERROR;

    if (core->usac_valid) {
        USDBG("usac: core decode route (aot=%d len=%zu)\n", core->cfg.aot,
              len);
        return usac_decode_frame(core, data, len, consumed, pcm, channels,
                                 sample_rate);
    }
    USDBG("core decode: aot=%d usac_valid=%d\n", core->cfg.aot,
          core->usac_valid);

    aac_bits_init(&b, data, len);
    frame = &core->work;
    ret = decode_raw_data_block(core, &b, frame);
    if (ret < 0)
        return ret;
    if (b.error)
        return AAC_CORE_NEED_MORE;

    *consumed = aac_bits_bytes_consumed(&b);

    /* A PCE inside the first raw_data_block defines the channel layout. */
    if (core->cfg.channel_config == 0 && core->cfg.pce_valid) {
        if (core->cfg.channels <= 0)
            return AAC_CORE_ERROR;
    }
    if (core->cfg.channels <= 0)
        return AAC_CORE_ERROR;

    /* Apply the per-element tools in the standard order: PNS, M/S, I/S,
     * TNS, filterbank. */
    for (e = 0; e < frame->element_count; e++) {
        int first = frame->element_first_channel[e];
        int count = frame->element_channels[e];
        int c0, c1;

        if (count == 2 && first + 1 < frame->channels) {
            c0 = first;
            c1 = first + 1;
            apply_pns_cpe(core, &frame->ch[c0], &frame->ch[c1],
                          core->cfg.frame_size);
            if (frame->ch[c0].ics.ms_mask_present)
                apply_ms(&frame->ch[c0], &frame->ch[c1],
                         core->cfg.frame_size);
            apply_is(&frame->ch[c0], &frame->ch[c1], core->cfg.frame_size);
            apply_tns(&frame->ch[c0], core->cfg.frame_size);
            apply_tns(&frame->ch[c1], core->cfg.frame_size);
            if (core->cfg.aot == 39) {
                synth_channel_eld(&core->state[c0], frame->ch[c0].spec,
                                  core->cfg.frame_size, frame->ch[c0].pcm);
                synth_channel_eld(&core->state[c1], frame->ch[c1].spec,
                                  core->cfg.frame_size, frame->ch[c1].pcm);
            } else {
                synth_channel(&core->state[c0], frame->ch[c0].spec,
                              frame->ch[c0].ics.window_sequence,
                              frame->ch[c0].ics.window_shape,
                              core->cfg.frame_size, frame->ch[c0].pcm,
                              core->cfg.aot == 23);
                synth_channel(&core->state[c1], frame->ch[c1].spec,
                              frame->ch[c1].ics.window_sequence,
                              frame->ch[c1].ics.window_shape,
                              core->cfg.frame_size, frame->ch[c1].pcm,
                              core->cfg.aot == 23);
            }
        } else {
            c0 = first;
            apply_pns_sce(core, &frame->ch[c0], core->cfg.frame_size);
            apply_tns(&frame->ch[c0], core->cfg.frame_size);
            if (core->cfg.aot == 39) {
                synth_channel_eld(&core->state[c0], frame->ch[c0].spec,
                                  core->cfg.frame_size, frame->ch[c0].pcm);
            } else {
                synth_channel(&core->state[c0], frame->ch[c0].spec,
                              frame->ch[c0].ics.window_sequence,
                              frame->ch[c0].ics.window_shape,
                              core->cfg.frame_size, frame->ch[c0].pcm,
                              core->cfg.aot == 23);
            }
        }
    }

    {
        int out_channels = core->cfg.channels;
        int core_F = core->cfg.frame_size;
        int use_sbr = core->cfg.has_sbr;
        int F = use_sbr ? core_F * 2 : core_F;
        int covered[AAC_CORE_MAX_CHANNELS];
        int i, o;

        if (core->cfg.has_ps && out_channels == 1)
            out_channels = 2; /* parametric stereo output */
        else if (out_channels > frame->channels)
            out_channels = frame->channels;
        if (out_channels <= 0)
            return AAC_CORE_ERROR;
        memset(covered, 0, sizeof(covered));

        if (use_sbr) {
            for (e = 0; e < frame->element_count; e++) {
                int first = frame->element_first_channel[e];
                int count = frame->element_channels[e];
                aac_sbr_decoder* sbr;
                const float* src[2];
                float* dst[2];
                int c;

                if (first < 0 || first + count > frame->channels ||
                    first >= AAC_CORE_MAX_CHANNELS)
                    continue;
                sbr = ensure_sbr(core, first, count);
                src[0] = frame->ch[first].pcm;
                src[1] = count > 1 ? frame->ch[first + 1].pcm : NULL;
                dst[0] = core->pcm_hi[first];
                dst[1] = count > 1 ? core->pcm_hi[first + 1] :
                         (core->cfg.has_ps && first + 1 < AAC_CORE_MAX_CHANNELS ?
                          core->pcm_hi[first + 1] : NULL);
                if (sbr != NULL &&
                    aac_sbr_apply(sbr, src, core_F, dst, F) == 0) {
                    for (c = 0; c < count; c++)
                        covered[first + c] = 1;
                    if (count == 1 && core->cfg.has_ps &&
                        first + 1 < AAC_CORE_MAX_CHANNELS)
                        covered[first + 1] = 1;
                }
            }
            /* Elements without SBR state (for example an LFE element) are
             * upsampled by a simple half-band interpolation. */
            for (e = 0; e < frame->channels; e++) {
                int c;
                if (covered[e])
                    continue;
                for (c = 0; c < core_F; c++) {
                    float v = frame->ch[e].pcm[c];
                    float n = (c + 1 < core_F) ? frame->ch[e].pcm[c + 1] : v;
                    core->pcm_hi[e][2 * c] = v;
                    core->pcm_hi[e][2 * c + 1] = 0.5f * (v + n);
                }
            }
        } else {
            for (e = 0; e < frame->channels; e++)
                memcpy(core->pcm_hi[e], frame->ch[e].pcm,
                       (size_t)core_F * sizeof(float));
        }

        memset(pcm, 0, (size_t)out_channels * (size_t)F * sizeof(float));
        for (o = 0; o < out_channels; o++) {
            int src;
            if (core->cfg.has_ps && core->cfg.channels == 1)
                src = o;
            else {
                src = core->mpeg_order ? o : core->cfg.output_order[o];
                if (src < 0 || src >= frame->channels)
                    src = o;
            }
            for (i = 0; i < F; i++)
                pcm[(size_t)i * out_channels + o] = core->pcm_hi[src][i];
        }
        *channels = out_channels;
        *sample_rate = aac_core_output_sample_rate(core);
    }

    core->frames++;
    return AAC_CORE_OK;
}
