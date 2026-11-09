/*
 * aacdec - Spectral Band Replication (HE-AAC) decoder.
 *
 * Clean-room implementation of the MPEG-4 SBR signal flow described in
 * ISO/IEC 14496-3 (subpart 4, SBR).  The QMF prototype filters, the SBR
 * Huffman code lengths and the SBR noise table are standardized numeric
 * constants (see aac_sbr_tables.h).  The algorithms (band table construction,
 * envelope/noise parsing, LPC based HF generation, envelope adjustment and
 * the QMF banks) are implemented from the specification.
 */
#include "aac_sbr.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "aac_sbr_tables.h"

#ifdef AAC_USAC_DEBUG
#include <stdio.h>
#define SBRDBG(...) do { fprintf(stderr, __VA_ARGS__); } while (0)
#else
#define SBRDBG(...)
#endif

#define SBR_PI          3.14159265358979323846
#define SBR_EPS         1e-20f
#define SBR_ENV_OFFSET  2       /* ENVELOPE_ADJUSTMENT_OFFSET */
#define SBR_NOISE_OFFSET 6      /* NOISE_FLOOR_OFFSET */

/* The SBR envelope scalefactors refer to the standard 64-band QMF signal
 * scale.  The QMF implementation here keeps the low band at unity gain, so a
 * fixed calibration factor is applied to the generated high band. */
#define SBR_HF_GAIN     0.015625f

enum
{
    SBR_FIXFIX = 0,
    SBR_FIXVAR = 1,
    SBR_VARFIX = 2,
    SBR_VARVAR = 3
};

/* ------------------------------------------------------------------------- */
/* Small helpers                                                             */
/* ------------------------------------------------------------------------- */



static int sbr_max_int(int a, int b)
{
    return a > b ? a : b;
}

static float sbr_exp2i(int x)
{
    if (x > 128)
        return 1.0e30f;
    if (x < -150)
        return 0.0f;
    return ldexpf(1.0f, x);
}

static void sbr_sort_u16(uint16_t* a, int n)
{
    int i, j;
    for (i = 1; i < n; i++) {
        uint16_t key = a[i];
        j = i - 1;
        while (j >= 0 && a[j] > key) {
            a[j + 1] = a[j];
            j--;
        }
        a[j + 1] = key;
    }
}

static void sbr_sort_i16(int16_t* a, int n)
{
    int i, j;
    for (i = 1; i < n; i++) {
        int16_t key = a[i];
        j = i - 1;
        while (j >= 0 && a[j] > key) {
            a[j + 1] = a[j];
            j--;
        }
        a[j + 1] = key;
    }
}

static int sbr_int16_in_table(const int16_t* table, int n, int value)
{
    int i;
    for (i = 0; i < n; i++)
        if (table[i] == value)
            return 1;
    return 0;
}

/* ------------------------------------------------------------------------- */
/* SBR Huffman decoding                                                      */
/* ------------------------------------------------------------------------- */

typedef struct
{
    int16_t child[2];
    int16_t symbol;
} sbr_huff_node;

static sbr_huff_node* sbr_trees[10];
static int sbr_tables_ready;

static int sbr_tree_build(int table)
{
    int symbols = (int)aac_sbr_huff_sizes[table];
    int offset = (int)aac_sbr_huff_offsets[table];
    const uint8_t(*tab)[2] = aac_sbr_huff_tab;
    int i;
    int code = 0;
    int prev_len = 0;
    int count = 1;
    int cap = symbols * 2 + 2;
    sbr_huff_node* nodes;

    for (i = 0; i < table; i++)
        tab += aac_sbr_huff_sizes[i];

    nodes = (sbr_huff_node*)malloc((size_t)cap * sizeof(*nodes));
    if (nodes == NULL)
        return -1;
    memset(nodes, 0, (size_t)cap * sizeof(*nodes));
    nodes[0].child[0] = -1;
    nodes[0].child[1] = -1;
    nodes[0].symbol = -32768;

    for (i = 0; i < symbols; i++) {
        int len = tab[i][1];
        int sym = (int)tab[i][0] + offset;
        int node = 0;
        int bit;

        if (len <= 0 || len > 20) {
            free(nodes);
            return -1;
        }
        if (i == 0)
            code = 0;
        else
            code = (code + 1) << (len - prev_len);
        prev_len = len;

        for (bit = len - 1; bit >= 0; bit--) {
            int b = (code >> bit) & 1;
            int next = nodes[node].child[b];
            if (next < 0) {
                if (count >= cap) {
                    free(nodes);
                    return -1;
                }
                next = count++;
                nodes[next].child[0] = -1;
                nodes[next].child[1] = -1;
                nodes[next].symbol = -32768;
                nodes[node].child[b] = (int16_t)next;
            }
            node = next;
        }
        if (nodes[node].symbol != -32768) {
            free(nodes);
            return -1;
        }
        nodes[node].symbol = (int16_t)sym;
    }

    sbr_trees[table] = nodes;
    return 0;
}

static int sbr_huffman_init(void)
{
    int i;
    if (sbr_tables_ready)
        return 0;
    for (i = 0; i < 10; i++) {
        if (sbr_tree_build(i) < 0) {
            int j;
            for (j = 0; j < i; j++) {
                free(sbr_trees[j]);
                sbr_trees[j] = NULL;
            }
            return -1;
        }
    }
    sbr_tables_ready = 1;
    return 0;
}

static int sbr_huff_decode(aac_bits* b, int table)
{
    int node = 0;
    sbr_huff_node* nodes = sbr_trees[table];

    if (nodes == NULL)
        return 0;
    while (node >= 0 && nodes[node].symbol == -32768) {
        int bit;
        if (aac_bits_left(b) == 0) {
            b->error = 1;
            return 0;
        }
        bit = (int)aac_bits_read(b, 1);
        node = nodes[node].child[bit];
    }
    if (node < 0)
        return 0;
    return nodes[node].symbol;
}

/* ------------------------------------------------------------------------- */
/* QMF modulation tables                                                     */
/* ------------------------------------------------------------------------- */

static float ana_cos_tab[32][320];
static float ana_sin_tab[32][320];
static float syn_cos_tab[640][64];
static float syn_sin_tab[640][64];
static int qmf_tables_ready;

static int sbr_qmf_tables_init(void)
{
    int k, n;
    if (qmf_tables_ready)
        return 0;
    for (k = 0; k < 32; k++) {
        for (n = 0; n < 320; n++) {
            double phase = SBR_PI / 32.0 * ((double)k + 0.5) *
                           ((double)n - 16.0);
            ana_cos_tab[k][n] = (float)((double)qmf_ds[n] * cos(phase));
            ana_sin_tab[k][n] = (float)((double)qmf_ds[n] * sin(phase));
        }
    }
    for (n = 0; n < 640; n++) {
        for (k = 0; k < 64; k++) {
            double phase = SBR_PI / 64.0 * ((double)k + 0.5) *
                           ((double)n - 32.0);
            syn_cos_tab[n][k] = (float)((double)qmf_h[n] * cos(phase));
            syn_sin_tab[n][k] = (float)((double)qmf_h[n] * sin(phase));
        }
    }
    qmf_tables_ready = 1;
    return 0;
}

/* ------------------------------------------------------------------------- */
/* Frequency band tables                                                     */
/* ------------------------------------------------------------------------- */

static void sbr_make_bands(int16_t* bands, int start, int stop, int num_bands)
{
    int k;
    int previous = start;
    float base = powf((float)stop / (float)start, 1.0f / (float)num_bands);
    float prod = (float)start;

    for (k = 0; k < num_bands - 1; k++) {
        int present;
        prod *= base;
        present = (int)lrintf(prod);
        bands[k] = (int16_t)(present - previous);
        previous = present;
    }
    bands[num_bands - 1] = (int16_t)(stop - previous);
}

static int sbr_check_n_master(int n_master, int xover_band)
{
    return n_master > 0 && xover_band < n_master;
}

static int sbr_make_f_master(aac_sbr_decoder* sbr)
{
    const int8_t* off;
    unsigned int temp, max_qmf_subbands;
    unsigned int start_min, stop_min;
    int16_t stop_dk[13];
    int k;
    const aac_sbr_header* hdr = &sbr->hdr;

    switch (sbr->out_rate) {
    case 16000:
        off = aac_sbr_offset[0];
        break;
    case 22050:
        off = aac_sbr_offset[1];
        break;
    case 24000:
        off = aac_sbr_offset[2];
        break;
    case 32000:
        off = aac_sbr_offset[3];
        break;
    case 44100:
    case 48000:
    case 64000:
        off = aac_sbr_offset[4];
        break;
    case 88200:
    case 96000:
    case 128000:
    case 176400:
    case 192000:
        off = aac_sbr_offset[5];
        break;
    default:
        return -1;
    }

    if (sbr->out_rate < 32000)
        temp = 3000;
    else if (sbr->out_rate < 64000)
        temp = 4000;
    else
        temp = 5000;

    start_min = ((temp << 7) + (sbr->out_rate >> 1)) / (unsigned)sbr->out_rate;
    stop_min = ((temp << 8) + (sbr->out_rate >> 1)) / (unsigned)sbr->out_rate;

    sbr->k[0] = (int)start_min + off[hdr->start_freq];

    if (hdr->stop_freq < 14) {
        sbr->k[2] = (int)stop_min;
        sbr_make_bands(stop_dk, (int)stop_min, 64, 13);
        sbr_sort_i16(stop_dk, 13);
        for (k = 0; k < hdr->stop_freq; k++)
            sbr->k[2] += stop_dk[k];
    } else if (hdr->stop_freq == 14) {
        sbr->k[2] = 2 * sbr->k[0];
    } else if (hdr->stop_freq == 15) {
        sbr->k[2] = 3 * sbr->k[0];
    } else {
        return -1;
    }
    if (sbr->k[2] > 64)
        sbr->k[2] = 64;

    if (sbr->out_rate <= 32000)
        max_qmf_subbands = 48;
    else if (sbr->out_rate == 44100)
        max_qmf_subbands = 35;
    else
        max_qmf_subbands = 32;

    if ((unsigned)(sbr->k[2] - sbr->k[0]) > max_qmf_subbands)
        return -1;

    if (!hdr->freq_scale) {
        int dk = hdr->alter_scale + 1;
        int k2diff;

        sbr->n_master = ((sbr->k[2] - sbr->k[0] + (dk & 2)) >> dk) << 1;
        if (!sbr_check_n_master(sbr->n_master, hdr->xover_band))
            return -1;
        for (k = 1; k <= sbr->n_master; k++)
            sbr->f_master[k] = (uint16_t)dk;
        k2diff = sbr->k[2] - sbr->k[0] - sbr->n_master * dk;
        if (k2diff < 0) {
            sbr->f_master[1]--;
            sbr->f_master[2] -= (uint16_t)(k2diff < -1);
        } else if (k2diff) {
            sbr->f_master[sbr->n_master]++;
        }
        sbr->f_master[0] = (uint16_t)sbr->k[0];
        for (k = 1; k <= sbr->n_master; k++)
            sbr->f_master[k] += sbr->f_master[k - 1];
    } else {
        int half_bands = 7 - hdr->freq_scale;
        int two_regions;
        int num_bands_0, num_bands_1 = 0;
        int vdk0_max, vdk1_min;
        int16_t vk0[49];
        int16_t vk1[49];
        float invwarp;

        if (49 * sbr->k[2] > 110 * sbr->k[0]) {
            two_regions = 1;
            sbr->k[1] = 2 * sbr->k[0];
        } else {
            two_regions = 0;
            sbr->k[1] = sbr->k[2];
        }

        num_bands_0 = (int)lrintf(half_bands *
                                  log2f((float)sbr->k[1] / (float)sbr->k[0])) *
                      2;
        if (num_bands_0 <= 0)
            return -1;

        vk0[0] = 0;
        sbr_make_bands(vk0 + 1, sbr->k[0], sbr->k[1], num_bands_0);
        sbr_sort_i16(vk0 + 1, num_bands_0);
        vdk0_max = vk0[num_bands_0];
        vk0[0] = (int16_t)sbr->k[0];
        for (k = 1; k <= num_bands_0; k++) {
            if (vk0[k] <= 0)
                return -1;
            vk0[k] = (int16_t)(vk0[k] + vk0[k - 1]);
        }

        if (two_regions) {
            invwarp = hdr->alter_scale ? 0.76923076923076923077f : 1.0f;
            num_bands_1 = (int)lrintf(half_bands * invwarp *
                                      log2f((float)sbr->k[2] /
                                            (float)sbr->k[1])) *
                          2;
            if (num_bands_1 <= 0)
                return -1;
            sbr_make_bands(vk1 + 1, sbr->k[1], sbr->k[2], num_bands_1);
            vdk1_min = 32767;
            for (k = 1; k <= num_bands_1; k++)
                if (vk1[k] < vdk1_min)
                    vdk1_min = vk1[k];
            if (vdk1_min < vdk0_max) {
                int change;
                sbr_sort_i16(vk1 + 1, num_bands_1);
                change = vdk0_max - vk1[1];
                if ((vk1[num_bands_1] - vk1[1]) >> 1 < change)
                    change = (vk1[num_bands_1] - vk1[1]) >> 1;
                vk1[1] = (int16_t)(vk1[1] + change);
                vk1[num_bands_1] = (int16_t)(vk1[num_bands_1] - change);
            }
            sbr_sort_i16(vk1 + 1, num_bands_1);
            vk1[0] = (int16_t)sbr->k[1];
            for (k = 1; k <= num_bands_1; k++) {
                if (vk1[k] <= 0)
                    return -1;
                vk1[k] = (int16_t)(vk1[k] + vk1[k - 1]);
            }
            sbr->n_master = num_bands_0 + num_bands_1;
            if (!sbr_check_n_master(sbr->n_master, hdr->xover_band))
                return -1;
            memcpy(sbr->f_master, vk0,
                   (size_t)(num_bands_0 + 1) * sizeof(sbr->f_master[0]));
            memcpy(sbr->f_master + num_bands_0 + 1, vk1 + 1,
                   (size_t)num_bands_1 * sizeof(sbr->f_master[0]));
        } else {
            sbr->n_master = num_bands_0;
            if (!sbr_check_n_master(sbr->n_master, hdr->xover_band))
                return -1;
            memcpy(sbr->f_master, vk0,
                   (size_t)(num_bands_0 + 1) * sizeof(sbr->f_master[0]));
        }
    }

    return 0;
}

static int sbr_hf_calc_npatches(aac_sbr_decoder* sbr)
{
    int i, k, last_k = -1, last_msb = -1, sb = 0;
    int msb = sbr->k[0];
    int usb = sbr->kx[1];
    int goal_sb = ((1000 << 11) + (sbr->out_rate >> 1)) / sbr->out_rate;

    sbr->num_patches = 0;
    if (goal_sb < sbr->kx[1] + sbr->m[1]) {
        for (k = 0; sbr->f_master[k] < goal_sb; k++)
            ;
    } else {
        k = sbr->n_master;
    }

    do {
        int odd = 0;
        if (k == last_k && msb == last_msb)
            return -1;
        last_k = k;
        last_msb = msb;
        for (i = k; i >= 0 &&
                        (i == k || sb > (sbr->k[0] - 1 + msb - odd)); i--) {
            sb = sbr->f_master[i];
            odd = (sb + sbr->k[0]) & 1;
        }
        if (i < 0)
            return -1;
        if (sbr->num_patches > 5)
            return -1;
        sbr->patch_num_subbands[sbr->num_patches] =
            (uint8_t)(sb - usb > 0 ? sb - usb : 0);
        sbr->patch_start_subband[sbr->num_patches] =
            (uint8_t)(sbr->k[0] - odd -
                      sbr->patch_num_subbands[sbr->num_patches]);
        if (sbr->patch_num_subbands[sbr->num_patches] > 0) {
            usb = sb;
            msb = sb;
            sbr->num_patches++;
        } else {
            msb = sbr->kx[1];
        }
        if (sbr->f_master[k] - sb < 3)
            k = sbr->n_master;
    } while (sb != sbr->kx[1] + sbr->m[1]);

    if (sbr->num_patches > 1 &&
        sbr->patch_num_subbands[sbr->num_patches - 1] < 3)
        sbr->num_patches--;
    return 0;
}

static void sbr_make_f_tablelim(aac_sbr_decoder* sbr)
{
    int k;

    if (sbr->bs_limiter_bands > 0) {
        static const float bands_warped[3] = {
            1.32715174233856803909f,
            1.18509277094158210129f,
            1.11987160404675912501f
        };
        float lim = bands_warped[sbr->bs_limiter_bands - 1];
        int16_t patch_borders[7];
        uint16_t* in = sbr->f_tablelim + 1;
        uint16_t* out = sbr->f_tablelim;

        patch_borders[0] = (int16_t)sbr->kx[1];
        for (k = 1; k <= sbr->num_patches; k++)
            patch_borders[k] =
                (int16_t)(patch_borders[k - 1] +
                          sbr->patch_num_subbands[k - 1]);

        memcpy(sbr->f_tablelim, sbr->f_tablelow,
               (size_t)(sbr->n[0] + 1) * sizeof(sbr->f_tablelow[0]));
        if (sbr->num_patches > 1)
            memcpy(sbr->f_tablelim + sbr->n[0] + 1, patch_borders + 1,
                   (size_t)(sbr->num_patches - 1) *
                       sizeof(patch_borders[0]));

        sbr_sort_u16(sbr->f_tablelim, sbr->num_patches + sbr->n[0]);

        sbr->n_lim = sbr->n[0] + sbr->num_patches - 1;
        while (out < sbr->f_tablelim + sbr->n_lim) {
            if ((float)*in >= (float)*out * lim) {
                *++out = *in++;
            } else if (*in == *out ||
                       !sbr_int16_in_table(patch_borders, sbr->num_patches,
                                           *in)) {
                in++;
                sbr->n_lim--;
            } else if (!sbr_int16_in_table(patch_borders, sbr->num_patches,
                                           *out)) {
                *out = *in++;
                sbr->n_lim--;
            } else {
                *++out = *in++;
            }
        }
    } else {
        sbr->f_tablelim[0] = sbr->f_tablelow[0];
        sbr->f_tablelim[1] = sbr->f_tablelow[sbr->n[0]];
        sbr->n_lim = 1;
    }
}

static int sbr_make_f_derived(aac_sbr_decoder* sbr)
{
    int k, temp;

    sbr->n[1] = sbr->n_master - sbr->hdr.xover_band;
    sbr->n[0] = (sbr->n[1] + 1) >> 1;
    if (sbr->n[1] <= 0)
        return -1;

    memcpy(sbr->f_tablehigh, &sbr->f_master[sbr->hdr.xover_band],
           (size_t)(sbr->n[1] + 1) * sizeof(sbr->f_master[0]));
    sbr->m[1] = sbr->f_tablehigh[sbr->n[1]] - sbr->f_tablehigh[0];
    sbr->kx[1] = sbr->f_tablehigh[0];

    if (sbr->kx[1] + sbr->m[1] > 64 || sbr->kx[1] > 32) {
        SBRDBG("sbr: f_derived kx=%d m=%d invalid (k0=%d k2=%d nm=%d xover=%d)\n",
               sbr->kx[1], sbr->m[1], sbr->k[0], sbr->k[2], sbr->n_master,
               sbr->hdr.xover_band);
        return -1;
    }

    sbr->f_tablelow[0] = sbr->f_tablehigh[0];
    temp = sbr->n[1] & 1;
    for (k = 1; k <= sbr->n[0]; k++)
        sbr->f_tablelow[k] = sbr->f_tablehigh[2 * k - temp];

    sbr->n_q = (int)lrintf(sbr->hdr.noise_bands *
                           log2f((float)sbr->k[2] / (float)sbr->kx[1]));
    if (sbr->n_q < 1)
        sbr->n_q = 1;
    if (sbr->n_q > 5) {
        SBRDBG("sbr: f_derived n_q=%d (k2=%d kx=%d nb=%d)\n", sbr->n_q,
               sbr->k[2], sbr->kx[1], sbr->hdr.noise_bands);
        return -1;
    }

    sbr->f_tablenoise[0] = sbr->f_tablelow[0];
    temp = 0;
    for (k = 1; k <= sbr->n_q; k++) {
        temp += (sbr->n[0] - temp) / (sbr->n_q + 1 - k);
        sbr->f_tablenoise[k] = sbr->f_tablelow[temp];
    }

    if (sbr_hf_calc_npatches(sbr) < 0)
        return -1;

    sbr_make_f_tablelim(sbr);

    sbr->data[0].f_indexnoise = 0;
    sbr->data[1].f_indexnoise = 0;
    sbr->data[0].f_indexsine = 0;
    sbr->data[1].f_indexsine = 0;
    return 0;
}

/* ------------------------------------------------------------------------- */
/* Bitstream parsing                                                         */
/* ------------------------------------------------------------------------- */

static int sbr_read_header(aac_sbr_decoder* sbr, aac_bits* b)
{
    aac_sbr_header old = sbr->hdr;
    int extra1, extra2;

    sbr->bs_amp_res_header = (int)aac_bits_read(b, 1);
    sbr->hdr.amp_res = sbr->bs_amp_res_header;
    sbr->hdr.start_freq = (int)aac_bits_read(b, 4);
    sbr->hdr.stop_freq = (int)aac_bits_read(b, 4);
    sbr->hdr.xover_band = (int)aac_bits_read(b, 3);
    (void)aac_bits_read(b, 2); /* bs_reserved */
    extra1 = (int)aac_bits_read(b, 1);
    extra2 = (int)aac_bits_read(b, 1);

    if (extra1) {
        sbr->hdr.freq_scale = (int)aac_bits_read(b, 2);
        sbr->hdr.alter_scale = (int)aac_bits_read(b, 1);
        sbr->hdr.noise_bands = (int)aac_bits_read(b, 2);
    } else {
        sbr->hdr.freq_scale = 2;
        sbr->hdr.alter_scale = 1;
        sbr->hdr.noise_bands = 2;
    }
    if (extra2) {
        sbr->hdr.limiter_bands = (int)aac_bits_read(b, 2);
        sbr->hdr.limiter_gains = (int)aac_bits_read(b, 2);
        sbr->hdr.interpol_freq = (int)aac_bits_read(b, 1);
        sbr->hdr.smoothing_mode = (int)aac_bits_read(b, 1);
    } else {
        sbr->hdr.limiter_bands = 2;
        sbr->hdr.limiter_gains = 2;
        sbr->hdr.interpol_freq = 1;
        sbr->hdr.smoothing_mode = 1;
    }

    if (memcmp(&old, &sbr->hdr, sizeof(old)) != 0) {
        sbr->reset = 1;
        sbr->data[0].f_indexnoise = 0;
        sbr->data[1].f_indexnoise = 0;
    }
    sbr->hdr_valid = 1;
    sbr->started = 1;
    return b->error ? -1 : 0;
}

/* USAC SBR header (ISO/IEC 23003-3): no bs_amp_res (carried in SbrInfo) and
 * bs_xover_band is fixed at 0 with two reserved bits. */
static int sbr_read_header_usac(aac_sbr_decoder* sbr, aac_bits* b)
{
    aac_sbr_header old = sbr->hdr;
    int extra1, extra2;

    sbr->hdr.start_freq = (int)aac_bits_read(b, 4);
    sbr->hdr.stop_freq = (int)aac_bits_read(b, 4);
    (void)aac_bits_read(b, 2); /* bs_reserved */
    extra1 = (int)aac_bits_read(b, 1);
    extra2 = (int)aac_bits_read(b, 1);

    if (extra1) {
        sbr->hdr.freq_scale = (int)aac_bits_read(b, 2);
        sbr->hdr.alter_scale = (int)aac_bits_read(b, 1);
        sbr->hdr.noise_bands = (int)aac_bits_read(b, 2);
    } else {
        sbr->hdr.freq_scale = 2;
        sbr->hdr.alter_scale = 1;
        sbr->hdr.noise_bands = 2;
    }
    if (extra2) {
        sbr->hdr.limiter_bands = (int)aac_bits_read(b, 2);
        sbr->hdr.limiter_gains = (int)aac_bits_read(b, 2);
        sbr->hdr.interpol_freq = (int)aac_bits_read(b, 1);
        sbr->hdr.smoothing_mode = (int)aac_bits_read(b, 1);
    } else {
        sbr->hdr.limiter_bands = 2;
        sbr->hdr.limiter_gains = 2;
        sbr->hdr.interpol_freq = 1;
        sbr->hdr.smoothing_mode = 1;
    }

    if (memcmp(&old, &sbr->hdr, sizeof(old)) != 0) {
        sbr->reset = 1;
        sbr->data[0].f_indexnoise = 0;
        sbr->data[1].f_indexnoise = 0;
    }
    sbr->hdr_valid = 1;
    sbr->started = 1;
    return b->error ? -1 : 0;
}

/* Apply the sbrConfig() defaults carried in UsacConfig. */
static void sbr_usac_default_header(aac_sbr_decoder* sbr,
                                    const aac_sbr_usac_defaults* dflt)
{
    sbr->hdr.start_freq = dflt->start_freq;
    sbr->hdr.stop_freq = dflt->stop_freq;
    sbr->hdr.freq_scale = dflt->freq_scale;
    sbr->hdr.alter_scale = dflt->alter_scale;
    sbr->hdr.noise_bands = dflt->noise_bands;
    sbr->hdr.limiter_bands = dflt->limiter_bands;
    sbr->hdr.limiter_gains = dflt->limiter_gains;
    sbr->hdr.interpol_freq = dflt->interpol_freq;
    sbr->hdr.smoothing_mode = dflt->smoothing_mode;
}

#ifdef AAC_USAC_DEBUG
#include <stdio.h>
#define SBRDBG(...) do { fprintf(stderr, __VA_ARGS__); } while (0)
#else
#define SBRDBG(...)
#endif

/* Rebuild the frequency tables after a header change. */
static int sbr_usac_update_tables(aac_sbr_decoder* sbr)
{
    if (sbr_make_f_master(sbr) < 0) {
        SBRDBG("sbr: f_master failed (start=%d stop=%d fs=%d as=%d nb=%d xover=%d)\n",
               sbr->hdr.start_freq, sbr->hdr.stop_freq, sbr->hdr.freq_scale,
               sbr->hdr.alter_scale, sbr->hdr.noise_bands, sbr->hdr.xover_band);
        sbr->initialized = 0;
        sbr->started = 0;
        return -1;
    }
    if (sbr_make_f_derived(sbr) < 0) {
        SBRDBG("sbr: f_derived failed\n");
        sbr->initialized = 0;
        sbr->started = 0;
        return -1;
    }
    sbr->initialized = 1;
    return 0;
}

static void sbr_dequant(aac_sbr_decoder* sbr);
static int sbr_read_grid(aac_sbr_decoder* sbr, aac_sbr_data* d, aac_bits* b);
static void sbr_read_dtdf(aac_sbr_decoder* sbr, aac_sbr_data* d, aac_bits* b,
                          int indep_flag);
static void sbr_read_invf(aac_sbr_decoder* sbr, aac_sbr_data* d, aac_bits* b);
static int sbr_read_envelope(aac_sbr_decoder* sbr, aac_sbr_data* d,
                             aac_bits* b, int coupling, int channel);
static int sbr_read_noise(aac_sbr_decoder* sbr, aac_sbr_data* d,
                          aac_bits* b, int coupling, int channel);
static void sbr_copy_grid(aac_sbr_data* dst, const aac_sbr_data* src);

int aac_sbr_parse_usac(aac_sbr_decoder* sbr, aac_bits* b, int sbr_ch,
                       int indep_flag, const aac_sbr_usac_defaults* dflt)
{
    aac_sbr_data* d0 = &sbr->data[0];
    aac_sbr_data* d1 = &sbr->data[1];
    int info_present = 1;
    int header_present = 1;
    int j;

    if (sbr == NULL || b == NULL || sbr_ch < 1 || sbr_ch > 2 || dflt == NULL)
        return -1;

    sbr->reset = 0;
    sbr->usac = 1;

    SBRDBG("sbr: entry bit=%zu bits=%02X %02X %02X %02X\n", b->bitpos,
           b->data[b->bitpos >> 3], b->data[(b->bitpos >> 3) + 1],
           b->data[(b->bitpos >> 3) + 2], b->data[(b->bitpos >> 3) + 3]);

    if (!indep_flag) {
        info_present = (int)aac_bits_read(b, 1);
        if (info_present)
            header_present = (int)aac_bits_read(b, 1);
        else
            header_present = 0;
    }

    if (info_present) {
        /* SbrInfo() */
        sbr->bs_amp_res_header = (int)aac_bits_read(b, 1);
        sbr->hdr.xover_band = (int)aac_bits_read(b, 4);
        (void)aac_bits_read(b, 1); /* bs_sbr_preprocessing */
        SBRDBG("sbr: info amp_res=%d xover=%d\n", sbr->bs_amp_res_header, sbr->hdr.xover_band);
    }

    if (header_present) {
        if (aac_bits_read(b, 1)) {
            aac_sbr_header old = sbr->hdr;
            int old_limiter_bands = sbr->bs_limiter_bands;
            SBRDBG("sbr: default header (start=%d stop=%d)\n",
                   dflt->start_freq, dflt->stop_freq);
            sbr_usac_default_header(sbr, dflt);
            if (memcmp(&old, &sbr->hdr, sizeof(old)) != 0) {
                sbr->reset = 1;
                sbr->data[0].f_indexnoise = 0;
                sbr->data[1].f_indexnoise = 0;
            }
            if (sbr->bs_limiter_bands != old_limiter_bands && !sbr->reset)
                sbr_make_f_tablelim(sbr);
        } else {
            if (sbr_read_header_usac(sbr, b) < 0)
                return -1;
        }
        sbr->hdr.amp_res = sbr->bs_amp_res_header;
        sbr->started = 1;
        SBRDBG("sbr: hdr start=%d stop=%d fs=%d as=%d nb=%d\n", sbr->hdr.start_freq, sbr->hdr.stop_freq, sbr->hdr.freq_scale, sbr->hdr.alter_scale, sbr->hdr.noise_bands);
        /* Failure to build the frequency tables must not break the frame
         * bitstream sync; the data is still parsed and the SBR tool simply
         * stays off (pure upsampling) like the reference decoder. */
        if (sbr_usac_update_tables(sbr) < 0) {
            SBRDBG("sbr: table update failed\n");
            sbr->started = 0;
        }
    }

    /* Save state from the previous frame. */
    d0->e_a[0] = d0->e_a[1];
    d1->e_a[0] = d1->e_a[1];
    sbr->kx[0] = sbr->kx[1];
    sbr->m[0] = sbr->m[1];

    SBRDBG("sbr: data sbr_ch=%d indep=%d bit=%zu\n", sbr_ch, indep_flag, b->bitpos);
    if (sbr_ch == 1) {
        if (sbr_read_grid(sbr, d0, b) < 0)
            return -1;
        sbr_read_dtdf(sbr, d0, b, indep_flag);
        sbr_read_invf(sbr, d0, b);
        if (sbr_read_envelope(sbr, d0, b, 0, 0) < 0)
            return -1;
        if (sbr_read_noise(sbr, d0, b, 0, 0) < 0)
            return -1;
        d0->add_harmonic_flag = (int)aac_bits_read(b, 1);
        if (d0->add_harmonic_flag)
            for (j = 0; j < sbr->n[1]; j++)
                d0->add_harmonic[j] = (int)aac_bits_read(b, 1);
        sbr->bs_coupling = 0;
    } else if (aac_bits_read(b, 1)) { /* bs_coupling */
        sbr->bs_coupling = 1;
        if (sbr_read_grid(sbr, d0, b) < 0)
            return -1;
        sbr_copy_grid(d1, d0);
        sbr_read_dtdf(sbr, d0, b, indep_flag);
        sbr_read_dtdf(sbr, d1, b, indep_flag);
        sbr_read_invf(sbr, d0, b);
        memcpy(d1->invf_mode[1], d1->invf_mode[0],
               sizeof(d1->invf_mode[0]));
        memcpy(d1->invf_mode[0], d0->invf_mode[0],
               sizeof(d1->invf_mode[0]));
        if (sbr_read_envelope(sbr, d0, b, 1, 0) < 0)
            return -1;
        if (sbr_read_noise(sbr, d0, b, 1, 0) < 0)
            return -1;
        if (sbr_read_envelope(sbr, d1, b, 1, 1) < 0)
            return -1;
        if (sbr_read_noise(sbr, d1, b, 1, 1) < 0)
            return -1;
        d0->add_harmonic_flag = (int)aac_bits_read(b, 1);
        if (d0->add_harmonic_flag)
            for (j = 0; j < sbr->n[1]; j++)
                d0->add_harmonic[j] = (int)aac_bits_read(b, 1);
        d1->add_harmonic_flag = (int)aac_bits_read(b, 1);
        if (d1->add_harmonic_flag)
            for (j = 0; j < sbr->n[1]; j++)
                d1->add_harmonic[j] = (int)aac_bits_read(b, 1);
    } else {
        sbr->bs_coupling = 0;
        if (sbr_read_grid(sbr, d0, b) < 0 ||
            sbr_read_grid(sbr, d1, b) < 0)
            return -1;
        sbr_read_dtdf(sbr, d0, b, indep_flag);
        sbr_read_dtdf(sbr, d1, b, indep_flag);
        sbr_read_invf(sbr, d0, b);
        sbr_read_invf(sbr, d1, b);
        if (sbr_read_envelope(sbr, d0, b, 0, 0) < 0)
            return -1;
        if (sbr_read_envelope(sbr, d1, b, 0, 1) < 0)
            return -1;
        if (sbr_read_noise(sbr, d0, b, 0, 0) < 0)
            return -1;
        if (sbr_read_noise(sbr, d1, b, 0, 1) < 0)
            return -1;
        d0->add_harmonic_flag = (int)aac_bits_read(b, 1);
        if (d0->add_harmonic_flag)
            for (j = 0; j < sbr->n[1]; j++)
                d0->add_harmonic[j] = (int)aac_bits_read(b, 1);
        d1->add_harmonic_flag = (int)aac_bits_read(b, 1);
        if (d1->add_harmonic_flag)
            for (j = 0; j < sbr->n[1]; j++)
                d1->add_harmonic[j] = (int)aac_bits_read(b, 1);
    }

    if (b->error) {
        SBRDBG("sbr: bitstream error at bit=%zu\n", b->bitpos);
        return -1;
    }

    if (!sbr->initialized) {
        /* No usable frequency tables (header not yet received): the frame
         * is upsampled without SBR. */
        sbr->active = 0;
        return 0;
    }

    sbr_dequant(sbr);
    sbr->active = 1;
    return 0;
}

/* ceil(log2(index + 1)) for index 0..5 */
static const int8_t sbr_ceil_log2[6] = {0, 1, 2, 2, 3, 3};

static int sbr_read_grid(aac_sbr_decoder* sbr, aac_sbr_data* d, aac_bits* b)
{
    int bs_pointer = 0;
    int abs_bord_trail = sbr->num_time_slots;
    int num_rel_lead, num_rel_trail;
    int bs_num_env_old = d->num_env;
    int bs_frame_class, bs_num_env;
    int i;

    d->freq_res[0] = d->freq_res[d->num_env];
    d->amp_res = sbr->bs_amp_res_header;
    d->t_env_num_env_old = d->t_env[d->num_env];

    bs_frame_class = (int)aac_bits_read(b, 2);
    switch (bs_frame_class) {
    case SBR_FIXFIX:
        bs_num_env = 1 << (int)aac_bits_read(b, 2);
        if (bs_num_env > (sbr->usac ? 8 : 5))
            return -1;
        d->num_env = bs_num_env;
        num_rel_lead = d->num_env - 1;
        if (d->num_env == 1)
            d->amp_res = 0;
        d->t_env[0] = 0;
        d->t_env[d->num_env] = abs_bord_trail;
        abs_bord_trail = (abs_bord_trail + (d->num_env >> 1)) / d->num_env;
        for (i = 0; i < num_rel_lead; i++)
            d->t_env[i + 1] = d->t_env[i] + abs_bord_trail;
        d->freq_res[1] = (int)aac_bits_read(b, 1);
        for (i = 1; i < d->num_env; i++)
            d->freq_res[i + 1] = d->freq_res[1];
        break;
    case SBR_FIXVAR:
        abs_bord_trail += (int)aac_bits_read(b, 2);
        num_rel_trail = (int)aac_bits_read(b, 2);
        d->num_env = num_rel_trail + 1;
        d->t_env[0] = 0;
        d->t_env[d->num_env] = abs_bord_trail;
        for (i = 0; i < num_rel_trail; i++)
            d->t_env[d->num_env - 1 - i] =
                d->t_env[d->num_env - i] - 2 * (int)aac_bits_read(b, 2) - 2;
        bs_pointer = (int)aac_bits_read(b, sbr_ceil_log2[d->num_env]);
        for (i = 0; i < d->num_env; i++)
            d->freq_res[d->num_env - i] = (int)aac_bits_read(b, 1);
        break;
    case SBR_VARFIX:
        d->t_env[0] = (int)aac_bits_read(b, 2);
        num_rel_lead = (int)aac_bits_read(b, 2);
        d->num_env = num_rel_lead + 1;
        d->t_env[d->num_env] = abs_bord_trail;
        for (i = 0; i < num_rel_lead; i++)
            d->t_env[i + 1] =
                d->t_env[i] + 2 * (int)aac_bits_read(b, 2) + 2;
        bs_pointer = (int)aac_bits_read(b, sbr_ceil_log2[d->num_env]);
        for (i = 1; i <= d->num_env; i++)
            d->freq_res[i] = (int)aac_bits_read(b, 1);
        break;
    case SBR_VARVAR:
        d->t_env[0] = (int)aac_bits_read(b, 2);
        abs_bord_trail += (int)aac_bits_read(b, 2);
        num_rel_lead = (int)aac_bits_read(b, 2);
        num_rel_trail = (int)aac_bits_read(b, 2);
        bs_num_env = num_rel_lead + num_rel_trail + 1;
        if (bs_num_env > (sbr->usac ? 8 : 5))
            return -1;
        d->num_env = bs_num_env;
        d->t_env[d->num_env] = abs_bord_trail;
        for (i = 0; i < num_rel_lead; i++)
            d->t_env[i + 1] =
                d->t_env[i] + 2 * (int)aac_bits_read(b, 2) + 2;
        for (i = 0; i < num_rel_trail; i++)
            d->t_env[d->num_env - 1 - i] =
                d->t_env[d->num_env - i] - 2 * (int)aac_bits_read(b, 2) - 2;
        bs_pointer = (int)aac_bits_read(b, sbr_ceil_log2[d->num_env]);
        for (i = 1; i <= d->num_env; i++)
            d->freq_res[i] = (int)aac_bits_read(b, 1);
        break;
    default:
        return -1;
    }
    d->frame_class = bs_frame_class;

    if (bs_pointer > d->num_env + 1)
        return -1;
    for (i = 1; i <= d->num_env; i++)
        if (d->t_env[i - 1] >= d->t_env[i])
            return -1;

    d->num_noise = (d->num_env > 1) + 1;
    d->t_q[0] = d->t_env[0];
    d->t_q[d->num_noise] = d->t_env[d->num_env];
    if (d->num_noise > 1) {
        int idx;
        if (d->frame_class == SBR_FIXFIX) {
            idx = d->num_env >> 1;
        } else if (d->frame_class & 1) {
            idx = d->num_env - sbr_max_int(bs_pointer - 1, 1);
        } else {
            if (!bs_pointer)
                idx = 1;
            else if (bs_pointer == 1)
                idx = d->num_env - 1;
            else
                idx = bs_pointer - 1;
        }
        d->t_q[1] = d->t_env[idx];
    }

    d->e_a[0] = -(d->e_a[1] != bs_num_env_old);
    d->e_a[1] = -1;
    if ((d->frame_class & 1) && bs_pointer)
        d->e_a[1] = d->num_env + 1 - bs_pointer;
    else if (d->frame_class == SBR_VARFIX && bs_pointer > 1)
        d->e_a[1] = bs_pointer - 1;


    return b->error ? -1 : 0;
}

static void sbr_copy_grid(aac_sbr_data* dst, const aac_sbr_data* src)
{
    dst->freq_res[0] = dst->freq_res[dst->num_env];
    dst->t_env_num_env_old = dst->t_env[dst->num_env];
    dst->e_a[0] = -(dst->e_a[1] != dst->num_env);

    memcpy(dst->freq_res + 1, src->freq_res + 1,
           sizeof(dst->freq_res) - sizeof(dst->freq_res[0]));
    memcpy(dst->t_env, src->t_env, sizeof(dst->t_env));
    memcpy(dst->t_q, src->t_q, sizeof(dst->t_q));
    dst->num_env = src->num_env;
    dst->amp_res = src->amp_res;
    dst->num_noise = src->num_noise;
    dst->frame_class = src->frame_class;
    dst->e_a[1] = src->e_a[1];
}

static void sbr_read_dtdf(aac_sbr_decoder* sbr, aac_sbr_data* d, aac_bits* b,
                          int indep_flag)
{
    int i;
    if (sbr->usac && indep_flag) {
        d->df_env[0] = 0;
        for (i = 1; i < d->num_env; i++)
            d->df_env[i] = (int)aac_bits_read(b, 1);
        d->df_noise[0] = 0;
        for (i = 1; i < d->num_noise; i++)
            d->df_noise[i] = (int)aac_bits_read(b, 1);
    } else {
        for (i = 0; i < d->num_env; i++)
            d->df_env[i] = (int)aac_bits_read(b, 1);
        for (i = 0; i < d->num_noise; i++)
            d->df_noise[i] = (int)aac_bits_read(b, 1);
    }
}

static void sbr_read_invf(aac_sbr_decoder* sbr, aac_sbr_data* d, aac_bits* b)
{
    int i;
    memcpy(d->invf_mode[1], d->invf_mode[0], 5 * sizeof(d->invf_mode[0][0]));
    for (i = 0; i < sbr->n_q; i++)
        d->invf_mode[0][i] = (int)aac_bits_read(b, 2);
}

static int sbr_read_envelope(aac_sbr_decoder* sbr, aac_sbr_data* d,
                             aac_bits* b, int coupling, int channel)
{
    int bits;
    int t_huff, f_huff;
    int delta = (channel == 1 && coupling) ? 2 : 1;
    int odd = sbr->n[1] & 1;
    int i, j;

    if (coupling && channel) {
        if (d->amp_res) {
            bits = 5;
            t_huff = 3;
            f_huff = 2;
        } else {
            bits = 6;
            t_huff = 2;
            f_huff = 3;
        }
    } else {
        if (d->amp_res) {
            bits = 6;
            t_huff = 4;
            f_huff = 5;
        } else {
            bits = 7;
            t_huff = 0;
            f_huff = 1;
        }
    }

    for (i = 0; i < d->num_env; i++) {
        int res = d->freq_res[i + 1];
        if (d->df_env[i]) {
            if (res == d->freq_res[i]) {
                for (j = 0; j < sbr->n[res]; j++)
                    d->env_facs_q[i + 1][j] =
                        d->env_facs_q[i][j] +
                        delta * sbr_huff_decode(b, t_huff);
            } else if (res) {
                for (j = 0; j < sbr->n[res]; j++) {
                    int k = (j + odd) >> 1;
                    d->env_facs_q[i + 1][j] =
                        d->env_facs_q[i][k] +
                        delta * sbr_huff_decode(b, t_huff);
                }
            } else {
                for (j = 0; j < sbr->n[res]; j++) {
                    int k = j ? 2 * j - odd : 0;
                    d->env_facs_q[i + 1][j] =
                        d->env_facs_q[i][k] +
                        delta * sbr_huff_decode(b, t_huff);
                }
            }
        } else {
            d->env_facs_q[i + 1][0] =
                delta * (int)aac_bits_read(b, bits);
            for (j = 1; j < sbr->n[res]; j++)
                d->env_facs_q[i + 1][j] =
                    d->env_facs_q[i + 1][j - 1] +
                    delta * sbr_huff_decode(b, f_huff);
        }
        for (j = 0; j < sbr->n[res]; j++) {
            if (d->env_facs_q[i + 1][j] < 0)
                d->env_facs_q[i + 1][j] = 0;
            if (d->env_facs_q[i + 1][j] > 127)
                d->env_facs_q[i + 1][j] = 127;
        }
        if (sbr->usac && sbr->inter_tes) {
            d->temp_shape[i] = (int)aac_bits_read(b, 1);
            if (d->temp_shape[i])
                d->temp_shape_mode[i] = (int)aac_bits_read(b, 2);
        }
        if (b->error)
            return -1;
    }
    memcpy(d->env_facs_q[0], d->env_facs_q[d->num_env],
           sizeof(d->env_facs_q[0]));
    return 0;
}

static int sbr_read_noise(aac_sbr_decoder* sbr, aac_sbr_data* d,
                          aac_bits* b, int coupling, int channel)
{
    int delta = (channel == 1 && coupling) ? 2 : 1;
    int t_huff, f_huff;
    int i, j;

    if (coupling && channel) {
        t_huff = 9;
        f_huff = 2;
    } else {
        t_huff = 8;
        f_huff = 5;
    }

    for (i = 0; i < d->num_noise; i++) {
        if (d->df_noise[i]) {
            for (j = 0; j < sbr->n_q; j++) {
                int sym = sbr_huff_decode(b, t_huff);
                d->noise_facs_q[i + 1][j] =
                    d->noise_facs_q[i][j] + delta * sym;
            }
        } else {
            d->noise_facs_q[i + 1][0] =
                delta * (int)aac_bits_read(b, 5);
            for (j = 1; j < sbr->n_q; j++)
                d->noise_facs_q[i + 1][j] =
                    d->noise_facs_q[i + 1][j - 1] +
                    delta * sbr_huff_decode(b, f_huff);
        }
        for (j = 0; j < sbr->n_q; j++) {
            if (d->noise_facs_q[i + 1][j] < 0)
                d->noise_facs_q[i + 1][j] = 0;
            if (d->noise_facs_q[i + 1][j] > 30)
                d->noise_facs_q[i + 1][j] = 30;
        }
        if (b->error)
            return -1;
    }
    memcpy(d->noise_facs_q[0], d->noise_facs_q[d->num_noise],
           sizeof(d->noise_facs_q[0]));
    return 0;
}

static void sbr_dequant(aac_sbr_decoder* sbr)
{
    int k, e;
    int ch;
    static const float exp2_tab[2] = {1.0f, 1.4142135623730951f};

    if (sbr->channels == 2 && sbr->bs_coupling) {
        int pan_offset = sbr->data[0].amp_res ? 12 : 24;
        for (e = 1; e <= sbr->data[0].num_env; e++) {
            int nenv = sbr->n[sbr->data[0].freq_res[e]];
            for (k = 0; k < nenv; k++) {
                float temp1, temp2, fac;
                if (sbr->data[0].amp_res) {
                    temp1 = sbr_exp2i(sbr->data[0].env_facs_q[e][k] + 7);
                    temp2 = sbr_exp2i(pan_offset -
                                      sbr->data[1].env_facs_q[e][k]);
                } else {
                    temp1 = sbr_exp2i((sbr->data[0].env_facs_q[e][k] >> 1) +
                                      7) *
                            exp2_tab[sbr->data[0].env_facs_q[e][k] & 1];
                    temp2 = sbr_exp2i((pan_offset -
                                       sbr->data[1].env_facs_q[e][k]) >> 1) *
                            exp2_tab[(pan_offset -
                                      sbr->data[1].env_facs_q[e][k]) & 1];
                }
                if (temp1 > 1.0e20f)
                    temp1 = 1.0f;
                fac = temp1 / (1.0f + temp2);
                sbr->data[0].env_facs[e][k] = fac;
                sbr->data[1].env_facs[e][k] = fac * temp2;
            }
        }
        for (e = 1; e <= sbr->data[0].num_noise; e++) {
            for (k = 0; k < sbr->n_q; k++) {
                float temp1 = sbr_exp2i(SBR_NOISE_OFFSET -
                                        sbr->data[0].noise_facs_q[e][k] + 1);
                float temp2 = sbr_exp2i(12 -
                                        sbr->data[1].noise_facs_q[e][k]);
                float fac = temp1 / (1.0f + temp2);
                sbr->data[0].noise_facs[e][k] = fac;
                sbr->data[1].noise_facs[e][k] = fac * temp2;
            }
        }
    } else {
        for (ch = 0; ch < (sbr->channels > 1 ? 2 : 1); ch++) {
            aac_sbr_data* d = &sbr->data[ch];
            for (e = 1; e <= d->num_env; e++) {
                int nenv = sbr->n[d->freq_res[e]];
                for (k = 0; k < nenv; k++) {
                    if (d->amp_res)
                        d->env_facs[e][k] =
                            sbr_exp2i(d->env_facs_q[e][k] + 6);
                    else
                        d->env_facs[e][k] =
                            sbr_exp2i((d->env_facs_q[e][k] >> 1) + 6) *
                            exp2_tab[d->env_facs_q[e][k] & 1];
                    if (d->env_facs[e][k] > 1.0e20f)
                        d->env_facs[e][k] = 1.0f;
                }
            }
            for (e = 1; e <= d->num_noise; e++)
                for (k = 0; k < sbr->n_q; k++)
                    d->noise_facs[e][k] =
                        sbr_exp2i(SBR_NOISE_OFFSET -
                                  d->noise_facs_q[e][k]);
        }
    }
}

static int sbr_read_single_channel(aac_sbr_decoder* sbr, aac_bits* b)
{
    aac_sbr_data* d = &sbr->data[0];
    int j;

    if (aac_bits_read(b, 1)) /* bs_data_extra */
        (void)aac_bits_read(b, 4);

    if (sbr_read_grid(sbr, d, b) < 0)
        return -1;
    sbr_read_dtdf(sbr, d, b, 0);
    sbr_read_invf(sbr, d, b);
    if (sbr_read_envelope(sbr, d, b, 0, 0) < 0)
        return -1;
    if (sbr_read_noise(sbr, d, b, 0, 0) < 0)
        return -1;

    d->add_harmonic_flag = (int)aac_bits_read(b, 1);
    if (d->add_harmonic_flag)
        for (j = 0; j < sbr->n[1]; j++)
            d->add_harmonic[j] = (int)aac_bits_read(b, 1);
    return b->error ? -1 : 0;
}

static int sbr_read_channel_pair(aac_sbr_decoder* sbr, aac_bits* b)
{
    aac_sbr_data* d0 = &sbr->data[0];
    aac_sbr_data* d1 = &sbr->data[1];
    int j;

    if (aac_bits_read(b, 1)) /* bs_data_extra */
        (void)aac_bits_read(b, 8);

    sbr->bs_coupling = (int)aac_bits_read(b, 1);
    if (sbr->bs_coupling) {
        if (sbr_read_grid(sbr, d0, b) < 0)
            return -1;
        sbr_copy_grid(d1, d0);
        sbr_read_dtdf(sbr, d0, b, 0);
        sbr_read_dtdf(sbr, d1, b, 0);
        sbr_read_invf(sbr, d0, b);
        memcpy(d1->invf_mode[1], d1->invf_mode[0],
               sizeof(d1->invf_mode[0]));
        memcpy(d1->invf_mode[0], d0->invf_mode[0],
               sizeof(d1->invf_mode[0]));
        if (sbr_read_envelope(sbr, d0, b, 1, 0) < 0)
            return -1;
        if (sbr_read_noise(sbr, d0, b, 1, 0) < 0)
            return -1;
        if (sbr_read_envelope(sbr, d1, b, 1, 1) < 0)
            return -1;
        if (sbr_read_noise(sbr, d1, b, 1, 1) < 0)
            return -1;
    } else {
        if (sbr_read_grid(sbr, d0, b) < 0 ||
            sbr_read_grid(sbr, d1, b) < 0)
            return -1;
        sbr_read_dtdf(sbr, d0, b, 0);
        sbr_read_dtdf(sbr, d1, b, 0);
        sbr_read_invf(sbr, d0, b);
        sbr_read_invf(sbr, d1, b);
        if (sbr_read_envelope(sbr, d0, b, 0, 0) < 0)
            return -1;
        if (sbr_read_envelope(sbr, d1, b, 0, 1) < 0)
            return -1;
        if (sbr_read_noise(sbr, d0, b, 0, 0) < 0)
            return -1;
        if (sbr_read_noise(sbr, d1, b, 0, 1) < 0)
            return -1;
    }

    d0->add_harmonic_flag = (int)aac_bits_read(b, 1);
    if (d0->add_harmonic_flag)
        for (j = 0; j < sbr->n[1]; j++)
            d0->add_harmonic[j] = (int)aac_bits_read(b, 1);
    d1->add_harmonic_flag = (int)aac_bits_read(b, 1);
    if (d1->add_harmonic_flag)
        for (j = 0; j < sbr->n[1]; j++)
            d1->add_harmonic[j] = (int)aac_bits_read(b, 1);
    return b->error ? -1 : 0;
}

/* Parse the bs_extended_data block.  Extension id 2 is Parametric Stereo. */
static int sbr_skip_extended_data(aac_sbr_decoder* sbr, aac_bits* b)
{
    int flag = (int)aac_bits_read(b, 1);
    if (!flag)
        return 0;
    {
        int num_bits_left = (int)aac_bits_read(b, 4);
        if (num_bits_left == 15)
            num_bits_left += (int)aac_bits_read(b, 8);
        num_bits_left <<= 3;
        if (num_bits_left > 7) {
            int id;
            num_bits_left -= 2;
            id = (int)aac_bits_read(b, 2);
            if (id == 2) {
                (void)aac_ps_parse(&sbr->ps, b, num_bits_left);
                sbr->ps_present = 1;
                num_bits_left = 0;
            }
            if (num_bits_left > 0)
                aac_bits_skip(b, (size_t)num_bits_left);
        }
    }
    return b->error ? -1 : 0;
}

/* ------------------------------------------------------------------------- */
/* QMF analysis / LF generation                                              */
/* ------------------------------------------------------------------------- */

static void sbr_qmf_analysis(aac_sbr_decoder* sbr, aac_sbr_data* d,
                             const float* input, int samples)
{
    float buf[320 + 1024];
    int slots = samples / 32;
    int l, k, n;

    if (slots > 32)
        slots = 32;

    memcpy(buf, d->ana_hist, 320 * sizeof(float));
    memcpy(buf + 320, input, (size_t)samples * sizeof(float));

    memset(d->W[d->Ypos], 0, sizeof(d->W[d->Ypos]));
    for (l = 0; l < slots; l++) {
        const float* base = buf + 320 + l * 32;
        for (k = 0; k < 32; k++) {
            float re = 0.0f, im = 0.0f;
            const float* ac = ana_cos_tab[k];
            const float* as = ana_sin_tab[k];
            for (n = 0; n < 320; n++) {
                float v = base[n - 319];
                re += v * ac[n];
                im -= v * as[n];
            }
            d->W[d->Ypos][l][k][0] = re * (1.0f / 32.0f);
            d->W[d->Ypos][l][k][1] = im * (1.0f / 32.0f);
        }
    }

    memcpy(d->ana_hist, buf + samples, 320 * sizeof(float));
    (void)sbr;
}

static void sbr_lf_gen(aac_sbr_decoder* sbr, aac_sbr_data* d)
{
    int k, i;
    int cur = d->Ypos;
    int prev = 1 - cur;
    int iframe = sbr->num_time_slots * 2;

    memset(d->X_low, 0, sizeof(d->X_low));
    for (k = 0; k < sbr->kx[1] && k < 32; k++)
        for (i = 8; i < iframe + 8; i++)
            d->X_low[k][i][0] = d->W[cur][i - 8][k][0],
            d->X_low[k][i][1] = d->W[cur][i - 8][k][1];
    for (k = 0; k < sbr->kx[0] && k < 32; k++)
        for (i = 0; i < 8; i++)
            d->X_low[k][i][0] = d->W[prev][iframe - 8 + i][k][0],
            d->X_low[k][i][1] = d->W[prev][iframe - 8 + i][k][1];
}

/* ------------------------------------------------------------------------- */
/* Inverse filtering, chirp and HF generation                                */
/* ------------------------------------------------------------------------- */

static void sbr_autocorrelate(const float x[40][2], float phi[3][2][2])
{
    float real_sum2 = x[0][0] * x[2][0] + x[0][1] * x[2][1];
    float imag_sum2 = x[0][0] * x[2][1] - x[0][1] * x[2][0];
    float real_sum1 = 0.0f, imag_sum1 = 0.0f, real_sum0 = 0.0f;
    int i;

    for (i = 1; i < 38; i++) {
        real_sum0 += x[i][0] * x[i][0] + x[i][1] * x[i][1];
        real_sum1 += x[i][0] * x[i + 1][0] + x[i][1] * x[i + 1][1];
        imag_sum1 += x[i][0] * x[i + 1][1] - x[i][1] * x[i + 1][0];
        real_sum2 += x[i][0] * x[i + 2][0] + x[i][1] * x[i + 2][1];
        imag_sum2 += x[i][0] * x[i + 2][1] - x[i][1] * x[i + 2][0];
    }
    phi[0][1][0] = real_sum2;
    phi[0][1][1] = imag_sum2;
    phi[2][1][0] = real_sum0 + x[0][0] * x[0][0] + x[0][1] * x[0][1];
    phi[1][0][0] = real_sum0 + x[38][0] * x[38][0] + x[38][1] * x[38][1];
    phi[1][1][0] = real_sum1 + x[0][0] * x[1][0] + x[0][1] * x[1][1];
    phi[1][1][1] = imag_sum1 + x[0][0] * x[1][1] - x[0][1] * x[1][0];
    phi[0][0][0] = real_sum1 + x[38][0] * x[39][0] + x[38][1] * x[39][1];
    phi[0][0][1] = imag_sum1 + x[38][0] * x[39][1] - x[38][1] * x[39][0];
    phi[2][0][0] = phi[2][0][1] = 0.0f;
}

static void sbr_hf_inverse_filter(aac_sbr_decoder* sbr, aac_sbr_data* d)
{
    int k;
    for (k = 0; k < sbr->k[0] && k < 32; k++) {
        float phi[3][2][2];
        float dk;
        sbr_autocorrelate(d->X_low[k], phi);

        dk = phi[2][1][0] * phi[1][0][0] -
             (phi[1][1][0] * phi[1][1][0] +
              phi[1][1][1] * phi[1][1][1]) / 1.000001f;
        if (dk == 0.0f) {
            sbr->alpha1[k][0] = 0.0f;
            sbr->alpha1[k][1] = 0.0f;
        } else {
            float tr = phi[0][0][0] * phi[1][1][0] -
                       phi[0][0][1] * phi[1][1][1] -
                       phi[0][1][0] * phi[1][0][0];
            float ti = phi[0][0][0] * phi[1][1][1] +
                       phi[0][0][1] * phi[1][1][0] -
                       phi[0][1][1] * phi[1][0][0];
            sbr->alpha1[k][0] = tr / dk;
            sbr->alpha1[k][1] = ti / dk;
        }
        if (phi[1][0][0] == 0.0f) {
            sbr->alpha0[k][0] = 0.0f;
            sbr->alpha0[k][1] = 0.0f;
        } else {
            float tr = phi[0][0][0] +
                       sbr->alpha1[k][0] * phi[1][1][0] +
                       sbr->alpha1[k][1] * phi[1][1][1];
            float ti = phi[0][0][1] +
                       sbr->alpha1[k][1] * phi[1][1][0] -
                       sbr->alpha1[k][0] * phi[1][1][1];
            sbr->alpha0[k][0] = -tr / phi[1][0][0];
            sbr->alpha0[k][1] = -ti / phi[1][0][0];
        }
        if (sbr->alpha1[k][0] * sbr->alpha1[k][0] +
                sbr->alpha1[k][1] * sbr->alpha1[k][1] >= 16.0f ||
            sbr->alpha0[k][0] * sbr->alpha0[k][0] +
                sbr->alpha0[k][1] * sbr->alpha0[k][1] >= 16.0f) {
            sbr->alpha1[k][0] = sbr->alpha1[k][1] = 0.0f;
            sbr->alpha0[k][0] = sbr->alpha0[k][1] = 0.0f;
        }
    }
}

static void sbr_chirp(aac_sbr_decoder* sbr, aac_sbr_data* d)
{
    static const float bw_tab[4] = {0.0f, 0.75f, 0.9f, 0.98f};
    int i;
    for (i = 0; i < sbr->n_q; i++) {
        float new_bw;
        if (d->invf_mode[0][i] + d->invf_mode[1][i] == 1)
            new_bw = 0.6f;
        else
            new_bw = bw_tab[d->invf_mode[0][i]];
        if (new_bw < d->bw_array[i])
            new_bw = 0.75f * new_bw + 0.25f * d->bw_array[i];
        else
            new_bw = 0.90625f * new_bw + 0.09375f * d->bw_array[i];
        d->bw_array[i] = new_bw < 0.015625f ? 0.0f : new_bw;
    }
}

static void sbr_hf_gen(aac_sbr_decoder* sbr, aac_sbr_data* d)
{
    int j, x;
    int g = 0;
    int k = sbr->kx[1];
    int start = 2 * d->t_env[0] + SBR_ENV_OFFSET;
    int end = 2 * d->t_env[d->num_env] + SBR_ENV_OFFSET;

    memset(d->X_high, 0, sizeof(d->X_high));
    for (j = 0; j < sbr->num_patches; j++) {
        for (x = 0; x < sbr->patch_num_subbands[j] && k < 64; x++, k++) {
            int p = sbr->patch_start_subband[j] + x;
            float a0r, a0i, a1r, a1i, bw, bw2;
            int t;

            while (g <= sbr->n_q && k >= sbr->f_tablenoise[g])
                g++;
            g--;
            if (g < 0)
                g = 0;
            if (p < 0 || p >= 32)
                continue;
            a0r = sbr->alpha0[p][0];
            a0i = sbr->alpha0[p][1];
            a1r = sbr->alpha1[p][0];
            a1i = sbr->alpha1[p][1];
            bw = d->bw_array[g];
            bw2 = bw * bw;

            for (t = start; t < end; t++) {
                float l0r = d->X_low[p][t][0];
                float l0i = d->X_low[p][t][1];
                float l1r = d->X_low[p][t - 1][0];
                float l1i = d->X_low[p][t - 1][1];
                float l2r = d->X_low[p][t - 2][0];
                float l2i = d->X_low[p][t - 2][1];

                d->X_high[k][t][0] =
                    l0r + (a0r * l1r - a0i * l1i) * bw +
                    (a1r * l2r - a1i * l2i) * bw2;
                d->X_high[k][t][1] =
                    l0i + (a0r * l1i + a0i * l1r) * bw +
                    (a1r * l2i + a1i * l2r) * bw2;
            }
        }
    }
}

/* ------------------------------------------------------------------------- */
/* Envelope mapping, estimation, gain calculation and assembly               */
/* ------------------------------------------------------------------------- */

static void sbr_mapping(aac_sbr_decoder* sbr, aac_sbr_data* d)
{
    int e, i, m, k;

    memset(&d->s_indexmapped[1][0], 0,
           7 * 48 * sizeof(d->s_indexmapped[0][0]));

    for (e = 0; e < d->num_env; e++) {
        int res = d->freq_res[e + 1];
        int ilim = sbr->n[res];
        const uint16_t* table = res ? sbr->f_tablehigh : sbr->f_tablelow;

        for (i = 0; i < ilim; i++)
            for (m = table[i]; m < table[i + 1]; m++)
                sbr->e_origmapped[e][m - sbr->kx[1]] =
                    d->env_facs[e + 1][i];

        k = (d->num_noise > 1) && (d->t_env[e] >= d->t_q[1]);
        for (i = 0; i < sbr->n_q; i++)
            for (m = sbr->f_tablenoise[i]; m < sbr->f_tablenoise[i + 1]; m++)
                sbr->q_mapped[e][m - sbr->kx[1]] =
                    d->noise_facs[k + 1][i];

        for (i = 0; i < sbr->n[1]; i++) {
            if (d->add_harmonic_flag) {
                int mid = (sbr->f_tablehigh[i] +
                           sbr->f_tablehigh[i + 1]) >> 1;
                d->s_indexmapped[e + 1][mid - sbr->kx[1]] =
                    d->add_harmonic[i] *
                    (e >= d->e_a[1] ||
                     d->s_indexmapped[0][mid - sbr->kx[1]] == 1);
            }
        }

        for (i = 0; i < ilim; i++) {
            int present = 0;
            for (m = table[i]; m < table[i + 1]; m++) {
                if (d->s_indexmapped[e + 1][m - sbr->kx[1]]) {
                    present = 1;
                    break;
                }
            }
            for (m = table[i]; m < table[i + 1]; m++)
                sbr->s_mapped[e][m - sbr->kx[1]] = (float)present;
        }
    }

    memcpy(d->s_indexmapped[0], d->s_indexmapped[d->num_env],
           sizeof(d->s_indexmapped[0]));
}

static void sbr_env_estimate(aac_sbr_decoder* sbr, aac_sbr_data* d)
{
    int e, m, k, p;
    int kx1 = sbr->kx[1];
    int m_max = sbr->m[1];

    if (sbr->bs_interpol_freq) {
        for (e = 0; e < d->num_env; e++) {
            float recip = 0.5f /
                          (float)(d->t_env[e + 1] - d->t_env[e]);
            int ilb = d->t_env[e] * 2 + SBR_ENV_OFFSET;
            int iub = d->t_env[e + 1] * 2 + SBR_ENV_OFFSET;
            if (ilb >= 40)
                return;
            for (m = 0; m < m_max; m++) {
                float sum = 0.0f;
                for (k = ilb; k < iub && k < 40; k++)
                    sum += d->X_high[m + kx1][k][0] *
                               d->X_high[m + kx1][k][0] +
                           d->X_high[m + kx1][k][1] *
                               d->X_high[m + kx1][k][1];
                sbr->e_curr[e][m] = sum * recip;
            }
        }
    } else {
        for (e = 0; e < d->num_env; e++) {
            int env_size = 2 * (d->t_env[e + 1] - d->t_env[e]);
            int ilb = d->t_env[e] * 2 + SBR_ENV_OFFSET;
            int iub = d->t_env[e + 1] * 2 + SBR_ENV_OFFSET;
            int res = d->freq_res[e + 1];
            const uint16_t* table = res ? sbr->f_tablehigh :
                                          sbr->f_tablelow;
            if (ilb >= 40)
                return;
            for (p = 0; p < sbr->n[res]; p++) {
                float sum = 0.0f;
                int den = env_size * (table[p + 1] - table[p]);
                if (den <= 0)
                    den = 1;
                for (k = table[p]; k < table[p + 1] && k < kx1 + m_max; k++) {
                    int i;
                    for (i = ilb; i < iub && i < 40; i++)
                        sum += d->X_high[k][i][0] * d->X_high[k][i][0] +
                               d->X_high[k][i][1] * d->X_high[k][i][1];
                }
                sum /= (float)den;
                for (k = table[p]; k < table[p + 1] && k < kx1 + m_max; k++)
                    sbr->e_curr[e][k - kx1] = sum;
            }
        }
    }
}

static void sbr_gain_calc(aac_sbr_decoder* sbr, aac_sbr_data* d)
{
    static const float limgain[4] = {
        0.70795f, 1.0f, 1.41254f, 10000000000.0f
    };
    int e, k, m;

    for (e = 0; e < d->num_env; e++) {
        int delta = !((e == d->e_a[1]) || (e == d->e_a[0]));
        for (k = 0; k < sbr->n_lim; k++) {
            int lo = sbr->f_tablelim[k] - sbr->kx[1];
            int hi = sbr->f_tablelim[k + 1] - sbr->kx[1];
            float gain_boost, gain_max;
            float sum[2] = {0.0f, 0.0f};

            if (lo < 0)
                lo = 0;
            if (hi > sbr->m[1])
                hi = sbr->m[1];
            for (m = lo; m < hi; m++) {
                float temp = sbr->e_origmapped[e][m] /
                             (1.0f + sbr->q_mapped[e][m]);
                sbr->q_m[e][m] = sqrtf(temp * sbr->q_mapped[e][m]);
                sbr->s_m[e][m] = sqrtf(temp *
                                       (float)d->s_indexmapped[e + 1][m]);
                if (!sbr->s_mapped[e][m]) {
                    sbr->gain[e][m] =
                        sqrtf(sbr->e_origmapped[e][m] /
                              ((1.0f + sbr->e_curr[e][m]) *
                               (1.0f + sbr->q_mapped[e][m] * delta)));
                } else {
                    sbr->gain[e][m] =
                        sqrtf(sbr->e_origmapped[e][m] *
                              sbr->q_mapped[e][m] /
                              ((1.0f + sbr->e_curr[e][m]) *
                               (1.0f + sbr->q_mapped[e][m])));
                }
                sbr->gain[e][m] += FLT_MIN;
            }
            for (m = lo; m < hi; m++) {
                sum[0] += sbr->e_origmapped[e][m];
                sum[1] += sbr->e_curr[e][m];
            }
            gain_max = limgain[sbr->bs_limiter_gains] *
                       sqrtf((FLT_EPSILON + sum[0]) /
                             (FLT_EPSILON + sum[1]));
            if (gain_max > 100000.0f)
                gain_max = 100000.0f;
            for (m = lo; m < hi; m++) {
                float q_m_max = sbr->q_m[e][m] * gain_max /
                                sbr->gain[e][m];
                if (sbr->q_m[e][m] > q_m_max)
                    sbr->q_m[e][m] = q_m_max;
                if (sbr->gain[e][m] > gain_max)
                    sbr->gain[e][m] = gain_max;
            }
            sum[0] = sum[1] = 0.0f;
            for (m = lo; m < hi; m++) {
                sum[0] += sbr->e_origmapped[e][m];
                sum[1] += sbr->e_curr[e][m] * sbr->gain[e][m] *
                              sbr->gain[e][m] +
                          sbr->s_m[e][m] * sbr->s_m[e][m] +
                          ((delta && !sbr->s_m[e][m]) ?
                               sbr->q_m[e][m] * sbr->q_m[e][m] : 0.0f);
            }
            gain_boost = sqrtf((FLT_EPSILON + sum[0]) /
                               (FLT_EPSILON + sum[1]));
            if (gain_boost > 1.584893192f)
                gain_boost = 1.584893192f;
            for (m = lo; m < hi; m++) {
                sbr->gain[e][m] *= gain_boost;
                sbr->q_m[e][m] *= gain_boost;
                sbr->s_m[e][m] *= gain_boost;
            }
        }
    }
}

static void sbr_hf_assemble(aac_sbr_decoder* sbr, aac_sbr_data* d,
                            float Y1[38][64][2])
{
    static const float h_smooth[5] = {
        0.33333333333333f,
        0.30150283239582f,
        0.21816949906249f,
        0.11516383427084f,
        0.03183050093751f
    };
    const int h_SL = 4 * !sbr->bs_smoothing_mode;
    const int kx = sbr->kx[1];
    const int m_max = sbr->m[1];
    float(*g_temp)[48] = d->g_temp;
    float(*q_temp)[48] = d->q_temp;
    int indexnoise = d->f_indexnoise;
    int indexsine = d->f_indexsine;
    int e, i, j, m;

    if (sbr->reset) {
        for (i = 0; i < h_SL; i++) {
            memcpy(g_temp[i + 2 * d->t_env[0]], sbr->gain[0],
                   (size_t)m_max * sizeof(sbr->gain[0][0]));
            memcpy(q_temp[i + 2 * d->t_env[0]], sbr->q_m[0],
                   (size_t)m_max * sizeof(sbr->q_m[0][0]));
        }
    } else if (h_SL) {
        for (i = 0; i < 4; i++) {
            memcpy(g_temp[i + 2 * d->t_env[0]],
                   g_temp[i + 2 * d->t_env_num_env_old],
                   sizeof(g_temp[0]));
            memcpy(q_temp[i + 2 * d->t_env[0]],
                   q_temp[i + 2 * d->t_env_num_env_old],
                   sizeof(q_temp[0]));
        }
    }

    for (e = 0; e < d->num_env; e++) {
        for (i = 2 * d->t_env[e]; i < 2 * d->t_env[e + 1]; i++) {
            memcpy(g_temp[h_SL + i], sbr->gain[e],
                   (size_t)m_max * sizeof(sbr->gain[0][0]));
            memcpy(q_temp[h_SL + i], sbr->q_m[e],
                   (size_t)m_max * sizeof(sbr->q_m[0][0]));
        }
    }

    for (e = 0; e < d->num_env; e++) {
        for (i = 2 * d->t_env[e]; i < 2 * d->t_env[e + 1] && i < 32; i++) {
            float g_filt[48];
            float q_filt[48];

            if (h_SL && e != d->e_a[0] && e != d->e_a[1]) {
                for (m = 0; m < m_max; m++) {
                    int idx1 = i + h_SL;
                    g_filt[m] = 0.0f;
                    q_filt[m] = 0.0f;
                    for (j = 0; j <= h_SL; j++) {
                        g_filt[m] += g_temp[idx1 - j][m] * h_smooth[j];
                        q_filt[m] += q_temp[idx1 - j][m] * h_smooth[j];
                    }
                }
            } else {
                memcpy(g_filt, g_temp[i + h_SL],
                       (size_t)m_max * sizeof(g_filt[0]));
                memcpy(q_filt, q_temp[i],
                       (size_t)m_max * sizeof(q_filt[0]));
            }

            for (m = 0; m < m_max; m++) {
                int bk = kx + m;
                if (bk < 64) {
                    Y1[i][bk][0] = d->X_high[bk][i + SBR_ENV_OFFSET][0] *
                                   g_filt[m];
                    Y1[i][bk][1] = d->X_high[bk][i + SBR_ENV_OFFSET][1] *
                                   g_filt[m];
                }
            }

            if (e != d->e_a[0] && e != d->e_a[1]) {
                float phi0, phi1;
                switch (indexsine & 3) {
                case 0: phi0 = 1.0f; phi1 = 0.0f; break;
                case 1: phi0 = 0.0f;
                        phi1 = (kx & 1) ? -1.0f : 1.0f; break;
                case 2: phi0 = -1.0f; phi1 = 0.0f; break;
                default: phi0 = 0.0f;
                         phi1 = (kx & 1) ? 1.0f : -1.0f; break;
                }
                for (m = 0; m < m_max; m++) {
                    int bk = kx + m;
                    indexnoise = (indexnoise + 1) & 0x1ff;
                    if (bk >= 64)
                        continue;
                    if (sbr->s_m[e][m] != 0.0f) {
                        Y1[i][bk][0] += sbr->s_m[e][m] * phi0;
                        Y1[i][bk][1] += sbr->s_m[e][m] * phi1;
                    } else {
                        Y1[i][bk][0] += q_filt[m] *
                                        aac_sbr_noise[indexnoise][0];
                        Y1[i][bk][1] += q_filt[m] *
                                        aac_sbr_noise[indexnoise][1];
                    }
                    phi1 = -phi1;
                }
            } else {
                int idx = indexsine & 1;
                int A = 1 - ((indexsine + (kx & 1)) & 2);
                int B = idx ? -A : A;
                for (m = 0; m + 1 < m_max; m += 2) {
                    int bk = kx + m;
                    if (bk < 64)
                        Y1[i][bk][idx] += sbr->s_m[e][m] * (float)A;
                    if (bk + 1 < 64)
                        Y1[i][bk + 1][idx] +=
                            sbr->s_m[e][m + 1] * (float)B;
                }
                if ((m_max & 1) && kx + m < 64)
                    Y1[i][kx + m][idx] += sbr->s_m[e][m] * (float)A;
                indexnoise = (indexnoise + m_max) & 0x1ff;
            }
            indexsine = (indexsine + 1) & 3;
        }
    }
    d->f_indexnoise = indexnoise;
    d->f_indexsine = indexsine;
}

/* ------------------------------------------------------------------------- */
/* 64-band QMF synthesis                                                     */
/* ------------------------------------------------------------------------- */

static void sbr_qmf_synthesis(aac_sbr_decoder* sbr, aac_sbr_data* d,
                              const float Xre[38][64],
                              const float Xim[38][64],
                              float* out, int out_samples)
{
    int slots = out_samples / 64;
    int i, n, k;

    if (slots > 32)
        slots = 32;

    for (i = 0; i < slots; i++) {
        for (n = 0; n < 640; n++) {
            float acc = 0.0f;
            const float* cs = syn_cos_tab[n];
            const float* sn = syn_sin_tab[n];
            for (k = 0; k < 64; k++)
                acc += Xre[i][k] * cs[k] - Xim[i][k] * sn[k];
            d->synth_buf[i * 64 + n] += acc;
        }
    }

    memcpy(out, d->synth_buf, (size_t)out_samples * sizeof(float));
    memmove(d->synth_buf, d->synth_buf + out_samples,
            640 * sizeof(d->synth_buf[0]));
    memset(d->synth_buf + 640, 0,
           (sizeof(d->synth_buf) / sizeof(d->synth_buf[0]) - 640) *
               sizeof(d->synth_buf[0]));
    (void)sbr;
}

/* ------------------------------------------------------------------------- */
/* Public API                                                                */
/* ------------------------------------------------------------------------- */

static void sbr_data_reset(aac_sbr_data* d)
{
    memset(d, 0, sizeof(*d));
    d->e_a[1] = -1;
    d->f_indexnoise = 0;
    d->f_indexsine = 0;
}

static void sbr_channel_reset(aac_sbr_data* d)
{
    sbr_data_reset(d);
}

int aac_sbr_open(aac_sbr_decoder* sbr, int channels, int core_rate,
                 int out_rate, int frame_size)
{
    int i;

    if (sbr == NULL || channels < 1 || channels > AAC_SBR_MAX_CHANNELS)
        return -1;
    if (sbr_huffman_init() < 0 || sbr_qmf_tables_init() < 0)
        return -1;

    memset(sbr, 0, sizeof(*sbr));
    sbr->channels = channels;
    sbr->core_rate = core_rate;
    sbr->out_rate = out_rate > 0 ? out_rate : 2 * core_rate;
    sbr->frame_size = frame_size;
    sbr->num_time_slots = (frame_size == 960) ? 15 : 16;
    sbr->kx[0] = 32;
    sbr->kx[1] = 32;
    sbr->m[0] = 0;
    sbr->m[1] = 0;
    sbr->n[0] = 0;
    sbr->n[1] = 0;
    sbr->n_q = 1;
    sbr->num_patches = 0;
    sbr->bs_limiter_bands = 2;
    sbr->bs_limiter_gains = 2;
    sbr->bs_interpol_freq = 1;
    sbr->bs_smoothing_mode = 1;
    sbr->bs_amp_res_header = 1;
    for (i = 0; i < channels; i++)
        sbr_channel_reset(&sbr->data[i]);
    aac_ps_reset(&sbr->ps);
    return 0;
}

void aac_sbr_close(aac_sbr_decoder* sbr)
{
    (void)sbr;
}

void aac_sbr_reset(aac_sbr_decoder* sbr)
{
    int ch, core_rate, out_rate, frame_size;
    if (sbr == NULL)
        return;
    ch = sbr->channels;
    core_rate = sbr->core_rate;
    out_rate = sbr->out_rate;
    frame_size = sbr->frame_size;
    (void)aac_sbr_open(sbr, ch, core_rate, out_rate, frame_size);
}

int aac_sbr_parse(aac_sbr_decoder* sbr, aac_bits* b, int len_bits,
                  int element_type, int crc)
{
    int start_bit;
    int header_flag;
    int i;
    int end_bit;

    if (sbr == NULL || b == NULL || sbr->channels <= 0)
        return -1;
    if (element_type != 0 && element_type != 1)
        return -1;

    start_bit = (int)b->bitpos;
    end_bit = start_bit + len_bits;

    sbr->active = 0;
    sbr->reset = 0;
    sbr->ps_present = 0;
    for (i = 0; i < sbr->channels; i++)

    if (crc)
        (void)aac_bits_read(b, 10);

    for (i = 0; i < sbr->channels; i++) {
        sbr->data[i].e_a[0] = sbr->data[i].e_a[1];
        sbr->kx[0] = sbr->kx[1];
        sbr->m[0] = sbr->m[1];
    }

    header_flag = (int)aac_bits_read(b, 1);
    if (header_flag) {
        if (sbr_read_header(sbr, b) < 0)
            return -1;
        if (sbr_make_f_master(sbr) < 0 ||
            sbr_make_f_derived(sbr) < 0) {
            sbr->initialized = 0;
            sbr->started = 0;
            return -1;
        }
        sbr->initialized = 1;
    }

    if (!sbr->started || !sbr->initialized) {
        aac_bits_skip(b, (size_t)(end_bit - (int)b->bitpos));
        return 0;
    }

    if (element_type == 0) {
        if (sbr_read_single_channel(sbr, b) < 0)
            return -1;
    } else {
        if (sbr_read_channel_pair(sbr, b) < 0)
            return -1;
    }

    if (sbr_skip_extended_data(sbr, b) < 0)
        return -1;

    sbr_dequant(sbr);
    for (i = 0; i < (element_type == 0 ? 1 : 2); i++)
        sbr->active = 1;

    if ((int)b->bitpos < end_bit)
        aac_bits_skip(b, (size_t)(end_bit - (int)b->bitpos));
    return b->error ? -1 : 0;
}

/* Build the 64-band spectrum for synthesis from X_low and Y. */
static void sbr_build_x(aac_sbr_decoder* sbr, aac_sbr_data* d,
                        float Xre[38][64], float Xim[38][64])
{
    int k, i;
    int kx1 = sbr->kx[1];
    int m1 = sbr->m[1];
    int iframe = sbr->num_time_slots * 2;

    memset(Xre, 0, 38 * 64 * sizeof(float));
    memset(Xim, 0, 38 * 64 * sizeof(float));

    for (k = 0; k < kx1 && k < 64; k++) {
        for (i = 0; i < 38; i++) {
            int li = i + SBR_ENV_OFFSET;
            if (li < 40) {
                Xre[i][k] = d->X_low[k][li][0];
                Xim[i][k] = d->X_low[k][li][1];
            }
        }
    }
    for (k = kx1; k < kx1 + m1 && k < 64; k++) {
        for (i = 0; i < iframe && i < 38; i++) {
            Xre[i][k] = d->Y[d->Ypos][i][k][0] * SBR_HF_GAIN;
            Xim[i][k] = d->Y[d->Ypos][i][k][1] * SBR_HF_GAIN;
        }
    }
}

int aac_sbr_apply(aac_sbr_decoder* sbr, const float* core_pcm[2],
                  int core_samples, float* out_pcm[2], int out_samples)
{
    int ch;
    float Xre[38][64];
    float Xim[38][64];

    if (sbr == NULL || core_pcm == NULL || out_pcm == NULL ||
        core_samples <= 0 || out_samples < core_samples * 2)
        return -1;

    for (ch = 0; ch < sbr->channels; ch++) {
        aac_sbr_data* d = &sbr->data[ch];
        int pos;

        sbr_qmf_analysis(sbr, d, core_pcm[ch], core_samples);
        sbr_lf_gen(sbr, d);
        pos = d->Ypos;

        if (sbr->active) {
            sbr_hf_inverse_filter(sbr, d);
            sbr_chirp(sbr, d);
            sbr_hf_gen(sbr, d);
            sbr_mapping(sbr, d);
            sbr_env_estimate(sbr, d);
            sbr_gain_calc(sbr, d);
            sbr_hf_assemble(sbr, d, d->Y[pos]);
        } else {
            memset(d->Y[pos], 0, sizeof(d->Y[pos]));
        }

        sbr_build_x(sbr, d, Xre, Xim);

        if (ch == 0 && sbr->channels == 1 && out_pcm[1] != NULL) {
            if (sbr->ps.start) {
                float Rre[38][64];
                float Rim[38][64];
                int i, j;
                aac_ps_apply(&sbr->ps, Xre, Xim, Rre, Rim,
                             sbr->kx[1] + sbr->m[1]);
                for (i = 0; i < 38; i++) {
                    for (j = 0; j < 64; j++) {
                        Rre[i][j] -= Xre[i][j];
                        Rim[i][j] -= Xim[i][j];
                    }
                }
                sbr_qmf_synthesis(sbr, &sbr->data[0], Xre, Xim,
                                  out_pcm[0], out_samples);
                sbr_qmf_synthesis(sbr, &sbr->data[1], Rre, Rim,
                                  out_pcm[1], out_samples);
                for (i = 0; i < out_samples; i++)
                    out_pcm[1][i] += out_pcm[0][i];
            } else {
                sbr_qmf_synthesis(sbr, &sbr->data[0], Xre, Xim,
                                  out_pcm[0], out_samples);
                memcpy(out_pcm[1], out_pcm[0],
                       (size_t)out_samples * sizeof(float));
            }
        } else {
            sbr_qmf_synthesis(sbr, d, Xre, Xim, out_pcm[ch], out_samples);
        }

        d->Ypos ^= 1;
    }

    sbr->active = 0;
    return 0;
}

int aac_sbr_is_active(const aac_sbr_decoder* sbr)
{
    if (sbr == NULL)
        return 0;
    return sbr->active;
}
