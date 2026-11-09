/*
 * aacdec - MPEG-4 Parametric Stereo (HE-AACv2).
 *
 * Clean-room implementation of the baseline 20-band PS tools: IID/ICC
 * parsing, the hybrid 8/4/2-band analysis and synthesis filterbank, the
 * transient-controlled all-pass decorrelator and the stereo mixing stage.
 * IPD/OPD data is parsed and then ignored (baseline PS).
 */
#include "aac_ps.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "aac_ps_tables.h"


#define PS_PI 3.14159265358979323846

/* ------------------------------------------------------------------------- */
/* Tables generated from the standard formulas                               */
/* ------------------------------------------------------------------------- */

static float ps_ha[46][8][4];
static float ps_qfract[50][3][2];
static float ps_phi[50][2];
static float ps_f20[8][8][2];
static float ps_f34_0_12[12][8][2];
static float ps_f34_1_8[8][8][2];
static float ps_f34_2_4[4][8][2];
static float ps_qfract34[50][3][2];
static float ps_phi34[50][2];
static float ps_pd_re[512];
static float ps_pd_im[512];
static int ps_tables_ready;

typedef struct
{
    int16_t child[2];
    int16_t symbol;
} ps_huff_node;

static ps_huff_node* ps_trees[10];
static int ps_huff_ready;

static int ps_tree_build(int table)
{
    int symbols = (int)aac_ps_huff_sizes[table];
    int offset = (int)aac_ps_huff_offsets[table];
    const uint8_t(*tab)[2] = aac_ps_huff_tab;
    int i;
    uint32_t code = 0;
    int count = 1;
    int cap = symbols * 2 + 2;
    ps_huff_node* nodes;

    for (i = 0; i < table; i++)
        tab += aac_ps_huff_sizes[i];

    nodes = (ps_huff_node*)malloc((size_t)cap * sizeof(*nodes));
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
        if (len < 32 && (code & ((1u << (32 - len)) - 1u)) != 0) {
            free(nodes);
            return -1;
        }

        for (bit = len - 1; bit >= 0; bit--) {
            int b = (int)((code >> (32 - len + bit)) & 1u);
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
        code += 1u << (32 - len);
    }

    ps_trees[table] = nodes;
    return 0;
}

static int ps_huff_init(void)
{
    int i;
    if (ps_huff_ready)
        return 0;
    for (i = 0; i < 10; i++) {
        if (ps_tree_build(i) < 0) {
            int j;
            for (j = 0; j < i; j++) {
                free(ps_trees[j]);
                ps_trees[j] = NULL;
            }
            return -1;
        }
    }
    ps_huff_ready = 1;
    return 0;
}

static int ps_huff_decode(aac_bits* b, int table)
{
    int node = 0;
    ps_huff_node* nodes = ps_trees[table];

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

static void ps_make_filters(float (*filter)[8][2], const float* proto, int bands)
{
    int q, n;
    for (q = 0; q < bands; q++) {
        for (n = 0; n < 7; n++) {
            double theta = 2.0 * PS_PI * (q + 0.5) * (n - 6) / bands;
            filter[q][n][0] = (float)(proto[n] * cos(theta));
            filter[q][n][1] = (float)(proto[n] * -sin(theta));
        }
    }
}

static int ps_tables_init(void)
{
    static const float iid_par_dequant[46] = {
        0.05623413251903f, 0.12589254117942f, 0.19952623149689f,
        0.31622776601684f, 0.44668359215096f, 0.63095734448019f,
        0.79432823472428f, 1.0f, 1.25892541179417f, 1.58489319246111f,
        2.23872113856834f, 3.16227766016838f, 5.01187233627272f,
        7.94328234724282f, 17.7827941003892f,
        0.00316227766017f, 0.00562341325190f, 0.01f, 0.01778279410039f,
        0.03162277660168f, 0.05623413251903f, 0.07943282347243f,
        0.11220184543020f, 0.15848931924611f, 0.22387211385683f,
        0.31622776601684f, 0.39810717055350f, 0.50118723362727f,
        0.63095734448019f, 0.79432823472428f, 1.0f, 1.25892541179417f,
        1.58489319246111f, 1.99526231496888f, 2.51188643150958f,
        3.16227766016838f, 4.46683592150963f, 6.30957344480193f,
        8.91250938133745f, 12.5892541179417f, 17.7827941003892f,
        31.6227766016838f, 56.2341325190349f, 100.0f,
        177.827941003892f, 316.227766016837f
    };
    static const double acos_icc[8] = {
        0.0, 0.35685527, 0.57133466, 0.92614472, 1.1943263,
        PS_PI / 2.0, 2.2006171, PS_PI
    };
    static const int8_t f_center_20[10] = {
        -3, -1, 1, 3, 5, 7, 10, 14, 18, 22
    };
    static const float fractional_delay_links[3] = {0.43f, 0.75f, 0.347f};
    static const float fractional_delay_gain = 0.39f;
    static const float g0_Q8[7] = {
        0.00746082949812f, 0.02270420949825f, 0.04546865930473f,
        0.07266113929591f, 0.09885108575264f, 0.11793710567217f, 0.125f
    };
    static const float g0_Q12[7] = {
        0.04081179924692f, 0.03812810994926f, 0.05144908135699f,
        0.06399831151592f, 0.07428313801106f, 0.08100347892914f, 0.08333333333333f
    };
    static const float g1_Q8[7] = {
        0.01565675600122f, 0.03752716391991f, 0.05417891378782f,
        0.08417044116767f, 0.10307344158036f, 0.12222452249753f, 0.125f
    };
    static const float g2_Q4[7] = {
        -0.05908211155639f, -0.04871498374946f, 0.0f, 0.07778723915851f,
        0.16486303567403f, 0.23279856662996f, 0.25f
    };
    static const int8_t f_center_34[32] = {
        2, 6, 10, 14, 18, 22, 26, 30,
        34, -10, -6, -2, 51, 57, 15, 21,
        27, 33, 39, 45, 54, 66, 78, 42,
        102, 66, 78, 90, 102, 114, 126, 90
    };
    static const float ipdopd_sin[8] = {
        0.0f, 0.70710678118655f, 1.0f, 0.70710678118655f,
        0.0f, -0.70710678118655f, -1.0f, -0.70710678118655f
    };
    static const float ipdopd_cos[8] = {
        1.0f, 0.70710678118655f, 0.0f, -0.70710678118655f,
        -1.0f, -0.70710678118655f, 0.0f, 0.70710678118655f
    };
    int iid, icc, k, m;

    if (ps_tables_ready)
        return 0;
    if (ps_huff_init() < 0)
        return -1;

    for (iid = 0; iid < 46; iid++) {
        double c = iid_par_dequant[iid];
        double c1 = sqrt(2.0) / sqrt(1.0 + c * c);
        double c2 = c * c1;
        for (icc = 0; icc < 8; icc++) {
            double alpha = 0.5 * acos_icc[icc];
            double beta = alpha * (c1 - c2) / sqrt(2.0);
            ps_ha[iid][icc][0] = (float)(c2 * cos(beta + alpha));
            ps_ha[iid][icc][1] = (float)(c1 * cos(beta - alpha));
            ps_ha[iid][icc][2] = (float)(c2 * sin(beta + alpha));
            ps_ha[iid][icc][3] = (float)(c1 * sin(beta - alpha));
        }
    }

    for (k = 0; k < 50; k++) {
        double fc = (k < 10) ? (double)f_center_20[k] * 0.125 :
                               (double)k - 6.5;
        double theta;
        for (m = 0; m < 3; m++) {
            theta = -PS_PI * fractional_delay_links[m] * fc;
            ps_qfract[k][m][0] = (float)cos(theta);
            ps_qfract[k][m][1] = (float)sin(theta);
        }
        theta = -PS_PI * fractional_delay_gain * fc;
        ps_phi[k][0] = (float)cos(theta);
        ps_phi[k][1] = (float)sin(theta);
    }

    ps_make_filters(ps_f20, g0_Q8, 8);
    ps_make_filters(ps_f34_0_12, g0_Q12, 12);
    ps_make_filters(ps_f34_1_8, g1_Q8, 8);
    ps_make_filters(ps_f34_2_4, g2_Q4, 4);

    for (k = 0; k < 50; k++) {
        double fc = (k < 32) ? (double)f_center_34[k] / 24.0 :
                               (double)k - 26.5;
        double theta;
        for (m = 0; m < 3; m++) {
            theta = -PS_PI * 0.43 * fc;
            if (m == 1) theta = -PS_PI * 0.75 * fc;
            if (m == 2) theta = -PS_PI * 0.347 * fc;
            ps_qfract34[k][m][0] = (float)cos(theta);
            ps_qfract34[k][m][1] = (float)sin(theta);
        }
        theta = -PS_PI * 0.39 * fc;
        ps_phi34[k][0] = (float)cos(theta);
        ps_phi34[k][1] = (float)sin(theta);
    }

    for (k = 0; k < 8; k++) {
        int j, l;
        for (j = 0; j < 8; j++) {
            for (l = 0; l < 8; l++) {
                double re = 0.25 * ipdopd_cos[k] + 0.5 * ipdopd_cos[j] +
                            ipdopd_cos[l];
                double im = 0.25 * ipdopd_sin[k] + 0.5 * ipdopd_sin[j] +
                            ipdopd_sin[l];
                double mag = 1.0 / sqrt(re * re + im * im);
                ps_pd_re[k * 64 + j * 8 + l] = (float)(re * mag);
                ps_pd_im[k * 64 + j * 8 + l] = (float)(im * mag);
            }
        }
    }
    ps_tables_ready = 1;
    return 0;
}

/* ------------------------------------------------------------------------- */
/* PS bitstream parsing                                                      */
/* ------------------------------------------------------------------------- */

static const int8_t ps_num_env_tab[2][4] = {
    {0, 1, 2, 4},
    {1, 2, 3, 4}
};

static const int8_t ps_nr_iidicc[6] = {10, 20, 34, 10, 20, 34};

static int ps_read_iid(aac_ps_decoder* ps, aac_bits* b, int e, int dt)
{
    int bnd, num = ps->iid_symbols;
    int table = dt ? (ps->iid_quant ? 1 : 3) : (ps->iid_quant ? 0 : 2);
    int prev = e ? e - 1 : ps->num_env_old - 1;
    int max = 7 + 8 * ps->iid_quant;
    int val = 0;
    if (prev < 0)
        prev = 0;
    for (bnd = 0; bnd < num; bnd++) {
        int sym = ps_huff_decode(b, table);
        if (dt)
            val = ps->iid_par[prev][bnd] + sym;
        else
            val += sym;
        ps->iid_par[e][bnd] = (int8_t)val;
        if (val > max || val < -max) {
            return -1;
        }
    }
    return b->error ? -1 : 0;
}

static int ps_read_icc(aac_ps_decoder* ps, aac_bits* b, int e, int dt)
{
    int bnd, num = ps->icc_symbols;
    int table = dt ? 5 : 4;
    int prev = e ? e - 1 : ps->num_env_old - 1;
    int val = 0;
    if (prev < 0)
        prev = 0;
    for (bnd = 0; bnd < num; bnd++) {
        if (dt)
            val = ps->icc_par[prev][bnd] + ps_huff_decode(b, table);
        else
            val += ps_huff_decode(b, table);
        ps->icc_par[e][bnd] = (int8_t)val;
        if (val > 7 || val < 0) {
            return -1;
        }
    }
    return b->error ? -1 : 0;
}

/* Read one IPD or OPD parameter set. */
static int ps_read_phase(aac_ps_decoder* ps, aac_bits* b,
                         int8_t par[AAC_PS_MAX_ENV][AAC_PS_MAX_PAR],
                         int table_df, int table_dt, int e, int dt)
{
    int bnd, num = ps->nr_ipdopd_par;
    int val = 0;
    int prev = e ? e - 1 : ps->num_env_old - 1;
    if (prev < 0)
        prev = 0;
    for (bnd = 0; bnd < num; bnd++) {
        int sym;
        if (dt)
            val = par[prev][bnd] + ps_huff_decode(b, table_dt);
        else {
            sym = ps_huff_decode(b, table_df);
            val += sym;
        }
        val &= 7;
        par[e][bnd] = (int8_t)val;
    }
    return b->error ? -1 : 0;
}

static int ps_read_ipdopd(aac_ps_decoder* ps, aac_bits* b)
{
    int e;
    for (e = 0; e < ps->num_env; e++) {
        int dt = (int)aac_bits_read(b, 1);
        if (ps_read_phase(ps, b, ps->ipd_par, 6, 7, e, dt) < 0)
            return -1;
        dt = (int)aac_bits_read(b, 1);
        if (ps_read_phase(ps, b, ps->opd_par, 8, 9, e, dt) < 0)
            return -1;
    }
    (void)aac_bits_read(b, 1); /* reserved_ps */
    return b->error ? -1 : 0;
}

void aac_ps_reset(aac_ps_decoder* ps)
{
    if (ps != NULL)
        memset(ps, 0, sizeof(*ps));
}

int aac_ps_parse(aac_ps_decoder* ps, aac_bits* b, int bits_left)
{
    int start_bit;
    int header;
    int e;

    if (ps == NULL || b == NULL)
        return 0;
    if (ps_tables_init() < 0) {
        return 0;
    }

    start_bit = (int)b->bitpos;
    header = (int)aac_bits_read(b, 1);
    if (header) {
        ps->enable_iid = (int)aac_bits_read(b, 1);
        if (ps->enable_iid) {
            int mode = (int)aac_bits_read(b, 3);
            if (mode > 5)
                goto fail;
            ps->iid_symbols = ps_nr_iidicc[mode];
            ps->freq_res_iid = mode > 2 ? mode - 3 : mode;
            ps->nr_iid_par = ps->freq_res_iid == 0 ? 20 : ps->iid_symbols;
            ps->nr_ipdopd_par = (mode == 0 || mode == 3) ? 5 :
                                (mode == 1 || mode == 4) ? 11 : 17;
            ps->iid_quant = mode > 2;
        }
        ps->enable_icc = (int)aac_bits_read(b, 1);
        if (ps->enable_icc) {
            int mode = (int)aac_bits_read(b, 3);
            if (mode > 5)
                goto fail;
            ps->icc_mode = mode;
            ps->icc_symbols = ps_nr_iidicc[mode];
            ps->freq_res_icc = mode > 2 ? mode - 3 : mode;
            ps->nr_icc_par = ps->freq_res_icc == 0 ? 20 : ps->icc_symbols;
        }
        ps->enable_ext = (int)aac_bits_read(b, 1);
    }

    ps->frame_class = (int)aac_bits_read(b, 1);
    ps->num_env_old = ps->num_env;
    ps->num_env = ps_num_env_tab[ps->frame_class][
        (int)aac_bits_read(b, 2)];
    ps->border_position[0] = -1;
    if (ps->frame_class) {
        for (e = 1; e <= ps->num_env; e++) {
            ps->border_position[e] = (int)aac_bits_read(b, 5);
            if (ps->border_position[e] < ps->border_position[e - 1])
                goto fail;
        }
    } else {
        for (e = 1; e <= ps->num_env; e++) {
            int shift = 0;
            int n = ps->num_env;
            while ((1 << shift) < n)
                shift++;
            ps->border_position[e] = (e * AAC_PS_QMF_SLOTS >> shift) - 1;
        }
    }

    if (ps->enable_iid) {
        for (e = 0; e < ps->num_env; e++) {
            int dt = (int)aac_bits_read(b, 1);
            if (ps_read_iid(ps, b, e, dt) < 0)
                goto fail;
        }
    } else {
        memset(ps->iid_par, 0, sizeof(ps->iid_par));
    }

    if (ps->enable_icc) {
        for (e = 0; e < ps->num_env; e++) {
            int dt = (int)aac_bits_read(b, 1);
            if (ps_read_icc(ps, b, e, dt) < 0)
                goto fail;
        }
    } else {
        memset(ps->icc_par, 0, sizeof(ps->icc_par));
    }

    if (ps->enable_ext) {
        int cnt = (int)aac_bits_read(b, 4);
        if (cnt == 15)
            cnt += (int)aac_bits_read(b, 8);
        cnt *= 8;
        if (cnt > 7) {
            int id = (int)aac_bits_read(b, 2);
            if (id == 0) {
                ps->enable_ipdopd = (int)aac_bits_read(b, 1);
                if (ps->enable_ipdopd) {
                    if (ps_read_ipdopd(ps, b) < 0)
                        goto fail;
                }
            }
            /* Consume any remaining extension/padding bits. */
            {
                int used = (int)b->bitpos;
                int end = start_bit + bits_left - (cnt - 2);
                if (end > used)
                    aac_bits_skip(b, (size_t)(end - used));
                else if (used > end)
                    b->bitpos = (size_t)end;
            }
        }
    }

    /* Create the final pseudo-envelope when needed. */
    if (ps->num_env == 0 ||
        ps->border_position[ps->num_env] < AAC_PS_QMF_SLOTS - 1) {
        int source = ps->num_env ? ps->num_env - 1 : ps->num_env_old - 1;
        if (source >= 0 && source != ps->num_env) {
            if (ps->enable_iid)
                memcpy(ps->iid_par[ps->num_env], ps->iid_par[source],
                       sizeof(ps->iid_par[0]));
            if (ps->enable_icc)
                memcpy(ps->icc_par[ps->num_env], ps->icc_par[source],
                       sizeof(ps->icc_par[0]));
        }
        ps->num_env++;
        ps->border_position[ps->num_env] = AAC_PS_QMF_SLOTS - 1;
    }

    ps->is34bands_old = ps->is34bands;
    ps->is34bands = (ps->enable_iid && ps->nr_iid_par == 34) ||
                    (ps->enable_icc && ps->nr_icc_par == 34);
    if (header)
        ps->start = 1;

    {
        int consumed = (int)b->bitpos - start_bit;
        if (consumed <= bits_left + 15) {
            if (consumed < bits_left)
                aac_bits_skip(b, (size_t)(bits_left - consumed));
            return consumed;
        }
    }

fail:
    ps->start = 0;
    memset(ps->iid_par, 0, sizeof(ps->iid_par));
    memset(ps->icc_par, 0, sizeof(ps->icc_par));
    aac_bits_skip(b, (size_t)bits_left);
    return bits_left;
}

/* ------------------------------------------------------------------------- */
/* Hybrid filterbank                                                         */
/* ------------------------------------------------------------------------- */

static void ps_hybrid_analysis(float out[8][2], const float in[][2],
                               const float (*filter)[8][2], int n)
{
    int i, j;
    for (i = 0; i < n; i++) {
        double sum_re = (double)filter[i][6][0] * in[6][0];
        double sum_im = (double)filter[i][6][0] * in[6][1];
        for (j = 0; j < 6; j++) {
            double inre0 = in[j][0] + in[12 - j][0];
            double inre1 = in[j][1] - in[12 - j][1];
            double inim0 = in[j][1] + in[12 - j][1];
            double inim1 = in[j][0] - in[12 - j][0];
            sum_re += (double)filter[i][j][0] * inre0 -
                      (double)filter[i][j][1] * inre1;
            sum_im += (double)filter[i][j][0] * inim0 +
                      (double)filter[i][j][1] * inim1;
        }
        out[i][0] = (float)sum_re;
        out[i][1] = (float)sum_im;
    }
}

static void ps_hybrid6(float in[][2], float out[][32][2], int len)
{
    int i;
    float temp[8][2];
    for (i = 0; i < len; i++, in++) {
        ps_hybrid_analysis(temp, (const float(*)[2])in, ps_f20, 8);
        out[0][i][0] = temp[6][0];
        out[0][i][1] = temp[6][1];
        out[1][i][0] = temp[7][0];
        out[1][i][1] = temp[7][1];
        out[2][i][0] = temp[0][0];
        out[2][i][1] = temp[0][1];
        out[3][i][0] = temp[1][0];
        out[3][i][1] = temp[1][1];
        out[4][i][0] = temp[2][0] + temp[5][0];
        out[4][i][1] = temp[2][1] + temp[5][1];
        out[5][i][0] = temp[3][0] + temp[4][0];
        out[5][i][1] = temp[3][1] + temp[4][1];
    }
}

static void ps_hybrid2(float in[][2], float out[][32][2], int len, int reverse)
{
    static const float g1[7] = {
        0.0f, 0.01899487526049f, 0.0f, -0.07293139167538f,
        0.0f, 0.30596630545168f, 0.5f
    };
    int i, j;
    for (i = 0; i < len; i++, in++) {
        double re_in = (double)g1[6] * in[6][0];
        double im_in = (double)g1[6] * in[6][1];
        double re_op = 0.0, im_op = 0.0;
        for (j = 0; j < 6; j += 2) {
            re_op += (double)g1[j + 1] * (in[j + 1][0] + in[12 - j - 1][0]);
            im_op += (double)g1[j + 1] * (in[j + 1][1] + in[12 - j - 1][1]);
        }
        out[reverse][i][0] = (float)(re_in + re_op);
        out[reverse][i][1] = (float)(im_in + im_op);
        out[!reverse][i][0] = (float)(re_in - re_op);
        out[!reverse][i][1] = (float)(im_in - im_op);
    }
}

static void ps_hybrid_analysis_all(aac_ps_decoder* ps,
                                   float Lre[38][64], float Lim[38][64],
                                   float out[91][32][2])
{
    int i, j;
    for (i = 0; i < 5; i++) {
        for (j = 0; j < 38; j++) {
            ps->in_buf[i][j + 6][0] = Lre[j][i];
            ps->in_buf[i][j + 6][1] = Lim[j][i];
        }
    }
    if (ps->is34bands) {
        for (j = 0; j < 32; j++) {
            float temp[12][2];
            int b;
            ps_hybrid_analysis(temp, ps->in_buf[0] + j, ps_f34_0_12, 12);
            for (b = 0; b < 12; b++) {
                out[b][j][0] = temp[b][0];
                out[b][j][1] = temp[b][1];
            }
            ps_hybrid_analysis(temp, ps->in_buf[1] + j, ps_f34_1_8, 8);
            for (b = 0; b < 8; b++) {
                out[12 + b][j][0] = temp[b][0];
                out[12 + b][j][1] = temp[b][1];
            }
            for (b = 0; b < 4; b++) {
                ps_hybrid_analysis(temp, ps->in_buf[2 + b] + j,
                                   ps_f34_2_4, 4);
                {
                    int c;
                    for (c = 0; c < 4; c++) {
                        out[20 + b * 4 + c][j][0] = temp[c][0];
                        out[20 + b * 4 + c][j][1] = temp[c][1];
                    }
                }
            }
            for (b = 5; b < 64; b++) {
                out[b + 27][j][0] = Lre[j][b];
                out[b + 27][j][1] = Lim[j][b];
            }
        }
    } else {
        ps_hybrid6(ps->in_buf[0], out, 32);
        ps_hybrid2(ps->in_buf[1], out + 6, 32, 1);
        ps_hybrid2(ps->in_buf[2], out + 8, 32, 0);
        for (i = 3; i < 64; i++) {
            for (j = 0; j < 32; j++) {
                out[i + 7][j][0] = Lre[j][i];
                out[i + 7][j][1] = Lim[j][i];
            }
        }
    }
    for (i = 0; i < 5; i++) {
        memcpy(ps->in_buf[i], ps->in_buf[i] + 32,
               6 * sizeof(ps->in_buf[i][0]));
    }
}

static void ps_hybrid_synthesis_all(aac_ps_decoder* ps,
                                    float Lre[38][64], float Lim[38][64],
                                    float Rre[38][64], float Rim[38][64])
{
    int n, i;
    if (ps->is34bands) {
        for (n = 0; n < 32; n++) {
            Lre[n][0] = 0.0f; Lim[n][0] = 0.0f;
            Rre[n][0] = 0.0f; Rim[n][0] = 0.0f;
            for (i = 0; i < 12; i++) {
                Lre[n][0] += ps->Lbuf[i][n][0];
                Lim[n][0] += ps->Lbuf[i][n][1];
                Rre[n][0] += ps->Rbuf[i][n][0];
                Rim[n][0] += ps->Rbuf[i][n][1];
            }
            Lre[n][1] = 0.0f; Lim[n][1] = 0.0f;
            Rre[n][1] = 0.0f; Rim[n][1] = 0.0f;
            for (i = 0; i < 8; i++) {
                Lre[n][1] += ps->Lbuf[12 + i][n][0];
                Lim[n][1] += ps->Lbuf[12 + i][n][1];
                Rre[n][1] += ps->Rbuf[12 + i][n][0];
                Rim[n][1] += ps->Rbuf[12 + i][n][1];
            }
            Lre[n][2] = Lim[n][2] = 0.0f;
            Lre[n][3] = Lim[n][3] = 0.0f;
            Lre[n][4] = Lim[n][4] = 0.0f;
            Rre[n][2] = Rim[n][2] = 0.0f;
            Rre[n][3] = Rim[n][3] = 0.0f;
            Rre[n][4] = Rim[n][4] = 0.0f;
            for (i = 0; i < 4; i++) {
                Lre[n][2] = Lre[n][2] + ps->Lbuf[20 + i][n][0];
                Lim[n][2] = Lim[n][2] + ps->Lbuf[20 + i][n][1];
                Rre[n][2] = Rre[n][2] + ps->Rbuf[20 + i][n][0];
                Rim[n][2] = Rim[n][2] + ps->Rbuf[20 + i][n][1];
                Lre[n][3] = Lre[n][3] + ps->Lbuf[24 + i][n][0];
                Lim[n][3] = Lim[n][3] + ps->Lbuf[24 + i][n][1];
                Rre[n][3] = Rre[n][3] + ps->Rbuf[24 + i][n][0];
                Rim[n][3] = Rim[n][3] + ps->Rbuf[24 + i][n][1];
                Lre[n][4] = Lre[n][4] + ps->Lbuf[28 + i][n][0];
                Lim[n][4] = Lim[n][4] + ps->Lbuf[28 + i][n][1];
                Rre[n][4] = Rre[n][4] + ps->Rbuf[28 + i][n][0];
                Rim[n][4] = Rim[n][4] + ps->Rbuf[28 + i][n][1];
            }
        }
        for (i = 5; i < 64; i++) {
            for (n = 0; n < 32; n++) {
                Lre[n][i] = ps->Lbuf[i + 27][n][0];
                Lim[n][i] = ps->Lbuf[i + 27][n][1];
                Rre[n][i] = ps->Rbuf[i + 27][n][0];
                Rim[n][i] = ps->Rbuf[i + 27][n][1];
            }
        }
    } else {
        for (n = 0; n < 32; n++) {
            Lre[n][0] = ps->Lbuf[0][n][0] + ps->Lbuf[1][n][0] +
                        ps->Lbuf[2][n][0] + ps->Lbuf[3][n][0] +
                        ps->Lbuf[4][n][0] + ps->Lbuf[5][n][0];
            Lim[n][0] = ps->Lbuf[0][n][1] + ps->Lbuf[1][n][1] +
                        ps->Lbuf[2][n][1] + ps->Lbuf[3][n][1] +
                        ps->Lbuf[4][n][1] + ps->Lbuf[5][n][1];
            Rre[n][0] = ps->Rbuf[0][n][0] + ps->Rbuf[1][n][0] +
                        ps->Rbuf[2][n][0] + ps->Rbuf[3][n][0] +
                        ps->Rbuf[4][n][0] + ps->Rbuf[5][n][0];
            Rim[n][0] = ps->Rbuf[0][n][1] + ps->Rbuf[1][n][1] +
                        ps->Rbuf[2][n][1] + ps->Rbuf[3][n][1] +
                        ps->Rbuf[4][n][1] + ps->Rbuf[5][n][1];

            Lre[n][1] = ps->Lbuf[6][n][0] + ps->Lbuf[7][n][0];
            Lim[n][1] = ps->Lbuf[6][n][1] + ps->Lbuf[7][n][1];
            Rre[n][1] = ps->Rbuf[6][n][0] + ps->Rbuf[7][n][0];
            Rim[n][1] = ps->Rbuf[6][n][1] + ps->Rbuf[7][n][1];

            Lre[n][2] = ps->Lbuf[8][n][0] + ps->Lbuf[9][n][0];
            Lim[n][2] = ps->Lbuf[8][n][1] + ps->Lbuf[9][n][1];
            Rre[n][2] = ps->Rbuf[8][n][0] + ps->Rbuf[9][n][0];
            Rim[n][2] = ps->Rbuf[8][n][1] + ps->Rbuf[9][n][1];
        }
        for (i = 3; i < 64; i++) {
            for (n = 0; n < 32; n++) {
                Lre[n][i] = ps->Lbuf[i + 7][n][0];
                Lim[n][i] = ps->Lbuf[i + 7][n][1];
                Rre[n][i] = ps->Rbuf[i + 7][n][0];
                Rim[n][i] = ps->Rbuf[i + 7][n][1];
            }
        }
    }
}

/* ------------------------------------------------------------------------- */
/* Decorrelation and stereo mixing                                           */
/* ------------------------------------------------------------------------- */

static void ps_decorrelate_slot(float out[][2], float delay[][2],
                                float ap_delay[3][37][2],
                                const float phi[2], const float qf[3][2],
                                const float* transient_gain, float g_decay,
                                int len)
{
    static const float a[3] = {
        0.65143905753106f, 0.56471812200776f, 0.48954165955695f
    };
    float ag[3];
    int m, n;

    for (m = 0; m < 3; m++)
        ag[m] = a[m] * g_decay;

    for (n = 0; n < len; n++) {
        float in_re = delay[n][0] * phi[0] - delay[n][1] * phi[1];
        float in_im = delay[n][0] * phi[1] + delay[n][1] * phi[0];
        for (m = 0; m < 3; m++) {
            float a_re = ag[m] * in_re;
            float a_im = ag[m] * in_im;
            float ld_re = ap_delay[m][n + 2 - m][0];
            float ld_im = ap_delay[m][n + 2 - m][1];
            float apd_re = in_re;
            float apd_im = in_im;
            in_re = ld_re * qf[m][0] - ld_im * qf[m][1] - a_re;
            in_im = ld_re * qf[m][1] + ld_im * qf[m][0] - a_im;
            ap_delay[m][n + 5][0] = apd_re + ag[m] * in_re;
            ap_delay[m][n + 5][1] = apd_im + ag[m] * in_im;
        }
        out[n][0] = transient_gain[n] * in_re;
        out[n][1] = transient_gain[n] * in_im;
    }
}

static void ps_decorrelate(aac_ps_decoder* ps, float out[91][32][2],
                           const float s[91][32][2])
{
    const int is34 = ps->is34bands;
    const int nr_bands = is34 ? 91 : 71;
    const int nr_par = is34 ? 34 : 20;
    const int nr_allpass = is34 ? 50 : 30;
    const int decay_cutoff = is34 ? 32 : 10;
    const int short_delay = is34 ? 62 : 42;
    const int8_t* k_to_i = is34 ? aac_ps_k_to_i_34 : aac_ps_k_to_i_20;
    const float (*phi)[2] = is34 ? ps_phi34 : ps_phi;
    const float (*qf)[3][2] = is34 ? ps_qfract34 : ps_qfract;
    float power[34][32];
    float transient_gain[34][32];
    int i, k, m, n;

    memset(power, 0, sizeof(power));
    for (k = 0; k < nr_bands; k++) {
        i = k_to_i[k];
        for (n = 0; n < 32; n++)
            power[i][n] += s[k][n][0] * s[k][n][0] +
                           s[k][n][1] * s[k][n][1];
    }
    for (i = 0; i < nr_par; i++) {
        for (n = 0; n < 32; n++) {
            float decayed = 0.76592833836465f * ps->peak_decay_nrg[i];
            float denom;
            ps->peak_decay_nrg[i] = decayed > power[i][n] ?
                                        decayed : power[i][n];
            ps->power_smooth[i] +=
                0.25f * (power[i][n] - ps->power_smooth[i]);
            ps->peak_decay_diff_smooth[i] +=
                0.25f * (ps->peak_decay_nrg[i] - power[i][n] -
                         ps->peak_decay_diff_smooth[i]);
            denom = 1.5f * ps->peak_decay_diff_smooth[i];
            transient_gain[i][n] = denom > ps->power_smooth[i] ?
                                       ps->power_smooth[i] / denom : 1.0f;
        }
    }

    for (k = 0; k < nr_allpass; k++) {
        float g_decay = 1.0f - 0.05f * (float)(k - decay_cutoff);
        if (g_decay < 0.0f)
            g_decay = 0.0f;
        if (g_decay > 1.0f)
            g_decay = 1.0f;
        memcpy(ps->delay[k], ps->delay[k] + 32,
               14 * sizeof(ps->delay[k][0]));
        memcpy(ps->delay[k] + 14, s[k], 32 * sizeof(ps->delay[k][0]));
        for (m = 0; m < 3; m++)
            memcpy(ps->ap_delay[k][m], ps->ap_delay[k][m] + 32,
                   5 * sizeof(ps->ap_delay[k][m][0]));
        ps_decorrelate_slot(out[k], ps->delay[k] + 12, ps->ap_delay[k],
                            phi[k], qf[k], transient_gain[k_to_i[k]],
                            g_decay, 32);
    }
    for (; k < short_delay; k++) {
        i = k_to_i[k];
        memcpy(ps->delay[k], ps->delay[k] + 32,
               14 * sizeof(ps->delay[k][0]));
        memcpy(ps->delay[k] + 14, s[k], 32 * sizeof(ps->delay[k][0]));
        for (n = 0; n < 32; n++) {
            out[k][n][0] = ps->delay[k][n][0] * transient_gain[i][n];
            out[k][n][1] = ps->delay[k][n][1] * transient_gain[i][n];
        }
    }
    for (; k < nr_bands; k++) {
        i = k_to_i[k];
        memcpy(ps->delay[k], ps->delay[k] + 32,
               14 * sizeof(ps->delay[k][0]));
        memcpy(ps->delay[k] + 14, s[k], 32 * sizeof(ps->delay[k][0]));
        for (n = 0; n < 32; n++) {
            out[k][n][0] = ps->delay[k][13 + n][0] * transient_gain[i][n];
            out[k][n][1] = ps->delay[k][13 + n][1] * transient_gain[i][n];
        }
    }
}

static void ps_map_10_to_20(int8_t* dst, const int8_t* src, int full)
{
    int b = full ? 9 : 4;
    if (!full)
        dst[10] = 0;
    for (; b >= 0; b--)
        dst[2 * b + 1] = dst[2 * b] = src[b];
}

static void ps_map_34_to_20(int8_t* dst, const int8_t* src, int full)
{
    dst[0] = (int8_t)((2 * src[0] + src[1]) / 3);
    dst[1] = (int8_t)((src[1] + 2 * src[2]) / 3);
    dst[2] = (int8_t)((2 * src[3] + src[4]) / 3);
    dst[3] = (int8_t)((src[4] + 2 * src[5]) / 3);
    dst[4] = (int8_t)((src[6] + src[7]) / 2);
    dst[5] = (int8_t)((src[8] + src[9]) / 2);
    dst[6] = src[10];
    dst[7] = src[11];
    dst[8] = (int8_t)((src[12] + src[13]) / 2);
    dst[9] = (int8_t)((src[14] + src[15]) / 2);
    dst[10] = src[16];
    if (full) {
        dst[11] = src[17];
        dst[12] = src[18];
        dst[13] = src[19];
        dst[14] = (int8_t)((src[20] + src[21]) / 2);
        dst[15] = (int8_t)((src[22] + src[23]) / 2);
        dst[16] = (int8_t)((src[24] + src[25]) / 2);
        dst[17] = (int8_t)((src[26] + src[27]) / 2);
        dst[18] = (int8_t)((src[28] + src[29] + src[30] + src[31]) / 4);
        dst[19] = (int8_t)((src[32] + src[33]) / 2);
    }
}

static void ps_map_10_to_34(int8_t* dst, const int8_t* src, int full)
{
    if (full) {
        dst[33] = dst[32] = dst[31] = dst[30] = dst[29] = dst[28] = src[9];
        dst[27] = dst[26] = dst[25] = dst[24] = src[8];
        dst[23] = dst[22] = dst[21] = dst[20] = src[7];
        dst[19] = dst[18] = src[6];
        dst[17] = dst[16] = src[5];
    } else {
        dst[16] = 0;
    }
    dst[15] = dst[14] = dst[13] = dst[12] = src[4];
    dst[11] = dst[10] = src[3];
    dst[9] = dst[8] = dst[7] = dst[6] = src[2];
    dst[5] = dst[4] = dst[3] = src[1];
    dst[2] = dst[1] = dst[0] = src[0];
}

static void ps_map_20_to_34(int8_t* dst, const int8_t* src, int full)
{
    if (full) {
        dst[33] = dst[32] = src[19];
        dst[31] = dst[30] = dst[29] = dst[28] = src[18];
        dst[27] = dst[26] = src[17];
        dst[25] = dst[24] = src[16];
        dst[23] = dst[22] = src[15];
        dst[21] = dst[20] = src[14];
        dst[19] = src[13];
        dst[18] = src[12];
        dst[17] = src[11];
    }
    dst[16] = src[10];
    dst[15] = dst[14] = src[9];
    dst[13] = dst[12] = src[8];
    dst[11] = src[7];
    dst[10] = src[6];
    dst[9] = dst[8] = src[5];
    dst[7] = dst[6] = src[4];
    dst[5] = src[3];
    dst[4] = (int8_t)((src[2] + src[3]) / 2);
    dst[3] = src[2];
    dst[2] = src[1];
    dst[1] = (int8_t)((src[0] + src[1]) / 2);
    dst[0] = src[0];
}

static void ps_map_val_34_to_20(float* p)
{
    p[0] = (2*p[0] + p[1]) * 0.33333333f;
    p[1] = (p[1] + 2*p[2]) * 0.33333333f;
    p[2] = (2*p[3] + p[4]) * 0.33333333f;
    p[3] = (p[4] + 2*p[5]) * 0.33333333f;
    p[4] = (p[6] + p[7]) * 0.5f;
    p[5] = (p[8] + p[9]) * 0.5f;
    p[6] = p[10];
    p[7] = p[11];
    p[8] = (p[12] + p[13]) * 0.5f;
    p[9] = (p[14] + p[15]) * 0.5f;
    p[10] = p[16];
    p[11] = p[17];
    p[12] = p[18];
    p[13] = p[19];
    p[14] = (p[20] + p[21]) * 0.5f;
    p[15] = (p[22] + p[23]) * 0.5f;
    p[16] = (p[24] + p[25]) * 0.5f;
    p[17] = (p[26] + p[27]) * 0.5f;
    p[18] = (p[28] + p[29] + p[30] + p[31]) * 0.25f;
    p[19] = (p[32] + p[33]) * 0.5f;
}

static void ps_map_val_20_to_34(float* p)
{
    p[33] = p[19];
    p[32] = p[19];
    p[31] = p[18];
    p[30] = p[18];
    p[29] = p[18];
    p[28] = p[18];
    p[27] = p[17];
    p[26] = p[17];
    p[25] = p[16];
    p[24] = p[16];
    p[23] = p[15];
    p[22] = p[15];
    p[21] = p[14];
    p[20] = p[14];
    p[19] = p[13];
    p[18] = p[12];
    p[17] = p[11];
    p[16] = p[10];
    p[15] = p[9];
    p[14] = p[9];
    p[13] = p[8];
    p[12] = p[8];
    p[11] = p[7];
    p[10] = p[6];
    p[9] = p[5];
    p[8] = p[5];
    p[7] = p[4];
    p[6] = p[4];
    p[5] = p[3];
    p[4] = (p[2] + p[3]) * 0.5f;
    p[3] = p[2];
    p[2] = p[1];
    p[1] = (p[0] + p[1]) * 0.5f;
}

static void ps_remap34(int8_t dst[5][34], const int8_t src[5][34], int num,
                       int num_env)
{
    int e;
    if (num == 34) {
        for (e=0;e<num_env;e++) memcpy(dst[e], src[e], 34);
    } else if (num == 20 || num == 11) {
        for (e=0;e<num_env;e++) ps_map_20_to_34(dst[e], src[e], 1);
    } else if (num == 10 || num == 5) {
        for (e=0;e<num_env;e++) ps_map_10_to_34(dst[e], src[e], 1);
    } else {
        for (e=0;e<num_env;e++) { memset(dst[e],0,34); memcpy(dst[e],src[e],17); }
    }
}

static void ps_remap20(int8_t dst[5][34], const int8_t src[5][34], int num,
                       int num_env)
{
    int e;
    if (num == 20) {
        for (e = 0; e < num_env; e++)
            memcpy(dst[e], src[e], 20);
    } else if (num == 10) {
        for (e = 0; e < num_env; e++)
            ps_map_10_to_20(dst[e], src[e], 1);
    } else if (num == 34) {
        for (e = 0; e < num_env; e++)
            ps_map_34_to_20(dst[e], src[e], 1);
    } else {
        for (e = 0; e < num_env; e++) {
            memset(dst[e], 0, 20);
            memcpy(dst[e], src[e], 16);
        }
    }
}

static void ps_stereo_interpolate(float l[][2], float r[][2],
                                  const float h0[4], const float h1[4],
                                  float width, int len)
{
    float h[4];
    float hs[4];
    int n, j;
    for (j = 0; j < 4; j++) {
        h[j] = h0[j];
        hs[j] = (h1[j] - h0[j]) * width;
    }
    for (n = 0; n < len; n++) {
        float lre = l[n][0], lim = l[n][1];
        float rre = r[n][0], rim = r[n][1];
        for (j = 0; j < 4; j++)
            h[j] += hs[j];
        l[n][0] = h[0] * lre + h[2] * rre;
        l[n][1] = h[0] * lim + h[2] * rim;
        r[n][0] = h[1] * lre + h[3] * rre;
        r[n][1] = h[1] * lim + h[3] * rim;
    }
}

static void ps_stereo_interpolate_ipdopd(float l[][2], float r[][2],
                                         const float h0[4], const float h0i[4],
                                         const float h1[4], const float h1i[4],
                                         float width, int len)
{
    float hr[4], hi[4], hsr[4], hsi[4];
    int n, j;
    for (j = 0; j < 4; j++) {
        hr[j] = h0[j];
        hi[j] = h0i[j];
        hsr[j] = (h1[j] - h0[j]) * width;
        hsi[j] = (h1i[j] - h0i[j]) * width;
    }
    for (n = 0; n < len; n++) {
        float lre = l[n][0], lim = l[n][1];
        float rre = r[n][0], rim = r[n][1];
        for (j = 0; j < 4; j++) {
            hr[j] += hsr[j];
            hi[j] += hsi[j];
        }
        l[n][0] = hr[0]*lre + hr[2]*rre - hi[0]*lim - hi[2]*rim;
        l[n][1] = hr[0]*lim + hr[2]*rim + hi[0]*lre + hi[2]*rre;
        r[n][0] = hr[1]*lre + hr[3]*rre - hi[1]*lim - hi[3]*rim;
        r[n][1] = hr[1]*lim + hr[3]*rim + hi[1]*lre + hi[3]*rre;
    }
}

static void ps_stereo_processing(aac_ps_decoder* ps)
{
    const int is34 = ps->is34bands;
    const int nbands = is34 ? 34 : 20;
    const int nr_bands = is34 ? 91 : 71;
    const int nr_ipd = is34 ? 17 : 11;
    const int8_t* k_to_i = is34 ? aac_ps_k_to_i_34 : aac_ps_k_to_i_20;
    int8_t iid_m[5][34], icc_m[5][34], ipd_m[5][34], opd_m[5][34];
    int e, b, k;

    if (is34) {
        ps_remap34(iid_m, ps->iid_par, ps->nr_iid_par, ps->num_env);
        ps_remap34(icc_m, ps->icc_par, ps->nr_icc_par, ps->num_env);
        if (ps->enable_ipdopd) {
            ps_remap34(ipd_m, ps->ipd_par, ps->nr_ipdopd_par, ps->num_env);
            ps_remap34(opd_m, ps->opd_par, ps->nr_ipdopd_par, ps->num_env);
        }
    } else {
        ps_remap20(iid_m, ps->iid_par, ps->nr_iid_par, ps->num_env);
        ps_remap20(icc_m, ps->icc_par, ps->nr_icc_par, ps->num_env);
        if (ps->enable_ipdopd) {
            ps_remap20(ipd_m, ps->ipd_par, ps->nr_ipdopd_par, ps->num_env);
            ps_remap20(opd_m, ps->opd_par, ps->nr_ipdopd_par, ps->num_env);
        }
    }

    if (ps->num_env_old) {
        memcpy(ps->H11[0][0], ps->H11[0][ps->num_env_old], 34*sizeof(float));
        memcpy(ps->H12[0][0], ps->H12[0][ps->num_env_old], 34*sizeof(float));
        memcpy(ps->H21[0][0], ps->H21[0][ps->num_env_old], 34*sizeof(float));
        memcpy(ps->H22[0][0], ps->H22[0][ps->num_env_old], 34*sizeof(float));
        memcpy(ps->H11[1][0], ps->H11[1][ps->num_env_old], 34*sizeof(float));
        memcpy(ps->H12[1][0], ps->H12[1][ps->num_env_old], 34*sizeof(float));
        memcpy(ps->H21[1][0], ps->H21[1][ps->num_env_old], 34*sizeof(float));
        memcpy(ps->H22[1][0], ps->H22[1][ps->num_env_old], 34*sizeof(float));
    }
    if (is34 && !ps->is34bands_old) {
        ps_map_val_20_to_34(ps->H11[0][0]);
        ps_map_val_20_to_34(ps->H12[0][0]);
        ps_map_val_20_to_34(ps->H21[0][0]);
        ps_map_val_20_to_34(ps->H22[0][0]);
        ps_map_val_20_to_34(ps->H11[1][0]);
        ps_map_val_20_to_34(ps->H12[1][0]);
        ps_map_val_20_to_34(ps->H21[1][0]);
        ps_map_val_20_to_34(ps->H22[1][0]);
        memset(ps->ipd_hist, 0, sizeof(ps->ipd_hist));
        memset(ps->opd_hist, 0, sizeof(ps->opd_hist));
    } else if (!is34 && ps->is34bands_old) {
        ps_map_val_34_to_20(ps->H11[0][0]);
        ps_map_val_34_to_20(ps->H12[0][0]);
        ps_map_val_34_to_20(ps->H21[0][0]);
        ps_map_val_34_to_20(ps->H22[0][0]);
        ps_map_val_34_to_20(ps->H11[1][0]);
        ps_map_val_34_to_20(ps->H12[1][0]);
        ps_map_val_34_to_20(ps->H21[1][0]);
        ps_map_val_34_to_20(ps->H22[1][0]);
        memset(ps->ipd_hist, 0, sizeof(ps->ipd_hist));
        memset(ps->opd_hist, 0, sizeof(ps->opd_hist));
    }

    for (e = 0; e < ps->num_env; e++) {
        for (b = 0; b < nbands; b++) {
            int iid = iid_m[e][b] + 7 + 23 * ps->iid_quant;
            int icc = icc_m[e][b] & 7;
            float h11, h12, h21, h22;
            if (iid < 0)
                iid = 0;
            if (iid > 45)
                iid = 45;
            h11 = ps_ha[iid][icc][0];
            h12 = ps_ha[iid][icc][1];
            h21 = ps_ha[iid][icc][2];
            h22 = ps_ha[iid][icc][3];
            ps->H11[0][e + 1][b] = h11;
            ps->H12[0][e + 1][b] = h12;
            ps->H21[0][e + 1][b] = h21;
            ps->H22[0][e + 1][b] = h22;
            if (ps->enable_ipdopd && b < nr_ipd) {
                int opd_idx = ps->opd_hist[b] * 8 + opd_m[e][b];
                int ipd_idx = ps->ipd_hist[b] * 8 + ipd_m[e][b];
                float opd_re = ps_pd_re[opd_idx & 0x1ff];
                float opd_im = ps_pd_im[opd_idx & 0x1ff];
                float ipd_re = ps_pd_re[ipd_idx & 0x1ff];
                float ipd_im = ps_pd_im[ipd_idx & 0x1ff];
                float ipd_adj_re = opd_re * ipd_re + opd_im * ipd_im;
                float ipd_adj_im = opd_im * ipd_re - opd_re * ipd_im;
                ps->opd_hist[b] = (int8_t)(opd_idx & 0x3f);
                ps->ipd_hist[b] = (int8_t)(ipd_idx & 0x3f);
                ps->H11[1][e + 1][b] = h11 * opd_im;
                ps->H12[1][e + 1][b] = h12 * ipd_adj_im;
                ps->H21[1][e + 1][b] = h21 * opd_im;
                ps->H22[1][e + 1][b] = h22 * ipd_adj_im;
                ps->H11[0][e + 1][b] = h11 * opd_re;
                ps->H12[0][e + 1][b] = h12 * ipd_adj_re;
                ps->H21[0][e + 1][b] = h21 * opd_re;
                ps->H22[0][e + 1][b] = h22 * ipd_adj_re;
            }
        }
        for (k = 0; k < nr_bands; k++) {
            int start = ps->border_position[e];
            int stop = ps->border_position[e + 1];
            int len = stop - start;
            float width = 1.0f / (float)(len > 0 ? len : 1);
            float h0[4], h1[4], h0i[4], h1i[4];
            b = k_to_i[k];
            h0[0] = ps->H11[0][e][b];
            h0[1] = ps->H12[0][e][b];
            h0[2] = ps->H21[0][e][b];
            h0[3] = ps->H22[0][e][b];
            h1[0] = ps->H11[0][e + 1][b];
            h1[1] = ps->H12[0][e + 1][b];
            h1[2] = ps->H21[0][e + 1][b];
            h1[3] = ps->H22[0][e + 1][b];
            if (ps->enable_ipdopd) {
                int neg = (is34 && k <= 13 && k >= 9) || (!is34 && k <= 1);
                h0i[0] = neg ? -ps->H11[1][e][b] : ps->H11[1][e][b];
                h0i[1] = neg ? -ps->H12[1][e][b] : ps->H12[1][e][b];
                h0i[2] = neg ? -ps->H21[1][e][b] : ps->H21[1][e][b];
                h0i[3] = neg ? -ps->H22[1][e][b] : ps->H22[1][e][b];
                h1i[0] = ps->H11[1][e + 1][b];
                h1i[1] = ps->H12[1][e + 1][b];
                h1i[2] = ps->H21[1][e + 1][b];
                h1i[3] = ps->H22[1][e + 1][b];
                if (len > 0)
                    ps_stereo_interpolate_ipdopd(&ps->Lbuf[k][1 + start],
                                                 &ps->Rbuf[k][1 + start],
                                                 h0, h0i, h1, h1i, width, len);
            } else if (len > 0) {
                ps_stereo_interpolate(&ps->Lbuf[k][1 + start],
                                      &ps->Rbuf[k][1 + start],
                                      h0, h1, width, len);
            }
        }
    }
}

void aac_ps_apply(aac_ps_decoder* ps, float Lre[38][64],
                  float Lim[38][64], float Rre[38][64],
                  float Rim[38][64], int top)
{
    if (ps == NULL || !ps->start || ps->num_env <= 0)
        return;

    top += 71 - 64;
    if (top < 0)
        top = 0;
    if (top < 91)
        memset(&ps->delay[top][0][0], 0,
               (size_t)(91 - top) * sizeof(ps->delay[0]));
    if (top < 50)
        memset(&ps->ap_delay[top][0][0][0], 0,
               (size_t)(50 - top) * sizeof(ps->ap_delay[0]));

    ps_hybrid_analysis_all(ps, Lre, Lim, ps->Lbuf);
    ps_decorrelate(ps, ps->Rbuf, (const float(*)[32][2])ps->Lbuf);
    ps_stereo_processing(ps);
    ps_hybrid_synthesis_all(ps, Lre, Lim, Rre, Rim);
}
