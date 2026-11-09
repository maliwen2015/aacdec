/*
 * ps_unit - unit checks for the PS hybrid filterbank and 34-band mixing.
 *
 * The implementation file is included directly so the internal hybrid and
 * mixing stages can be exercised without a bitstream.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "aac_ps.c"

static void fill_input(float re[38][64], float im[38][64])
{
    int j, k;
    for (j = 0; j < 38; j++) {
        for (k = 0; k < 64; k++) {
            re[j][k] = (float)sin(j * 0.7 + k * 0.3) * 1000.0f;
            im[j][k] = (float)cos(j * 0.4 + k * 0.2) * 1000.0f;
        }
    }
}

int main(void)
{
    aac_ps_decoder ps;
    float re[38][64], im[38][64], rr[38][64], ri[38][64];
    float save_re[38][64], save_im[38][64];
    double sig = 0.0, err = 0.0, el = 0.0, er = 0.0;
    int j, k, bad = 0;

    memset(&ps, 0, sizeof(ps));
    if (ps_tables_init() < 0) {
        fprintf(stderr, "table init failed\n");
        return 1;
    }

    /* 34-band hybrid analysis/synthesis must reconstruct the QMF input. */
    ps.is34bands = 1;
    ps.start = 1;
    ps.num_env = 1;
    ps.num_env_old = 1;
    ps.border_position[0] = -1;
    ps.border_position[1] = 31;
    fill_input(re, im);
    memcpy(save_re, re, sizeof(save_re));
    memcpy(save_im, im, sizeof(save_im));
    ps_hybrid_analysis_all(&ps, re, im, ps.Lbuf);
    ps_hybrid_synthesis_all(&ps, re, im, rr, ri);
    for (j = 0; j < 32; j++) {
        for (k = 0; k < 64; k++) {
            double dr = re[j][k] - save_re[j][k];
            double di = im[j][k] - save_im[j][k];
            sig += save_re[j][k] * save_re[j][k] +
                   save_im[j][k] * save_im[j][k];
            err += dr * dr + di * di;
        }
    }
    if (10.0 * log10(sig / (err + 1e-30)) < 80.0) {
        fprintf(stderr, "34-band hybrid transparency failed: %g dB\n",
                10.0 * log10(sig / (err + 1e-30)));
        return 1;
    }

    /* 34-band IID mixing: after two frames the steady-state separation must
     * be strong. */
    memset(&ps, 0, sizeof(ps));
    ps.is34bands = 1;
    ps.start = 1;
    ps.num_env = 1;
    ps.num_env_old = 1;
    ps.border_position[0] = -1;
    ps.border_position[1] = 31;
    ps.nr_iid_par = 34;
    ps.nr_icc_par = 34;
    ps.iid_quant = 1;
    for (j = 0; j < 34; j++) {
        ps.iid_par[0][j] = 15;
        ps.icc_par[0][j] = 0;
    }
    for (j = 0; j < 34; j++) {
        ps.H11[0][0][j] = 1.0f;
        ps.H12[0][0][j] = 1.0f;
        ps.H11[0][1][j] = 1.0f;
        ps.H12[0][1][j] = 1.0f;
    }
    fill_input(re, im);
    aac_ps_apply(&ps, re, im, rr, ri, 64);
    aac_ps_apply(&ps, re, im, rr, ri, 64);
    for (j = 0; j < 32; j++) {
        for (k = 0; k < 64; k++) {
            el += re[j][k] * re[j][k] + im[j][k] * im[j][k];
            er += rr[j][k] * rr[j][k] + ri[j][k] * ri[j][k];
            if (!isfinite(re[j][k]) || !isfinite(rr[j][k]) ||
                !isfinite(ri[j][k]))
                bad = 1;
        }
    }
    if (bad || 10.0 * log10(el / (er + 1e-30)) < 40.0) {
        fprintf(stderr, "34-band IID separation failed: %g dB finite=%d\n",
                10.0 * log10(el / (er + 1e-30)), !bad);
        return 1;
    }

    /* IPD/OPD path should stay finite. */
    memset(&ps, 0, sizeof(ps));
    ps.is34bands = 1;
    ps.start = 1;
    ps.num_env = 1;
    ps.num_env_old = 1;
    ps.border_position[0] = -1;
    ps.border_position[1] = 31;
    ps.nr_iid_par = 34;
    ps.nr_icc_par = 34;
    ps.iid_quant = 1;
    ps.enable_ipdopd = 1;
    ps.nr_ipdopd_par = 17;
    for (j = 0; j < 34; j++) {
        ps.iid_par[0][j] = 3;
        ps.icc_par[0][j] = 2;
    }
    for (j = 0; j < 17; j++) {
        ps.ipd_par[0][j] = (int8_t)(j & 7);
        ps.opd_par[0][j] = (int8_t)((j + 3) & 7);
    }
    fill_input(re, im);
    aac_ps_apply(&ps, re, im, rr, ri, 64);
    aac_ps_apply(&ps, re, im, rr, ri, 64);
    for (j = 0; j < 32; j++) {
        for (k = 0; k < 64; k++) {
            if (!isfinite(re[j][k]) || !isfinite(rr[j][k]) ||
                !isfinite(ri[j][k])) {
                fprintf(stderr, "IPD/OPD produced non-finite data\n");
                return 1;
            }
        }
    }

    printf("ok: PS 34-band hybrid, IID mixing and IPD/OPD\n");
    return 0;
}
