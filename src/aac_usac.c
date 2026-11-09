/*
 * aacdec - USAC (xHE-AAC) AudioSpecificConfig parser.
 *
 * Syntax follows ISO/IEC 23003-3 (USAC): AudioSpecificConfig with
 * audioObjectType 42 carries UsacConfig() directly; see aac_usac.h.
 */
#include "aac_usac.h"

#include <string.h>

#define USAC_OK        0
#define USAC_NEED_MORE (-1)
#define USAC_ERROR     (-2)

#ifdef AAC_USAC_DEBUG
#include <stdio.h>
#define DBG(...) do { fprintf(stderr, __VA_ARGS__); } while (0)
#else
#define DBG(...)
#endif

/* ISO/IEC 23003-3, usacSamplingFrequency table (32 entries). */
static const int aac_usac_sample_rates[32] = {
    96000, 88200, 64000, 48000, 44100, 32000, 24000, 22050,
    16000, 12000, 11025,  8000,  7350,    -1,    -1, 57600,
    51200, 40000, 38400, 34150, 28800, 25600, 20000, 19200,
    17075, 14400, 12800,  9600,    -1,    -1,    -1,    -1,
};

/* EscapedValue(nValueBits, nEscapeValueBits, nEscValBits), 23003-3 1.6. */
static uint32_t read_escaped_value(aac_bits* b, int nb1, int nb2, int nb3)
{
    uint32_t val, val2;

    val = aac_bits_read(b, nb1);
    if (b->error)
        return 0;
    if (val < (uint32_t)((1 << nb1) - 1))
        return val;
    val2 = aac_bits_read(b, nb2);
    val += val2;
    if (b->error)
        return 0;
    if (nb3 && val2 == (uint32_t)((1 << nb2) - 1))
        val += aac_bits_read(b, nb3);
    return val;
}

/* loudnessInfo(), ISO/IEC 23003-4 Table 58.  The measurements are parsed but
 * not exported (loudness normalization is not applied). */
static int parse_loudness_info(aac_bits* b, int v1)
{
    int nb_meas, i;

    (void)aac_bits_read(b, 6);              /* bsDrcSetId */
    if (v1)
        (void)aac_bits_read(b, 6);          /* bsEqSetId */
    (void)aac_bits_read(b, 7);              /* bsDownmixId */
    if (aac_bits_read(b, 1))                /* bsSamplePeakLevelPresent */
        (void)aac_bits_read(b, 12);
    if (aac_bits_read(b, 1)) {              /* bsTruePeakLevelPresent */
        (void)aac_bits_read(b, 12);
        (void)aac_bits_read(b, 4);
        (void)aac_bits_read(b, 2);
    }
    nb_meas = (int)aac_bits_read(b, 4);     /* bsNumMeasurements */
    for (i = 0; i < nb_meas; i++) {
        int method = (int)aac_bits_read(b, 4);
        int width = (method == 7) ? 5 : (method == 8) ? 2 : 8;
        (void)aac_bits_read(b, width);      /* bsMethodValue */
        (void)aac_bits_read(b, 4);          /* bsMeasurementStatus */
        (void)aac_bits_read(b, 2);          /* bsMeasurementReliability */
    }
    return b->error ? USAC_ERROR : USAC_OK;
}

/* loudnessInfoSet(), ISO/IEC 23003-4 Table 57. */
static int parse_loudness_set(aac_bits* b)
{
    int nb_album, nb_info, i;

    nb_album = (int)aac_bits_read(b, 6);
    nb_info = (int)aac_bits_read(b, 6);
    for (i = 0; i < nb_album; i++)
        if (parse_loudness_info(b, 0) != USAC_OK)
            return USAC_ERROR;
    for (i = 0; i < nb_info; i++)
        if (parse_loudness_info(b, 0) != USAC_OK)
            return USAC_ERROR;
    if (b->error)
        return USAC_ERROR;

    if (aac_bits_read(b, 1)) {              /* loudnessInfoSetExtPresent */
        for (;;) {
            int type = (int)aac_bits_read(b, 4);
            int size_bits;
            uint32_t bit_size;
            size_t start;

            if (b->error || type == 0)      /* UNIDRCLOUDEXT_TERM */
                break;
            size_bits = (int)aac_bits_read(b, 4) + 4;
            bit_size = aac_bits_read(b, size_bits) + 1;
            if (b->error)
                return USAC_ERROR;
            start = b->bitpos;
            if (type == 1) {                /* UNIDRCLOUDEXT_EQ */
                int n;
                nb_album = (int)aac_bits_read(b, 6);
                nb_info = (int)aac_bits_read(b, 6);
                for (n = 0; n < nb_album; n++)
                    if (parse_loudness_info(b, 1) != USAC_OK)
                        return USAC_ERROR;
                for (n = 0; n < nb_info; n++)
                    if (parse_loudness_info(b, 1) != USAC_OK)
                        return USAC_ERROR;
            }
            if (b->error)
                return USAC_ERROR;
            if (b->bitpos + (bit_size - (b->bitpos - start)) > b->size * 8u ||
                bit_size < (b->bitpos - start))
                return USAC_ERROR;
            b->bitpos = start + bit_size;
        }
    }
    DBG("usac: loudness set consumed %zu bits\n", b->bitpos - start);
    return b->error ? USAC_ERROR : USAC_OK;
}

static int parse_sbr_config(aac_bits* b, aac_usac_elem_config* e)
{
    int harmonic, intertes, pvc;
    int extra1, extra2;

    harmonic = (int)aac_bits_read(b, 1);
    intertes = (int)aac_bits_read(b, 1);
    pvc = (int)aac_bits_read(b, 1);
    if (harmonic || intertes || pvc) {
        /* eSBR tools are not supported. */
        return USAC_ERROR;
    }
    e->sbr_start_freq = (int)aac_bits_read(b, 4);
    e->sbr_stop_freq = (int)aac_bits_read(b, 4);
    DBG("usac: sbr cfg harmonic=%d intertes=%d pvc=%d start=%d stop=%d bit=%zu\n",
        harmonic, intertes, pvc, e->sbr_start_freq, e->sbr_stop_freq,
        b->bitpos);
    extra1 = (int)aac_bits_read(b, 1);
    extra2 = (int)aac_bits_read(b, 1);
    e->sbr_freq_scale = 2;
    e->sbr_alter_scale = 1;
    e->sbr_noise_bands = 2;
    if (extra1) {
        e->sbr_freq_scale = (int)aac_bits_read(b, 2);
        e->sbr_alter_scale = (int)aac_bits_read(b, 1);
        e->sbr_noise_bands = (int)aac_bits_read(b, 2);
    }
    e->sbr_limiter_bands = 2;
    e->sbr_limiter_gains = 2;
    e->sbr_interpol_freq = 1;
    e->sbr_smoothing_mode = 1;
    if (extra2) {
        e->sbr_limiter_bands = (int)aac_bits_read(b, 2);
        e->sbr_limiter_gains = (int)aac_bits_read(b, 2);
        e->sbr_interpol_freq = (int)aac_bits_read(b, 1);
        e->sbr_smoothing_mode = (int)aac_bits_read(b, 1);
    }
    return b->error ? USAC_NEED_MORE : USAC_OK;
}

static int parse_element_config(aac_bits* b, aac_usac_elem_config* e,
                                int sbr_ratio)
{
    e->tw_mdct = 0;
    e->noise_fill = 0;
    e->stereo_config_index = 0;
    e->sbr_ratio = 0;
    e->ext_type = -1;
    e->ext_default_len = 0;
    e->ext_payload_frag = 0;
    e->sbr_start_freq = 0;
    e->sbr_stop_freq = 0;
    e->sbr_freq_scale = 2;
    e->sbr_alter_scale = 1;
    e->sbr_noise_bands = 2;
    e->sbr_limiter_bands = 2;
    e->sbr_limiter_gains = 2;
    e->sbr_interpol_freq = 1;
    e->sbr_smoothing_mode = 1;

    switch (e->type) {
    case AAC_USAC_ID_SCE:
    case AAC_USAC_ID_LFE:
        e->tw_mdct = (int)aac_bits_read(b, 1);
        e->noise_fill = (int)aac_bits_read(b, 1);
        e->sbr_ratio = sbr_ratio;
        DBG("usac: SCE/LFE tw=%d nf=%d\n", e->tw_mdct, e->noise_fill);
        if (sbr_ratio > 0 && parse_sbr_config(b, e) < 0)
            return USAC_NEED_MORE;
        if (e->tw_mdct) {
            DBG("usac: tw_mdct unsupported\n");
            return USAC_ERROR;
        }
        break;

    case AAC_USAC_ID_CPE:
        e->tw_mdct = (int)aac_bits_read(b, 1);
        e->noise_fill = (int)aac_bits_read(b, 1);
        e->sbr_ratio = sbr_ratio;
        DBG("usac: CPE tw=%d nf=%d\n", e->tw_mdct, e->noise_fill);
        if (sbr_ratio > 0 && parse_sbr_config(b, e) < 0)
            return USAC_NEED_MORE;
        if (e->tw_mdct) {
            DBG("usac: tw_mdct unsupported (CPE)\n");
            return USAC_ERROR;
        }
        if (sbr_ratio > 0)
            e->stereo_config_index = (int)aac_bits_read(b, 2);
        DBG("usac: CPE stereo_config_index=%d\n", e->stereo_config_index);
        /* MPS (stereo_config_index != 0) is not supported. */
        if (e->stereo_config_index) {
            DBG("usac: MPS unsupported\n");
            return USAC_ERROR;
        }
        break;

    case AAC_USAC_ID_EXT: {
        int ext_len;
        e->ext_type = (int)read_escaped_value(b, 4, 8, 16);
        ext_len = (int)read_escaped_value(b, 4, 8, 16);
        if (b->error)
            return USAC_NEED_MORE;
        DBG("usac: EXT type=%d len=%d\n", e->ext_type, ext_len);
        if (aac_bits_read(b, 1)) /* usacExtElementDefaultLengthPresent */
            e->ext_default_len =
                (int)read_escaped_value(b, 8, 16, 0) + 1;
        e->ext_payload_frag = (int)aac_bits_read(b, 1);
        DBG("usac: EXT default_len=%d frag=%d\n", e->ext_default_len,
            e->ext_payload_frag);
        /* The config payload is byte aligned to ext_config_len bytes. */
        if (ext_len > 0) {
            size_t start = b->bitpos;
            /* uniDrcConfig() parsing is not implemented; skip its payload
             * so the stream remains parseable. */
            b->bitpos = start + (size_t)ext_len * 8u;
            if (b->bitpos > b->size * 8u) {
                b->bitpos = b->size * 8u;
                return USAC_NEED_MORE;
            }
        }
        break;
    }

    default:
        return USAC_ERROR;
    }
    return USAC_OK;
}

/* UsacConfig() body; the caller has already positioned the stream past the
 * audioObjectType. */
static int parse_usac_config_body(aac_bits* b, aac_usac_config* c)
{
    int freq_idx;
    int sbr_ratio;
    int i;
    int sbr_present = 0;

    freq_idx = (int)aac_bits_read(b, 5);
    DBG("usac: bit=%zu freq_idx=%d\n", b->bitpos, freq_idx);
    if (b->error)
        return USAC_NEED_MORE;
    if (freq_idx == 0x1f) {
        c->sample_rate = (int)aac_bits_read(b, 24);
    } else {
        c->sample_rate = aac_usac_sample_rates[freq_idx];
    }
    if (c->sample_rate <= 0) {
        DBG("usac: bad sample rate idx %d\n", freq_idx);
        return USAC_ERROR;
    }

    c->core_sbr_frame_len_idx = (int)aac_bits_read(b, 3);
    DBG("usac: bit=%zu frame_len_idx=%d\n", b->bitpos,
        c->core_sbr_frame_len_idx);
    if (c->core_sbr_frame_len_idx > 5) {
        DBG("usac: bad frame_len_idx\n");
        return USAC_ERROR;
    }
    c->core_frame_len =
        (c->core_sbr_frame_len_idx == 0 || c->core_sbr_frame_len_idx == 2)
            ? 768 : 1024;

    switch (c->core_sbr_frame_len_idx) {
    case 2: sbr_ratio = 2; break;   /* 8/3 */
    case 3: sbr_ratio = 3; break;   /* 2x  */
    case 4: sbr_ratio = 1; break;   /* 4/3 */
    default: sbr_ratio = 0; break;
    }
    c->sbr_ratio = sbr_ratio;

    c->channel_config_index = (int)aac_bits_read(b, 5);
    DBG("usac: bit=%zu cc_idx=%d\n", b->bitpos, c->channel_config_index);
    if (b->error)
        return USAC_NEED_MORE;
    if (c->channel_config_index == 0) {
        /* UsacChannelConfig(): output channels with speaker positions. */
        int num_out = (int)read_escaped_value(b, 5, 8, 16);
        DBG("usac: custom ch num_out=%d\n", num_out);
        if (b->error || num_out <= 0 || num_out > 8)
            return USAC_ERROR;
        for (i = 0; i < num_out; i++)
            (void)aac_bits_read(b, 5);  /* bsOutputChannelPos */
        c->channels = num_out;
    } else if (c->channel_config_index <= 7) {
        static const int chan_of_conf[8] = { 0, 1, 2, 3, 4, 5, 6, 8 };
        c->channels = chan_of_conf[c->channel_config_index];
    } else {
        DBG("usac: unsupported channel_config_index %d\n",
            c->channel_config_index);
        /* 22.2 / 6.1(back) / 7.1(wide) etc. are not mapped yet. */
        return USAC_ERROR;
    }

    /* UsacDecoderConfig */
    c->nb_elements = (int)read_escaped_value(b, 4, 8, 16) + 1;
    DBG("usac: bit=%zu nb_elems=%d\n", b->bitpos, c->nb_elements);
    if (b->error || c->nb_elements > AAC_USAC_MAX_ELEMENTS)
        return USAC_ERROR;

    {
        int audio_channels = 0;
        int audio_elems = 0;
        for (i = 0; i < c->nb_elements; i++) {
            int ret;
            c->elems[i].type = (int)aac_bits_read(b, 2);
            DBG("usac: elem %d type=%d bit=%zu\n", i, c->elems[i].type,
                b->bitpos);
            ret = parse_element_config(b, &c->elems[i], sbr_ratio);
            if (ret < 0) {
                DBG("usac: elem %d config parse failed (%d)\n", i, ret);
                return ret;
            }
            if (c->elems[i].type != AAC_USAC_ID_EXT) {
                int ch = (c->elems[i].type == AAC_USAC_ID_CPE) ? 2 : 1;
                if (audio_channels + ch > c->channels)
                    return USAC_ERROR;
                audio_channels += ch;
                audio_elems++;
            }
            if (c->elems[i].sbr_ratio > 0 &&
                c->elems[i].type != AAC_USAC_ID_EXT)
                sbr_present = 1;
        }
        if (c->channel_config_index != 0 && audio_elems == 0)
            return USAC_ERROR;
    }
    c->sbr_present = sbr_present;

    if (aac_bits_read(b, 1)) { /* usacConfigExtensionPresent */
        int nb_ext = (int)read_escaped_value(b, 2, 4, 8) + 1;
        DBG("usac: config extensions=%d bit=%zu\n", nb_ext, b->bitpos);
        if (b->error)
            return USAC_NEED_MORE;
        for (i = 0; i < nb_ext; i++) {
            int type = (int)read_escaped_value(b, 4, 8, 16);
            int len = (int)read_escaped_value(b, 4, 8, 16);
            size_t payload_start = b->bitpos;
            DBG("usac: ext %d type=%d len=%d bit=%zu\n", i, type, len,
                b->bitpos);
            if (b->error)
                return USAC_NEED_MORE;
            if (type == 2) {
                /* loudnessInfoSet(): parsed per syntax; the payload is
                 * declared in whole bytes and padded to the byte boundary. */
                if (parse_loudness_set(b) != USAC_OK)
                    return USAC_ERROR;
                if ((size_t)len * 8u >= b->bitpos - payload_start)
                    b->bitpos = payload_start + (size_t)len * 8u;
                else
                    return USAC_ERROR;
            } else if (type == 7) {
                (void)aac_bits_read(b, 16); /* stream identifier */
                if ((size_t)len * 8u >= b->bitpos - payload_start)
                    b->bitpos = payload_start + (size_t)len * 8u;
            } else {
                /* Fill and unknown payloads are byte aligned. */
                b->bitpos += (size_t)len * 8u;
                if (b->bitpos > b->size * 8u) {
                    b->bitpos = b->size * 8u;
                    return USAC_NEED_MORE;
                }
            }
        }
    }

    return USAC_OK;
}

static int parse_usac_config(aac_bits* b, aac_usac_config* c)
{
    return parse_usac_config_body(b, c);
}

int aac_usac_parse_asc(const uint8_t* data, size_t len, aac_usac_asc* asc)
{
    aac_bits b;
    int aot;
    int ret;

    if (data == NULL || len < 2 || asc == NULL)
        return -1;
    memset(asc, 0, sizeof(*asc));

    aac_bits_init(&b, data, len);
    aot = (int)aac_bits_read(&b, 5);
    if (aot == 31)
        aot = 32 + (int)aac_bits_read(&b, 6);
    if (b.error || aot != 42)
        return -1;

    /* AudioSpecificConfig() common header (ISO/IEC 14496-3 1.A.2.3):
     * samplingFrequencyIndex(4) + channelConfiguration(4) precede the
     * UsacConfig() specific config for AOT 42. */
    (void)aac_bits_read(&b, 4); /* samplingFrequencyIndex */
    (void)aac_bits_read(&b, 4); /* channelConfiguration */
    if (b.error)
        return -1;

    ret = parse_usac_config(&b, &asc->config);
    asc->config_bits = b.bitpos;
    if (ret != USAC_OK)
        return -1;

    aac_bits_align(&b);
    asc->asc_end_bit = (int)b.bitpos;
    return 0;
}

int aac_usac_parse_asc_bits_state(aac_bits* b, aac_usac_asc* asc)
{
    int ret;

    if (b == NULL || asc == NULL)
        return -1;
    memset(asc, 0, sizeof(*asc));

    ret = parse_usac_config_body(b, &asc->config);
    asc->config_bits = b->bitpos;
    if (ret != USAC_OK)
        return -1;

    /* No byte alignment here: within LATM/LOAS StreamMuxConfig the
     * frameLengthType field follows the AudioSpecificConfig immediately
     * without alignment.  asc_end_bit is the exact bit position after the
     * UsacConfig(). */
    asc->asc_end_bit = (int)b->bitpos;
    return 0;
}
