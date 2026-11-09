/*
 * Minimal AAC decoder demo.
 *
 *   aacdec_example [options] input output
 *
 *   -t, --transport auto|adts|adif|latm|raw|mcp0|mcp1
 *   -c, --config FILE       AudioSpecificConfig / StreamMuxConfig
 *   -f, --format s16|s16p|f32|f32p
 *   -o, --order wav|mpeg
 *   -s, --chunk N           input chunk size in bytes (default 65536)
 *
 * Interleaved S16 output is written as a WAV file.  The other PCM formats are
 * written as raw little-endian data.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "aacdec.h"

typedef struct
{
    FILE* file;
    int is_wav;
    int channels;
    int sample_rate;
    int bits;
    uint32_t data_bytes;
    int header_written;
} output_file;

static void put16(FILE* f, unsigned v)
{
    fputc((int)(v & 0xFF), f);
    fputc((int)((v >> 8) & 0xFF), f);
}

static void put32(FILE* f, uint32_t v)
{
    put16(f, v & 0xFFFF);
    put16(f, (v >> 16) & 0xFFFF);
}

static int write_wav_header(output_file* o)
{
    long pos;

    if (!o->is_wav || o->header_written)
        return 0;

    pos = ftell(o->file);
    if (pos < 0 || fseek(o->file, 0, SEEK_SET) != 0)
        return -1;

    fwrite("RIFF", 1, 4, o->file);
    put32(o->file, 36);
    fwrite("WAVE", 1, 4, o->file);
    fwrite("fmt ", 1, 4, o->file);
    put32(o->file, 16);
    put16(o->file, 1);                    /* PCM */
    put16(o->file, (unsigned)o->channels);
    put32(o->file, (uint32_t)o->sample_rate);
    put32(o->file, (uint32_t)(o->sample_rate * o->channels * o->bits / 8));
    put16(o->file, (unsigned)(o->channels * o->bits / 8));
    put16(o->file, (unsigned)o->bits);
    fwrite("data", 1, 4, o->file);
    put32(o->file, 0);

    o->header_written = 1;
    return (fseek(o->file, pos, SEEK_SET) == 0) ? 0 : -1;
}

static void finish_wav(output_file* o)
{
    if (o->is_wav && o->header_written) {
        if (fseek(o->file, 4, SEEK_SET) == 0)
            put32(o->file, 36 + o->data_bytes);
        if (fseek(o->file, 40, SEEK_SET) == 0)
            put32(o->file, o->data_bytes);
        fflush(o->file);
    }
}

static int write_frame(output_file* o, const aac_frame_t* frame)
{
    int p;

    if (frame->format == AAC_OUTPUT_S16 && !frame->is_planar) {
        size_t bytes = (size_t)frame->frame_size *
                       (size_t)frame->channels * 2;
        if (o->is_wav) {
            o->channels = frame->channels;
            o->sample_rate = frame->sample_rate;
            o->bits = 16;
            if (write_wav_header(o) < 0)
                return -1;
        }
        if (fwrite(frame->planes[0], 1, bytes, o->file) != bytes)
            return -1;
        o->data_bytes += (uint32_t)bytes;
        return 0;
    }

    for (p = 0; p < frame->num_planes; p++) {
        size_t bytes = (size_t)frame->frame_size *
                       (size_t)frame->bytes_per_sample;
        if (fwrite(frame->planes[p], 1, bytes, o->file) != bytes)
            return -1;
        o->data_bytes += (uint32_t)bytes;
    }
    return 0;
}

static int parse_transport(const char* s)
{
    if (s == NULL || strcmp(s, "auto") == 0)
        return AAC_TRANSPORT_AUTO;
    if (strcmp(s, "adts") == 0)
        return AAC_TRANSPORT_ADTS;
    if (strcmp(s, "adif") == 0)
        return AAC_TRANSPORT_ADIF;
    if (strcmp(s, "latm") == 0 || strcmp(s, "loas") == 0)
        return AAC_TRANSPORT_LATM;
    if (strcmp(s, "raw") == 0)
        return AAC_TRANSPORT_RAW;
    if (strcmp(s, "mcp0") == 0)
        return AAC_TRANSPORT_LATM_MCP0;
    if (strcmp(s, "mcp1") == 0)
        return AAC_TRANSPORT_LATM_MCP1;
    return -1;
}

static int parse_format(const char* s)
{
    if (s == NULL || strcmp(s, "s16") == 0)
        return AAC_OUTPUT_S16;
    if (strcmp(s, "s16p") == 0)
        return AAC_OUTPUT_S16_PLANAR;
    if (strcmp(s, "f32") == 0)
        return AAC_OUTPUT_F32;
    if (strcmp(s, "f32p") == 0)
        return AAC_OUTPUT_F32_PLANAR;
    return -1;
}

static uint8_t* read_file(const char* name, int* size)
{
    FILE* f;
    long n;
    uint8_t* data;

    *size = 0;
    f = fopen(name, "rb");
    if (f == NULL)
        return NULL;
    if (fseek(f, 0, SEEK_END) != 0) {
        fclose(f);
        return NULL;
    }
    n = ftell(f);
    if (n < 0 || n > 0x7FFFFFFF) {
        fclose(f);
        return NULL;
    }
    rewind(f);
    data = (uint8_t*)malloc((size_t)n);
    if (data == NULL) {
        fclose(f);
        return NULL;
    }
    if (n > 0 && fread(data, 1, (size_t)n, f) != (size_t)n) {
        free(data);
        fclose(f);
        return NULL;
    }
    fclose(f);
    *size = (int)n;
    return data;
}

static void print_info(const aac_stream_info* info)
{
    int i;
    fprintf(stderr,
            "stream: %d Hz, %d samples/frame, %d channel(s), AOT %d",
            info->sample_rate, info->frame_size, info->channels, info->aot);
    if (info->ext_aot != 0)
        fprintf(stderr, " (extension AOT %d @ %d Hz)", info->ext_aot,
                info->ext_sample_rate);
    fprintf(stderr, ", config %d", info->channel_config);
    if (info->bitrate > 0)
        fprintf(stderr, ", %d bit/s", info->bitrate);
    fprintf(stderr, "\n");
    if (getenv("AACDEC_PRINT_CHANNELS")) {
        fprintf(stderr, "  channel map:");
        for (i = 0; i < info->channels && i < AAC_DEC_MAX_CHANNELS; i++)
            fprintf(stderr, " %02x.%d", info->channel_types[i],
                    info->channel_indices[i]);
        fprintf(stderr, "\n");
    }
}

static void usage(const char* argv0)
{
    fprintf(stderr,
            "usage: %s [options] input output\n"
            "  -t, --transport auto|adts|adif|latm|raw|mcp0|mcp1\n"
            "  -c, --config FILE       ASC/SMC for raw and packetized LATM\n"
            "  -f, --format s16|s16p|f32|f32p\n"
            "  -o, --order wav|mpeg\n"
            "  -s, --chunk N           input chunk size (default 65536)\n",
            argv0);
}

int main(int argc, char** argv)
{
    const char* input = NULL;
    const char* output = NULL;
    const char* config_file = NULL;
    int transport = AAC_TRANSPORT_AUTO;
    int format = AAC_OUTPUT_S16;
    int order = AAC_CHANNEL_ORDER_WAV;
    size_t chunk = 65536;
    uint8_t* config = NULL;
    int config_len = 0;
    FILE* in = NULL;
    output_file out;
    aac_param_t param;
    aac_packet_t pkt;
    aac_stream_info info;
    aac_frame_t frame;
    void* dec = NULL;
    uint8_t* buf = NULL;
    int64_t pts = 0;
    int frames = 0;
    int headers = 0;
    int eof = 0;
    int ret;
    int i;

    for (i = 1; i < argc; i++) {
        const char* a = argv[i];
        if (a[0] == '-' && a[1] != '\0') {
            const char* val = NULL;
            if (a[1] == '-') {
                const char* name = a + 2;
                const char* eq = strchr(name, '=');
                size_t n = eq ? (size_t)(eq - name) : strlen(name);
                if (strncmp(name, "transport", n) == 0 && n >= 3) {
                    if (eq)
                        val = eq + 1;
                    else if (i + 1 < argc)
                        val = argv[++i];
                    transport = parse_transport(val);
                } else if (strncmp(name, "config", n) == 0 && n >= 3) {
                    if (eq)
                        val = eq + 1;
                    else if (i + 1 < argc)
                        val = argv[++i];
                    config_file = val;
                } else if (strncmp(name, "format", n) == 0 && n >= 2) {
                    if (eq)
                        val = eq + 1;
                    else if (i + 1 < argc)
                        val = argv[++i];
                    format = parse_format(val);
                } else if (strncmp(name, "order", n) == 0 && n >= 2) {
                    if (eq)
                        val = eq + 1;
                    else if (i + 1 < argc)
                        val = argv[++i];
                    if (val && strcmp(val, "mpeg") == 0)
                        order = AAC_CHANNEL_ORDER_MPEG;
                    else
                        order = AAC_CHANNEL_ORDER_WAV;
                } else if (strncmp(name, "chunk", n) == 0 && n >= 2) {
                    if (eq)
                        val = eq + 1;
                    else if (i + 1 < argc)
                        val = argv[++i];
                    if (val && atoi(val) > 0)
                        chunk = (size_t)atoi(val);
                } else if (strcmp(name, "help") == 0) {
                    usage(argv[0]);
                    return 2;
                } else {
                    usage(argv[0]);
                    return 2;
                }
            } else {
                switch (a[1]) {
                case 't':
                    if (i + 1 < argc)
                        transport = parse_transport(argv[++i]);
                    break;
                case 'c':
                    if (i + 1 < argc)
                        config_file = argv[++i];
                    break;
                case 'f':
                    if (i + 1 < argc)
                        format = parse_format(argv[++i]);
                    break;
                case 'o':
                    if (i + 1 < argc && strcmp(argv[i + 1], "mpeg") == 0) {
                        order = AAC_CHANNEL_ORDER_MPEG;
                        i++;
                    } else if (i + 1 < argc) {
                        order = AAC_CHANNEL_ORDER_WAV;
                        i++;
                    }
                    break;
                case 's':
                    if (i + 1 < argc && atoi(argv[i + 1]) > 0) {
                        chunk = (size_t)atoi(argv[++i]);
                    }
                    break;
                case 'h':
                    usage(argv[0]);
                    return 2;
                default:
                    usage(argv[0]);
                    return 2;
                }
            }
            if (transport < 0 || format < 0) {
                fprintf(stderr, "invalid transport or format\n");
                return 2;
            }
        } else if (input == NULL) {
            input = a;
        } else if (output == NULL) {
            output = a;
        } else {
            usage(argv[0]);
            return 2;
        }
    }

    if (input == NULL || output == NULL) {
        usage(argv[0]);
        return 2;
    }

    if (config_file != NULL) {
        config = read_file(config_file, &config_len);
        if (config == NULL || config_len <= 0) {
            fprintf(stderr, "cannot read config file %s\n", config_file);
            return 2;
        }
    }

    in = fopen(input, "rb");
    if (in == NULL) {
        perror(input);
        free(config);
        return 2;
    }

    memset(&out, 0, sizeof(out));
    out.is_wav = (format == AAC_OUTPUT_S16 && !getenv("AACDEC_RAW"));
    out.file = fopen(output, "wb");
    if (out.file == NULL) {
        perror(output);
        fclose(in);
        free(config);
        return 2;
    }
    if (out.is_wav) {
        out.bits = 16;
        if (fwrite("\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
                   "\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"
                   "\0\0\0\0\0\0\0\0\0\0\0\0", 1, 44, out.file) != 44) {
            fprintf(stderr, "write failed\n");
            fclose(out.file);
            fclose(in);
            free(config);
            return 2;
        }
    }

    memset(&param, 0, sizeof(param));
    param.transport = transport;
    param.output_format = format;
    param.config = config;
    param.config_len = config_len;
    param.channel_order = order;

    dec = aac_dec_open(&param);
    if (dec == NULL) {
        fprintf(stderr, "aac_dec_open failed\n");
        fclose(out.file);
        fclose(in);
        free(config);
        return 2;
    }

    buf = (uint8_t*)malloc(chunk);
    if (buf == NULL) {
        fprintf(stderr, "out of memory\n");
        aac_dec_close(dec);
        fclose(out.file);
        fclose(in);
        free(config);
        return 2;
    }

    fprintf(stderr, "%s\n", aac_dec_version());

    while (!eof) {
        size_t n = fread(buf, 1, chunk, in);
        if (n == 0) {
            eof = 1;
        } else {
            memset(&pkt, 0, sizeof(pkt));
            pkt.data = buf;
            pkt.len = (int)n;
            pkt.pts = pts;
            pkt.dts = pts;
            pts += 1024;
            if (aac_dec_send_packet(dec, &pkt) < 0)
                fprintf(stderr, "send_packet: %s\n",
                        aac_dec_strerror(AAC_DEC_ERROR));
        }

        for (;;) {
            memset(&info, 0, sizeof(info));
            memset(&frame, 0, sizeof(frame));
            ret = aac_dec_recv_frame(dec, &info, &frame);
            if (ret == AAC_DEC_GOT_HEADER) {
                headers++;
                print_info(&info);
                continue;
            }
            if (ret == AAC_DEC_GOT_FRAME) {
                if (out.is_wav) {
                    out.channels = frame.channels;
                    out.sample_rate = frame.sample_rate;
                }
                if (write_frame(&out, &frame) < 0) {
                    fprintf(stderr, "write failed\n");
                    aac_dec_frame_unref(dec, &frame);
                    goto done;
                }
                frames++;
                aac_dec_frame_unref(dec, &frame);
                continue;
            }
            if (ret == AAC_DEC_DEFAULT)
                break;
            if (ret == AAC_DEC_END)
                break;
            fprintf(stderr, "decode: %s\n", aac_dec_strerror(ret));
            break;
        }
    }

    for (;;) {
        memset(&info, 0, sizeof(info));
        memset(&frame, 0, sizeof(frame));
        ret = aac_dec_flush(dec, &info, &frame);
        if (ret == AAC_DEC_GOT_HEADER) {
            headers++;
            print_info(&info);
            continue;
        }
        if (ret == AAC_DEC_GOT_FRAME) {
            if (out.is_wav) {
                out.channels = frame.channels;
                out.sample_rate = frame.sample_rate;
            }
            if (write_frame(&out, &frame) < 0) {
                fprintf(stderr, "write failed\n");
                aac_dec_frame_unref(dec, &frame);
                break;
            }
            frames++;
            aac_dec_frame_unref(dec, &frame);
            continue;
        }
        if (ret == AAC_DEC_END)
            break;
        if (ret < 0) {
            fprintf(stderr, "flush: %s\n", aac_dec_strerror(ret));
            break;
        }
    }

done:
    finish_wav(&out);
    aac_dec_close(dec);
    fclose(out.file);
    fclose(in);
    free(buf);
    free(config);
    fprintf(stderr, "decoded %d frames, %d configuration update(s)\n", frames,
            headers);
    return 0;
}
