﻿# aacdec

一个从零实现的 MPEG-2 / MPEG-4 **AAC-LC / HE-AAC / AAC-LD / AAC-ELD**
音频解码器。接口模式参考
[`mpeg2dec`](../mpeg2dec) / `h264dec`：打开解码器 → 发送码流包 → 接收引用计数
的 PCM 帧 → flush → 关闭。

算法参考 ISO/IEC 13818-7（MPEG-2 AAC）和 ISO/IEC 14496-3（MPEG-4 AAC）标准，
并结合 Fraunhofer FDK AAC 与 FFmpeg 的解码流程进行独立实现；**没有直接复制
FDK 的解码器源码**。Huffman 码字、scalefactor band 偏移等数值表是标准中规定的
常量表。

---

## 已实现的功能

### 传输格式

| 格式 | 说明 |
| --- | --- |
| ADTS | `.aac` 同步流，支持 MPEG-2 / MPEG-4 ADTS，单/多 raw_data_block，头部 CRC |
| ADIF | 标准 ADIF 头 + PCE，随后是连续的 raw_data_block |
| LATM / LOAS | LOAS 同步流，带内 StreamMuxConfig / AudioSpecificConfig；AudioMuxVersion 0/1，多 sub-frame |
| RAW | 裸 access unit，配合 AudioSpecificConfig（ASC）或 StreamMuxConfig（SMC） |
| packetized LATM | `LATM_MCP0`（带外 SMC）与 `LATM_MCP1`（带内配置） |

`AAC_TRANSPORT_AUTO` 会自动探测 ADTS、ADIF 和 LOAS；提供带外配置时默认为 RAW。

### AAC-LC 工具

* MPEG-2 / MPEG-4 AAC-LC（AOT 2）、ER AAC-LC（AOT 17）、LFE 元素。
* 1024 与 960 samples/frame 两种帧长。
* ER AAC-LC 使用隐式元素类型与 `epConfig = 0`；ER protection/CRC 未实现。
* 通道配置 1–7（mono、stereo、3.0、4.0、5.0、5.1、7.1），以及
  `channelConfiguration = 0` 的 PCE 通道映射。
* 所有标准采样率：96k、88.2k、64k、48k、44.1k、32k、24k、22.05k、16k、
  12k、11.025k、8k、7.35k。
* ICS：`ONLY_LONG`、`LONG_START`、`EIGHT_SHORT`、`LONG_STOP`。
* sine 与 KBD 窗（长窗 alpha=4，短窗 alpha=6）。
* 光谱 Huffman 码本 1–11（含 escape）、scalefactor Huffman 码本。
* 逆量化 `|q|^(4/3) * 2^((sf-100)/4)`，短窗 scalefactor grouping。
* Pulse data、TNS（反射系数 → LPC → AR 滤波）、PNS（含 CPE 相关噪声）。
* CPE 的 M/S stereo、intensity stereo（码本 14/15）。
* IMDCT：DCT-IV 映射 + 预计算余弦表；窗、重叠相加。
* 输出 S16/F32、交错/平面四种 PCM 格式。
* 包时间戳（`pts`/`dts`/`opaque`）按 access unit 传递；支持任意大小输入包
  （ADTS/ADIF/LOAS 可 1 字节一包，RAW 要求一个完整 AU 一包）。
* MPEG-4 DRC（FIL extension type 11）：解析 `dynamic_range_info()`，并通过
  `drc_program_ref_level` 导出 program reference level。
* 错误隐藏：损坏的 access unit 默认由上一帧的淡出副本替代；调用方检测到
  丢包时也可以调用 `aac_dec_conceal_frame()` 生成替代帧，替代帧的
  `frame->error` 为 1。

### HE-AAC / SBR

* 从 `raw_data_block` 的 `FIL` 元素解析 `EXT_SBR_DATA` / `EXT_SBR_DATA_CRC`
  （类型 13/14），也支持 ASC 中显式 AOT 5 / 29。
* SBR 头、频率表（master / low / high / noise / limiter tables）、grid、
  envelope、noise floor、inverse filtering、chirp、HF patching、envelope
  adjustment 和 32→64 band QMF 分析/合成。
* 输出采样率为 `ext_sample_rate`（通常为 core 的 2 倍），帧长为 2048 / 1920。
* HE-AACv2（AOT 29）解析 PS 扩展，并使用 20-band baseline 或 34-band
  高分辨率 PS 合成输出立体声：IID/ICC、hybrid filterbank、all-pass
  decorrelator 和 mixing 均已启用。20-band PS 输出与 FDK 的逐样本结果
  不是 bit-exact（典型通道相关系数约 0.975，左右分离度约 25–30 dB，主要
  受 QMF 混叠限制）。
* 支持 34-band 高分辨率 hybrid 模式，以及 IPD/OPD 相位参数。
* 20-band 与 34-band 模式之间切换时会对混合系数做重映射。
* SBR 的 QMF 使用标准原型滤波器数值表但采用直接调制实现，因此与 FFmpeg /
  FDK 的逐样本结果不是 bit-exact（典型波形相关系数约 0.999，误差主要来自
  QMF 混叠）。

### AAC-LD / AAC-ELD

* AAC-LD（AOT 23）与 AAC-ELD（AOT 39），支持 512 / 480 samples/frame。
* 通过 LATM/LOAS、RAW + AudioSpecificConfig 或 packetized LATM 输入配置；
  ADTS / ADIF 不承载这两种对象类型。
* AAC-LD：128 样本低重叠窗与低延迟合成滤波器组，仅长帧；支持
  48 / 44.1 / 32 / 24 / 22.05 kHz。
* AAC-ELD：3N 低延迟窗与 3N 合成滤波器组，隐式元素标签、
  `common_window = 1`，仅长帧；支持 48 / 44.1 / 32 / 24 / 22.05 kHz。
* AAC-LD 的 ER 语法中 TNS 数据位于 gain control 标志之后；ELD 按 ELD
  语法解析 TNS 标志。
* 测试流与 FDK 参考 PCM 的稳态波形 SNR 约 88 dB（容差 4 LSB）。

### USAC (xHE-AAC, AOT 42)

* 解析 USAC `AudioSpecificConfig`：`usacSamplingFrequency`、
  `coreSbrFrameLengthIndex`（core 帧长 768/1024）、
  `channelConfigurationIndex`、`UsacDecoderConfig` 的元素列表
  （SCE / CPE / LFE / EXT）、`sbrConfig` 默认头参数以及
  `usacConfigExtension`（stream identifier、loudnessInfoSet）。
* 实现 ISO/IEC 23003-3 的上下文自适应算术解码器
  （`scalefactorband_data`）：标准 CDF 表、2-tuple 上下文映射、
  逃逸/MSB/LSB 分层解码与跨帧上下文记忆。解码结果与 FDK-AAC 参考
  逐符号、逐位一致（包括长逃逸串与 8 个短窗的上下文切换）。
* 实现 `UsacFrame`：`usacIndependencyFlag`、元素循环、
  `UsacExtElement`（filler / UniDrc / elevation / audio preroll 载荷跳过）、
  CPE 公共窗 / M-S 掩码 / 公共 TNS、`tns_data`、`noiseFillingData`、
  `facData`、USAC `sbrData` 与 2x SBR 上采样。
* 实现复杂立体声预测（`ms_mask_mode == 3`）：`cplx_pred_data` 的
  频域/时域 DPCM（`delta_code_time`）、7 抽头 MDST 插值、上下混合，
  以及跨帧的 alpha 系数与 MDST 历史。
* 测试流 `tests/data/usac_stereo.loas`（48 kHz 立体声，12 个 LOAS 帧，
  由 FDK xHE-AAC 编码器生成）可完整解码，无错误帧；输出与 FDK 参考
  解码器对齐后波形相关系数约 0.999。
* **限制**：768 帧长、8/3 与 4/3 SBR 比率、MPS212 立体声、
  LPD（ACELP/TCX）与时间扭曲尚未实现。

### 尚未实现

* AAC-LD / AAC-ELD 的 LD-SBR、ELDv2 / LD-MPS 扩展。
* ER protection（`epConfig > 0`）、ER CRC；AAC-LD / AAC-ELD 仅支持
  `epConfig = 0`。
* USAC 的 768 帧长、8/3 与 4/3 SBR 比率、MPS212、LPD、时间扭曲
  （见上节）。
* AAC Main / SSR / LTP、CCE（coupling channel element）。
* MPEG-4 DRC 仅解析并导出 program reference level，不应用增益；
  MPEG-D DRC / uniDrc、downmix、limiter、loudness metadata 未实现。
* 多声道 > 8、MPEG Surround / MPS。
* LATM 多 program、AudioMuxVersionA；length-delimited ASC 已支持。
* ADTS 多 RDB 的头部 CRC 会校验；单 RDB、每 RDB 的 CRC 字段和 LATM
  配置 CRC 目前只跳过，不做校验。
* AAC-LD / AAC-ELD 仅支持 48 / 44.1 / 32 / 24 / 22.05 kHz 和长帧；
  更低采样率（16 kHz 及以下）没有 SWB 表，无法解码有效帧。

---

## 构建

### Make（MinGW-w64 / MSYS2 / Linux / macOS）

```sh
mingw32-make -j8        # 或 make -j8
mingw32-make check
```

生成：

* `libaacdec.a`
* `aacdec_example.exe` / `aacdec_example`

### CMake

```sh
cmake -S . -B build-cmake -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release
cmake --build build-cmake -j8
ctest --test-dir build-cmake --output-on-failure
```

CMake 选项：`AACDEC_BUILD_EXAMPLE`、`AACDEC_BUILD_TESTS`、
`AACDEC_ENABLE_LTO`。

---

## 使用示例

```c
#include "aacdec.h"

aac_param_t param;
memset(&param, 0, sizeof(param));
param.transport     = AAC_TRANSPORT_AUTO; /* ADTS/ADIF/LOAS 自动探测 */
param.output_format = AAC_OUTPUT_S16;

void *dec = aac_dec_open(&param);

aac_packet_t pkt;
pkt.data = data;
pkt.len  = size;
pkt.pts  = pts;
aac_dec_send_packet(dec, &pkt);

aac_stream_info info;
aac_frame_t frame;
int ret;
while ((ret = aac_dec_recv_frame(dec, &info, &frame)) != AAC_DEC_DEFAULT) {
    if (ret == AAC_DEC_GOT_HEADER)
        continue;                       /* info 描述新的流配置 */
    if (ret == AAC_DEC_GOT_FRAME) {
        /* frame.planes[0..num_planes-1], frame.strides[],
           frame.frame_size, frame.channels, frame.sample_rate */
        aac_dec_frame_unref(dec, &frame);
        continue;
    }
    break;                              /* 错误 */
}

/* 流结束 */
for (;;) {
    ret = aac_dec_flush(dec, &info, &frame);
    if (ret == AAC_DEC_GOT_HEADER)
        continue;
    if (ret == AAC_DEC_GOT_FRAME) {
        aac_dec_frame_unref(dec, &frame);
        continue;
    }
    if (ret == AAC_DEC_END || ret == AAC_DEC_DEFAULT)
        break;
}

aac_dec_close(dec);
```

### 输出布局

* `AAC_OUTPUT_S16`（默认）：`num_planes = 1`，`planes[0]` 为交错
  `int16_t`，`strides[0] = channels * 2`。
* `AAC_OUTPUT_S16_PLANAR`：`num_planes = channels`，每个 plane 是一个声道的
  `int16_t`，`strides[c] = 2`。
* `AAC_OUTPUT_F32` / `AAC_OUTPUT_F32_PLANAR`：同样的布局，样本为 `float`，
  满刻度 ±1.0。

`channel_order` 默认为 `AAC_CHANNEL_ORDER_WAV`（5.1 顺序为
L R C LFE Ls Rs）；设为 `AAC_CHANNEL_ORDER_MPEG` 时输出规范 MPEG 顺序
（5.1 为 C L R Ls Rs LFE）。

### 错误隐藏与 DRC

损坏的 access unit 默认会被上一帧的淡出副本替代，`frame->error` 为 1；
连续丢失最多 4 帧，之后输出静音，收到正常帧后自动恢复。若传输层能检测到
丢包，可以在 `recv` 返回 `AAC_DEC_DEFAULT` 后主动插入替代帧：

```c
int ret = aac_dec_conceal_frame(dec, &info, &frame);
if (ret == AAC_DEC_GOT_FRAME) {
    /* frame.error == 1 */
    aac_dec_frame_unref(dec, &frame);
}
```

关闭自动隐藏：

```c
param.flags = AAC_DEC_FLAG_DISABLE_CONCEALMENT;
```

MPEG-4 DRC 的解析结果通过 `info.drc_program_ref_level` 导出，没有 DRC 时为
-1。当前只解析和导出，不把 DRC 增益应用到 PCM。

### 带外配置

RAW / packetized LATM 需要配置：

```c
static const uint8_t asc[] = { 0x12, 0x10 };  /* 44.1 kHz, stereo, AAC-LC */
aac_param_t p = { 0 };
p.transport = AAC_TRANSPORT_RAW;
p.config = asc;
p.config_len = sizeof(asc);
void *dec = aac_dec_open(&p);
```

也可以在开始解码前调用 `aac_dec_set_config()`。

### 命令行示例

```sh
aacdec_example input.aac output.wav
aacdec_example -t adts -s 1 input.aac output.wav     # 1 字节输入包
aacdec_example -t raw -c stereo.asc input.aac raw.pcm
aacdec_example -f f32 -o mpeg input.latm output.pcm
aacdec_example -h
```

`-f s16|s16p|f32|f32p` 选择输出格式，`-o wav|mpeg` 选择声道顺序，
`-s N` 选择输入块大小。S16 交错输出写 WAV（设置环境变量 `AACDEC_RAW=1`
可写裸 PCM）。

---

## 测试

`make check` / `ctest` 会运行：

| 测试 | 内容 |
| --- | --- |
| `api_check` | 用 1/7/4096 字节三类输入包解码，比较 PCM 字节完全一致；检查 header/frame/flush、S16/F32/planar 输出 |
| `compare_check` | 与 FFmpeg 9.0.2 的 AAC 解码结果逐样本比较，容差 32 LSB |
| `he_check` | 检查 HE-AAC 输出为 44100 Hz / 2048 samples / 2 ch，并验证 SBR 帧可解码 |
| `usac_check` | 解码 xHE-AAC（AOT 42）LOAS 流，检查 48000 Hz / 2048 samples / 2 ch，并验证所有 USAC 帧无错误 |
| `ps_check` | 检查 HE-AACv2 输出存在真实左右差异（非 dual-mono），验证 PS 混合生效 |
| `ps_unit` | 验证 34-band hybrid 透明性、34-band IID 分离度和 IPD/OPD 路径 |
| `ld_check` | 解码 AAC-LD LOAS 流，与 FDK 参考 PCM 逐样本比较（SNR ≥ 60 dB） |
| `eld_check` | 解码 AAC-ELD LOAS 流，与 FDK 参考 PCM 逐样本比较（SNR ≥ 60 dB） |
| `drc_check` | 合成 FIL/DRC 载荷，验证 `dynamic_range_info()` 解析与 stream info 导出 |
| `conceal_check` | 验证损坏帧/声明丢帧的隐藏、淡出与恢复，以及禁用标志 |
| `er_check` | 合成 ER AAC-LC（AOT 17）ASC 与 raw_data_block，验证 ER 路径 |
| `crc_check` | 破坏 ADTS 多 RDB 头 CRC，验证检测、隐藏与后续恢复 |
| LATM/CRC transport 测试 | 验证 AudioMuxVersion 1、多 sub-frame LOAS，以及单/多 RDB ADTS CRC 流的解码结果 |
| transport 测试 | 同一组 AAC access unit 的 ADTS、LATM/LOAS、ADIF、RAW+ASC 解码结果互相一致 |

`tests/data/` 带有很小的测试流：44.1 kHz stereo（ADTS/LATM/ADIF/RAW+ASC）、
48 kHz 5.1 ADTS、HE-AAC stereo 和 HE-AACv2 stereo、AAC-LD / AAC-ELD
LOAS 流、48 kHz stereo xHE-AAC（USAC）LOAS 流，以及用 FDK-AAC 3.0 生成的
AudioMuxVersion 1 / 双 sub-frame LOAS 和 ADTS CRC 流。可以用 FFmpeg
重新生成 ADTS/LATM：

```sh
ffmpeg -f lavfi -i "sine=frequency=440:sample_rate=44100:duration=0.2" \
       -ac 2 -c:a aac -b:a 96k -f adts tests/data/tiny_stereo.aac
ffmpeg -i tests/data/tiny_stereo.aac -f s16le tests/data/tiny_stereo.s16
```

在 44.1 kHz 立体声上，本解码器与 FFmpeg 的 AAC 解码器逐样本最大误差通常为
1–7 LSB（浮点舍入差异）。

---

## 实现结构

```
aacdec.h                公共 API
src/aacdec.c            传输解复用、LATM 多子帧、ADTS CRC、错误隐藏、PCM 格式转换、公共 API
src/aac_bits.h          位读取器
src/aac_tables.[ch]     采样率、SFB 表、窗表、2^(e/4)/|q|^(4/3) 表
src/aac_tables_data.h   标准 Huffman 码字/码长、SFB 偏移、TNS 系数表
src/aac_ld_tables.h     AAC-LD 512/480 SWB 与 TNS 表
src/aac_eld_tables.h    AAC-ELD 窗表
src/aac_eld_ld_tables.h AAC-ELD 低延迟合成窗表
src/aac_huffman.[ch]    运行时构建 Huffman 解码树
src/aac_core.[ch]       raw_data_block 语法、ICS、逆量化、PNS/MS/IS/TNS、IMDCT、DRC
src/aac_sbr.[ch]        HE-AAC SBR：FIL 扩展、频率表、QMF、HF 生成与包络调整
src/aac_ps.[ch]         HE-AACv2 PS 解析、20/34-band hybrid、decorrelator/mixing、IPD/OPD
src/aac_sbr_tables.h    标准 SBR QMF/Huffman/noise 数值表
src/aac_ps_tables.h     标准 PS Huffman 表与 hybrid 映射表
tests/                  api_check、compare_check、测试流
example_decoder.c       命令行示例
```

### 解码流程

1. 传输层从 ADTS/ADIF/LOAS/RAW 中提取 access unit。
2. `aac_core` 解析 `raw_data_block`：SCE/CPE/LFE/PCE/DSE/FIL/END。
3. ICS 解析：窗序列、section data、scalefactor、pulse、TNS、spectral data。
4. Huffman 解码 → 有符号量化系数；escape 序列在符号位之后读取。
5. 逆量化与 scalefactor 缩放。
6. PNS → M/S → intensity stereo → TNS。
7. IMDCT（DCT-IV 映射）→ 加窗 → 与上一帧 overlap 相加。AAC-LD 使用
   低重叠窗，AAC-ELD 使用 3N 低延迟合成滤波器组（不经过常规 overlap-add）。
8. 若存在 SBR：32-band QMF 分析 → HF 生成/包络调整 → 64-band QMF 合成，
   输出 2 倍采样率；无 SBR 的元素仍做 2 倍上采样。
9. 按 WAV 或 MPEG 顺序重排输出声道，转换成请求的 PCM 格式。

### 性能

当前 IMDCT 使用直接 DCT-IV（O(N²)，N=1024/960；短窗 N=128/120），
并通过预计算余弦表避免三角函数调用。SBR 的 QMF 也使用预计算调制表的直接
实现（32×320 分析、64×640 合成）。在开发机上 AAC-LC 44.1 kHz 立体声约
10 倍实时；HE-AAC 因 QMF 计算量更大，在开发机上仍在 10 倍实时以上。
后续可替换为 FFT 加速的 IMDCT/QMF。

---

## 许可

本项目源码采用 MIT License（见 `LICENSE`）。Huffman 码字、scalefactor
band 偏移、窗/TNS 系数等数值表来自 ISO/IEC 13818-7 与 ISO/IEC 14496-3
标准；它们是标准规定的常量，不是从 FDK 解码器复制的实现代码。
