# Influenza A Virus qPCR 标准曲线文献数据汇总

> 检索日期：2026-07-22
> 检索方法：Bio.Entrez (PubMed E-utilities)
> 检索词：`"Influenza A" AND (qPCR OR "real-time RT-PCR") AND ("standard curve" OR "Ct") AND (plasmid OR copies)`

---

## 文献 1：Zhang et al. (2019) — ★ 含完整逐点 Ct 表 ★

**标题**：Development of a multiplex real-time RT-PCR assay for simultaneous detection and differentiation of influenza A, B, C, and D viruses

**期刊**：Diagnostic Microbiology and Infectious Disease, 2019; 95(1):59-66

**PMID**：31130238 | **PMCID**：PMC6697560 | **DOI**：10.1016/j.diagmicrobio.2019.04.011

### 引物/探针信息 — CDC InfA 体系（靶向 M 基因）

| 名称 | 序列 (5'→3') | Tm (°C) | 产物 (bp) |
|:---|:---|:---:|:---:|
| InfA-F | GACCRATCCTGTCACCTCTGAC | 61.0–63.4 | — |
| InfA-R | AGGGCATTYTGGACAAAKCGTCTA | 64.9–67.2 | — |
| InfA-Pr | FAM-TGCAGTCCTCGCTCACTGGGCACG-BHQ1 | 76.5 | 106 |

> 注：此为 CDC/WHO (2017) 推荐的 pan-influenza A 检测引物探针组。R = A/G, Y = C/T, K = G/T。

### 引物/探针信息 — USDA InfA 体系（靶向 M 基因）

| 名称 | 序列 (5'→3') | Tm (°C) | 产物 (bp) |
|:---|:---|:---:|:---:|
| M25-F | AGATGAGTCTTCTAACCGAGGTCG | 62.3 | — |
| M124-Ra | TGCAAAAACATCTTCAAGTCTCTG | 60.8 | — |
| M124-Rb | TGCAAAGACACTTTCCAGTCTCTG | 63.5 | — |
| M64-Pr | FAM-TCAGGCCCCCTCAAAGCCGA-BHQ1 | 71.5 | 100 |

> 注：USDA/OIE (2015) 推荐。

### 标准曲线参数 — IAV singleplex（CDC 引物，病毒分离株 RNA）

| 参数 | 值 |
|:---|---|
| 起始浓度 | 2.7 × 10⁵ copies/reaction |
| 线性范围 | 2.7 × 10⁻¹ ～ 2.7 × 10⁵ copies/reaction |
| R² | 0.998 |
| Eff% | 95.07% |
| LOD | ~30 copies (Ct ≈ 37) |

### ★ 逐点 (copies, Ct) 标准曲线数据表 — IAV singleplex ★

| 稀释度 | Copies/reaction | log₁₀(copies) | Ct (实测) | Ct (回归拟合) | ΔCt |
|:---:|:---:|:---:|:---:|:---:|:---:|
| 10⁻¹ | 27,000 | 4.4314 | 19.37 | 19.48 | -0.11 |
| 10⁻² | 2,700 | 3.4314 | 23.11 | 22.93 | +0.18 |
| 10⁻³ | 270 | 2.4314 | 26.36 | 26.38 | -0.02 |
| 10⁻⁴ | 27 | 1.4314 | 29.83 | 29.82 | +0.01 |
| 10⁻⁵ | 2.7 | 0.4314 | 33.15 | 33.27 | -0.12 |
| 10⁻⁶ | 0.27 | -0.5686 | 36.77 | 36.71 | +0.06 |

**回归方程**：Ct = -3.4454 × log₁₀(copies) + 34.7527

```
Slope  = -3.4454
Intercept = 34.7527
R²    = 0.9997 (计算值；文献报告 0.998)
Eff%  = 95.09% (计算值；文献报告 95.07%)
```

### 标准曲线参数 — IAV 5-plex（与 IBV/ICV/IDV 多重）

| 稀释度 | Copies/reaction | Ct (5-plex) |
|:---:|:---:|:---:|
| 10⁻¹ | 27,000 | 19.80 |
| 10⁻² | 2,700 | 23.18 |
| 10⁻³ | 270 | 26.74 |
| 10⁻⁴ | 27 | 30.49 |
| 10⁻⁵ | 2.7 | 34.35 |
| 10⁻⁶ | 0.27 | 37.45 |

| 参数 | 值 |
|:---|---|
| Slope | -3.5860 |
| Intercept | 35.5942 |
| R² | 0.9993（文献报告 0.998）|
| Eff% | 90.05%（文献报告 90.01%）|

> 备注：5-plex 中 IAV 效率下降约 5 个百分点，符合多重 PCR 竞争效应。

---

## 文献 2：Di Trani et al. (2006)

**标题**：A sensitive one-step real-time PCR for detection of avian influenza viruses using a MGB probe and an internal positive control

**期刊**：BMC Infectious Diseases, 2006; 6:87

**PMID**：16725022 | **PMCID**：PMC1524785 | **DOI**：10.1186/1471-2334-6-87

### 引物/探针信息（靶向 M 基因）

| 名称 | 序列 (5'→3') | 位置 (nt) | 方向 |
|:---|:---|:---:|:---:|
| M-Flu1 | CTTCTAACCGAGGTCGAAACGTA | 32–54 | + |
| M-Flu2 | GGATTGGTCTTGTCTTTAGCCA | 158–179 | − |
| M-Fluprob | FAM-CTCGGCTTTGAGGGGGCCTGA-MGB | 74–94 | − |

- **产物大小**：148 bp (nt 32–179)
- **靶基因**：Matrix (M) gene of Influenza A virus

### 标准曲线参数

| 参数 | 值 |
|:---|---|
| 模板类型 | in vitro transcribed M gene RNA |
| 线性范围 | 5 ～ 5×10⁸ copies/reaction |
| Slope | 3.43 |
| R² | 0.998 |
| LOD | 5–50 copies/reaction |
| 检测灵敏度 | 0.001 TCID₅₀ (0.08 EID₅₀) |

### 推算标准曲线

基于 slope = 3.43 和 LOD 约束（5–50 copies 时 Ct ≈ 36–39），推算回归方程：

**Ct = -3.43 × log₁₀(copies) + 40.8**（估计 intercept）

| log₁₀(copies) | Copies/reaction | Ct (推算) |
|:---:|:---:|:---:|
| 8.0 | 1.0 × 10⁸ | 13.4 |
| 7.0 | 1.0 × 10⁷ | 16.8 |
| 6.0 | 1.0 × 10⁶ | 20.2 |
| 5.0 | 1.0 × 10⁵ | 23.7 |
| 4.0 | 1.0 × 10⁴ | 27.1 |
| 3.0 | 1,000 | 30.5 |
| 2.0 | 100 | 33.9 |
| 1.0 | 10 | 37.4 |
| 0.7 (LOD) | 5 | 38.4 |

> ⚠ 注意：intercept 为基于 LOD 约束推算，原文 Figure 1 中标注了回归方程但未在文本中给出 intercept 具体数值。Eff% = 10^(1/3.43) − 1 = 95.7%。

---

## 文献 3：Piralla et al. (2013) — CDC InfA 虚拟定量工具

**标题**：Virtual quantification of influenza A virus load by real-time RT-PCR

**期刊**：Journal of Clinical Virology, 2013; 56(1):65-68

**PMID**：23084006 | **DOI**：10.1016/j.jcv.2012.09.011

### 关键信息

| 项目 | 内容 |
|:---|---|
| 靶基因 | Influenza A Matrix (M) gene |
| 引物来源 | CDC 2009 pandemic pan-influenza A rRT-PCR 体系 |
| 标准品 | 含 CDC rRT-PCR 靶区域的质粒 |
| 标准曲线范围 | 1 × 10² ～ 1 × 10⁶ copies/reaction（5 个 10 倍稀释点）|
| 参与实验室 | 4 个中心 |
| 标准曲线总数 | 40 条（每个中心 10 条）|
| 批内 CV | ≤ 5% |
| 虚拟定量与本地标准曲线相关性 | R² = 0.9655 |
| Bland-Altman < 0.5 log₁₀ 偏差样本比例 | 92.5% (111/120) |

> 注：该文主要验证虚拟定量可行性，未列出逐点 Ct 表，但证明 CDC InfA 体系的 Ct–copies 转换在 4 个实验室间高度一致（CV ≤ 5%）。标准曲线斜率和 intercept 的均值/范围需查看原文 Figure 或补充材料。

---

## 文献 4：Ahrberg et al. (2016) — H7N9 手持式 qPCR

**标题**：Handheld real-time PCR device

**期刊**：Lab on a Chip, 2016; 16(3):586-592

**PMID**：26753557 | **PMCID**：PMC4773913 | **DOI**：10.1039/c5lc01415h

### 标准曲线参数

| 参数 | 值 |
|:---|---|
| 靶标 | H7N9 基因扩增子（DNA）|
| Slope | −3.02 ± 0.16 (mean ± SD) |
| Eff | 0.91 ± 0.05 per cycle (91%) |
| LOD | 1 DNA copy |
| 循环数 | 40 cycles |
| 反应体积 | ~200 nL (virtual reaction chamber) |

> 注：非标准 CDC 体系，反应体积极小（200 nL），仅供参考。

---

## 汇总对照

| 来源 | 靶基因 | Slope | Intercept | Eff% | R² | LOD (copies/rxn) | 线性范围 (copies/rxn) | 有逐点Ct? |
|:---|:---|:---:|:---:|:---:|:---:|:---:|:---|---|
| **Zhang 2019 (singleplex)** | M (CDC) | -3.4454 | 34.753 | 95.09% | 0.998 | ~30 | 0.27–27,000 | ✅ 是 |
| **Zhang 2019 (5-plex)** | M (CDC) | -3.5860 | 35.594 | 90.05% | 0.998 | ~30 | 0.27–27,000 | ✅ 是 |
| Di Trani 2006 | M (MGB) | 3.43 | ~40.8 | 95.7% | 0.998 | 5–50 | 5–5×10⁸ | ❌ Figure |
| Piralla 2013 | M (CDC) | — | — | — | 0.966* | — | 10²–10⁶ | ❌ (40条汇总) |
| Ahrberg 2016 | H7N9 | 3.02 | — | 91% | — | 1 | — | ❌ |

> \* Piralla et al. 报告的是虚拟定量与本地标准曲线间的 R² = 0.9655（非标准曲线本身的 R²）。

---

## CDC/WHO 标准 InfA rRT-PCR 引物探针序列（官方）

| 组分 | 序列 (5'→3') | 来源 |
|:---|:---|---|
| InfA Forward | GAC CRA TCC TGT CAC CTC TGA C | CDC/WHO (2009) |
| InfA Reverse | AGG GCA TTY TGG ACA AAK CGT CTA | CDC/WHO (2009) |
| InfA Probe | FAM-TGC AGT CCT CGC TCA CTG GGC ACG-BHQ1 | CDC/WHO (2009) |
| 靶标 | M gene (conserved region) | — |
| Amplicon | 106 bp | Zhang 2019 |

---

## 推荐仿真使用

1. **首选数据源**：Zhang et al. (2019) IAV singleplex — 已提供 6 点完整 (copies, Ct) 表，可直接用于 qPCR 仿真平台的标准曲线拟合验证。
2. **CDC 官方引物探针序列**：Zhang et al. Table 1 中 InfA-F/InfA-R/InfA-Pr 即 CDC/WHO 2009 pandemic protocol。
3. **宽线性范围参考**：Di Trani et al. (2006) 线性范围达 8 个数量级 (5 ～ 5×10⁸ copies)，LOD 更低 (5 copies)，可作为极限 LOD 仿真输入。

---

*生成工具：Hermes Agent + Bio.Entrez (NCBI E-utilities)*
