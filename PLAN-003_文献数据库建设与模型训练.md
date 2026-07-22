# PLAN-003: 三种方法文献数据库建设与模型训练

## §1 目标

PLAN-002 的教训：数据不够就跑交叉验证，LAMP 仅 2 个来源 10 个点，留出法形同虚设，V3 泛化失败。根因是"边检索边验证"的策略导致数据不足时硬跑验证。

**本 PLAN 策略转向：数据先行，验证后行。**

分两阶段：
1. **Phase 1 — 建库**：为 qPCR、RPA、LAMP 三种方法各自建立结构化文献数据库，每篇文献提取引物/探针序列 + 全部验证数据（Ct/Tt、扩增曲线、LOD、干扰实验、线性范围），数据量达到统计要求后才进入 Phase 2
2. **Phase 2 — 验证与训练**：基于数据库做模型验证、交叉验证、参数训练

**Phase 1 不做完不启动 Phase 2。**

## §2 背景

### 2.1 PLAN-002 暴露的问题

| 问题 | 根因 | 后果 |
|:---|:---|:---|
| LAMP Tt 模型 V3 泛化失败（RMSE 差 8.7min） | 仅 2 个文献来源，交叉验证无统计意义 | 无法证明模型泛化能力 |
| qPCR 仅 2/4 病原体验证 | Dengue/HBV/MTB 全文 Ct 数据未提取 | 病原体覆盖不足 |
| Tt 文献 14+8 PMID 已列清单但未筛选 | 逐篇筛选含 Tt 表格的文献未推进 | 数据停留在 PMID 清单层面 |
| Behrmann 2020 和 "lit2" 数据从 probit/估算得来 | 非真实逐点 Tt，数据质量存疑 | 交叉验证输入不可靠 |

### 2.2 当前已有数据（不丢弃，纳入新数据库）

| 已有数据 | 路径 | 状态 |
|:---|:---|:---|
| CDC N1 qPCR 黄金对照 (3 套标准曲线) | `validation/CDC_N1_golden_standard.py` | ✅ 可用 |
| FluA qPCR 黄金对照 (2 套) | `validation/expanded_golden_standards.py` | ✅ 可用 |
| RPA Tt 数据 (MHV 5 梯度×3 重复) | `rpa_standard_curve_Tt_data_compilation.md` | ✅ 可用 |
| LAMP Tt 数据 (SARS-CoV-2 5 梯度) | `lamp_rpa_kinetics_research_report.md` | ✅ 可用 |
| RPA 14 + LAMP 8 PMID 检索清单 | `rpa_lamp_Tt_search_results.md` | ⏸ 待逐篇处理 |
| Dengue 6 + HBV 3 + MTB 8 PMID | `qPCR_*_Ct_data.md` | ⏸ 待全文提取 |

## §3 方法

### 3.1 文献筛选标准

**纳入标准（全部满足）**：
1. 期刊影响因子 ≥3.0，或被引用 ≥50 次，或为方法学奠基文献（如 Notomi 2000, Piepenburg 2006）
2. **必须含引物/探针序列**（F/R primer + probe，非仅"商业试剂盒黑盒"）
3. 含至少 1 项以下数据：标准曲线参数 (slope/intercept/Eff/R²)、逐点 (copies→Ct/Tt)、LOD、干扰实验、线性范围
4. 有 GenBank 号或可追溯的靶序列

**排除标准**：
- 仅做定性检测（无定量数据）
- 引物序列未公开（仅写"使用 XX 试剂盒"）
- 纯综述无原始数据（但可提取参考文献清单）

**质量分级**：
| 等级 | 标准 | 数据完整度 |
|:---:|:---|:---|
| ⭐⭐⭐ A 级 | 引物+探针序列 + 逐点 Ct/Tt 表 + 标准曲线参数 + LOD + 干扰实验 | 完整 |
| ⭐⭐ B 级 | 引物+探针序列 + 标准曲线参数 (slope/intercept/Eff/R²) + LOD | 较完整 |
| ⭐ C 级 | 引物+探针序列 + LOD 或部分 Ct/Tt 数据 | 基本可用 |

### 3.2 数据库字段设计（三种方法统一框架 + 方法特异字段）

**通用字段（三种方法共有）**：

```
record_id          : QR-L-001 (qPCR/RPA/LAMP-文献序号)
pmid               : 数字
doi                : 字符串
journal            : 字符串
impact_factor      : 数字
first_author       : 姓 年 (如 Vogels 2020)
year               : 数字
pathogen           : 病原体名称
target_gene        : 靶基因
genbank_accession  : GenBank 号
primer_F           : 正向引物序列 5'→3'
primer_R           : 反向引物序列 5'→3'
probe              : 探针序列 (如有)
probe_type         : TaqMan / exo / nfo / fpg / SYBR / Beacon / 侧流
product_size       : 产物大小 bp
reaction_conditions: 酶/温度/Mg²⁺/dNTP/体系体积/仪器
slope              : 标准曲线斜率
intercept          : 标准曲线截距
efficiency_pct     : 扩增效率 %
r_squared          : R²
lod_copies         : LOD (copies/reaction)
linear_range       : 线性范围 (low, high)
interference_test  : 干扰实验结果 (特异性/交叉反应/基质效应)
data_points        : [(copies, Ct 或 Tt), ...]  逐点数据
quality_grade      : A/B/C
notes              : 备注
```

**方法特异字段**：

| 方法 | 特异字段 |
|:---|:---|
| qPCR | `master_mix`(试剂盒), `instrument`(仪器), `fluor_channel`(荧光通道), `baseline_method`(基线校正) |
| RPA | `recombinase_source`(TwistDx/国产), `mg_concentration`(MgAc mM), `probe_modification`(3'封闭方式), `amplicon_optimal`(是否落在 100-200bp 最优区间) |
| LAMP | `primer_set`(F3/B3/FIP/BIP/Loop F/Loop B 全序列), `calcein`(钙黄绿素/SYBR/酚红), `loop_primer`(是否含 Loop 引物), `pyrophosphate`(焦磷酸镁比浊) |

### 3.3 数据库存储格式

```
data/
├── database/
│   ├── qpcr_literature_database.json     ← qPCR 文献数据库 (JSON, 结构化)
│   ├── rpa_literature_database.json      ← RPA 文献数据库
│   ├── lamp_literature_database.json     ← LAMP 文献数据库
│   ├── qpcr_literature_database.md       ← Markdown 可读版 (人类审阅用)
│   ├── rpa_literature_database.md
│   └── lamp_literature_database.md
```

JSON 用于程序读取，Markdown 用于人类审阅和三级文件体系登记。

### 3.4 检索策略

**PubMed E-utilities 检索词**（每方法分 3 轮）：

| 轮次 | 策略 | 目标 |
|:---:|:---|:---|
| R1 | 病原体名 + 方法名 + "standard curve" OR "limit of detection" OR "linearity" | 定量验证文献 |
| R2 | 病原体名 + 方法名 + "primer" OR "probe" + "sequence" | 含引物序列文献 |
| R3 | 方法名 + "interference" OR "specificity" OR "cross-reactivity" OR "matrix effect" | 干扰实验文献 |

**病原体覆盖目标**（每方法至少覆盖以下，GC% 从 38% 到 65%）：

| # | 病原体 | 靶基因 | GC% 参考 | 备注 |
|:---:|:---|:---|---:|:---|
| 1 | SARS-CoV-2 | N / RdRp | ~48% | 已有 CDC N1 |
| 2 | Influenza A | M / HA | ~40% | 已有 Zhang 2019 |
| 3 | Dengue virus | NS1 / 3'UTR | ~47% | 待提取 |
| 4 | HBV | X / S | ~55% | 待提取 |
| 5 | M. tuberculosis | IS6110 / IS1081 | ~63% | 待提取 |
| 6 | HIV-1 | gag / LTR | ~45% | 新增 |
| 7 | ZIKV | E / NS5 | ~51% | 新增 |
| 8 | E. coli O157 | rfbE / stx2 | ~51% | 新增（细菌 GC 接近中位） |

### 3.5 文献全文获取方法

1. **PMC Open Access**: `efetch -db pmc -id XXXX -format xml`（优先，全文免费）
2. **PubMed Abstract**: `efetch -db pubmed -id XXXX -format abstract`（退而求其次）
3. **期刊官网 PDF**: 部分文献需从期刊网站获取（如 NAR, Clinical Chemistry 有开放获取）
4. **补充材料**: 引物序列和逐点 Ct 表常在 Supplementary Table 中
5. 记录获取状态：`full_text` / `abstract_only` / `supplementary_only`

## §4 步骤

### Phase 1 — 建库（步骤 1-8）

| 步骤 | 内容 | 工具 | 交付物 | 预估工时 |
|:---:|:---|:---|:---|:---|
| 1 | 创建 `data/database/` 目录结构 + 数据库 JSON schema 模板 | 手动 | 目录 + schema | 0.5h |
| 2 | qPCR 文献检索 R1-R3（8 病原体 × 3 轮 = 24 次检索） | Biopython Entrez | `qpcr_search_raw.json` | 2h |
| 3 | qPCR 文献逐篇筛选 + 全文获取 + 数据提取 | Python 脚本 + 手动 | `qpcr_literature_database.json` (≥15 篇) | 6h |
| 4 | RPA 文献检索 R1-R3（8 病原体 × 3 轮） | Biopython Entrez | `rpa_search_raw.json` | 2h |
| 5 | RPA 文献逐篇筛选 + 全文获取 + 数据提取 | Python 脚本 + 手动 | `rpa_literature_database.json` (≥10 篇) | 5h |
| 6 | LAMP 文献检索 R1-R3（8 病原体 × 3 轮） | Biopython Entrez | `lamp_search_raw.json` | 2h |
| 7 | LAMP 文献逐篇筛选 + 全文获取 + 数据提取 | Python 脚本 + 手动 | `lamp_literature_database.json` (≥10 篇) | 5h |
| 8 | 数据库质量审核：字段完整率、A 级文献占比、病原体覆盖率、数据点总数 | Python 脚本 | `REPORT_PLAN-003_数据库质量审核.md` | 1h |

**Phase 1 完成判定**（不达标不进入 Phase 2）：

| 指标 | 目标 | 说明 |
|:---|:---|:---|
| qPCR A 级文献数 | ≥8 篇 | 含逐点 Ct 表 |
| RPA A 级文献数 | ≥6 篇 | 含逐点 Tt 表 |
| LAMP A 级文献数 | ≥6 篇 | 含逐点 Tt 表 |
| 三方法各自病原体覆盖 | ≥6 种 | GC% 覆盖 38-65% |
| 三方法各自数据点总数 | ≥50 个 | (copies, Ct/Tt) 对 |
| 引物探针序列完整率 | ≥90% | A 级文献必须 100% |
| 干扰实验数据文献数 | ≥5 篇/方法 | 特异性/交叉反应/基质 |

### Phase 2 — 验证与训练（步骤 9-15，Phase 1 达标后启动）

| 步骤 | 内容 | 工具 | 交付物 |
|:---:|:---|:---|:---|
| 9 | qPCR 多病原体 Ct 对照（用数据库 A 级文献逐点对照） | Python | `REPORT_PLAN-003_qPCR验证.md` |
| 10 | Tt 模型交叉验证（用数据库 A 级文献留出法，训练 60% / 测试 40%） | Python | `REPORT_PLAN-003_Tt交叉验证.md` |
| 11 | Tt 模型参数训练（按反应体系分组：酶浓度/Mg²⁺/温度） | Python/scipy | 更新 `core/amplification.py` |
| 12 | 扩增曲线形状验证（如数据库含 raw fluorescence 数据） | Python | `REPORT_PLAN-003_曲线形状.md` |
| 13 | Poisson 低拷贝模型验证（用数据库 LOD/probit 数据） | Python | `REPORT_PLAN-003_Poisson验证.md` |
| 14 | 干扰实验仿真验证（用数据库特异性/基质数据） | Python | `REPORT_PLAN-003_干扰验证.md` |
| 15 | 汇总 REPORT + skill 更新 | 手动 | `REPORT_PLAN-003_总结.md` |

## §5 验收标准

### Phase 1 验收（建库）

| # | 验收项 | 目标 |
|:---:|:---|:---|
| P1-1 | qPCR 文献数 | ≥15 篇（A 级 ≥8） |
| P1-2 | RPA 文献数 | ≥10 篇（A 级 ≥6） |
| P1-3 | LAMP 文献数 | ≥10 篇（A 级 ≥6） |
| P1-4 | 病原体覆盖 | 每方法 ≥6 种，GC% 覆盖 38-65% |
| P1-5 | 数据点总数 | 每方法 ≥50 个 (copies, Ct/Tt) |
| P1-6 | 引物序列完整率 | ≥90% |
| P1-7 | 干扰实验文献 | 每方法 ≥5 篇 |

### Phase 2 验收（验证与训练）

| # | 验收项 | 目标 |
|:---:|:---|:---|
| P2-1 | qPCR ΔCt | ≥6 病原体，各 ≤2.0，avg ≤1.0 |
| P2-2 | Tt 测试集 RMSE | RPA ≤2.0min, LAMP ≤2.0min |
| P2-3 | Tt 模型泛化 | 训练/测试 RMSE 差 ≤1.0min |
| P2-4 | LOD 偏差 | 与文献偏差 ≤0.5 个对数级 |
| P2-5 | 干扰实验预测 | 特异性 ≥90% 与文献一致 |

## §6 数据存放

```
data/database/
├── qpcr_literature_database.json      ← R-0013
├── rpa_literature_database.json       ← R-0014
├── lamp_literature_database.json      ← R-0015
├── qpcr_literature_database.md        ← R-0016
├── rpa_literature_database.md         ← R-0017
└── lamp_literature_database.md        ← R-0018

_shared/调研文献/真实数据/
├── qpcr_search_raw.json               ← 检索原始结果
├── rpa_search_raw.json
└── lamp_search_raw.json

_shared/验证报告/
├── REPORT_PLAN-003_数据库质量审核.md   ← R-0019
├── REPORT_PLAN-003_qPCR验证.md
├── REPORT_PLAN-003_Tt交叉验证.md
├── REPORT_PLAN-003_曲线形状.md
├── REPORT_PLAN-003_Poisson验证.md
├── REPORT_PLAN-003_干扰验证.md
└── REPORT_PLAN-003_总结.md

_shared/流水账/
└── 2026-07-XX_PLAN-003_文献数据库建设.md
```

## §7 记录规则

- 每步完成后立即在流水账追加事件记录
- 新文档从 R-0013 开始编号
- 同步更新 INDEX.md
- JSON 数据库文件登记 R-NNNN（因含结构化文献数据，视为文档）
- 代码脚本不登记 R-NNNN

## §8 交付物

1. **三个结构化文献数据库**（JSON + Markdown 双格式）
2. **数据库质量审核报告**（Phase 1 达标判定）
3. **验证与训练报告集**（Phase 2，6 个子 REPORT + 1 个总结）
4. **更新后的 Tt 模型参数**（按反应体系分组训练）
5. **skill 更新**：nucleic-acid-assay-design-validation

---

## 附：执行者须知（Bot A）

### 项目路径
- 项目根：`~/Documents/16人工智能AI/projects/qPCR_引物探针设计仿真平台/`
- 全局共享：`~/Documents/16人工智能AI/projects/_shared/`

### 运行环境
- Python 3.13 (`/usr/local/bin/python3`)，hermes venv
- 依赖：primer3-py 2.3.0, Biopython 1.87, numpy/scipy/pandas, plotly
- Streamlit: http://localhost:8501

### 文件体系（三级文件体系，强制执行）
1. 先查 `_shared/PROGRAM-QPCR_项目文件管理程序.md` → 确认规则
2. 执行中同步更新 `_shared/INDEX.md`（R-NNNN）+ 流水账
3. 产出 REPORT 到 `_shared/验证报告/`
4. R-NNNN 从 R-0013 开始

### 已有数据（纳入数据库，不重复检索）
- `validation/CDC_N1_golden_standard.py` — CDC N1 qPCR 3 套标准曲线（Vogels 2020, Stevens 2020）
- `validation/expanded_golden_standards.py` — FluA 2 套（Zhang 2019 PMID:31130238）+ Di Trani 2006
- `_shared/调研文献/真实数据/rpa_standard_curve_Tt_data_compilation.md` — RPA 3 篇含 Tt 数据 + 20 篇文献清单
- `_shared/调研文献/真实数据/rpa_lamp_Tt_search_results.md` — RPA 14 + LAMP 8 PMID 待逐篇筛选
- `_shared/调研文献/真实数据/qPCR_Dengue_Ct_data.md` — 6 PMID 待提取
- `_shared/调研文献/真实数据/qPCR_HBV_Ct_data.md` — 3 PMID 待提取
- `_shared/调研文献/真实数据/qPCR_MTB_Ct_data.md` — 8 PMID 待提取

### 关键要求
1. **Phase 1 不做完不启动 Phase 2**——数据量达标前不跑验证
2. **文献必须含引物/探针序列**——黑盒试剂盒文献不纳入
3. **优先 A 级文献**——含逐点 Ct/Tt 表格的文献价值最高
4. **数据库 JSON + Markdown 双格式**——JSON 给程序读，Markdown 给人审阅
5. **逐篇提取全文**——不只是摘要，引物序列和 Ct 表常在正文/补充材料中
6. **记录反应条件**——酶/温度/Mg²⁺/体系/仪器，用于后续按体系分组训练

### 文献检索技术要点
- PubMed E-utilities: `esearch + efetch`（Biopython Bio.Entrez）
- PMC 全文: `efetch db=pmc`（优先用于提取引物序列和 Ct 表）
- OpenAlex API 可辅助筛选高引用文献
- 期刊官网 PDF：Clinical Chemistry, NAR, J Clin Microbiol 等有开放获取
- 补充材料：引物序列和逐点 Ct 表常在 Supplementary Table 1/2 中

### 数据库 JSON Schema 示例

```json
{
  "record_id": "QR-L-001",
  "method": "qPCR",
  "pmid": 32546853,
  "doi": "10.1038/s41598-020-68602-2",
  "journal": "Scientific Reports",
  "impact_factor": 4.6,
  "first_author": "Vogels 2020",
  "pathogen": "SARS-CoV-2",
  "target_gene": "N",
  "genbank_accession": "MN908947.1",
  "primer_F": "GACCCCAAAATCAGCGAAAT",
  "primer_R": "TCTGGTTACTGCCAGTTGAATCTG",
  "probe": "ACCCCGCATTACGTTTGGTGGACC",
  "probe_type": "TaqMan",
  "product_size": 93,
  "reaction_conditions": {
    "master_mix": "TaqPath 1-Step RT-qPCR",
    "instrument": "ABI 7500 Fast",
    "mg_conc": 3.0,
    "dntp_conc": 0.8,
    "primer_conc_nM": 200,
    "reaction_volume_uL": 20
  },
  "standard_curve": {
    "slope": -3.32,
    "intercept": 38.5,
    "efficiency_pct": 100.0,
    "r_squared": 0.998
  },
  "lod_copies": 5,
  "linear_range": [5, 5000000],
  "interference_test": {
    "specificity": "100% (no cross-react with SARS/MERS/OC43/229E)",
    "matrix_effect": "not tested",
    "inhibitors": "not tested"
  },
  "data_points": [
    {"copies": 1, "ct": 38.5},
    {"copies": 5, "ct": 36.1},
    {"copies": 10, "ct": 35.0},
    {"copies": 100, "ct": 31.7},
    {"copies": 1000, "ct": 28.4},
    {"copies": 10000, "ct": 25.1},
    {"copies": 100000, "ct": 21.7},
    {"copies": 1000000, "ct": 18.4}
  ],
  "quality_grade": "A",
  "full_text_source": "PMC7314152",
  "notes": "CDC N1 validation, multi-center"
}
```

### 病原体覆盖与 GC% 参考

| 病原体 | 靶基因 | GC% | GenBank |
|:---|:---|---:|:---|
| SARS-CoV-2 | N | ~48% | MN908947.1 |
| Influenza A | M | ~40% | CY230218.1 |
| Dengue 2 | NS1 | ~47% | M19197.1 |
| HBV | X | ~55% | V01460.1 |
| M. tuberculosis | IS6110 | ~63% | X17348.1 |
| HIV-1 | gag | ~45% | 检索确定 |
| ZIKV | E | ~51% | 检索确定 |
| E. coli O157 | rfbE | ~51% | 检索确定 |
