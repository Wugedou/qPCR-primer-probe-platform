# PLAN-002: 验证扩展与模型泛化

## §1 目标

当前仿真平台已完成概念验证（qPCR vs CDC N1 ΔCt=0.88, LAMP/RPA Tt RMSE<0.6min），但验证数据存在结构性缺口：qPCR 仅 1 个病原体验证、Tt 模型仅 2 篇文献 10 个点拟合（自证嫌疑）、无真实扩增曲线形状验证、低拷贝端无 Poisson 模型。

本 PLAN 目标：**将验证从"概念验证"升级到"可指导实验"级别**，证明模型在多病原体、多反应体系、低拷贝端具有泛化能力。

## §2 背景

### 2.1 当前验证数据清单

| 模块 | 验证数据 | 数据量 | 来源 |
|:---|:---|:---|:---|
| qPCR Ct 对照 | CDC N1 (SARS-CoV-2 N) | 3 条标准曲线 × 6-8 点 | Vogels 2020 PMID:32546853; Stevens 2020 PMID:32073526 |
| RPA Tt 回测 | MHV RT-RPA | 5 梯度 × 3 重复 | Wang 2022 PMID:36532497 |
| LAMP Tt 回测 | SARS-CoV-2 RT-LAMP 中位值 | 5 梯度 | PMID:34006453 |
| LAMP/RPA 设计成功率 | 5 个 NCBI 病原体序列 | 5 病原体 × (LAMP 3 组 + RPA 5 组) | MN908947, CY230218, M19197, V01460, X17348 |

### 2.2 结构性缺口（按严重度排序）

| # | 缺口 | 严重度 | 说明 |
|:---:|:---|:---:|:---|
| G1 | qPCR 仅 1 个病原体验证 | 🔴 高 | 无法证明 Ct 模型对不同 GC%/模板长度/二级结构泛化 |
| G2 | Tt 模型仅 2 篇文献 10 个点 | 🔴 高 | 斜率/截距从 2 篇拟合再回测 2 篇，RMSE 0.5min 有自证嫌疑 |
| G3 | 无真实扩增曲线回放 | 🟡 中 | 仅比对 Tt 终点值，S 曲线形状（陡度/平台/baseline）未验证 |
| G4 | 低拷贝 Poisson 缺失 | 🟡 中 | 1-10 copies 用确定性 S 曲线，应加随机过程 |
| G5 | LAMP 一手规范未读 | 🟢 低 | Notomi 2000/Tomita 2008 原文 2 次调研超时 |

## §3 方法

### 3.1 qPCR 多病原体 Ct 对照扩展（对应 G1）

**目标**：从文献提取 4 个额外病原体的 qPCR 标准曲线参数（slope/intercept/Eff/R²/LOD/实测 Ct），与仿真平台逐点对照。

| 病原体 | 靶基因 | 拟检索文献方向 | GC% 覆盖 |
|:---|:---|:---|---:|
| Influenza A | HA/M1 | CDC Flu RT-PCR Panel, WHO PCR protocol | ~40% |
| Dengue virus | NS1/3'UTR | CDC Arbovirus RT-PCR, Balm 2012 | ~47% |
| HBV | X/S | WHO HBV qPCR standard, Candotti 2018 | ~55% |
| M. tuberculosis | IS6110 | CDC TB qPCR, Helb 2010 | ~63% |

**数据提取要求**：
- 每个病原体至少 1 篇含逐拷贝数 Ct 表格的文献（非仅有 slope/intercept）
- 提取字段：primer/probe 序列（若有）、product size、slope、intercept、Eff%、R²、LOD、线性范围、逐点 (copies, Ct)
- 保存格式：扩展 `validation/CDC_N1_golden_standard.py`，新增 4 个 `GoldenCurve` 对象

**验证方式**：
1. 对每个病原体靶序列运行平台引物设计 → 检查设计是否与文献引物吻合（Tm/产物大小/位置）
2. 对每个病原体运行扩增仿真（6 梯度 copies）→ 与文献 Ct 逐点对照
3. 计算 ΔCt、RMSE、效率偏差

### 3.2 Tt 模型交叉验证（对应 G2）

**目标**：收集 5-10 篇独立文献的标准曲线数据，按"训练集/测试集"分割做交叉验证。

**分两步**：

**Step A — 数据收集（目标 8-12 篇）**

RPA 标准曲线文献（已有 3 篇含 Tt 数据，需补充 5-7 篇）：

| # | 病原体 | 来源 | 当前状态 |
|:---:|:---|:---|:---|
| 1 | MHV | Wang 2022 PMID:36532497 | ✅ 已有 5 梯度×3 重复 |
| 2 | SHIV | Chen 2020 PMID:31655105 | ✅ 已有但仅 1 点 |
| 3 | 铜绿假单胞菌 | Jin 2018 PMID:29690742 | ⚠ 仅检出率无 Tt |
| 4 | SARS-CoV-2 | Behrmann 2020 PMID:32384153 | 需提取 probit → Tt |
| 5 | SARS-CoV-2 | Lau 2021 PMID:33406112 | 需提取 |
| 6 | ARV | Ma 2020 Front Vet Sci | 需提取 |
| 7 | HIV-1 | Lillis 2016 PMID:26821087 | 需提取 |
| 8-12 | 新检索 | SARS-CoV-2/流感/登革/寨卡 RPA | 需检索补充 |

LAMP 标准曲线文献（已有 1 篇 5 梯度，需补充 4-6 篇）：

| # | 病原体 | 来源 | 当前状态 |
|:---:|:---|:---|:---|
| 1 | SARS-CoV-2 | PMID:34006453 | ✅ 已有 5 梯度 |
| 2 | SARS-CoV-2 | 其他 RT-LAMP 文献 | 需检索 |
| 3 | MTB | IS6110 LAMP 文献 | 需检索 |
| 4 | 流感 | IA LAMP 文献 | 需检索 |
| 5 | 登革 | DV LAMP 文献 | 需检索 |
| 6-8 | 其他 | 寨卡/HBV/HIV LAMP | 需检索 |

**数据提取要求**：
- 优先含逐拷贝数 Tt/Ct 表格的文献（非仅有 LOD）
- 记录反应条件：酶/温度/Mg²⁺/反应体积/探针类型
- 保存到 `_shared/调研文献/真实数据/` 下新建 `lamp_standard_curve_compilation.md` 和扩展 `rpa_standard_curve_Tt_data_compilation.md`

**Step B — 交叉验证**

1. 将所有文献标准曲线按"留出法"分组：
   - 训练集：60% 文献（含当前已用的 2 篇）
   - 测试集：40% 文献（留出，不参与拟合）
2. 用训练集重新拟合 Tt = a × log10(copies) + b（分别 RPA/LAMP）
3. 用测试集计算 RMSE/平均偏差
4. 对比：当前 v2 模型（2 篇拟合）vs 交叉验证模型（多文献拟合）的参数差异
5. 如测试集 RMSE 显著高于训练集 → 说明过拟合，需参数化体系变量

### 3.3 真实扩增曲线回放（对应 G3）

**目标**：获取等温扩增原始荧光曲线数据（raw fluorescence vs time），验证 S 曲线形状而非仅 Tt 终点。

**数据源**：
1. **qpcR R 包**：内置 26 个数据集（但为 qPCR 数据，可验证 qPCR S 曲线形状）
2. **PCRedux GitHub**：含扩增曲线特征提取数据（RDML 格式）
3. **RDML 官网**：rdml.org 示例文件
4. **文献补充材料**：部分 LAMP/RPA 论文含 Supplementary Data 的 raw fluorescence

**验证方式**：
1. 从 qpcR/PCRedux 提取 qPCR 荧光曲线 → 拟合平台的 S 曲线模型 → 比较 baseline/plateau/steepness
2. 搜索含 LAMP/RPA raw fluorescence 的文献补充材料 → 同上
3. 如找不到等温 raw data，至少完成 qPCR S 曲线形状验证 + 记录等温数据的不可获得性

### 3.4 低拷贝 Poisson 模型（对应 G4）

**目标**：1-10 copies 区间加入随机过程模型，替代当前确定性 S 曲线。

**方法**：
1. 实现 Poisson 采样：每次仿真从 Poisson(λ=copies) 抽取实际起始分子数
2. 对 1/5/10 copies 各运行 N=1000 次蒙特卡洛仿真
3. 统计检出概率（Ct < 40 或 Tt < max_time 的比例）
4. 与文献报告的 LOD probit 数据对照（Behrmann 2020 SARS-CoV-2 probit 已有）
5. 输出：检出概率-拷贝数曲线 + 95% LOD 置信区间

### 3.5 LAMP 一手规范补充（对应 G5，低优先级）

**目标**：获取 Notomi 2000 NAR 原文和 Tomita 2008 IJMM 原文的关键设计约束。

**方法**：
1. 从 PubMed Central 全文获取 Notomi 2000 (PMID: 10807277)
2. 从期刊获取 Tomita 2008 IJMM
3. 提取：F1c-F2 Tm 差的实际描述（是否严格要求 ≥2°C）、Loop 引物设计规范、稳定性判据
4. 与当前 v2 算法的降级策略对照（当前 1.5-2.0°C 降级是否合理）

## §4 步骤

| 步骤 | 内容 | 工具 | 交付物 | 预估工时 |
|:---:|:---|:---|:---|:---|
| 1 | 建立 PLAN + 文件体系登记 | 手动 | 本 PLAN, INDEX 更新, 流水账 | 0.5h |
| 2 | 并行文献检索：4 个 qPCR 病原体标准曲线 | delegate_task ×4 或 PubMed API | `调研文献/真实数据/qPCR_multi_pathogen_Ct_data.md` | 2-4h |
| 3 | 并行文献检索：补充 RPA/LAMP Tt 文献 (各 5-7 篇) | delegate_task ×2 | 扩展 rpa_*.md, 新建 lamp_*.md | 2-4h |
| 4 | 汇总数据，提取逐点 (copies, Ct/Tt) 表 | Python 脚本 | `validation/expanded_golden_standards.py` | 1h |
| 5 | qPCR 多病原体验证：设计+仿真+对照 | Python 脚本 | `validation/multi_pathogen_qpcr_validate.py` | 2h |
| 6 | Tt 模型交叉验证：训练/测试分割+重拟合+对比 | Python 脚本 | `validation/tt_cross_validation.py` | 2h |
| 7 | 真实扩增曲线回放（qpcR/PCRedux） | R 或 Python | `validation/curve_shape_validate.py` | 2h |
| 8 | 低拷贝 Poisson 蒙特卡洛模型 | Python 脚本 | `core/poisson_monte_carlo.py` + `validation/poisson_validate.py` | 2h |
| 9 | LAMP 一手规范补充（PMC 全文） | 浏览器/PubMed | `调研文献/LAMP设计规范_Notomi_Tomita原文.md` | 1h |
| 10 | 汇总出 REPORT | 手动 | `REPORT_PLAN-002_验证扩展总结.md` | 1h |
| 11 | 更新 INDEX + 流水账 + AGENTS.md 进度 | 手动 | 文件更新 | 0.5h |
| 12 | 沉淀更新 skill: nucleic-acid-assay-design-validation | skill_manage | skill 更新 | 0.5h |

**并行策略**：步骤 2+3 可并行（4+2=6 个子任务），步骤 5/6/7/8 在数据就绪后可部分并行。

## §5 验收标准

| # | 验收项 | 目标 | 判定方式 |
|:---:|:---|:---|:---|
| V1 | qPCR 多病原体对照 | ≥4 个病原体，每个 ΔCt ≤2.0，平均 ΔCt ≤1.5 | 逐点对照表 |
| V2 | Tt 模型交叉验证 | 测试集 RMSE：RPA ≤2.0min, LAMP ≤2.0min | 留出法交叉验证 |
| V3 | Tt 模型泛化判定 | 训练集与测试集 RMSE 差值 ≤1.0min | 若 >1.0 则标记过拟合 |
| V4 | S 曲线形状验证 | qPCR 至少 3 条真实曲线拟合 R²≥0.95 | 曲线对照图 |
| V5 | Poisson 检出概率 | 95% LOD 与文献 probit 偏差 ≤0.5 个对数级 | probit 曲线对照 |
| V6 | LAMP 规范闭环 | Notomi 2000/Tomita 2008 关键约束已提取 | 文献摘要 |
| V7 | 文献数据量 | RPA Tt 文献 ≥8 篇, LAMP Tt 文献 ≥5 篇 | 数据汇总文件 |

## §6 数据存放

```
_shared/
├── 调研文献/
│   ├── 真实数据/
│   │   ├── rpa_standard_curve_Tt_data_compilation.md    ← 扩展（已有 3 篇 → ≥8 篇）
│   │   ├── lamp_standard_curve_compilation.md           ← 新建（≥5 篇）
│   │   └── qPCR_multi_pathogen_Ct_data.md               ← 新建（≥4 病原体）
│   └── LAMP设计规范_Notomi_Tomita原文.md                ← 新建
├── 验证报告/
│   ├── REPORT_PLAN-002_qPCR多病原体验证.md
│   ├── REPORT_PLAN-002_Tt交叉验证.md
│   ├── REPORT_PLAN-002_S曲线形状验证.md
│   ├── REPORT_PLAN-002_Poisson低拷贝验证.md
│   └── REPORT_PLAN-002_验证扩展总结.md
└── 流水账/
    └── 2026-07-XX_PLAN-002_验证扩展.md

validation/
├── expanded_golden_standards.py       ← 扩展黄金对照（4+ 病原体）
├── multi_pathogen_qpcr_validate.py   ← qPCR 多病原体验证脚本
├── tt_cross_validation.py            ← Tt 交叉验证脚本
├── curve_shape_validate.py           ← S 曲线形状验证
└── poisson_validate.py               ← Poisson 蒙特卡洛验证

core/
└── poisson_monte_carlo.py            ← Poisson 随机过程模型
```

## §7 记录规则

- 每步完成后立即在 `_shared/流水账/` 追加事件记录
- 每个新文档分配 R-NNNN 编号（从 R-0013 开始）
- 同步更新 `_shared/INDEX.md` 仓库清单
- SOP 修订（若有）另存新版本
- 代码文件不登记 R-NNNN，但需在流水账中记录

## §8 交付物

1. **数据层**：≥4 病原体 qPCR 标准曲线 + ≥8 篇 RPA Tt 文献 + ≥5 篇 LAMP Tt 文献
2. **代码层**：5 个验证脚本 + 1 个 Poisson 模型模块
3. **报告层**：5 个子 REPORT + 1 个专题总结 REPORT
4. **判定层**：V1-V7 验收标准逐项判定表
5. **skill 更新**：nucleic-acid-assay-design-validation skill 同步更新验证数据和方法

---

## 附：执行者须知（Bot A）

### 项目路径
- 项目根：`~/Documents/16人工智能AI/projects/qPCR_引物探针设计仿真平台/`
- 全局共享：`~/Documents/16人工智能AI/projects/_shared/`

### 运行环境
- Python 3.13 (`/usr/local/bin/python3`)，已安装于 hermes venv
- 依赖：primer3-py 2.3.0, Biopython 1.87, numpy/scipy/pandas, plotly
- Streamlit: http://localhost:8501

### 文件体系（三级文件体系，强制执行）
1. **先查** `_shared/PROGRAM-QPCR_项目文件管理程序.md` → 确认规则
2. **执行中** 同步更新 `_shared/INDEX.md`（R-NNNN 登记）+ 流水账
3. **产出** REPORT 到 `_shared/验证报告/`
4. R-NNNN 从 R-0013 开始编号
5. 代码文件(.py)不登记 R-NNNN，但流水账记录

### 关键代码文件（需读取理解后再修改）
- `core/amplification.py` — Tt 模型当前参数（RPA: midpoint=13.4-1.55*log10(N0); LAMP: midpoint=40.77-4.60*log10(N0)）
- `core/lamp_design.py` — LAMP 设计 v2（F1c 优先 + 动态降级）
- `core/rpa_design.py` — RPA 设计 v2（产物收紧 + homopolymer 检查）
- `validation/CDC_N1_golden_standard.py` — 当前黄金对照（需扩展）
- `validation/real_pathogen_sequences.json` — 5 个真实病原体序列

### 已有数据文件（不要重复检索）
- `_shared/调研文献/真实数据/rpa_standard_curve_Tt_data_compilation.md` — 已有 3 篇 RPA Tt 数据（MHV/SHIV/铜绿假单胞菌）+ 20 篇文献清单
- `_shared/调研文献/lamp_rpa_kinetics_research_report.md` — 21 篇等温动力学文献
- `_shared/调研文献/RPA引物探针设计规范_5篇核心文献.md` — RPA 设计规范
- `_shared/调研文献/LAMP设计规范_文献调研.md` — LAMP 设计规范

### 文献检索方法
- PubMed E-utilities: `esearch + efetch`（Biopython Bio.Entrez）
- OpenAlex API: `https://api.openalex.org/works?search=...`
- PMC 全文: `efetch db=pmc`
- 优先提取含**逐拷贝数 Ct/Tt 表格**的文献（非仅有 slope/intercept 或 LOD）
- 记录反应条件：酶/温度/Mg²⁺/反应体积/探针类型

### 注意事项
1. **不要修改已有 REPORT**（R-0004 到 R-0012），只新增
2. **Tt 模型交叉验证是核心**——当前 v2 的 RMSE 0.5min 是用拟合数据回测自身，必须用留出法证明泛化
3. **qPCR 多病原体对照**——如果某病原体找不到含逐点 Ct 的文献，记录缺口并降级为仅用 slope/intercept 对照
4. **Poisson 模型**——Behrmann 2020 (PMID:32384153) 已有 probit 数据，可作为低拷贝验证基准
5. 完成后更新 `AGENTS.md` 的进度表
