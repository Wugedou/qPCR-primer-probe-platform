# PLAN-LAMP-RPA: LAMP/RPA 引物设计与仿真验证改进专题

## §1 目标

对 LAMP 和 RPA 两个等温扩增方法的设计算法和仿真模型进行系统性文献调研 + 真实数据回测 + 算法改进, 使其在真实病原体序列上的设计质量和仿真精度达到可指导实验的水平。

## §2 背景

当前状态(2026-07-22):
- LAMP 设计: 用 SARS-CoV-2 N 基因 560bp 片段验证, Notomi 2000 规范 42/42 项通过
- RPA 设计: 用 419bp 片段验证, TwistDX 规范 38/40 项通过
- 等温仿真: Tt 值落在文献报告范围
- **不足**: 仅用一个序列验证; Tt 模型参数固定未参数化; 无真实实验数据回测; LAMP 的 Loop 引物未生成

## §3 方法

### 3.1 文献调研(3 个子方向, 并行)

| 方向 | 关键问题 | 目标文献 |
|:---|:---|:---|
| LAMP 设计规范 | Tm 约束、Loop 引物、稳定性判据 | Notomi 2000 NAR, Tomita 2008 IJMM, Sahin 2011 BMC Bioinformatics |
| RPA 设计规范 | 长引物优化、exo/nfo probe、3' 端约束 | Piepenburg 2006 NAR, Li 2019, TwistDX manual |
| 等温扩增动力学 | Tt-拷贝数关系、效率模型、LOD | Tomita 2008, Iwamoto 2003, Dart 2016 SARS-CoV-2 |

### 3.2 真实数据收集

- LAMP: 公开文献中报告的标准曲线(Tt vs 拷贝数)
- RPA: TwistDX 官方验证数据 + 论文报告的标准曲线
- 已确认数据源: qpcR R 包内置数据(26 个数据集), PCRedux GitHub, RDML 官网示例

### 3.3 算法改进方向(基于调研结果确定优先级)

1. **LAMP Loop 引物生成** — 当前未生成, 需在 F1-B1c 间扫描满足 Tm/GC 的子序列
2. **LAMP 稳定性评分** — 引入综合质量分(类似 primer3 penalty), 替代当前简单的 Tm 平衡分
3. **LAMP GC 极端序列处理** — 低 GC(<40%) 或高 GC(>60%) 序列的 Tm 约束放宽策略
4. **RPA 探针 3' 端阻断** — exo-probe 3' 端需 C3-spacer 或磷酸化, 当前未标注
5. **RPA 引物 5' 端加性** — TwistDX 建议 5' 加 10-20nt 无关序列提高重组酶负载
6. **Tt 模型参数化** — 当前固定斜率/截距, 需按反应体系(Mg²⁺/酶浓度/温度)参数化
7. **Tt-拷贝数关系** — 当前用 log10 线性, 文献报告 LAMP/RPA 可能非线性(低拷贝端有平台)

### 3.4 验证策略

- 多序列验证: 至少 5 个不同病原体基因(SARS-CoV-2 N, Influenza HA, Dengue NS1, Zika E, HBV X)
- 多文献基准: 每个病原体至少 2 篇文献报告 Tt/Ct
- 量化指标: 设计成功率、Tt 误差、LOD 预测误差

## §4 步骤

| 步骤 | 内容 | 工具 | 交付物 |
|:---:|:---|:---|:---|
| 1 | 建立专题目录和文件体系 | 手动 | 本 PLAN + _shared/ |
| 2 | 并行文献调研(3 个子任务) | delegate_task | 调研文献/*.md |
| 3 | 并行真实数据源收集(3 个子任务) | delegate_task | 调研文献/真实数据源*.md |
| 4 | 汇总调研结果, 确定算法改进优先级 | 分析 | 改进路线图.md |
| 5 | 实施算法改进(LAMP + RPA) | CC / 手动 | core/lamp_design_v2.py, rpa_design_v2.py |
| 6 | 多序列多文献验证 | Python 脚本 | validation/multi_pathogen_validate.py |
| 7 | 出专题 REPORT | 手动 | REPORT_PLAN-LAMP-RPA_*.md |
| 8 | 沉淀为 skill | skill_manage | lamp-rpa-design skill |

## §5 验收标准

- LAMP 设计成功率: ≥80%(5 个病原体序列中至少 4 个成功)
- RPA 设计成功率: ≥80%
- Tt 预测误差: 与文献报告值平均偏差 ≤3 min
- LOD 预测误差: 与文献报告值偏差 ≤1 个数量级
- 每个改进点有文献依据 + 数据支撑

## §6 数据存放

- 调研文献: `qPCR项目/_shared/调研文献/`
- 真实数据: `qPCR项目/_shared/调研文献/真实数据/`
- 改进代码: `core/lamp_design.py`, `core/rpa_design.py`(直接更新)
- 验证脚本: `validation/multi_pathogen_validate.py`
- 专题报告: `qPCR项目/_shared/验证报告/REPORT_PLAN-LAMP-RPA_*.md`

## §7 记录规则

- 每步完成后立即在 `_shared/流水账/` 追加事件记录
- 算法改进后在 `_shared/仓库清单.md` 登记
- SOP 修订另存新版本

## §8 交付物

1. LAMP/RPA 设计算法改进版(代码)
2. 多病原体验证报告(REPORT)
3. 调研文献摘要集
4. 可复用 skill: lamp-rpa-design
