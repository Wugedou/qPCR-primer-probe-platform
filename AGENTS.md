# qPCR 引物探针设计仿真平台 — 项目工作手册

> 入口文档：项目目标、进度、数据位置、文件体系

## 一、项目目标

构建端到端核酸检测试剂设计平台：输入核酸序列 / GenBank 号 / 病原信息 → AI 辅助设计引物探针 → 仿真验证（不同酶、反应体系、扩增效率、LOD、线性） → 自动优化循环 → 部分替代真实实验。

目标用户：中国 CDC 团队，用于快速开发特定微生物核酸检测试剂。

## 二、当前进度（2026-07-22）

| 模块 | 状态 | 验证状态 |
|:---|:---|:---:|
| qPCR 引物设计 (TaqMan/SYBR/Beacon) | ✓ 完成 | ✓ CDC N1 黄金对照通过 |
| qPCR 扩增曲线仿真 | ✓ 完成 | ✓ 平均 ΔCt=0.88 |
| qPCR 热力学多体系参数化 | ✓ 完成 | ⬜ 待验证 |
| qPCR 优化循环 | ✓ 完成 | ⬜ 待验证 |
| RPA 引物设计 + exo-probe | ✓ 完成 | 🟡 基础验证通过, 需专题深化 |
| LAMP 4 引物设计 | ✓ 完成 | 🟡 基础验证通过, 需专题深化 |
| RPA 等温扩增仿真 | ✓ 完成 | 🟡 Tt 落在文献范围, 需专题深化 |
| LAMP 等温扩增仿真 | ✓ 完成 | 🟡 Tt 落在文献范围, 需专题深化 |
| 综合报告 | ✓ 完成 | — |

**当前专题**: [PLAN-LAMP-RPA 验证改进专题](PLAN-LAMP-RPA_验证改进专题.md)

## 三、文件体系

```
qPCR_引物探针设计仿真平台/
├── AGENTS.md                      ← 本文件
├── _shared/
│   ├── INDEX.md                   ← 项目专属文件索引
│   ├── PROGRAM-QPCR_项目文件管理程序.md
│   ├── SOP/                       ← 项目专属 SOP
│   ├── 流水账/                    ← 按日期
│   ├── 调研文献/                  ← 文献摘要、数据源
│   └── 验证报告/                  ← 验证 REPORT
├── core/                          ← 核心代码
├── validation/                    ← 验证脚本
├── data/                          ← 试剂盒预设等
├── CLAUDE.md                      ← 原始规格
└── CLAUDE_V2.md                   ← V2 修复规格
```

## 四、关键路径

- 项目根: `~/Documents/16人工智能AI/projects/qPCR_引物探针设计仿真平台/`
- 全局共享: `~/Documents/16人工智能AI/projects/_shared/`
- 验证报告(全局): `~/Documents/16人工智能AI/projects/_shared/验证报告/`
- 调研报告(全局): `~/Documents/16人工智能AI/projects/_shared/调研报告/`

## 五、运行环境

- Python 3.13 (`/usr/local/bin/python3`)
- Streamlit: http://localhost:8501
- 依赖: primer3-py 2.3.0, Biopython 1.87, numpy/scipy/pandas, plotly
- 已安装于 hermes venv

## 六、工作纪律

遵循全局 AGENTS.md (`~/Documents/16人工智能AI/projects/_shared/AGENTS.md`) 的三级文件体系：
1. 专项工作查 Program File → 按 SOP 执行 → 同步登记仓库清单 + 流水账 → 出 REPORT
2. 轻量任务只记流水账
3. SOP 修订另存新版本, 原版标 SUPERSEDED, 变更记入流水账
