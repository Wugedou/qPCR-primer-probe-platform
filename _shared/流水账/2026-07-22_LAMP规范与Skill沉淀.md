# 2026-07-22 LAMP 设计规范 + Skill 沉淀

## 工作项
1. **LAMP 设计规范调研** — 第 3 次自动化调研超时 (50 次 API 调用)
   - 决定: 不再重试, 基于已有调研结果 + 领域知识综合编写
   - 产出: `_shared/调研文献/LAMP设计规范_文献调研.md` (R-0012)

2. **LAMP 设计规范文档** (R-0012)
   - 核心文献: Notomi 2000, Nagamine 2002, Tomita 2008, Sahin 2011
   - 规范: F3/B3 Tm 60-65, F1c/B1c 65-70, F1c-F2 差 ≥2°C (降级 1.5°C)
   - 动力学: T_d = Sa/(ξ·Se), Tt = -4.60*log10(N0) + 39.2

3. **Skill 沉淀** (PLAN 步骤 8 完成)
   - 创建 skill: `lamp-rpa-primer-design` (molecular-diagnostics 类别)
   - 包含: LAMP/RPA 设计算法 + Tt 仿真 + 多病原体验证 + 文献依据
   - 路径: `~/.hermes/profiles/antigen-expression/skills/molecular-diagnostics/lamp-rpa-primer-design/`

## 专题完成状态

| PLAN 步骤 | 状态 |
|:---|:---:|
| 1. 建立专题目录和文件体系 | ✓ |
| 2. 并行文献调研 | ✓ (LAMP 综合, RPA 完成, 动力学完成) |
| 3. 并行真实数据源收集 | ✓ (RPA 完成, LAMP 间接获得) |
| 4. 汇总调研结果, 确定改进优先级 | ✓ |
| 5. 实施算法改进 | ✓ (LAMP + RPA + Tt 模型) |
| 6. 多序列多文献验证 | ✓ (5 病原体, 2 文献基准) |
| 7. 出专题 REPORT | ✓ (R-0010) |
| 8. 沉淀为 skill | ✓ |

## 专题最终成果

| 维度 | v1 基线 | v2 改进后 |
|:---|:---|:---|
| LAMP 设计成功率 | 40% (2/5) | 100% (5/5) |
| RPA 设计成功率 | 92% | 100% (5/5) |
| RPA Tt RMSE | 2.50 min | 0.51 min |
| LAMP Tt RMSE | 1.67 min | 0.56 min |
| 文献依据 | — | 26 篇核心文献 |

## 数据/文件
- `_shared/调研文献/LAMP设计规范_文献调研.md` (R-0012)
- Skill: `lamp-rpa-primer-design` (molecular-diagnostics)
- `_shared/INDEX.md` (更新, 共 12 个 R-NNNN)

专题 8 个步骤全部完成。
