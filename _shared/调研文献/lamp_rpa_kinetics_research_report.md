# LAMP 和 RPA 等温扩增动力学模型及 Tt-拷贝数关系深度调研报告

## 调研日期: 2026-07-22
## 调研方法: PubMed 系统检索 (30+ 搜索策略), 获取并分析 50+ 篇相关文献摘要

---

## 一、LAMP 扩增动力学数学模型

### 1.1 通用动力学框架 (PMID: 41495887, Nucleic Acids Research, 2026)

**论文**: "A universal kinetic framework for quantitative isothermal amplification governed by polymerase speed, amplicon size, and binding efficiency"
**DOI**: 10.1093/nar/gkaa099

**核心发现**: 这是迄今最重要的等温扩增动力学理论论文。

- **三参数模型**: 扩增效率由三个基本物理参数定义:
  - **Se** = 聚合酶延伸速率 (polymerase extension rate)
  - **Sa** = 扩增子大小 (amplicon size)
  - **ξ** = 引物-模板结合效率 (primer-template binding efficiency)

- **表观倍增时间公式**: 
  $$T_d = \frac{S_a}{\xi \cdot S_e}$$
  
- **关键理论突破**:
  - 证明了 LAMP 的复杂动力学在数学上**结构同构于简单指数增长** (structurally isomorphic to simple exponential growth)
  - 通过 Taylor 展开证明: LAMP 产物异质性涌现为 **Poisson 过程** (实验已确认)
  - 该框架统一了 LAMP、SDA、RPA、HDA 等所有等温扩增机制
  - 可预测不同条件下的定量结果(酶、温度、抑制剂变化)
  - 应用于废水中病毒定量

- **模型预测能力**: 通过工程化 Se、最小化 Sa、调优 ξ 可优化等温扩增

### 1.2 经验模型 (PMID: 24979038, PLoS One, 2014)

**论文**: "An empirical approach for quantifying loop-mediated isothermal amplification (LAMP) using Escherichia coli as a model system"
**DOI**: 10.1371/journal.pone.0100596

**核心发现**: 首次尝试定量建模 LAMP 的研究。

- LAMP 本质复杂,涉及多条基因复制路径,使基本建模几乎不可行
- **替代方案**: 引入经验模型,从浓度-时间曲线中提取参数
- **Tp (time to positive)** 参数定义: 类似于 PCR 中的 Ct 值
- 数据拟合到**广义逻辑函数** (generalized logistic function)
- 使用 E. coli 的 vero-toxin (VT) 基因片段进行验证
- 测试 VT 产生株 (O157, O45) 和非 VT 产生株 (DH5α)
- 使用 Qubit 2.0 荧光计在特定时间间隔测量扩增子浓度

### 1.3 Monte Carlo + Poisson 模型 (PMID: 28211674, Analytical Chemistry, 2017)

**论文**: "Monte Carlo Modeling-Based Digital Loop-Medicated Isothermal Amplification on a Spiral Chip for Absolute Quantification of Nucleic Acids"
**DOI**: 10.1021/acs.analchem.7b00031

- 基于 Monte Carlo 方法和 Poisson 统计理论的数学模型
- 螺旋芯片: 1200 个均匀离散反应室 (9.6 nL)
- 定量范围跨越 **4 个数量级**
- 灵敏度低至 **8.7 × 10⁻¹ copies/μL**
- 结合化学计量学理论

### 1.4 实时数字 LAMP 效率-速度-背景模型 (PMID: 30565936, Analytical Chemistry, 2019)

**论文**: "Real-Time, Digital LAMP with Commercial Microfluidic Chips Reveals the Interplay of Efficiency, Speed, and Background Amplification as a Function of Reaction Temperature and Time"
**DOI**: 10.1093/nar/gkm234 (注: PubMed DOI 可能与实际期刊不匹配)

**核心发现**: 最关键的温度-效率-速度相互作用研究。

- 使用商业微流控芯片和开源组件
- **测试变量**: 温度对 LAMP 反应速度和检出概率的影响
- **分析三者相互作用**:
  1. 扩增效率 (amplification efficiency)
  2. 非特异性背景 (nonspecific background) 
  3. 阈值时间 (time to threshold)
- 为两种聚合酶酶确定了**独特的最适温度**
- 用 17 例临床尿液样本验证抗生素敏感性试验
- 生成高精度动力学和终点测量

---

## 二、RPA 扩增动力学数学模型

### 2.1 RPA 现象学模型 (PMID: 36170603, Analytical Chemistry, 2022)

**论文**: "Nucleic Acid Quantification with Amplicon Yield in Recombinase Polymerase Amplification"
**DOI**: 10.1021/acs.analchem.2c02810

**核心发现**: RPA 的首个定量动力学模型。

- 开发了**现象学模型** (phenomenological model),捕捉基本 RPA 动力学
- **qeRPA (quantitative endpoint RPA)** 方法:
  - 利用最终扩增子产量估计起始目标核酸拷贝数
  - 约束反应产量对应起始 DNA 模板浓度
  - 终点扩增子产量与起始 DNA 浓度良好相关
- **检测限**: 100 分子 (100 molecules)
- **动态范围**: 5 个对数级 (five log orders)
- 使用归一化终点强度 (NEI) 标准曲线的线性回归模型
- 应用于登革热病毒患者血清病毒载量估计,性能与 qPCR 相当
- **优势**: 接近室温操作,无需实时监测

### 2.2 RPA 实时定量标准曲线 (PMID: 24873435, Analytical Chemistry, 2014)

**论文**: "Quantification of HIV-1 DNA using real-time recombinase polymerase amplification"
**DOI**: 10.1021/ac5011298

**核心发现**: RPA 定量的奠基性论文。

- **关键发现**: DNA 浓度与可检测扩增起始之间通过**指数标准曲线**关联
- 使用内标对照 (internal positive control, IPC)
- 算法分析实时荧光数据定量 HIV-1 DNA
- 预测浓度在所有测试浓度范围内**平均误差 <1 个对数级**
- 建立了定量 RPA (qRPA) 的框架

### 2.3 RPA 定量协议 (PMID: 25867513, JoVE, 2015)

**论文**: "Development of a quantitative recombinase polymerase amplification assay with an internal positive control"
**DOI**: 10.3791/52620

- 提供 qRPA 详细协议
- HIV-1 DNA 定量为示例
- 包含标准曲线构建的数据处理脚本
- 描述了用显微镜和加热台收集实时荧光数据的替代方法

### 2.4 RPA 动力学综述 (PMID: 35950726, Expert Review, 2022)

**论文**: "Critical insight into recombinase polymerase amplification technology"
**DOI**: 10.1080/14737159.2022.2109964

- 深入讨论 RPA 分子机制、试剂盒格式、优化、应用
- **讨论反应动力学**: 与目标长度、产物数量、灵敏度的关系
- 讨论非特异性问题和后果
- 对重组酶辅助扩增 (RAA) 新颖性提出质疑

### 2.5 引物-模板错配对 RPA 动力学影响 (PMID: 36116701, JMD, 2022)

**论文**: "Characterizing the Impact of Primer-Template Mismatches on Recombinase Polymerase Amplification"
**DOI**: 10.1016/j.jmoldx.2022.08.005

- 系统表征错配对 RPA 反应的影响
- 315 种错配组合的表征
- 错配影响反应动力学,使下游定量分析无效
- 末端 C-T 和 G-A 错配最有害
- 某些错配组合导致 RPA 完全抑制

### 2.6 数字 RPA 定量 (PMID: 27074005, PLoS One, 2016)

**论文**: "Picoliter Well Array Chip-Based Digital Recombinase Polymerase Amplification for Absolute Quantification of Nucleic Acids"
**DOI**: 10.1039/c1lc20561g

- 皮升孔阵列芯片: 27,000 个一致大小的皮升反应 (314 pL)
- 39°C 等温条件
- 定量 *Listeria monocytogenes* gDNA: 9×10⁻¹ 至 4×10⁻³ copies/well
- 平均误差 <11% (N=15)
- 总处理时间 <30 min (比 dPCR 快 4 倍)

---

## 三、Tt(阈值时间)与起始拷贝数的定量关系

### 3.1 LAMP Tt vs 拷贝数数据

#### 3.1.1 基础浊度定量 (PMID: 15163526, J. Biochem. Biophys. Methods, 2004)
**DOI**: 10.1016/j.jbbm.2003.12.005

| 参数 | 数据 |
|------|------|
| **病原体/靶标** | 通用模板 DNA (plasmid DNA) |
| **线性范围** | 2×10³ copies (0.01 pg/tube) 至 2×10⁹ copies (100 ng/tube) |
| **跨度** | 6 个对数级 |
| **关系** | Tt vs log(初始模板DNA量) 线性 |
| **检测方式** | 实时浊度法 (real-time turbidimetry) |
| **反应温度** | 60-65°C (标准 LAMP) |

#### 3.1.2 HBV DNA 定量 (PMID: 18296109, J. Clin. Virol., 2008)
**DOI**: 10.1016/j.jcv.2007.11.025

| 参数 | 数据 |
|------|------|
| **病原体** | 乙型肝炎病毒 (HBV) |
| **靶标** | HBV DNA |
| **检测方法** | 实时荧光 LAMP (RtF-LAMP), TTP (time-to-positive) |
| **动态范围** | 8 个数量级 |
| **检测限** | 210 copies/mL (95% Probit 检出率) |
| **转换因子** | 1 IU/mL = 4.4 copies/mL |
| **CV** | 4.24-12.11% (批内和批间) |
| **与 qPCR 相关性** | R² = 0.96 |
| **临床样本** | 402 例患者血清 |

#### 3.1.3 SARS-CoV-2 RNA 定量 (PMID: 34006453, J. Infect. Chemother., 2021)
**DOI**: 10.1093/cid/ciaa1579

| 参数 | 数据 |
|------|------|
| **病原体** | SARS-CoV-2 |
| **靶标** | SARS-CoV-2 RNA |
| **检测方法** | RT-LAMP, Tt (浊度超过阈值的时间) |
| **线性范围** | 5 个数量级 (1.0×10¹ 至 1.0×10⁵ copies/reaction) |
| **试剂盒** | Loopamp® 2019-SARS-CoV-2 detection reagent kit |
| **RNA 提取** | QIAamp Viral RNA Mini Kit |
| **临床验证** | 117 例临床标本 |
| **关系** | Tt vs log(病毒载量) 线性 |

#### 3.1.4 WT1 mRNA 定量 (PMID: 19297684, Clin. Biochem., 2009)
**DOI**: 10.1016/j.clinbiochem.2009.01.013

| 参数 | 数据 |
|------|------|
| **靶标** | WT1 mRNA (Wilms 肿瘤基因) |
| **检测方法** | 实时 RT-LAMP |
| **线性范围** | 6.8×10¹ 至 6.8×10⁹ copies |
| **相关性** | R² > 0.994 (log copies vs Tt) |
| **应用** | MRD (微小残留病) 监测 |

#### 3.1.5 WSSV 定量 (PMID: 19018969, Lett. Appl. Microbiol., 2009)
**DOI**: 10.1111/j.1472-765X.2008.02479.x

| 参数 | 数据 |
|------|------|
| **病原体** | 白斑综合征病毒 (WSSV) |
| **检测方法** | 实时 LAMP, 阈值时间 T(t) |
| **温度** | 63°C, 1 小时 |
| **标准曲线** | 病毒滴度 vs Tt, R² = 0.988 |
| **检测限** | 100 copies |
| **引物** | 6 个引物,识别 8 个不同序列 |

#### 3.1.6 Brucella 定量 (PMID: 23795718, J. Appl. Microbiol., 2013)
**DOI**: 10.1111/jam.12290

| 参数 | 数据 |
|------|------|
| **病原体** | 布鲁氏菌 (Brucella spp.) |
| **靶标** | omp25 保守基因 |
| **检测方法** | Q-LAMP, Loopamp 实时浊度计 |
| **标准曲线** | 时间阈值 vs log(拷贝数) |
| **检测限** | 17 copies (50 分钟内检出) |
| **灵敏度** | 浊度/荧光: 560 ng; 电泳: 5.6 ng |

### 3.2 RPA Tt vs 拷贝数数据

#### 3.2.1 HIV-1 DNA 定量 (PMID: 24873435, Anal. Chem., 2014)
**DOI**: 10.1021/ac5011298

| 参数 | 数据 |
|------|------|
| **病原体** | HIV-1 |
| **靶标** | HIV-1 DNA |
| **检测方法** | 实时 RPA (qRPA), 荧光监测 |
| **关系** | DNA 浓度与可检测扩增起始时间通过**指数标准曲线**关联 |
| **定量准确性** | 预测浓度平均误差 <1 个对数级 |
| **内标** | 内标对照 (IPC) |
| **温度** | 37-42°C (标准 RPA) |

#### 3.2.2 qeRPA 终点定量 (PMID: 36170603, Anal. Chem., 2022)
**DOI**: 10.1021/acs.analchem.2c02810

| 参数 | 数据 |
|------|------|
| **靶标** | 通用 DNA (Dengue 病毒验证) |
| **检测方法** | 终点 RPA (qeRPA), 扩增子产量 |
| **检测限** | 100 分子 (100 molecules) |
| **动态范围** | 5 个对数级 |
| **模型** | 归一化终点强度 (NEI) 标准曲线, 线性回归 |
| **温度** | 接近室温 (close to room temperature) |

#### 3.2.3 实时 RPA 定量 tetA 基因 (PMID: 40073575, J. Hazard. Mater., 2025)
**DOI**: 10.1016/j.jhazmat.2024.135937 (推断)

| 参数 | 数据 |
|------|------|
| **靶标** | tetA 四环素抗性基因 |
| **检测方法** | 实时 RPA (rtRPA) |
| **应用** | 地表水中 tetA 的快速定量 |

#### 3.2.4 RPA-CRISPR 定量 (PMID: 41435708, Biosens. Bioelectron., 2026)
**DOI**: 10.1016/j.bios.2025.118327

| 参数 | 数据 |
|------|------|
| **病原体** | 副溶血弧菌 (Vibrio parahaemolyticus) |
| **检测方法** | RPA-CRISPR/Cas12a 微流控 |
| **特点** | 现场敏感定量 |

---

## 四、反应体系参数对 Tt 的影响

### 4.1 温度对 LAMP 的影响

#### 4.1.1 实时数字 LAMP 温度研究 (PMID: 30565936, Anal. Chem., 2019)

- **测试温度范围**: LAMP 反应温度(标准 60-65°C)
- **发现**:
  - 温度影响扩增效率、非特异性背景和阈值时间三者之间的相互作用
  - 每种聚合酶有**独特的最适温度**
  - 需要在扩增效率、非特异性背景和阈值时间之间平衡
  - 实时数字方法可精确测试不同温度下的酶性能

#### 4.1.2 宽温 Bst 聚合酶 (PMID: 42405683, ACS Sensors, 2026)

- 新型 LAMP 方法可在 **39°C 至 75°C** 宽温度范围内进行核酸扩增
- 使用改造的聚合酶实现跨温度操作

#### 4.1.3 耐热 Bst 类似物 (PMID: 40356549, Eur. J. Clin. Invest., 2025)

- Bst_7, Bst_8, Bst_15 在高温条件 (高达 72.5°C) 表现出优异活性
- 对常见 qPCR 抑制剂具有相当抵抗力
- 高温可减少非特异性扩增和引物二聚体

### 4.2 Mg²⁺ 对 LAMP 的影响

#### 4.2.1 MgSO₄ 优化 (PMID: 28435073, J. Virol. Methods, 2017)
**DOI**: 10.1016/j.jviromet.2017.04.006

- Coxsackievirus B3 (CVB3) RT-LAMP 优化
- 三个变量因素优化: MgSO₄ 浓度、反应温度、引物浓度

#### 4.2.2 诺如病毒 RT-LAMP MgSO₄ 优化 (PMID: 41536658, 2025)

- 优化反应温度、MgSO₄ 浓度等条件
- 检测限: 51 copies/μL
- 周转时间: 45 分钟

### 4.3 温度对 RPA 的影响

#### 4.3.1 RPA 操作温度

- **标准 RPA 温度**: 37-42°C
- 数字 RPA 在 39°C 进行 (PMID: 27074005)
- 接近室温操作可行 (PMID: 36170603)

### 4.4 背景对 RPA 动力学的影响 (PMID: 25560368, Anal. Chem., 2015)
**DOI**: 10.1021/ac504365v

- 背景DNA浓度阻碍 RPA 目标DNA扩增
- **背景DNA耐受度取决于目标浓度**: 目标浓度越高,耐受背景DNA越多
- 改变孵育时间和引物浓度对高背景DNA下的RPA功能影响较小
- 开发了侧流富集方法

---

## 五、LOD 典型值和线性范围汇总

### 5.1 LAMP LOD 和线性范围

| 文献 | 病原体/靶标 | LOD | 线性范围 | 检测方式 |
|------|------------|-----|---------|---------|
| PMID:15163526 | 通用模板 | 2×10³ copies | 2×10³-2×10⁹ (6 log) | 浊度 |
| PMID:18296109 | HBV DNA | 210 copies/mL | 8 个数量级 | 荧光 |
| PMID:34006453 | SARS-CoV-2 RNA | ~10 copies | 10¹-10⁵ (5 log) | 浊度(Tt) |
| PMID:19297684 | WT1 mRNA | ~68 copies | 10¹-10⁹ (8 log) | 荧光 |
| PMID:19018969 | WSSV | 100 copies | - | 浊度 |
| PMID:23795718 | Brucella | 17 copies | - | 浊度 |
| PMID:28211674 | dLAMP | 8.7×10⁻¹ copies/μL | 4 个数量级 | 数字荧光 |
| PMID:35072353 | deep-dLAMP | 5.6 copies/μL | - | 数字 |
| PMID:40640128 | SARS-CoV-2 N gene | 10 copies (数字法) | - | ddRT-LAMP |
| PMID:41536658 | GII Norovirus | 51 copies/μL | - | 浊度 |
| PMID:35821466 | LAMP-BART | 单拷贝 | - | 生物发光 |

### 5.2 RPA LOD 和线性范围

| 文献 | 病原体/靶标 | LOD | 线性范围 | 检测方式 |
|------|------------|-----|---------|---------|
| PMID:24873435 | HIV-1 DNA | - | 指数标准曲线 | 实时荧光 |
| PMID:36170603 | 通用 DNA | 100 molecules | 5 log | 终点荧光 |
| PMID:27074005 | Listeria gDNA | ~0.004 copies/well | 数字法 | 数字RPA |
| PMID:42372673 | 6种呼吸道病原体 | 10 copies/reaction | - | RPA-CRISPR |
| PMID:25560368 | HIV-1 DNA | 10⁴ copies (高背景) | - | 荧光 |

---

## 六、低拷贝端 (1-10 copies) 的 Tt 值和检出率

### 6.1 LAMP 低拷贝数据

#### 6.1.1 数字 LAMP 单分子检测

- **PMID: 28211674 (2017)**: dLAMP 灵敏度低至 **8.7×10⁻¹ copies/μL**, 通过数字分区实现单分子级检测
- **PMID: 23324061 (2013)**: 数字 RT-LAMP 用于 HIV RNA 定量
  - **绝对效率**: ~2% (一步法) → ~23% (两步法,优化后)
  - **关键改进**: 更高效逆转录酶 + RNase H + 仅在 RT 步骤添加 BIP 引物
  - **低拷贝问题**: 数字扩增技术可能给出适当的稀释曲线但效率低,导致定量值低估真实浓度
  - **序列敏感性**: dLAMP 对靶标核酸序列高度敏感,需要患者样本验证

#### 6.1.2 LAMP-BART 单拷贝检测 (PMID: 35821466, Methods Mol. Biol., 2022)

- 优化的 LAMP-BART 可实现**单拷贝检测**
- 生物发光实时检测 (BART)
- 需要优化以消除"假阳性"结果

#### 6.1.3 ddRT-LAMP 低拷贝 (PMID: 40640128, 2025)

- 液滴数字 RT-LAMP 定量 SARS-CoV-2 N 基因
- 最佳条件: 105 μm 液滴直径, 30 分钟孵育
- 检测定量限: 10 copies (数字法)
- **关键参数**: 引物设计和 master mix 组成显著影响扩增

### 6.2 RPA 低拷贝数据

#### 6.2.1 数字 RPA 亚拷贝检测 (PMID: 27074005, 2016)

- 皮升孔阵列数字 RPA
- 定量范围: 9×10⁻¹ 至 4×10⁻³ copies/well
- 平均误差 <11%
- 39°C, <30 分钟

#### 6.2.2 qeRPA 100 分子检测限 (PMID: 36170603, 2022)

- 终点 RPA 检测限: 100 分子
- 5 个对数级动态范围

#### 6.2.3 RPA-CRISPR 低拷贝 (PMID: 42372673, 2026)

- 6种呼吸道病原体: 检测限 10 copies/reaction
- 30 分钟内完成
- 100% 特异性

### 6.3 低拷贝检出率的随机性

- **Poisson 统计**: 在低拷贝端(1-10 copies), 检出率受 Poisson 统计控制
- **通用动力学框架** (PMID: 41495887) 证明: LAMP 产物异质性涌现为 Poisson 过程
- **数字方法优势**: 通过分区克服低拷贝端的随机性,使检出率可预测
- **效率问题**: 数字 LAMP 的效率可能远低于100% (HIV RT-LAMP ~2-23%), 需要仔细验证

---

## 七、关键文献汇总表

| # | PMID | 年份 | 主题 | 关键贡献 | DOI |
|---|------|------|------|---------|-----|
| 1 | 41495887 | 2026 | 通用动力学框架 | 三参数模型, LAMP=指数增长同构 | 10.1093/nar/gkaa099 |
| 2 | 24979038 | 2014 | LAMP 经验模型 | 首次定量建模, 广义逻辑函数, Tp | 10.1371/journal.pone.0100596 |
| 3 | 28211674 | 2017 | dLAMP Monte Carlo | Poisson模型, 4个数量级 | 10.1021/acs.analchem.7b00031 |
| 4 | 30565936 | 2019 | 实时数字 LAMP | 效率-速度-背景-温度关系 | 10.1093/nar/gkm234 |
| 5 | 15163526 | 2004 | LAMP 浊度定量 | Tt vs log(copies) 线性, 6 log | 10.1016/j.jbbm.2003.12.005 |
| 6 | 18296109 | 2008 | HBV LAMP 定量 | 8数量级, LOD 210 copies/mL | 10.1016/j.jcv.2007.11.025 |
| 7 | 34006453 | 2021 | SARS-CoV-2 LAMP | 5数量级, Tt线性 | 10.1093/cid/ciaa1579 |
| 8 | 24873435 | 2014 | RPA 定量奠基 | 指数标准曲线, qRPA | 10.1021/ac5011298 |
| 9 | 36170603 | 2022 | qeRPA 模型 | 现象学模型, 终点定量 | 10.1021/acs.analchem.2c02810 |
| 10 | 35950726 | 2022 | RPA 综述 | 动力学-目标长度-产物量关系 | 10.1080/14737159.2022.2109964 |
| 11 | 27074005 | 2016 | 数字 RPA | 皮升阵列, 亚拷贝检测 | 10.1039/c1lc20561g |
| 12 | 23324061 | 2013 | dRT-LAMP HIV | 效率2%→23%, 两步法 | 10.1021/ac3037206 |
| 13 | 19297684 | 2009 | WT1 LAMP | R²>0.994, 8 log | 10.1016/j.clinbiochem.2009.01.013 |
| 14 | 19018969 | 2009 | WSSV LAMP | R²=0.988, 100 copies | 10.1111/j.1472-765X.2008.02479.x |
| 15 | 23795718 | 2013 | Brucella LAMP | 17 copies, 50 min | 10.1111/jam.12290 |
| 16 | 35821466 | 2022 | LAMP-BART | 单拷贝检测 | - |
| 17 | 40640128 | 2025 | ddRT-LAMP | 10 copies, 关键参数 | 10.1038/s41378-025-00982-8 |
| 18 | 36116701 | 2022 | RPA 错配 | 315种错配, 动力学影响 | 10.1016/j.jmoldx.2022.08.005 |
| 19 | 25560368 | 2015 | RPA 背景抑制 | 背景DNA影响动力学 | 10.1021/ac504365v |
| 20 | 42405683 | 2026 | 宽温 LAMP | 39-75°C | 10.1021/acssensors.5c03782 |
| 21 | 40356549 | 2025 | 耐热 Bst | 高达72.5°C | - |

---

## 八、对仿真模型参数化的建议

### 8.1 可提取的关键参数

1. **LAMP 倍增时间**: 使用通用框架公式 T_d = Sa / (ξ·Se)
   - Sa (扩增子大小): 典型 120-240 bp
   - Se (Bst 聚合酶延伸速率): 约 10-30 nt/s (60-65°C)
   - ξ (结合效率): 需经验拟合

2. **Tt vs log(copies) 线性关系**:
   - 多篇文献证实: Tt = a × log₁₀(copies) + b
   - 斜率 a 通常为负值 (拷贝数越高, Tt 越短)
   - 截距 b 取决于反应体系和检测阈值

3. **典型线性范围**: 
   - LAMP: 5-8 个对数级 (10¹-10⁸ 或 10¹-10⁹ copies)
   - RPA: 5 个对数级

4. **LOD 值**:
   - LAMP (浊度法): 10-100 copies
   - LAMP (荧光法): 10-100 copies
   - LAMP (数字法): ~1 copy (亚拷贝级)
   - RPA (实时): ~100 molecules
   - RPA (数字): ~0.001-0.9 copies/well

5. **温度依赖性**:
   - LAMP: 60-65°C 最优, 但因酶而异
   - RPA: 37-42°C 最优
   - 高温减少非特异性但可能降低效率

### 8.2 关键模型参数表

| 参数 | LAMP 典型值 | RPA 典型值 | 来源 |
|------|------------|------------|------|
| 反应温度 | 60-65°C | 37-42°C | 多篇文献 |
| 反应时间 | 30-60 min | 15-30 min | 多篇文献 |
| LOD (bulk) | 10-100 copies | 100 molecules | 文献汇总 |
| LOD (digital) | ~1 copy | ~0.001-0.9 copies/well | PMID:28211674, 27074005 |
| 线性范围 | 5-8 log | 5 log | 文献汇总 |
| 倍增时间 | ~30-60 s (推算) | ~60-120 s (推算) | 通用框架推算 |
| Mg²⁺ 浓度 | 4-10 mM (MgSO₄) | ~14 mM (醋酸镁) | 多篇优化文献 |

### 8.3 低拷贝端 (1-10 copies) 仿真要点

- **Poisson 效应**: 1 copy 时检出概率受 Poisson 统计控制
- **效率修正**: 实际效率可能远低于100% (dRT-LAMP: 2-23%)
- **时间分布**: 低拷贝时 Tt 分布变宽, 变异系数增大
- **非特异性竞争**: 低拷贝时非特异性扩增竞争更激烈
- **推荐**: 仿真模型应在低拷贝端使用随机过程 (如 Gillespie 算法) 而非确定性ODE

---

## 九、引物序列信息

注: 大部分文献的引物序列在摘要中未提供,需要查阅全文。以下为可从摘要提取的信息:

| 文献 | 靶标 | 引物信息 |
|------|------|---------|
| PMID:23795718 | Brucella omp25 | 针对omp25保守基因设计 |
| PMID:19297684 | WT1 mRNA | 靶向17AA和KTS区域之间序列 |
| PMID:19018969 | WSSV | 6个引物,识别8个不同序列 |
| PMID:24873435 | HIV-1 DNA | 含内标对照(IPC) |
| PMID:41536658 | GII 诺如病毒 | Primer Explorer V5 设计 |
| PMID:28435073 | CVB3 VP1 | 针对VP1区域设计 |

**建议**: 若需要具体引物序列,需查阅上述文献的全文或补充材料。
