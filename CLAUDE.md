# qPCR 引物探针设计 + 仿真验证平台

## 项目目标
构建一个 Web 可视化平台：给定核酸序列（或 GenBank 检索号 + 病原信息）→ 自动设计 qPCR 引物探针 → 模拟仿真扩增过程 → 给出优化建议 → 形成优化循环。

## 技术栈
- Python 3.13 (/usr/local/bin/python3)
- Streamlit（Web 界面 + 交互可视化）
- primer3-py 2.3.0（引物设计 + 热力学计算）
- biopython 1.87（序列处理 + GenBank 检索 + BLAST）
- numpy/scipy/pandas（数值计算）
- matplotlib/plotly（扩增曲线可视化）
- scikit-learn（效率预测模型，可选）

## 环境验证
```
primer3-py: 2.3.0 ✓ (import primer3; primer3.__version__)
biopython: 1.87 ✓
numpy/scipy/pandas/matplotlib ✓
```
安装 Streamlit: pip install streamlit plotly

## 功能模块规格

### 模块1: 序列输入界面 (sidebar)
- 文本框：粘贴核酸序列（FASTA 或纯序列）
- GenBank 号检索：输入 accession number → Biopython Entrez 自动下载
- 病原信息输入：病原名称、基因名称、靶标区域描述
- 参数设置区：
  - 检测类型：TaqMan / SYBR Green / Molecular Beacon
  - 产物大小范围（默认 70-150bp）
  - 引物 Tm 范围（默认 58-62°C）
  - 探针 Tm 范围（默认 68-72°C）

### 模块2: 引物探针设计
- 调用 primer3.design_primers() 设计多组候选
- 每组输出：正向/反向引物序列、Tm、GC%、探针序列、产物大小、penalty
- 可选：指定靶标区域（SEQUENCE_TARGET）
- 可选：排除区域（SEQUENCE_EXCLUDED_REGION）
- 结果表格展示，按 penalty 排序

### 模块3: 热力学仿真验证（核心）
- 使用 primer3.calc_tm/hairpin/homodimer/heterodimer
- 参数化不同反应体系（5种预设试剂盒 + 自定义）：
  - ABI TaqMan Fast (Mg2+=3mM, dNTP=0.8mM, primer=200nM)
  - QIAGEN QuantiFast (Mg2+=2.5mM, dNTP=0.6mM, primer=300nM)
  - Takara Premix Taq HS (Mg2+=4mM, dNTP=1.0mM, primer=100nM)
  - Roche LC480 (Mg2+=3.5mM, dNTP=0.7mM, primer=250nM)
  - 自定义体系（滑块调整 Mg2+/dNTP/引物浓度）
- 输出：不同体系下的 Tm 变化表、发夹/二聚体评估表
- 可视化：不同体系 Tm 对比柱状图

### 模块4: 扩增曲线仿真（核心）
- 指数扩增模型: N(C) = N0 * (1+E)^C
- 参数：
  - 起始模板量（拷贝数，10^1 ~ 10^8 梯度，或自定义）
  - 扩增效率（从模型3预测，或手动输入 50-100%）
  - 循环数（默认 40）
  - 荧光阈值（默认 0.1，可调）
  - 噪声水平（模拟实际实验波动）
- 输出：
  - 多条扩增曲线（不同模板量的 S 型曲线）
  - Ct 值表
  - 标准曲线（log浓度 vs Ct），斜率/R²/计算效率
  - LOD 估计（95% 检出概率的最低拷贝数）
  - 线性范围

### 模块5: 优化循环（核心）
- 当仿真结果不理想时（效率<90%、LOD太高、二级结构问题），自动给出建议：
  - "引物 Tm 差异过大，建议调整引物长度"
  - "探针发夹 Tm 过高（53°C），建议缩短或调整 GC 含量"
  - "Mg2+ 浓度偏低导致 Tm 偏低，建议使用 Takara 体系"
  - "检测下限偏高，建议增加引物浓度或优化探针"
- 根据建议自动调整参数重新设计和仿真
- 展示优化前后的对比
- 优化维度：
  1. 引物序列（换一组候选）
  2. 引物/探针长度
  3. 反应体系（Mg2+/dNTP/引物浓度）
  4. 标记体系（TaqMan → 分子信标 → PNA 探针等）
  5. 探针体系（LNA/修饰碱基/MGB 等）

### 模块6: 综合报告
- 最优引物探针方案汇总
- 理化参数表（MW、pI、GC、Tm）
- 不同体系下的性能预测对比
- 推荐反应体系
- 预期 LOD 和线性范围
- 可导出为文本报告

## 文件结构
```
qPCR_引物探针设计仿真平台/
├── app.py                 # Streamlit 主应用
├── core/
│   ├── __init__.py
│   ├── sequence_input.py  # 序列获取（手动/GenBank）
│   ├── primer_design.py   # 引物探针设计引擎
│   ├── thermo_sim.py       # 热力学仿真引擎
│   ├── amplification.py   # 扩增曲线模拟引擎
│   ├── optimizer.py        # 优化循环引擎
│   └── report.py           # 报告生成
├── data/
│   └── master_mix_presets.json  # 试剂盒预设参数
├── CLAUDE.md
└── requirements.txt
```

## 代码规范
- 所有中文注释
- 函数有 docstring
- 类型标注
- 每个模块可独立测试
- Streamlit 界面用中文标签
