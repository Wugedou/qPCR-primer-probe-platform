# qPCR/RPA/LAMP 引物探针设计 + 仿真验证平台 — 修订规格

## 需要修复的问题

### 问题1: RPA 和 LAMP 设计缺失
当前只有 qPCR (TaqMan/SYBR) 设计。需要新增：

- **RPA (重组酶聚合酶扩增)**:
  - 引物长度: 30-35 bp (比 PCR 长很多)
  - 产物大小: 100-500 bp (比 qPCR 大)
  - 不需要探针(基础 RPA)，但 exo-probe/nfo-probe 可选
  - 正向/反向引物 Tm 无严格限制(等温扩增)
  - 避免引物 3' 端发夹和二聚体
  - 设计引擎: 用 primer3 但参数调整(长引物、大产物)
  
- **LAMP (环介导等温扩增)**:
  - 需要 4 条核心引物: FIP, BIP, F3, B3
  - FIP = F1c + TTTT linker + F2 (长度 40-45bp)
  - BIP = B1c + TTTT linker + B2 (长度 40-45bp)
  - F3/B3: 18-22bp 外侧引物
  - 可选: Loop F/B, Stem primers
  - 产物: 梯状条带(100-1000bp+)
  - 设计引擎: 需要识别 6 个区域(F3-F2-F1-B1c-B2c-B3c)
  - 可用 primer3 设计 F3/B3，FIP/BIP 需要特殊组装逻辑

### 问题2: 热力学和优化循环改参数看不出差异

**根因已确认**: `core/thermo_sim.py` 中 `calc_hairpin`/`calc_homodimer`/`calc_heterodimer` 没有传入盐浓度参数 (mv_conc/dv_conc/dntp_conc/dna_conc)，导致不同体系的二级结构结果完全相同。

**修复方案**:
1. `thermo_sim.py`: `_hairpin`/`_homodimer`/`_heterodimer` 函数传入 `mv_conc`/`dv_conc`/`dntp_conc`/`dna_conc`
2. `optimizer.py`: 优化循环需要真正改变参数后产生不同的仿真结果
   - 改变 Mg2+ 浓度后重算 Tm/发夹/二聚体，体现差异
   - 换引物候选后重算所有热力学，体现差异
   - 优化前后的数值对比表要有实际变化

## 具体修改要求

### 修改 core/thermo_sim.py
```python
# 当前(错误):
def _hairpin(seq):
    res = primer3.calc_hairpin(seq)
    
# 修复后:
def _hairpin(seq, tm_kw):
    res = primer3.calc_hairpin(seq, **tm_kw)
    # tm_kw 包含 mv_conc, dv_conc, dntp_conc, dna_conc
```
同理修复 _homodimer 和 _heterodimer。

### 新增 core/rpa_design.py
RPA 引物设计引擎:
- `DesignParameters` 增加 detection_type="rpa"
- 引物长度 30-35bp
- 产物 100-500bp
- 不需要 Tm 严格约束(等温 37-42°C)
- 检查 3' 端发夹/二聚体(对 RPA 极重要)
- 可选 exo-probe 设计(46-52bp, 含 THF/abhCITE 修饰位点)

### 新增 core/lamp_design.py
LAMP 引物设计引擎:
- 识别靶序列上的 6 个区域: F3-F2-F1-(loop)-B1c-B2c-B3c
- 设计 F3/B3: 用 primer3, 18-22bp
- 组装 FIP = F1c(反补) + TTTT + F2
- 组装 BIP = B1c(反补) + TTTT + B2
- 可选 Loop F/B 引物
- Tm 约束: FIP/BIP 各区段 Tm 60-65°C, F3/B3 Tm 55-60°C

### 修改 app.py
- sidebar 检测类型下拉框增加 "RPA" 和 "LAMP" 选项
- Tab 1 根据 detection_type 调用不同的设计引擎
- Tab 2 热力学仿真: 对 RPA 显示引物二聚体/发夹(不需要 Tm 对比), 对 LAMP 显示 FIP/BIP 组装图
- Tab 3 扩增仿真: 
  - qPCR: S 型荧光曲线 + Ct + 标准曲线 (已有)
  - RPA: 时间-荧光曲线(等温, 5-20分钟), 不同模板量
  - LAMP: 时间-荧光曲线(等温, 15-60分钟), 可选 turbidity/fluorescence
- Tab 4 优化循环: 修复后体现参数变化差异

### 修改 core/optimizer.py
- 优化循环每次迭代后,用新参数重算热力学,展示前后数值差异
- 对比表要显示: 旧Tm → 新Tm, 旧发夹dG → 新发夹dG, 旧效率 → 新效率

## 文件结构(更新后)
```
qPCR_引物探针设计仿真平台/
├── app.py                  # Streamlit 主应用(更新)
├── core/
│   ├── __init__.py
│   ├── sequence_input.py   # 不变
│   ├── primer_design.py    # 重构: qPCR 设计(不变) + 路由
│   ├── rpa_design.py       # 新增: RPA 引物设计
│   ├── lamp_design.py      # 新增: LAMP 引物设计
│   ├── thermo_sim.py       # 修复: 传盐参数 + 适配 RPA/LAMP
│   ├── amplification.py    # 更新: 增加 RPA/LAMP 扩增模型
│   ├── optimizer.py        # 修复: 体现参数变化差异
│   └── report.py           # 更新: 支持 RPA/LAMP 报告
├── data/
│   └── master_mix_presets.json
├── CLAUDE.md
└── requirements.txt
```

## 技术约束
- Python 3.13 (/usr/local/bin/python3)
- 已安装: primer3-py 2.3.0, biopython 1.87, streamlit 1.59.1, plotly 6.6.0
- 不要改 width="stretch" (已修复的 use_container_width 替换)
- 所有 UI 中文标签
- 不要问问题,直接实现
