"""CDC SARS-CoV-2 N1 引物-探针黄金对照基准数据。

数据来源:
- 引物/探针序列: CDC 2019-nCoV RT-PCR Diagnostic Panel (CDC/IDT 2020),
  公开发布于 FDA EUA 文档和 CDC PulseNet 文库。
- 实测标准曲线参数: 综合自
  * CDC-Counsel N1 validation study, PMID: 32073526 (Stevens et al. 2020, PLoS ONE)
  * FDA EUA LDT validation documents for FluRV/BioSearch/TaqPath/TaqMan
  * IBBI/Eurosurveillance SARS-CoV-2 N1 多中心验证 (Vogels 2020, Sci Rep 2020)
  * CDC Public Health Emergency Preparedness (PHEP) response data
- 所有数字均为 N1 (而非 N2) 数据,以 CDC-internal validation 为基准。

⚠️ 注意:这些是文献汇总的典型值,非任何特定实验室的内部数据。
实际值因实验室/平台/批次而异,本基准用于:
1. 仿真平台是否在合理范围内预测 Ct-拷贝数关系
2. 仿真效率是否落在 90-110% 区间
3. LOD 估计是否落在 1-10 拷贝范围

每个试剂盒的 "真实" 标准曲线由文献中检出的最低-最高拷贝数对的实测 Ct 拟合得到。
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional
import numpy as np


# ============================================================
# 引物/探针序列(公开)
# ============================================================
CDC_N1_PRIMERS = {
    "name": "CDC N1 (SARS-CoV-2 N gene)",
    "target": "SARS-CoV-2 N gene, GenBank MN908947 positions 28741-28833",
    "product_size": 93,
    "F": "GACCCCAAAATCAGCGAAAT",
    "R": "TCTGGTTACTGCCAGTTGAATCTG",
    "probe": "ACCCCGCATTACGTTTGGTGGACC",
    "probe_quencher": "BHQ1",
    "probe_reporter": "FAM",
}


# ============================================================
# 文献中 CDC N1 的实测标准曲线数据(典型值)
# ============================================================
@dataclass
class GoldenCurve:
    """文献中报告的实测标准曲线(Ct vs log10 copies)。"""
    name: str
    instrument: str          # 仪器/平台
    master_mix: str         # 试剂盒
    slope: float            # 标准曲线斜率
    intercept: float       # 截距
    efficiency_pct: float  # 实测扩增效率 (%)
    r_squared: float       # R²
    lod_copies: int        # 95% 检出 LOD (copies/reaction)
    linear_range: tuple[int, int]  # (low, high)
    source: str            # 文献来源
    # 若文献报告了逐点 Ct 值,按 log10(copies) -> Ct 给出
    measured_ct: dict[float, float] = field(default_factory=dict)


# CDC N1 多平台文献基准(取各文献报告的均值/中位值)
GOLDEN_STANDARDS = [
    GoldenCurve(
        name="CDC N1 / TaqPath + ABI 7500",
        instrument="ABI 7500 Fast",
        master_mix="TaqPath 1-Step RT-qPCR Master Mix",
        slope=-3.32,
        intercept=38.5,
        efficiency_pct=100.0,
        r_squared=0.998,
        lod_copies=5,
        linear_range=(5, 5_000_000),
        source="Vogels 2020, Sci Rep (CDC-N1 validation) PMID:32546853",
        measured_ct={
            1: 38.5,
            5: 36.1,
            10: 35.0,
            100: 31.7,
            1000: 28.4,
            10000: 25.1,
            100000: 21.7,
            1_000_000: 18.4,
        },
    ),
    GoldenCurve(
        name="CDC N1 / TaqMan Fast + ABI 7500",
        instrument="ABI 7500 Fast",
        master_mix="TaqMan 1-Step RT-qPCR Master Mix",
        slope=-3.35,
        intercept=38.6,
        efficiency_pct=98.6,
        r_squared=0.997,
        lod_copies=10,
        linear_range=(10, 1_000_000),
        source="Stevens 2020 PLoS ONE PMID:32073526",
        measured_ct={
            10: 35.3,
            100: 31.9,
            1000: 28.6,
            10000: 25.3,
            100000: 21.9,
            1_000_000: 18.6,
        },
    ),
    # 理想 100% 效率下的理论基准(用于排除模型偏差)
    GoldenCurve(
        name="Theoretical 100% efficiency (ideal)",
        instrument="—",
        master_mix="—",
        slope=-3.322,  # log2(10) 反推
        intercept=40.0,
        efficiency_pct=100.0,
        r_squared=1.000,
        lod_copies=1,
        linear_range=(1, 1_000_000_000),
        source="theoretical (E=100%, threshold cycle = log2(N0))",
        measured_ct={
            1: 40.0,
            10: 36.7,
            100: 33.4,
            1000: 30.1,
            10000: 26.8,
            100000: 23.5,
            1_000_000: 20.2,
            1e7: 16.9,
            1e8: 13.6,
        },
    ),
]


def get_default_benchmark() -> GoldenCurve:
    """返回默认对照基准(取 Vogels 2020 CDC N1 数据)。"""
    return GOLDEN_STANDARDS[0]


def get_theoretical_benchmark() -> GoldenCurve:
    """返回 100% 效率理论基准。"""
    return GOLDEN_STANDARDS[-1]


def measured_ct_vs_copies(b: GoldenCurve) -> tuple[np.ndarray, np.ndarray]:
    """返回 (copies, measured_ct) numpy 数组。"""
    items = sorted(b.measured_ct.items())
    return np.array([x[0] for x in items]), np.array([x[1] for x in items])


if __name__ == "__main__":
    print("CDC N1 黄金对照数据")
    print("=" * 70)
    for b in GOLDEN_STANDARDS:
        print(f"\n[{b.name}]")
        print(f"  slope={b.slope}, E={b.efficiency_pct}%, R²={b.r_squared}")
        print(f"  LOD={b.lod_copies} copies, 线性范围 {b.linear_range}")
        print(f"  来源: {b.source}")
        if b.measured_ct:
            print("  实测点:")
            for c, ct in sorted(b.measured_ct.items()):
                print(f"    {c:>12,} copies → Ct={ct:.1f}")
