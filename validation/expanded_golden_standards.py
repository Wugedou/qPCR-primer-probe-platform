"""扩展黄金对照：FluA (Influenza A) qPCR 标准曲线数据

数据来源：
- Zhang et al. 2019, Diagn Microbiol Infect Dis, PMID:31130238
- Di Trani et al. 2006, BMC Infect Dis, PMID:16725022
- CDC/WHO 2009 pandemic InfA rRT-PCR protocol
"""

from dataclasses import dataclass, field

@dataclass
class GoldenCurve:
    name: str
    instrument: str
    master_mix: str
    slope: float
    intercept: float
    efficiency_pct: float
    r_squared: float
    lod_copies: int
    linear_range: tuple
    source: str
    measured_ct: dict = field(default_factory=dict)
    # 引物探针信息（若有）
    primer_F: str = ""
    primer_R: str = ""
    probe: str = ""
    product_size: int = 0

# FluA 引物序列
FLUA_PRIMERS = {
    "name": "CDC InfA (Influenza A M gene)",
    "F": "GACCRATCCTGTCACCTCTGAC",
    "R": "AGGGCATTYTGGACAAAKCGTCTA",
    "probe": "FAM-TGCAGTCCTCGCTCACTGGGCACG-BHQ1",
    "product_size": 106,
}

# ============================================================
# FluA 黄金对照曲线
# ============================================================
FLUA_GOLDEN_STANDARDS = [
    GoldenCurve(
        name="FluA / Zhang 2019 IAV singleplex (CDC InfA primers)",
        instrument="ABI QuantStudio",
        master_mix="AgPath-ID One-Step RT-PCR Kit",
        slope=-3.4454,
        intercept=34.7527,
        efficiency_pct=95.07,
        r_squared=0.998,
        lod_copies=30,
        linear_range=(0.27, 27000),
        source="Zhang 2019 PMID:31130238, singleplex",
        primer_F=FLUA_PRIMERS["F"],
        primer_R=FLUA_PRIMERS["R"],
        probe=FLUA_PRIMERS["probe"],
        product_size=FLUA_PRIMERS["product_size"],
        measured_ct={
            27_000: 19.37,
            2_700: 23.11,
            270: 26.36,
            27: 29.83,
            2.7: 33.15,
            0.27: 36.77,
        },
    ),
    GoldenCurve(
        name="FluA / Zhang 2019 IAV 5-plex (multiplex)",
        instrument="ABI QuantStudio",
        master_mix="AgPath-ID One-Step RT-PCR Kit",
        slope=-3.5860,
        intercept=35.5942,
        efficiency_pct=90.05,
        r_squared=0.998,
        lod_copies=30,
        linear_range=(0.27, 27000),
        source="Zhang 2019 PMID:31130238, 5-plex",
        measured_ct={
            27_000: 19.80,
            2_700: 23.18,
            270: 26.74,
            27: 30.49,
            2.7: 34.35,
            0.27: 37.45,
        },
    ),
    GoldenCurve(
        name="FluA / Di Trani 2006 (M gene, MGB probe)",
        instrument="ABI 7700",
        master_mix="QuantiTect Probe RT-PCR Kit",
        slope=-3.43,
        intercept=40.8,  # estimated from LOD constraint
        efficiency_pct=95.7,
        r_squared=0.998,
        lod_copies=5,
        linear_range=(5, 5e8),
        source="Di Trani 2006 PMID:16725022, M gene",
        # No full Ct table — estimated from slope
        measured_ct={},
    ),
]

# ============================================================
# 当前验证覆盖汇总
# ============================================================
VALIDATION_COVERAGE = {
    "CDC_N1 (SARS-CoV-2)": {
        "datasets": 2,
        "has_ct_table": True,
        "ΔCt_range": "0.1–1.5",
        "status": "✅ validated (PLAN-001)",
    },
    "FluA (Influenza A)": {
        "datasets": 2,
        "has_ct_table": True,
        "ΔCt_range": "待验证",
        "status": "🔍 PLAN-002 进行中",
    },
    "Dengue virus": {
        "datasets": 0,
        "has_ct_table": False,
        "status": "⏸ 文献已检索(6 PMIDs)，逐点Ct待全文提取",
    },
    "HBV": {
        "datasets": 0,
        "has_ct_table": False,
        "status": "⏸ 文献已检索(3 PMIDs)，逐点Ct待全文提取",
    },
    "MTB": {
        "datasets": 0,
        "has_ct_table": False,
        "status": "⏸ 文献已检索(8 PMIDs)，逐点Ct待全文提取",
    },
}
