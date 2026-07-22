"""qPCR/RPA/LAMP 引物探针设计仿真平台核心模块。

包含序列获取、qPCR/RPA/LAMP 引物设计、热力学仿真、扩增曲线仿真、
优化循环和报告生成。
"""
from core.sequence_input import fetch_sequence_manual, fetch_sequence_genbank
from core.primer_design import design_primer_pairs, design_single_pair
from core.thermo_sim import evaluate_thermodynamics, compare_master_mixes
from core.amplification import simulate_amplification, compute_ct, fit_standard_curve, estimate_lod
from core.amplification import simulate_rpa, estimate_rpa_tt, simulate_lamp, estimate_lamp_tt
from core.amplification import RPAIsothermalParams, LAMPIsothermalParams
from core.rpa_design import design_rpa_primers, RPADesignParameters, RPAPrimerPair, RPAProbe
from core.lamp_design import design_lamp_primers, LAMPDesignParameters, LAMPPrimerSet
from core.optimizer import OptimizationEngine, OptimizationResult
from core.report import generate_report

__all__ = [
    "fetch_sequence_manual",
    "fetch_sequence_genbank",
    "design_primer_pairs",
    "design_single_pair",
    "evaluate_thermodynamics",
    "compare_master_mixes",
    "simulate_amplification",
    "simulate_rpa",
    "estimate_rpa_tt",
    "simulate_lamp",
    "estimate_lamp_tt",
    "RPAIsothermalParams",
    "LAMPIsothermalParams",
    "design_rpa_primers",
    "RPADesignParameters",
    "RPAPrimerPair",
    "RPAProbe",
    "design_lamp_primers",
    "LAMPDesignParameters",
    "LAMPPrimerSet",
    "compute_ct",
    "fit_standard_curve",
    "estimate_lod",
    "OptimizationEngine",
    "OptimizationResult",
    "generate_report",
]
