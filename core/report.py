"""综合报告生成模块:汇总最优方案、可导出文本。"""
from __future__ import annotations

from datetime import datetime
from typing import Optional

from core.amplification import StandardCurveResult
from core.optimizer import OptimizationResult, OptimizationSuggestion
from core.primer_design import PrimerPairResult
from core.thermo_sim import ThermoResult


def _physchem(seq: str) -> dict:
    """基础理化参数估算(长度、GC%、MW、A260)。"""
    if not seq:
        return {"length": 0, "gc": 0.0, "mw": 0.0, "a260": 0.0}
    length = len(seq)
    gc = (seq.count("G") + seq.count("C")) / length * 100.0
    # 平均碱基分子量(dNMP):A=313.21, T=304.20, C=289.18, G=329.21 - 18 (脱水)
    mw_per_base = {
        "A": 313.21 - 18.0,
        "T": 304.20 - 18.0,
        "C": 289.18 - 18.0,
        "G": 329.21 - 18.0,
    }
    mw = sum(mw_per_base.get(b, 300.0) for b in seq.upper())
    # A260 估算(单链):1 μg ≈ 33 nmol × 1000 bp / length
    a260_per_ug = 33.0 * 1000.0 / length  # pmol/μg 简化
    return {
        "length": length,
        "gc": gc,
        "mw": mw,
        "a260_per_ug": a260_per_ug,
    }


def generate_report(
    sequence: str,
    organism: str,
    gene: str,
    detection_type: str,
    best_pair: PrimerPairResult,
    best_thermo: ThermoResult,
    best_system_name: str,
    std_curve: Optional[StandardCurveResult] = None,
    lod_info: Optional[dict] = None,
    suggestions: Optional[list[OptimizationSuggestion]] = None,
) -> str:
    """生成可导出的文本综合报告。

    Args:
        sequence: 模板序列。
        organism: 病原名称。
        gene: 基因名。
        detection_type: 检测类型。
        best_pair: 最优引物对。
        best_thermo: 最优体系下的热力学评估。
        best_system_name: 推荐体系名称。
        std_curve: 标准曲线结果(可空)。
        lod_info: LOD 信息(可空)。
        suggestions: 优化建议(可空)。

    Returns:
        文本报告字符串。
    """
    lines: list[str] = []
    sep = "=" * 70

    lines.append(sep)
    lines.append("        qPCR 引物探针设计方案报告")
    lines.append(sep)
    lines.append(f"生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append("")
    lines.append("【目标信息】")
    lines.append(f"  病原体: {organism or '-'}")
    lines.append(f"  靶基因: {gene or '-'}")
    lines.append(f"  检测类型: {detection_type}")
    lines.append(f"  模板序列长度: {len(sequence)} bp")
    lines.append("")

    lines.append("【最优引物探针】")
    lines.append(f"  正向引物(5'→3'): {best_pair.forward_seq}")
    lines.append(f"    长度: {len(best_pair.forward_seq)} nt    "
                 f"Tm: {best_pair.forward_tm:.2f}°C    GC: {best_pair.forward_gc:.1f}%")
    pc = _physchem(best_pair.forward_seq)
    lines.append(f"    MW: {pc['mw']:.0f} g/mol")
    lines.append(f"  反向引物(5'→3'): {best_pair.reverse_seq}")
    lines.append(f"    长度: {len(best_pair.reverse_seq)} nt    "
                 f"Tm: {best_pair.reverse_tm:.2f}°C    GC: {best_pair.reverse_gc:.1f}%")
    pc = _physchem(best_pair.reverse_seq)
    lines.append(f"    MW: {pc['mw']:.0f} g/mol")
    if best_pair.probe_seq:
        lines.append(f"  探针(5'→3'): {best_pair.probe_seq}")
        lines.append(f"    长度: {len(best_pair.probe_seq)} nt    "
                     f"Tm: {best_pair.probe_tm:.2f}°C    GC: {best_pair.probe_gc:.1f}%")
        pc = _physchem(best_pair.probe_seq)
        lines.append(f"    MW: {pc['mw']:.0f} g/mol")
    lines.append(f"  产物大小: {best_pair.product_size} bp")
    lines.append(f"  primer3 惩罚值: {best_pair.penalty:.4f}")
    lines.append("")

    lines.append("【推荐反应体系】")
    lines.append(f"  体系: {best_system_name}")
    lines.append(f"  评估备注: {best_thermo.notes}")
    lines.append(f"  引物 Tm 差异: {best_thermo.tm_difference:.2f}°C")
    if best_pair.probe_seq:
        lines.append(f"  探针发夹 Tm: {best_thermo.probe_hairpin_tm:.2f}°C")
    lines.append(f"  异源二聚体 ΔG: {best_thermo.heterodimer_dg:.2f} kcal/mol")
    lines.append("")

    if std_curve is not None:
        lines.append("【预期扩增性能】")
        lines.append(f"  斜率: {std_curve.slope:.3f}")
        lines.append(f"  截距: {std_curve.intercept:.3f}")
        lines.append(f"  R²: {std_curve.r_squared:.4f}")
        lines.append(f"  扩增效率: {std_curve.efficiency_pct:.2f}%")
        lines.append(f"  线性范围: log10 = {std_curve.linear_range[0]:.2f} ~ "
                     f"{std_curve.linear_range[1]:.2f} "
                     f"(约 {10**std_curve.linear_range[0]:.0f} ~ "
                     f"{10**std_curve.linear_range[1]:.0f} 拷贝)")
        lines.append("")

    if lod_info is not None and lod_info.get("lod_copies") is not None:
        lines.append("【检测下限 (LOD)】")
        lines.append(f"  估计 LOD: {lod_info['lod_copies']:.1f} 拷贝 "
                     f"(log10 = {lod_info['lod_log10']:.2f})")
        lines.append(f"  检出概率: {lod_info['detection_probability']*100:.0f}%")
        lines.append(f"  备注: {lod_info['note']}")
        lines.append("")

    if suggestions:
        lines.append("【优化建议】")
        for i, s in enumerate(suggestions, 1):
            lines.append(f"  {i}. [{s.priority.upper()}] {s.dimension}: {s.issue}")
            lines.append(f"     → 建议: {s.action}")
        lines.append("")

    lines.append(sep)
    lines.append("              报告结束")
    lines.append(sep)
    return "\n".join(lines)


def report_from_optimization(
    sequence: str,
    organism: str,
    gene: str,
    detection_type: str,
    opt_result: OptimizationResult,
) -> str:
    """从 OptimizationResult 直接生成报告(使用最优 after_pair,若无则用 before)。"""
    pair = opt_result.after_pair or opt_result.before_pair
    thermo = opt_result.after_thermo or opt_result.before_thermo
    std = opt_result.after_std or opt_result.before_std
    return generate_report(
        sequence=sequence,
        organism=organism,
        gene=gene,
        detection_type=detection_type,
        best_pair=pair,
        best_thermo=thermo,
        best_system_name=opt_result.selected_system,
        std_curve=std if isinstance(std, StandardCurveResult) else None,
        suggestions=opt_result.suggestions,
    )