"""热力学仿真模块:评估引物/探针在不同反应体系下的 Tm、二聚体、发夹。

基于 primer3.calc_tm / calc_hairpin / calc_homodimer / calc_heterodimer。
支持 5 种预设试剂盒 + 自定义体系。
"""
from __future__ import annotations

import json
from dataclasses import dataclass, asdict, field
from pathlib import Path
from typing import Optional

import primer3


@dataclass
class ReactionSystem:
    """反应体系参数。

    Attributes:
        name: 体系名称。
        mg_conc_mM: Mg2+ 浓度(mM)。
        dntp_conc_mM: dNTP 浓度(mM)。
        primer_conc_nM: 引物浓度(nM)。
        probe_conc_nM: 探针浓度(nM)。
        dna_conc_ng: 模板 DNA 浓度(ng/uL)。
        monovalent_conc_mM: 一价阳离子浓度(mM),primer3 盐校正使用。
    """

    name: str = "自定义体系"
    mg_conc_mM: float = 3.0
    dntp_conc_mM: float = 0.8
    primer_conc_nM: float = 200.0
    probe_conc_nM: float = 100.0
    dna_conc_ng: float = 50.0
    monovalent_conc_mM: float = 50.0


@dataclass
class ThermoResult:
    """单组引物探针在指定体系下的热力学评估结果。

    Attributes:
        system_name: 体系名称。
        forward_tm: 正向引物 Tm。
        reverse_tm: 反向引物 Tm。
        probe_tm: 探针 Tm。
        forward_hairpin_dg: 正向引物发夹 ΔG (kcal/mol)。
        reverse_hairpin_dg: 反向引物发夹 ΔG。
        probe_hairpin_dg: 探针发夹 ΔG。
        forward_homodimer_dg: 正向引物同源二聚体 ΔG。
        reverse_homodimer_dg: 反向引物同源二聚体 ΔG。
        heterodimer_dg: 正向/反向异源二聚体 ΔG。
        probe_homodimer_dg: 探针同源二聚体 ΔG。
        forward_hairpin_tm: 正向引物发夹 Tm(°C)。
        reverse_hairpin_tm: 反向引物发夹 Tm(°C)。
        probe_hairpin_tm: 探针发夹 Tm(°C)。
        tm_difference: 引物 Tm 差(绝对值)。
        acceptable: 是否满足基本可用条件。
        notes: 备注。
    """

    system_name: str = ""
    forward_tm: float = 0.0
    reverse_tm: float = 0.0
    probe_tm: float = 0.0
    forward_hairpin_dg: float = 0.0
    reverse_hairpin_dg: float = 0.0
    probe_hairpin_dg: float = 0.0
    forward_homodimer_dg: float = 0.0
    reverse_homodimer_dg: float = 0.0
    heterodimer_dg: float = 0.0
    probe_homodimer_dg: float = 0.0
    forward_hairpin_tm: float = 0.0
    reverse_hairpin_tm: float = 0.0
    probe_hairpin_tm: float = 0.0
    tm_difference: float = 0.0
    acceptable: bool = True
    notes: str = ""

    def to_dict(self) -> dict:
        return asdict(self)


def _primer3_tm_kwargs(sys: ReactionSystem) -> dict:
    """构造 primer3.calc_tm 的盐校正参数。

    primer3 接受 mv_conc (monovalent cations mM) 和 dv_conc (divalent cations mM)。
    dNTP 也会螯合 Mg2+,primer3 自动从 dv_conc 中扣除。
    """
    return {
        "mv_conc": sys.monovalent_conc_mM,
        "dv_conc": sys.mg_conc_mM,
        "dntp_conc": sys.dntp_conc_mM,
        "dna_conc": sys.dna_conc_ng * 1e-3,  # 近似转换为 μM 量级;primer3 接受任意正值
        "tm_method": "santalucia",
        "salt_corrections_method": "santalucia",
    }


def _primer3_struct_kwargs(sys: ReactionSystem) -> dict:
    """构造 primer3 发夹/二聚体结构计算所需的盐参数。

    primer3.calc_hairpin/calc_homodimer/calc_heterodimer 接受
    mv_conc (mM) / dv_conc (mM) / dntp_conc (mM) / dna_conc (nM)
    等参数;不传这些参数时,所有 master mix 体系会得到相同的默认值。
    这里使用与 calc_tm 一致的盐修正假设。
    """
    return {
        "mv_conc": sys.monovalent_conc_mM,
        "dv_conc": sys.mg_conc_mM,
        "dntp_conc": sys.dntp_conc_mM,
        "dna_conc": sys.dna_conc_ng * 1e-3,
    }


def evaluate_thermodynamics(
    forward: str,
    reverse: str,
    probe: str,
    system: ReactionSystem,
) -> ThermoResult:
    """在指定反应体系下评估引物/探针的热力学性质。

    Args:
        forward: 正向引物。
        reverse: 反向引物。
        probe: 探针(可空字符串)。
        system: 反应体系参数。

    Returns:
        ThermoResult 对象。
    """
    tm_kw = _primer3_tm_kwargs(system)
    try:
        f_tm = float(primer3.calc_tm(forward, **tm_kw))
        r_tm = float(primer3.calc_tm(reverse, **tm_kw))
    except Exception:
        # 后备:无盐校正
        f_tm = float(primer3.calc_tm(forward))
        r_tm = float(primer3.calc_tm(reverse))

    p_tm = 0.0
    if probe:
        try:
            p_tm = float(primer3.calc_tm(probe, **tm_kw))
        except Exception:
            p_tm = float(primer3.calc_tm(probe))

    # 二聚体/发夹使用 primer3 calc_homodimer / heterodimer / hairpin
    # 必须传入盐浓度参数(mv_conc/dv_conc/dntp_conc/dna_conc),
    # 否则不同 master mix 体系会得到完全相同的默认值结果。
    tm_kw_struct = _primer3_struct_kwargs(system)

    def _hairpin(seq: str) -> tuple[float, float]:
        if not seq:
            return 0.0, 0.0
        try:
            res = primer3.calc_hairpin(seq, **tm_kw_struct)
            return float(res.dg / 1000.0), float(res.tm)  # dg cal/mol -> kcal/mol
        except Exception:
            return 0.0, 0.0

    def _homodimer(seq: str) -> float:
        if not seq:
            return 0.0
        try:
            res = primer3.calc_homodimer(seq, **tm_kw_struct)
            return float(res.dg / 1000.0)
        except Exception:
            return 0.0

    def _heterodimer(s1: str, s2: str) -> float:
        if not s1 or not s2:
            return 0.0
        try:
            res = primer3.calc_heterodimer(s1, s2, **tm_kw_struct)
            return float(res.dg / 1000.0)
        except Exception:
            return 0.0

    f_hp_dg, f_hp_tm = _hairpin(forward)
    r_hp_dg, r_hp_tm = _hairpin(reverse)
    p_hp_dg, p_hp_tm = _hairpin(probe) if probe else (0.0, 0.0)

    f_homo_dg = _homodimer(forward)
    r_homo_dg = _homodimer(reverse)
    p_homo_dg = _homodimer(probe) if probe else 0.0
    het_dg = _heterodimer(forward, reverse)

    tm_diff = abs(f_tm - r_tm)

    # 评估可用性
    notes_list: list[str] = []
    acceptable = True
    if tm_diff > 2.0:
        notes_list.append(f"引物 Tm 差异较大({tm_diff:.1f}°C)")
        acceptable = False
    if p_tm > 0 and p_tm < f_tm:
        notes_list.append("探针 Tm 低于引物 Tm,可能导致探针先解链")
        acceptable = False
    if p_hp_tm > 40 and probe:
        notes_list.append(f"探针发夹 Tm 偏高({p_hp_tm:.1f}°C)")
        acceptable = False
    if het_dg < -9.0:
        notes_list.append(f"引物间二聚体 ΔG 过低({het_dg:.1f} kcal/mol)")
        acceptable = False

    return ThermoResult(
        system_name=system.name,
        forward_tm=f_tm,
        reverse_tm=r_tm,
        probe_tm=p_tm,
        forward_hairpin_dg=f_hp_dg,
        reverse_hairpin_dg=r_hp_dg,
        probe_hairpin_dg=p_hp_dg,
        forward_homodimer_dg=f_homo_dg,
        reverse_homodimer_dg=r_homo_dg,
        heterodimer_dg=het_dg,
        probe_homodimer_dg=p_homo_dg,
        forward_hairpin_tm=f_hp_tm,
        reverse_hairpin_tm=r_hp_tm,
        probe_hairpin_tm=p_hp_tm,
        tm_difference=tm_diff,
        acceptable=acceptable,
        notes="; ".join(notes_list) if notes_list else "通过基础评估",
    )


def _preset_to_system(name: str, preset: dict) -> ReactionSystem:
    """从预设 JSON 字典转换为 ReactionSystem。"""
    return ReactionSystem(
        name=preset.get("name", name),
        mg_conc_mM=float(preset.get("mg_conc_mM", 3.0)),
        dntp_conc_mM=float(preset.get("dntp_conc_mM", 0.8)),
        primer_conc_nM=float(preset.get("primer_conc_nM", 200)),
        probe_conc_nM=float(preset.get("probe_conc_nM", 100)),
        dna_conc_ng=float(preset.get("dna_conc_ng", 50)),
        monovalent_conc_mM=50.0,
    )


def load_preset_systems(path: Optional[str] = None) -> dict[str, ReactionSystem]:
    """加载所有预设体系为 ReactionSystem 字典。

    Args:
        path: JSON 路径;None 时使用默认。

    Returns:
        {key: ReactionSystem} 字典。
    """
    if path is None:
        path = str(Path(__file__).resolve().parent.parent / "data" / "master_mix_presets.json")
    with open(path, "r", encoding="utf-8") as f:
        raw = json.load(f)
    return {k: _preset_to_system(k, v) for k, v in raw.items()}


def compare_master_mixes(
    forward: str,
    reverse: str,
    probe: str,
    systems: dict[str, ReactionSystem],
) -> list[ThermoResult]:
    """在多个反应体系下批量评估同一组引物探针。

    Args:
        forward: 正向引物。
        reverse: 反向引物。
        probe: 探针。
        systems: 体系字典。

    Returns:
        ThermoResult 列表。
    """
    results: list[ThermoResult] = []
    for sys in systems.values():
        results.append(evaluate_thermodynamics(forward, reverse, probe, sys))
    return results


def recommend_best_system(
    results: list[ThermoResult],
    target_primer_tm: float = 60.0,
    target_probe_tm: float = 70.0,
) -> Optional[ThermoResult]:
    """根据评估结果推荐最佳体系。

    评分依据:Tm 接近目标、引物 Tm 差小、二聚体 ΔG 弱、发夹 Tm 低。

    Args:
        results: ThermoResult 列表。
        target_primer_tm: 引物目标 Tm。
        target_probe_tm: 探针目标 Tm。

    Returns:
        得分最低(Tm 偏差最小)的 ThermoResult;无结果时返回 None。
    """
    if not results:
        return None
    best: Optional[ThermoResult] = None
    best_score = float("inf")
    for r in results:
        score = (
            abs(r.forward_tm - target_primer_tm)
            + abs(r.reverse_tm - target_primer_tm)
            + abs(r.probe_tm - target_probe_tm) * 0.5
            + r.tm_difference * 0.5
            + (r.heterodimer_dg if r.heterodimer_dg < -7 else 0)
        )
        if r.acceptable and score < best_score:
            best_score = score
            best = r
    return best