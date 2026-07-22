"""优化循环引擎:根据仿真结果自动评估并生成优化建议。

支持多维度优化:引物序列、引物/探针长度、反应体系、标记体系、探针修饰。
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

from core.amplification import (
    AmplificationParams,
    CurveSample,
    fit_standard_curve,
    simulate_amplification,
)
from core.primer_design import DesignParameters, PrimerPairResult, design_primer_pairs
from core.thermo_sim import (
    ReactionSystem,
    ThermoResult,
    compare_master_mixes,
    evaluate_thermodynamics,
    load_preset_systems,
    recommend_best_system,
)


@dataclass
class OptimizationSuggestion:
    """单条优化建议。

    Attributes:
        dimension: 优化维度(引物序列/反应体系/标记体系/探针修饰/参数调整)。
        priority: 优先级(high/medium/low)。
        issue: 当前问题描述。
        action: 建议采取的动作。
    """

    dimension: str
    priority: str
    issue: str
    action: str


@dataclass
class OptimizationResult:
    """一轮优化的完整结果。

    Attributes:
        before_pair: 优化前引物对。
        after_pair: 优化后引物对(若无可用为 None);若仅调整体系,可能与 before_pair 相同。
        before_thermo: 优化前热力学评估(在 before_system 下)。
        after_thermo: 优化后热力学评估(在 after_system 下;若 None 则未执行)。
        before_system: 优化前体系名称。
        after_system: 优化后体系名称(若未调整则为 None)。
        before_curves: 优化前扩增曲线。
        after_curves: 优化后扩增曲线。
        before_std: 优化前标准曲线。
        after_std: 优化后标准曲线。
        suggestions: 建议列表。
        selected_system: 推荐的反应体系名称。
        summary: 优化总结文字。
        changed_primer: 是否重新设计了引物。
        changed_system: 是否切换了反应体系。
    """

    before_pair: PrimerPairResult
    after_pair: Optional[PrimerPairResult]
    before_thermo: ThermoResult
    after_thermo: Optional[ThermoResult]
    before_system: str = ""
    after_system: Optional[str] = None
    before_curves: list[CurveSample] = field(default_factory=list)
    after_curves: list[CurveSample] = field(default_factory=list)
    before_std: Optional[object] = None
    after_std: Optional[object] = None
    suggestions: list[OptimizationSuggestion] = field(default_factory=list)
    selected_system: str = ""
    summary: str = ""
    changed_primer: bool = False
    changed_system: bool = False


def _generate_suggestions(
    pair: PrimerPairResult,
    thermo: ThermoResult,
    std_efficiency: Optional[float],
    std_r2: Optional[float],
    lod_log: Optional[float],
) -> list[OptimizationSuggestion]:
    """根据当前状态生成优化建议。"""
    suggestions: list[OptimizationSuggestion] = []

    # 1. 引物 Tm 差异
    if thermo.tm_difference > 2.0:
        suggestions.append(OptimizationSuggestion(
            dimension="引物序列",
            priority="high",
            issue=f"引物 Tm 差异 {thermo.tm_difference:.1f}°C,可能影响扩增效率",
            action="重新设计,选择 Tm 更接近的引物对;或调整引物长度使 Tm 对齐",
        ))
    elif thermo.tm_difference > 1.0:
        suggestions.append(OptimizationSuggestion(
            dimension="引物序列",
            priority="medium",
            issue=f"引物 Tm 差异 {thermo.tm_difference:.1f}°C 偏大",
            action="尝试使用其他候选对或微调引物长度",
        ))

    # 2. 探针 Tm 异常
    if pair.probe_seq:
        if thermo.probe_tm < thermo.forward_tm:
            suggestions.append(OptimizationSuggestion(
                dimension="探针修饰",
                priority="high",
                issue=f"探针 Tm({thermo.probe_tm:.1f}°C)低于引物 Tm({thermo.forward_tm:.1f}°C)",
                action="加长探针 2~4 nt 或在探针中引入 LNA 碱基以提高 Tm",
            ))
        if thermo.probe_hairpin_tm > 40.0:
            suggestions.append(OptimizationSuggestion(
                dimension="探针修饰",
                priority="medium",
                issue=f"探针发夹 Tm {thermo.probe_hairpin_tm:.1f}°C 偏高,可能影响杂交效率",
                action="缩短探针长度或调整 GC 含量;可考虑 MGB 探针",
            ))

    # 3. 二聚体
    if thermo.heterodimer_dg < -9.0:
        suggestions.append(OptimizationSuggestion(
            dimension="引物序列",
            priority="high",
            issue=f"引物间异源二聚体 ΔG 较低({thermo.heterodimer_dg:.1f} kcal/mol)",
            action="重新设计引物对,避免 3' 端互补",
        ))

    # 4. 扩增效率
    if std_efficiency is not None and std_efficiency < 90.0:
        suggestions.append(OptimizationSuggestion(
            dimension="反应体系",
            priority="high",
            issue=f"预测扩增效率 {std_efficiency:.1f}% 偏低(<90%)",
            action="尝试 Mg2+ 浓度更高的体系(如 Takara 4 mM);或检查引物二聚体",
        ))
    elif std_efficiency is not None and std_efficiency < 95.0:
        suggestions.append(OptimizationSuggestion(
            dimension="反应体系",
            priority="medium",
            issue=f"预测扩增效率 {std_efficiency:.1f}% 可优化",
            action="微调 Mg2+ / 引物浓度,或选择更优试剂盒",
        ))

    # 5. 标准曲线 R²
    if std_r2 is not None and std_r2 < 0.98:
        suggestions.append(OptimizationSuggestion(
            dimension="参数调整",
            priority="medium",
            issue=f"标准曲线 R²={std_r2:.3f} 低于 0.98",
            action="检查最低/最高浓度点;增加中间梯度或扩大线性范围",
        ))

    # 6. LOD
    if lod_log is not None and lod_log > 2.0:
        suggestions.append(OptimizationSuggestion(
            dimension="反应体系",
            priority="medium",
            issue=f"预测 LOD≈10^{lod_log:.1f} 偏高",
            action="增加引物/探针浓度;改用高灵敏度试剂盒(如 Takara HS)",
        ))

    # 7. GC 含量
    if pair.forward_gc > 65 or pair.reverse_gc > 65:
        suggestions.append(OptimizationSuggestion(
            dimension="引物序列",
            priority="low",
            issue=f"引物 GC% 偏高(F:{pair.forward_gc:.0f}%, R:{pair.reverse_gc:.0f}%)",
            action="考虑替换候选引物对以降低 GC%",
        ))
    if pair.forward_gc < 35 or pair.reverse_gc < 35:
        suggestions.append(OptimizationSuggestion(
            dimension="引物序列",
            priority="low",
            issue=f"引物 GC% 偏低(F:{pair.forward_gc:.0f}%, R:{pair.reverse_gc:.0f}%)",
            action="考虑替换候选引物对以提高 GC%",
        ))

    # 如果没有任何问题
    if not suggestions:
        suggestions.append(OptimizationSuggestion(
            dimension="无",
            priority="low",
            issue="当前引物探针组合表现良好",
            action="无需调整,可直接进行实验验证",
        ))

    return suggestions


def _apply_to_design_params(
    base: DesignParameters,
    suggestions: list[OptimizationSuggestion],
) -> DesignParameters:
    """根据建议调整设计参数,生成下一轮设计输入。"""
    new = DesignParameters(**{k: getattr(base, k) for k in base.__dataclass_fields__})

    for sug in suggestions:
        if sug.dimension == "引物序列":
            # 放宽 Tm 范围,允许更长/更短的引物
            new.primer_tm_min = max(56.0, base.primer_tm_min - 1.0)
            new.primer_tm_max = min(64.0, base.primer_tm_max + 1.0)
            new.num_return = min(20, base.num_return + 5)
        elif sug.dimension == "反应体系":
            new.product_size_min = max(50, base.product_size_min - 10)
            new.product_size_max = min(200, base.product_size_max + 10)
        elif sug.dimension == "探针修饰":
            if "Tm" in sug.issue or "Tm" in sug.action:
                # 提高探针 Tm 范围
                new.probe_tm_min = max(65.0, base.probe_tm_min - 1.0)
                new.probe_tm_max = min(75.0, base.probe_tm_max + 1.0)

    return new


class OptimizationEngine:
    """优化循环引擎:评估 → 建议 → 重新设计 → 再评估。"""

    def __init__(self, sequence: str, system: ReactionSystem) -> None:
        """初始化引擎。

        Args:
            sequence: 模板序列。
            system: 初始反应体系。
        """
        self.sequence = sequence
        self.system = system

    def evaluate(
        self,
        pair: PrimerPairResult,
        copy_numbers: list[float],
        amp_params: AmplificationParams,
    ) -> dict:
        """评估当前引物对的热力学与扩增表现。

        Args:
            pair: 待评估的引物对。
            copy_numbers: 拷贝数梯度。
            amp_params: 扩增参数。

        Returns:
            包含 thermo / curves / std / lod 的字典。
        """
        thermo = evaluate_thermodynamics(pair.forward_seq, pair.reverse_seq, pair.probe_seq, self.system)
        curves = simulate_amplification(copy_numbers, amp_params)

        std = None
        try:
            std = fit_standard_curve(curves)
        except ValueError:
            std = None

        from core.amplification import estimate_lod

        lod_info = estimate_lod(curves)

        std_eff = std.efficiency_pct if std else None
        std_r2 = std.r_squared if std else None
        lod_log = lod_info.get("lod_log10")

        suggestions = _generate_suggestions(pair, thermo, std_eff, std_r2, lod_log)

        return {
            "thermo": thermo,
            "curves": curves,
            "std": std,
            "lod": lod_info,
            "suggestions": suggestions,
        }

    def optimize(
        self,
        initial_pair: PrimerPairResult,
        design_params: DesignParameters,
        copy_numbers: list[float],
        amp_params: Optional[AmplificationParams] = None,
        max_iterations: int = 2,
    ) -> OptimizationResult:
        """执行完整优化循环。

        Args:
            initial_pair: 初始引物对。
            design_params: 初始设计参数。
            copy_numbers: 拷贝数梯度。
            amp_params: 扩增参数。
            max_iterations: 最大迭代轮数。

        Returns:
            OptimizationResult 对象。
        """
        amp_params = amp_params or AmplificationParams()

        # 评估当前
        before_eval = self.evaluate(initial_pair, copy_numbers, amp_params)
        before_thermo = before_eval["thermo"]
        before_curves = before_eval["curves"]
        before_std = before_eval["std"]
        before_lod = before_eval["lod"]

        suggestions = before_eval["suggestions"]
        before_system_name = self.system.name
        before_summary = (
            f"初始[{before_system_name}]: 引物 Tm 差 {before_thermo.tm_difference:.1f}°C, "
            f"F_Tm={before_thermo.forward_tm:.1f}°C, "
            f"F_发夹dG={before_thermo.forward_hairpin_dg:.2f}, "
            f"标准曲线 R²={before_std.r_squared:.3f}" if before_std else
            f"初始[{before_system_name}]: Tm 差 {before_thermo.tm_difference:.1f}°C"
        )

        # 1. 推荐反应体系(基于初始引物)
        presets = load_preset_systems()
        thermo_all = compare_master_mixes(
            initial_pair.forward_seq,
            initial_pair.reverse_seq,
            initial_pair.probe_seq,
            presets,
        )
        best_sys = recommend_best_system(thermo_all)
        selected_system = best_sys.system_name if best_sys else self.system.name

        after_pair: Optional[PrimerPairResult] = None
        after_thermo: Optional[ThermoResult] = None
        after_curves: list[CurveSample] = []
        after_std = None
        after_system_name: Optional[str] = None
        changed_primer = False
        changed_system = False
        after_summary = "初始引物对表现良好,无需优化"

        # 2. 决定优化路径
        #    体系建议 → 换 master mix 后重新评估(展示真实数值差异)
        #    引物序列建议 → 重新设计引物对,再用选中的体系评估
        sys_suggestion = [s for s in suggestions if s.dimension == "反应体系"]
        primer_suggestion = [s for s in suggestions if s.dimension == "引物序列"]
        probe_suggestion = [s for s in suggestions if s.dimension == "探针修饰"]

        target_system: Optional[ReactionSystem] = None
        target_pair: Optional[PrimerPairResult] = None

        # 优先响应体系建议:选最佳体系(若当前不是最佳)
        if sys_suggestion and best_sys is not None and best_sys.system_name != self.system.name:
            target_system = best_sys
        # 其次,若建议提到提高 Mg2+,尝试 Mg2+ 最高的 preset
        if target_system is None and sys_suggestion:
            # 找 Mg2+ 最高的 preset
            highest_mg = max(presets.values(), key=lambda s: s.mg_conc_mM)
            if highest_mg.mg_conc_mM > self.system.mg_conc_mM:
                target_system = highest_mg

        # 引物序列/探针修饰建议 → 重新设计
        if primer_suggestion or probe_suggestion:
            new_params = _apply_to_design_params(design_params, suggestions)
            try:
                new_pairs = design_primer_pairs(self.sequence, new_params)
                if new_pairs:
                    initial_f = initial_pair.forward_seq
                    candidates = [p for p in new_pairs if p.forward_seq != initial_f]
                    if not candidates:
                        candidates = new_pairs
                    target_pair = candidates[0]
            except (RuntimeError, ValueError):
                target_pair = None

        # 3. 实际执行评估(必须产生可对比的 after_* 数据)
        eval_system = target_system if target_system is not None else self.system
        eval_pair = target_pair if target_pair is not None else initial_pair

        # 只有当 target_* 真的与当前不同时,才构造 after_*
        system_changed = target_system is not None and target_system.name != self.system.name
        primer_changed = target_pair is not None and target_pair.forward_seq != initial_pair.forward_seq

        if system_changed or primer_changed:
            try:
                new_engine = OptimizationEngine(self.sequence, eval_system)
                after_eval = new_engine.evaluate(eval_pair, copy_numbers, amp_params)
                after_pair = eval_pair
                after_thermo = after_eval["thermo"]
                after_curves = after_eval["curves"]
                after_std = after_eval["std"]
                after_system_name = eval_system.name
                changed_primer = primer_changed
                changed_system = system_changed

                # 综合得分
                def _score(p: PrimerPairResult, t: ThermoResult, std_obj) -> float:
                    score = t.tm_difference + abs(t.forward_tm - 60) + abs(t.reverse_tm - 60)
                    if std_obj:
                        score += (100 - std_obj.efficiency_pct) * 0.1
                        score += (1 - std_obj.r_squared) * 20
                    return score

                before_score = _score(initial_pair, before_thermo, before_std)
                after_score = _score(after_pair, after_thermo, after_std) if after_pair else float("inf")

                # 描述对比(突出实际数值变化)
                chg_bits = []
                if changed_system:
                    chg_bits.append(f"体系 {before_system_name}→{after_system_name}")
                if changed_primer:
                    chg_bits.append("重新设计引物对")
                after_summary = (
                    f"优化后[{after_system_name}]({'; '.join(chg_bits)}): "
                    f"F_Tm {before_thermo.forward_tm:.1f}→{after_thermo.forward_tm:.1f}°C, "
                    f"R_Tm {before_thermo.reverse_tm:.1f}→{after_thermo.reverse_tm:.1f}°C, "
                    f"Tm差 {before_thermo.tm_difference:.1f}→{after_thermo.tm_difference:.1f}°C, "
                    f"F_发夹dG {before_thermo.forward_hairpin_dg:.2f}→{after_thermo.forward_hairpin_dg:.2f}, "
                    f"异二聚体dG {before_thermo.heterodimer_dg:.2f}→{after_thermo.heterodimer_dg:.2f}"
                )
                if after_std is not None:
                    after_summary += f", 效率 {before_std.efficiency_pct:.1f}→{after_std.efficiency_pct:.1f}%, R² {before_std.r_squared:.3f}→{after_std.r_squared:.3f}"

                if after_score >= before_score:
                    # 未改善 — 但仍保留 after_*(用 None 表示未采用)
                    after_pair = None
                    after_thermo = None
                    after_curves = []
                    after_std = None
                    after_system_name = None
                    changed_primer = False
                    changed_system = False
                    after_summary = (
                        f"尝试调整(体系→{eval_system.name},引物重新设计)后综合评分"
                        f" {before_score:.2f} vs {after_score:.2f},未改善,保留初始方案"
                    )
            except (RuntimeError, ValueError) as exc:
                after_summary = f"优化评估失败: {exc},保留初始引物对"

        return OptimizationResult(
            before_pair=initial_pair,
            after_pair=after_pair,
            before_thermo=before_thermo,
            after_thermo=after_thermo,
            before_system=before_system_name,
            after_system=after_system_name,
            before_curves=before_curves,
            after_curves=after_curves,
            before_std=before_std,
            after_std=after_std,
            suggestions=suggestions,
            selected_system=selected_system,
            summary=before_summary + " || " + after_summary,
            changed_primer=changed_primer,
            changed_system=changed_system,
        )