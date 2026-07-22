"""仿真平台与真实文献数据的对照验证。

验证目标:
1. 对 CDC N1 (SARS-CoV-2 N gene) 设计引物探针
2. 用 amplification.py 仿真 10^1 ~ 10^6 拷贝梯度
3. 把仿真 Ct 与文献报告值对照(典型 CDC N1 数据)
4. 评估:
   - Ct-拷贝数关系是否落在真实范围
   - 扩增效率是否落在 90-110%
   - LOD 是否落在 1-10 拷贝
   - 标准曲线 R² 是否 ≥ 0.99
   - 低拷贝(1-10 copies)Ct 是否合理(35-38)
"""
from __future__ import annotations
import os, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
PROJ_ROOT = HERE.parent
sys.path.insert(0, str(PROJ_ROOT))

# 先检查关键依赖可用
try:
    import primer3
except ImportError:
    print("[ERR] primer3-py 未安装,请先 pip install primer3-py")
    sys.exit(1)

try:
    import bio                                                          # noqa: F401
except ImportError:
    print("[WARN] Biopython 未安装,GenBank 检索将不可用")

from core.amplification import (
    AmplificationParams,
    simulate_amplification,
    fit_standard_curve,
    estimate_lod,
)
from validation.CDC_N1_golden_standard import (
    CDC_N1_PRIMERS,
    GOLDEN_STANDARDS,
    get_default_benchmark,
    get_theoretical_benchmark,
    measured_ct_vs_copies,
)


def sep(title: str = ""):
    print()
    print("=" * 72)
    if title:
        print(title)
        print("=" * 72)


def run_design_check():
    """Step 1: 用 primer3 设计 CDC N1 区域引物,与黄金对照序列核对。"""
    sep("Step 1: 引物设计与黄金对照核对")
    print("靶标: SARS-CoV-2 N 基因, CDC N1 区域")
    print(f"  黄金对照 F: {CDC_N1_PRIMERS['F']} ({len(CDC_N1_PRIMERS['F'])} nt)")
    print(f"  黄金对照 R: {CDC_N1_PRIMERS['R']} ({len(CDC_N1_PRIMERS['R'])} nt)")
    print(f"  黄金对照 Probe: {CDC_N1_PRIMERS['probe']} ({len(CDC_N1_PRIMERS['probe'])} nt)")

    # 用 primer3 直接计算黄金对照的 Tm(无需 thermo_sim 自定义类型)
    p = CDC_N1_PRIMERS
    # 计算 Tm(理想体系)
    tm_f = primer3.calc_tm(p["F"], mv_conc=50, dv_conc=3.0, dntp_conc=0.8, dna_conc=200)
    tm_r = primer3.calc_tm(p["R"], mv_conc=50, dv_conc=3.0, dntp_conc=0.8, dna_conc=200)
    tm_p = primer3.calc_tm(p["probe"], mv_conc=50, dv_conc=3.0, dntp_conc=0.8, dna_conc=200)
    print(f"\n  primer3 计算 Tm (ABI 体系: Mg²⁺3.0/dNTP0.8/primer200nM):")
    print(f"    F:  {tm_f:.1f}°C  (CDC 文献 ~57-58°C)")
    print(f"    R:  {tm_r:.1f}°C  (CDC 文献 ~58-59°C)")
    print(f"    Pr: {tm_p:.1f}°C  (CDC 文献 ~68-69°C)")
    return p


def run_amplification_simulation(copies: list[float]) -> dict:
    """Step 2: 用 amplification.py 仿真梯度扩增,返回结果字典。"""
    sep("Step 2: 仿真扩增曲线(Ct vs copies)")
    print(f"  拷贝数梯度: {copies}")
    params = AmplificationParams(efficiency=0.95, threshold=0.1)
    samples = simulate_amplification(copies, params)
    for s in samples:
        ct_str = f"{s.ct:.2f}" if s.ct is not None else "Undetermined"
        print(f"    {s.start_copies:>12,.0f} copies → Ct = {ct_str}")
    return {"samples": samples, "params": params}


def run_standard_curve_fit(samples):
    """Step 3: 拟合仿真标准曲线。"""
    sep("Step 3: 仿真标准曲线拟合")
    sc = fit_standard_curve(samples)
    print(f"  slope     = {sc.slope:.4f}  (理想 -3.32)")
    print(f"  intercept = {sc.intercept:.4f}")
    print(f"  R²        = {sc.r_squared:.4f}  (应 ≥0.99)")
    print(f"  E         = {sc.efficiency_pct:.2f}%  (应 90-110%)")
    print(f"  线性范围  = {sc.linear_range}")
    return sc


def run_lod_estimate(samples):
    """Step 4: LOD 估计。"""
    sep("Step 4: LOD 估计")
    lod = estimate_lod(samples, detection_probability=0.95, n_replicates=3)
    print(f"  95% LOD  = {lod.get('lod_copies'):.2f} copies")
    print(f"  说明: {lod.get('note')}")
    return lod


def compare_with_golden(our_cts: dict, benchmark_name: str = None):
    """Step 5: 与文献基准逐点对照。"""
    sep(f"Step 5: 仿真 vs 文献基准对照 ({benchmark_name or 'default'})")
    benchmark = next(
        (b for b in GOLDEN_STANDARDS if b.name == benchmark_name),
        get_default_benchmark(),
    )
    measured = benchmark.measured_ct

    print(f"{'copies':>12} | {'仿真 Ct':>10} | {'文献 Ct':>10} | {'ΔCt':>8} | {'判定':>8}")
    print("-" * 60)
    total_dev = 0.0
    n = 0
    for c, ct_sim in sorted(our_cts.items()):
        if c in measured:
            ct_lit = measured[c]
            dev = ct_sim - ct_lit
            verdict = "✓" if abs(dev) <= 1.0 else ("⚠" if abs(dev) <= 2.0 else "✗")
            print(f"{c:>12,.0f} | {ct_sim:>10.2f} | {ct_lit:>10.2f} | {dev:>+8.2f} | {verdict:>8}")
            total_dev += abs(dev)
            n += 1
    if n > 0:
        print(f"\n  平均 |ΔCt| = {total_dev/n:.2f}  (判定: ≤1.0 优秀, ≤2.0 可接受, >2.0 需修模)")
    return {"avg_abs_dev": total_dev / max(n, 1), "n_points": n}


def main():
    print("# qPCR 仿真平台 vs CDC N1 真实文献数据对照验证")
    print(f"# 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")

    # Step 1: 引物核对
    p = run_design_check()

    # Step 2: 仿真扩增(覆盖文献拷贝数范围)
    copies = [10, 100, 1000, 10000, 100000, 1_000_000]
    sim = run_amplification_simulation(copies)
    samples = sim["samples"]

    # Step 3: 标准曲线
    sc = run_standard_curve_fit(samples)

    # Step 4: LOD
    lod = run_lod_estimate(samples)

    # Step 5: 与文献逐点对照
    our_cts = {s.start_copies: s.ct for s in samples if s.ct is not None}
    dev_real = compare_with_golden(our_cts, "CDC N1 / TaqPath + ABI 7500")
    dev_theory = compare_with_golden(our_cts, "Theoretical 100% efficiency (ideal)")

    # 总判定
    sep("总判定")
    # 注:"理论 100% 基准"假设从 cycle 0 起峰,但真实 qPCR
    # 前 3-5 个循环荧光基线未稳,阈值检测必然滞后 ~2-3 个循环。
    # 因此理论基准用 |ΔCt| ≤ 3.0,文献基准用 |ΔCt| ≤ 1.5。
    pass_criteria = [
        ("扩增效率 90-110%",       90.0 <= sc.efficiency_pct <= 110.0),
        ("R² ≥ 0.99",              sc.r_squared >= 0.99),
        ("LOD ≤ 10 copies",        lod.get("lod_copies", 1e9) <= 10.0),
        ("与文献 |ΔCt| ≤ 1.5",     dev_real["avg_abs_dev"] <= 1.5),
        ("与理论 |ΔCt| ≤ 3.0 (qPCR 基线滞后容差)", dev_theory["avg_abs_dev"] <= 3.0),
    ]
    for label, ok in pass_criteria:
        print(f"  [{'✓' if ok else '✗'}] {label}")
    all_pass = all(ok for _, ok in pass_criteria)
    print(f"\n  总结论: {'通过 — 仿真模型在合理范围内' if all_pass else '未通过 — 需调整模型'}")
    print()
    return 0 if all_pass else 1


if __name__ == "__main__":
    sys.exit(main())
