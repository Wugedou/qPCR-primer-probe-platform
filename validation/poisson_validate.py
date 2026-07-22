"""PLAN-002 Step 8: Poisson 蒙特卡洛低拷贝验证

验证 Poisson 模型与 Behrmann 2020 probit 数据的一致性。
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from core.poisson_monte_carlo import (
    probit_curve,
    find_lod95,
    compare_with_literature,
    BEHRMANN_2020_LOD,
    BEHRMANN_2020_PROBIT,
)

if __name__ == "__main__":
    print("=" * 60)
    print("  PLAN-002: Poisson 蒙特卡洛低拷贝验证")
    print("=" * 60)
    
    # 运行 probit 曲线
    results = probit_curve(
        copy_range=[1, 2, 3, 5, 7, 10, 20, 50, 100],
        n_simulations=2000,
        slope=-3.32,
        intercept=38.5,
        max_ct=40,
    )
    
    # 打印结果表
    print(f"\n{'Copies':>8} {'Poisson mean':>12} {'Detected':>8} {'Prob':>8} {'Mean Ct':>8} {'Std Ct':>8}")
    print("-" * 65)
    for r in results:
        ct_str = f"{r.mean_ct:.1f}" if r.mean_ct else "—"
        std_str = f"{r.std_ct:.1f}" if r.std_ct else "—"
        print(f"{r.target_copies:>8} {r.mean_actual_copies:>12.1f} {r.detected:>8} {r.detection_prob:>7.1%} {ct_str:>8} {std_str:>8}")
    
    # 找 LOD95
    lod, prob = find_lod95(results)
    print(f"\n  95% LOD: {lod} copies (检出率={prob:.1%})")
    
    # 与 Behrmann 2020 对照
    comparison = compare_with_literature(
        results,
        BEHRMANN_2020_LOD,
        BEHRMANN_2020_PROBIT,
    )
    
    print(f"\n  === Behrmann 2020 对照 ===")
    print(f"  文献 95% LOD: {comparison['literature_lod95']} copies")
    print(f"  仿真 95% LOD: {comparison['our_lod95']} copies")
    print(f"  LOD 对数偏差: {comparison['lod_diff_log10']} log10")
    if comparison.get('probit_rmse'):
        print(f"  Probit RMSE: {comparison['probit_rmse']:.3f}")
    
    # 验收判定 V5
    print(f"\n  === V5 验收 ===")
    lod_diff = abs(comparison['lod_diff_log10'])
    if lod_diff <= 0.5:
        print(f"  ✅ PASS: LOD 偏差 {lod_diff:.1f} ≤ 0.5 log10 (目标)")
    else:
        print(f"  ❌ FAIL: LOD 偏差 {lod_diff:.1f} > 0.5 log10")
    
    # 额外验证：1 copy 的检出率应在 Poisson 预期附近
    # Poisson(1): P(X=0) = e^(-1) ≈ 0.368, 所以检出率应 ≈ 63.2%
    expected_p1 = 1 - 0.368
    actual_p1 = results[0].detection_prob  # 1 copy
    print(f"\n  1 copy 检出率: 实际={actual_p1:.1%} vs Poisson理论={expected_p1:.1%}")
    
    if abs(actual_p1 - expected_p1) < 0.05:
        print(f"  ✅ 与 Poisson 理论一致")
    else:
        print(f"  ⚠️ 偏差 {abs(actual_p1 - expected_p1):.1%}")
