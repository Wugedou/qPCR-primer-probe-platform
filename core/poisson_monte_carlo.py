"""Poisson 蒙特卡洛低拷贝模型 — PLAN-002 Step 8

在 1-10 copies 区间用随机过程替代确定性 S 曲线。
每次仿真从 Poisson(λ=copies) 抽取实际起始分子数。
"""

import random
import math
from dataclasses import dataclass
from typing import List, Tuple, Optional, Dict
import statistics

# Poisson sampler
def _poisson(lam):
    """Poisson distributed random number (optimized for small lambda)"""
    if lam == 0:
        return 0
    # Knuth's algorithm for Poisson
    L = math.exp(-lam)
    k = 0
    p = 1.0
    while p > L:
        k += 1
        p *= random.random()
    return k - 1

@dataclass
class PoissonResult:
    """Poisson 仿真结果"""
    target_copies: int        # 目标拷贝数
    n_simulations: int         # 仿真次数
    detected: int              # 检出次数 (Ct < max_ct)
    detection_prob: float      # 检出概率
    mean_ct: Optional[float]   # 检出样本的平均 Ct
    std_ct: Optional[float]    # 检出样本的 Ct 标准差
    mean_actual_copies: float  # Poisson 抽取后的平均实际拷贝数

def poisson_monte_carlo(
    target_copies: int,
    max_ct: int = 40,
    threshold: float = 0.05,
    n_simulations: int = 1000,
    slope: float = -3.32,
    intercept: float = 38.5,
    efficiency: float = 1.0,
    seed: int = 42,
) -> PoissonResult:
    """运行 Poisson 蒙特卡洛低拷贝仿真

    Args:
        target_copies: 目标拷贝数（期望值）
        max_ct: 最大 Ct 值（超过此值视为未检出）
        threshold: 荧光阈值
        n_simulations: 蒙特卡洛仿真次数
        slope: 标准曲线斜率
        intercept: 标准曲线截距
        efficiency: 扩增效率 (0-1)
        seed: 随机种子

    Returns:
        PoissonResult 含检出概率和统计
    """
    random.seed(seed)
    
    actual_copies_list = []
    ct_values = []
    detected_count = 0
    
    for _ in range(n_simulations):
        # Poisson 抽取实际起始分子数
        actual = _poisson(target_copies)
        actual_copies_list.append(actual)
        
        if actual == 0:
            # 没有分子 → 未检出
            continue
        
        # 确定性 Ct 计算（仿真平台的 S 曲线模型）
        # Ct = slope * log10(actual) + intercept, 按效率校正
        log_copies = math.log10(actual)
        ct = slope * log_copies + intercept
        
        # 效率校正：效率 < 1 时 Ct 后移
        if efficiency < 1.0:
            ct += (1.0 / efficiency - 1.0) * (slope * log_copies)
        
        if ct < max_ct:
            detected_count += 1
            ct_values.append(ct)
    
    mean_actual = statistics.mean(actual_copies_list) if actual_copies_list else 0
    
    return PoissonResult(
        target_copies=target_copies,
        n_simulations=n_simulations,
        detected=detected_count,
        detection_prob=detected_count / n_simulations,
        mean_ct=statistics.mean(ct_values) if ct_values else None,
        std_ct=statistics.stdev(ct_values) if len(ct_values) > 1 else None,
        mean_actual_copies=mean_actual,
    )

def probit_curve(
    copy_range: List[int] = [1, 2, 3, 5, 10, 20, 50, 100],
    n_simulations: int = 1000,
    **kwargs,
) -> List[PoissonResult]:
    """生成检出概率-拷贝数曲线

    Returns:
        PoissonResult 列表，可用于拟合 probit 曲线
    """
    results = []
    for copies in copy_range:
        r = poisson_monte_carlo(
            target_copies=copies,
            n_simulations=n_simulations,
            **kwargs,
        )
        results.append(r)
    return results

def find_lod95(results: List[PoissonResult]) -> Tuple[Optional[int], float]:
    """从 probit 结果中找 95% 检出 LOD

    Returns:
        (LOD95_copies, detection_prob)
    """
    # 插值找第一个 ≥95% 的
    for r in results:
        if r.detection_prob >= 0.95:
            return r.target_copies, r.detection_prob
    return None, 0.0

def compare_with_literature(
    our_results: List[PoissonResult],
    literature_lod: int,
    literature_probit_data: Optional[Dict[int, float]] = None,
) -> dict:
    """与文献 probit 数据对照

    Args:
        our_results: 我们的 Poisson 仿真结果
        literature_lod: 文献报告的 95% LOD
        literature_probit_data: 文献的 {copies: detection_prob} 数据（可选）

    Returns:
        对照结果字典
    """
    our_lod, our_prob = find_lod95(our_results)
    
    comparison = {
        "our_lod95": our_lod,
        "literature_lod95": literature_lod,
        "lod_diff_log10": round(
            math.log10(our_lod / literature_lod) if our_lod and literature_lod else 0, 2
        ),
    }
    
    if literature_probit_data:
        rmse = 0
        n = 0
        for copies, lit_prob in literature_probit_data.items():
            # 找最接近的仿真结果
            for r in our_results:
                if r.target_copies == copies:
                    rmse += (r.detection_prob - lit_prob) ** 2
                    n += 1
                    break
        comparison["probit_rmse"] = math.sqrt(rmse / n) if n else None
    
    return comparison


# ============================================================
# 与 Behrmann 2020 对照（PMID:32384153, SARS-CoV-2 RT-RPA probit）
# ============================================================
BEHRMANN_2020_LOD = 7  # copies/reaction (文献 95% LOD)
BEHRMANN_2020_PROBIT = {
    # 近似的 probit 数据（从文献 Figure 估算）
    1: 0.05,
    2: 0.15,
    3: 0.35,
    5: 0.65,
    7: 0.95,
    10: 0.99,
}
