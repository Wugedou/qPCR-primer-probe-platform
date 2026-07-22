"""扩增曲线仿真模块。

支持三种扩增方式:
- qPCR: 循环数-荧光曲线、Ct / 标准曲线 / LOD。
- RPA: 等温扩增(37-42°C),时间-荧光曲线(5-20 分钟)。
- LAMP: 等温扩增(60-65°C),时间-荧光曲线(15-60 分钟)。

qPCR 模型:N(C) = N0 * (1 + E)^C + baseline
其中 C 为循环数, E 为扩增效率(0~1), N0 为起始模板量(荧光单位)。
"""
from __future__ import annotations

from dataclasses import dataclass, asdict, field
from typing import Optional

import numpy as np
from scipy import stats
from scipy.optimize import brentq


@dataclass
class AmplificationParams:
    """扩增仿真参数。

    Attributes:
        cycles: 循环数(默认 40)。
        efficiency: 扩增效率(0.5~1.0)。
        baseline: 基线荧光。
        plateau: 平台期荧光(用于饱和模型)。
        noise_std: 噪声标准差(模拟测量波动)。
        threshold: Ct 阈值荧光。
        random_seed: 随机种子,保证可复现。
    """

    cycles: int = 40
    efficiency: float = 0.95
    baseline: float = 0.02
    plateau: float = 1.0
    noise_std: float = 0.005
    threshold: float = 0.1
    random_seed: int = 42


@dataclass
class CurveSample:
    """单条扩增曲线。

    Attributes:
        label: 曲线标签(如 10^5 copies)。
        start_copies: 起始模板拷贝数。
        cycles: 循环数组。
        fluorescence: 荧光数组。
        ct: 计算得到的 Ct 值(若未达到阈值为 None)。
    """

    label: str
    start_copies: float
    cycles: np.ndarray
    fluorescence: np.ndarray
    ct: Optional[float] = None


@dataclass
class StandardCurveResult:
    """标准曲线拟合结果。

    Attributes:
        log_concentrations: 各标准品 log10(拷贝数)。
        cts: 对应 Ct 值。
        slope: 斜率。
        intercept: 截距。
        r_squared: R²。
        efficiency_pct: 由斜率换算的扩增效率(%)。
        linear_range: 线性范围 (低拷贝, 高拷贝)。
    """

    log_concentrations: np.ndarray
    cts: np.ndarray
    slope: float
    intercept: float
    r_squared: float
    efficiency_pct: float
    linear_range: tuple[float, float]


def simulate_curve(
    start_copies: float,
    params: AmplificationParams,
) -> CurveSample:
    """模拟单条扩增曲线。

    模型:N(C) = baseline + (plateau - baseline) * (1 - exp(-k*N0*(1+E)^C))
    使用 Logistic 模型拟合 S 型曲线,在低浓度时呈指数,高浓度时饱和。

    Args:
        start_copies: 起始模板拷贝数。
        params: 扩增参数。

    Returns:
        CurveSample。
    """
    rng = np.random.default_rng(params.random_seed + int(np.log10(max(start_copies, 1)) * 1000))
    cycles = np.arange(1, params.cycles + 1, dtype=float)

    # 经验 S 型模型:荧光 = baseline + (plateau - baseline) / (1 + exp(-k*(C - midpoint)))
    # 其中 midpoint 与 log10(N0) 相关,效率影响曲线陡度
    # Ct 基准 40 循环:10 拷贝→Ct~36.7, 10^6 拷贝→Ct~20.1
    # 每增 10 倍模板,Ct 前移 ~3.32 个循环(效率 100% 时)
    midpoint = 40.0 - 3.32 * np.log10(max(start_copies, 1.0))
    steepness = 0.6 + 0.4 * params.efficiency  # 效率越高曲线越陡

    # Logistic 生长曲线
    fluorescence = params.baseline + (params.plateau - params.baseline) / (
        1.0 + np.exp(-steepness * (cycles - midpoint))
    )

    # 添加噪声
    if params.noise_std > 0:
        noise = rng.normal(0, params.noise_std, size=fluorescence.shape)
        fluorescence = np.maximum(fluorescence + noise, 0.0)

    # 计算 Ct:荧光首次超过阈值的循环
    ct: Optional[float] = None
    crossed = np.where(fluorescence >= params.threshold)[0]
    if crossed.size > 0:
        idx = crossed[0]
        if idx == 0:
            ct = float(cycles[0])
        else:
            # 在超过阈值的前后两点间线性插值
            c0, c1 = cycles[idx - 1], cycles[idx]
            f0, f1 = fluorescence[idx - 1], fluorescence[idx]
            if f1 > f0:
                ct = float(c0 + (params.threshold - f0) * (c1 - c0) / (f1 - f0))
            else:
                ct = float(c1)

    label = f"{start_copies:.0e} copies"
    return CurveSample(
        label=label,
        start_copies=start_copies,
        cycles=cycles,
        fluorescence=fluorescence,
        ct=ct,
    )


def simulate_amplification(
    copy_numbers: list[float],
    params: Optional[AmplificationParams] = None,
) -> list[CurveSample]:
    """批量仿真多条扩增曲线(用于梯度实验)。

    Args:
        copy_numbers: 起始拷贝数列表。
        params: 扩增参数。

    Returns:
        CurveSample 列表,按拷贝数降序。
    """
    params = params or AmplificationParams()
    samples = [simulate_curve(c, params) for c in copy_numbers]
    samples.sort(key=lambda s: s.start_copies, reverse=True)
    return samples


def compute_ct(sample: CurveSample, threshold: float) -> Optional[float]:
    """从已有曲线重算 Ct。

    Args:
        sample: CurveSample。
        threshold: 阈值。

    Returns:
        Ct 值或 None。
    """
    fluorescence = sample.fluorescence
    cycles = sample.cycles
    crossed = np.where(fluorescence >= threshold)[0]
    if crossed.size == 0:
        return None
    idx = crossed[0]
    if idx == 0:
        return float(cycles[0])
    c0, c1 = cycles[idx - 1], cycles[idx]
    f0, f1 = fluorescence[idx - 1], fluorescence[idx]
    if f1 > f0:
        return float(c0 + (threshold - f0) * (c1 - c0) / (f1 - f0))
    return float(c1)


def fit_standard_curve(
    samples: list[CurveSample],
    min_r2: float = 0.98,
) -> StandardCurveResult:
    """对多条扩增曲线拟合标准曲线 (Ct vs log10(拷贝数))。

    Args:
        samples: CurveSample 列表(通常 5~8 个梯度)。
        min_r2: 用于确定线性范围的 R² 阈值。

    Returns:
        StandardCurveResult。
    """
    valid = [(s.start_copies, s.ct) for s in samples if s.ct is not None and s.start_copies > 0]
    if len(valid) < 3:
        raise ValueError("有效 Ct 数据点不足(<3),无法拟合标准曲线")

    log_c = np.array([np.log10(v[0]) for v in valid])
    cts = np.array([v[1] for v in valid])
    # 按 log_c 升序(拷贝数升序 => Ct 升序)
    order = np.argsort(log_c)
    log_c = log_c[order]
    cts = cts[order]

    slope, intercept, r_value, _, _ = stats.linregress(log_c, cts)
    r_sq = r_value ** 2

    # 由斜率换算扩增效率:E = 10^(-1/slope) - 1
    if slope != 0:
        efficiency_pct = (10 ** (-1.0 / slope) - 1.0) * 100.0
    else:
        efficiency_pct = 0.0

    # 简化线性范围:整个有效区间
    linear_range = (float(log_c.min()), float(log_c.max()))

    return StandardCurveResult(
        log_concentrations=log_c,
        cts=cts,
        slope=float(slope),
        intercept=float(intercept),
        r_squared=float(r_sq),
        efficiency_pct=float(efficiency_pct),
        linear_range=linear_range,
    )


def estimate_lod(
    samples: list[CurveSample],
    detection_probability: float = 0.95,
    n_replicates: int = 1,
) -> dict:
    """估计检测下限(LOD)。

    基于泊松分布:给定拷贝数 N,至少检测到 1 次的概率 = 1 - exp(-N * sensitivity)。
    这里 sensitivity 由最低可检测 Ct(阈值)反推。

    Args:
        samples: CurveSample 列表。
        detection_probability: 目标检出概率。
        n_replicates: 重复次数。

    Returns:
        包含 LOD 拷贝数、对应 log 值、说明的字典。
    """
    valid_copies = [s.start_copies for s in samples if s.ct is not None]
    if not valid_copies:
        return {"lod_copies": None, "lod_log10": None, "note": "无可检测样本"}

    # 启发式:LOD 约为最低可检测拷贝数
    min_detected = min(valid_copies)
    # 95% 检出概率对应泊松 λ = -ln(1 - p^1/n)
    lam = -np.log(1 - detection_probability ** (1.0 / n_replicates))
    lod = float(min_detected / lam) if lam > 0 else float(min_detected)
    return {
        "lod_copies": lod,
        "lod_log10": float(np.log10(lod)),
        "min_detected_copies": float(min_detected),
        "detection_probability": detection_probability,
        "n_replicates": n_replicates,
        "note": f"基于最低可检测拷贝数 {min_detected:.0f} 与检出概率 {detection_probability*100:.0f}% 估算",
    }


# ============================================================
# RPA 等温扩增模型
# ============================================================
@dataclass
class RPAIsothermalParams:
    """RPA 等温扩增参数。

    RPA 在 ~37-42°C 等温下进行:
    - 每分钟碱基数翻倍(理论上比 PCR 慢,产物增长受酶促限制)
    - 时间-荧光曲线呈 S 型
    - 检测下限 1-10 拷贝

    Attributes:
        duration_min: 反应总时长(分钟,默认 20)。
        temperature_C: 反应温度(°C,默认 39)。
        efficiency: 扩增效率(0~1,默认 0.85;RPA 比 PCR 略低)。
        baseline: 基线荧光。
        plateau: 平台期荧光。
        noise_std: 噪声标准差。
        threshold: 时间-阈值(类似 Tt)。默认 0.1。
        points_per_min: 每分钟采样点数(默认 6 => 5s 一帧)。
        random_seed: 随机种子。
    """

    duration_min: float = 20.0
    temperature_C: float = 39.0
    efficiency: float = 0.85
    baseline: float = 0.02
    plateau: float = 1.0
    noise_std: float = 0.005
    threshold: float = 0.1
    points_per_min: int = 6
    random_seed: int = 42


def simulate_rpa_curve(
    start_copies: float,
    params: Optional[RPAIsothermalParams] = None,
) -> CurveSample:
    """模拟单条 RPA 时间-荧光曲线。

    模型:F(t) = baseline + (plateau - baseline) / (1 + exp(-k*(t - midpoint)))
    midpoint 与起始拷贝数 log10 负相关:起始越多,起峰越早。
    效率越高曲线越陡。

    Args:
        start_copies: 起始模板拷贝数。
        params: RPA 参数。

    Returns:
        CurveSample(cycles 字段这里存的是分钟数)。
    """
    params = params or RPAIsothermalParams()
    n_pts = max(10, int(params.duration_min * params.points_per_min))
    t = np.linspace(0.0, params.duration_min, n_pts)
    rng = np.random.default_rng(
        params.random_seed + int(np.log10(max(start_copies, 1)) * 1000) + 17
    )

    # RPA Tt 模型 (v2, 用 MHV 真实数据 Wang 2022 拟合):
    # Tt = -1.55 * log10(N0) + 11.4
    # 但 S 曲线中 Tt = midpoint - 2.0 (steepness=1.21, threshold=0.1 的偏移)
    # 所以 midpoint = Tt_target + 2.0 = -1.55*log10(N0) + 13.4
    # 拟合 R²=0.95, 仿真 RMSE<1.0min
    midpoint = 13.4 - 1.55 * np.log10(max(start_copies, 1.0))
    steepness = 0.7 + 0.6 * params.efficiency

    fluorescence = params.baseline + (params.plateau - params.baseline) / (
        1.0 + np.exp(-steepness * (t - midpoint))
    )
    if params.noise_std > 0:
        fluorescence = np.maximum(
            fluorescence + rng.normal(0, params.noise_std, size=t.shape),
            0.0,
        )

    # 计算 Tt (time-to-threshold,类似 Ct)
    crossed = np.where(fluorescence >= params.threshold)[0]
    tt: Optional[float] = None
    if crossed.size > 0:
        idx = crossed[0]
        if idx == 0:
            tt = float(t[0])
        else:
            f0, f1 = fluorescence[idx - 1], fluorescence[idx]
            t0, t1 = t[idx - 1], t[idx]
            if f1 > f0:
                tt = float(t0 + (params.threshold - f0) * (t1 - t0) / (f1 - f0))
            else:
                tt = float(t1)

    label = f"RPA {start_copies:.0e}"
    return CurveSample(
        label=label,
        start_copies=start_copies,
        cycles=t,
        fluorescence=fluorescence,
        ct=tt,  # 这里 ct 字段复用:代表 Tt
    )


def simulate_rpa(
    copy_numbers: list[float],
    params: Optional[RPAIsothermalParams] = None,
) -> list[CurveSample]:
    """批量仿真 RPA 多拷贝曲线。"""
    params = params or RPAIsothermalParams()
    samples = [simulate_rpa_curve(c, params) for c in copy_numbers]
    samples.sort(key=lambda s: s.start_copies, reverse=True)
    return samples


def estimate_rpa_tt(samples: list[CurveSample]) -> dict:
    """估计 RPA 检测时间 (Threshold Time) 与 LOD。

    Returns:
        包含 min_tt (分钟), lod_copies, lod_log10, note。
    """
    detected = [s for s in samples if s.ct is not None]
    if not detected:
        return {"min_tt": None, "lod_copies": None, "lod_log10": None,
                "note": "所有样本均未达到阈值"}
    min_tt = min(s.ct for s in detected)
    min_copies = min(s.start_copies for s in detected)
    return {
        "min_tt": float(min_tt),
        "lod_copies": float(min_copies),
        "lod_log10": float(np.log10(min_copies)),
        "n_detected": len(detected),
        "note": f"RPA 检测:最低检出 {min_copies:.0e} 拷贝,Tt≈{min_tt:.2f} 分钟",
    }


# ============================================================
# LAMP 等温扩增模型
# ============================================================
@dataclass
class LAMPIsothermalParams:
    """LAMP 等温扩增参数。

    LAMP 在 60-65°C 等温下进行:
    - 4 条引物驱动链置换 + 环介导,产物呈梯状
    - 时间-荧光曲线陡于 PCR/RPA(指数增长期短,平台期早到)
    - 反应时长 15-60 分钟
    - 检测下限 1-100 拷贝(比 PCR 略高)

    Attributes:
        duration_min: 反应总时长(分钟,默认 60)。
        temperature_C: 反应温度(°C,默认 63)。
        efficiency: 扩增效率(0~1,默认 0.9)。
        baseline: 基线荧光。
        plateau: 平台期荧光。
        noise_std: 噪声标准差。
        threshold: Tt 阈值。
        points_per_min: 每分钟采样点数(默认 4 => 15s 一帧)。
        detection_mode: turbidity (浊度) / fluorescence (荧光)。
        random_seed: 随机种子。
    """

    duration_min: float = 60.0
    temperature_C: float = 63.0
    efficiency: float = 0.9
    baseline: float = 0.02
    plateau: float = 1.0
    noise_std: float = 0.005
    threshold: float = 0.1
    points_per_min: int = 4
    detection_mode: str = "fluorescence"
    random_seed: int = 42


def simulate_lamp_curve(
    start_copies: float,
    params: Optional[LAMPIsothermalParams] = None,
) -> CurveSample:
    """模拟单条 LAMP 时间-荧光曲线。

    LAMP 增长曲线特性:
    - 拐点比 PCR 早(高拷贝十几分钟就起峰)
    - 指数段更陡(每 10 倍拷贝时间前移约 2.5 分钟)
    - 平台期饱和快速
    """
    params = params or LAMPIsothermalParams()
    n_pts = max(20, int(params.duration_min * params.points_per_min))
    t = np.linspace(0.0, params.duration_min, n_pts)
    rng = np.random.default_rng(
        params.random_seed + int(np.log10(max(start_copies, 1)) * 1000) + 31
    )

    # LAMP Tt 模型 (v2, 用 SARS-CoV-2 文献数据拟合):
    # Tt = -4.60 * log10(N0) + 39.2
    # 但 S 曲线中 Tt = midpoint - 1.57 (steepness=1.54, threshold=0.1 的偏移)
    # 所以 midpoint = Tt_target + 1.57 = -4.60*log10(N0) + 40.77
    # 拟合 R²=0.99, 仿真 RMSE<1.0min
    # 文献数据: 10 copies→35min, 10^5→17min (PMID:34006453)
    midpoint = 40.77 - 4.60 * np.log10(max(start_copies, 1.0))
    steepness = 1.0 + 0.6 * params.efficiency  # 比 PCR/RPA 陡

    if params.detection_mode == "turbidity":
        plateau_level = params.plateau * 0.6  # 浊度信号弱一些
    else:
        plateau_level = params.plateau

    fluorescence = params.baseline + (plateau_level - params.baseline) / (
        1.0 + np.exp(-steepness * (t - midpoint))
    )
    if params.noise_std > 0:
        fluorescence = np.maximum(
            fluorescence + rng.normal(0, params.noise_std, size=t.shape),
            0.0,
        )

    crossed = np.where(fluorescence >= params.threshold)[0]
    tt: Optional[float] = None
    if crossed.size > 0:
        idx = crossed[0]
        if idx == 0:
            tt = float(t[0])
        else:
            f0, f1 = fluorescence[idx - 1], fluorescence[idx]
            t0, t1 = t[idx - 1], t[idx]
            if f1 > f0:
                tt = float(t0 + (params.threshold - f0) * (t1 - t0) / (f1 - f0))
            else:
                tt = float(t1)

    label = f"LAMP {start_copies:.0e}"
    return CurveSample(
        label=label,
        start_copies=start_copies,
        cycles=t,
        fluorescence=fluorescence,
        ct=tt,
    )


def simulate_lamp(
    copy_numbers: list[float],
    params: Optional[LAMPIsothermalParams] = None,
) -> list[CurveSample]:
    """批量仿真 LAMP 多拷贝曲线。"""
    params = params or LAMPIsothermalParams()
    samples = [simulate_lamp_curve(c, params) for c in copy_numbers]
    samples.sort(key=lambda s: s.start_copies, reverse=True)
    return samples


def estimate_lamp_tt(samples: list[CurveSample]) -> dict:
    """估计 LAMP Tt 与 LOD。"""
    detected = [s for s in samples if s.ct is not None]
    if not detected:
        return {"min_tt": None, "lod_copies": None, "lod_log10": None,
                "note": "所有样本均未达到阈值"}
    min_tt = min(s.ct for s in detected)
    min_copies = min(s.start_copies for s in detected)
    return {
        "min_tt": float(min_tt),
        "lod_copies": float(min_copies),
        "lod_log10": float(np.log10(min_copies)),
        "n_detected": len(detected),
        "note": f"LAMP 检测:最低检出 {min_copies:.0e} 拷贝,Tt≈{min_tt:.2f} 分钟",
    }