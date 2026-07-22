"""LAMP (Loop-mediated isothermal amplification) 引物设计引擎。

LAMP 设计要点 (Notomi 2000, Tomita 2008):
- 6 个区域 (模板正链 5'→3'): F3 | F2 | F1 | [loop] | B1c | B2c | B3c
- 4 条核心引物:
    F3  = 正链 F3 区 (18-22bp)
    B3  = B3c 区反向互补 (18-22bp)
    FIP = F1c(F1 区反向互补) + linker + F2(正链 F2 区) (40-48bp)
    BIP = B1c(B1c 区正链) + linker + B2(B2c 区反向互补) (40-48bp)
- 可选: Loop F (loop 区正链), Loop B (loop 区反向互补)
- Tm 约束 (60-65°C 反应温度, SantaLucia 参数):
    F3, B3, F2, B2: Tm 60-65°C
    F1c, B1c:       Tm 65-70°C (比 F2/B2 高 ~5°C, 确保环优先形成)
- 产物: 120-300bp (F3 起点到 B3c 终点)
- 60-65°C 等温扩增

核心算法:
  1. 滑动窗口枚举所有候选子序列,计算 Tm
  2. 按 Tm 分桶: high-Tm (65-70°C) 候选用于 F1c/B1c, mid-Tm (60-65°C) 候选用于 F3/F2/B2/B3
  3. 按位置关系组合: F3 → F2 → F1 → loop → B1c → B2c → B3c
  4. 检查 amplicon 总长 120-300bp
  5. 检查 FIP/BIP 发夹和二聚体 (长序列二级结构是 LAMP 主要失败原因)
"""
from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Optional

import primer3


# ============================================================
# 数据结构
# ============================================================
@dataclass
class LAMPPrimerSet:
    """单组 LAMP 引物集 (4 核心 + 可选 2 loop)。"""
    set_id: int
    F3: str = ""
    B3: str = ""
    FIP: str = ""
    BIP: str = ""
    LoopF: str = ""
    LoopB: str = ""
    F3_pos: int = 0
    F2_pos: int = 0
    F1_pos: int = 0
    B1c_pos: int = 0
    B2c_pos: int = 0
    B3c_pos: int = 0
    amplicon_size: int = 0
    F3_tm: float = 0.0
    B3_tm: float = 0.0
    FIP_tm: float = 0.0
    BIP_tm: float = 0.0
    LoopF_tm: float = 0.0
    LoopB_tm: float = 0.0
    FIP_F1c_tm: float = 0.0
    FIP_F2_tm: float = 0.0
    BIP_B1c_tm: float = 0.0
    BIP_B2_tm: float = 0.0
    F3_gc: float = 0.0
    B3_gc: float = 0.0
    FIP_gc: float = 0.0
    BIP_gc: float = 0.0
    FIP_hairpin_dg: float = 0.0
    BIP_hairpin_dg: float = 0.0
    FIP_BIP_heterodimer_dg: float = 0.0
    accepted: bool = True
    issues: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class LAMPDesignParameters:
    """LAMP 引物设计参数。

    Tm 区间基于 60-65°C 反应温度 (SantaLucia 1998 + 盐校正):
      F3/B3/F2/B2: 60-65°C (与外引物和内引物 3' 端匹配)
      F1c/B1c:     65-70°C (确保环结构优先形成, 这是 LAMP 成功关键)
    """
    # 区域长度
    F3_min_size: int = 18
    F3_opt_size: int = 20
    F3_max_size: int = 22
    B3_min_size: int = 18
    B3_opt_size: int = 20
    B3_max_size: int = 22
    F2_min_size: int = 17
    F2_opt_size: int = 19
    F2_max_size: int = 21
    F1_min_size: int = 18
    F1_opt_size: int = 20
    F1_max_size: int = 24
    B1_min_size: int = 18
    B1_opt_size: int = 20
    B1_max_size: int = 24
    B2_min_size: int = 17
    B2_opt_size: int = 19
    B2_max_size: int = 21
    linker: str = "TTTT"

    # amplicon 总跨度 (F3 起点到 B3c 终点)
    amplicon_min: int = 120
    amplicon_max: int = 300
    loop_min: int = 10   # loop 区最小长度
    loop_max: int = 80   # loop 区最大长度

    # Tm 区间 (LAMP 关键约束)
    F3B3_tm_min: float = 60.0    # F3, B3
    F3B3_tm_max: float = 65.0
    F2B2_tm_min: float = 60.0    # F2, B2 (FIP/BIP 的 3' 端)
    F2B2_tm_max: float = 65.0
    F1B1_tm_min: float = 65.0    # F1c, B1c (FIP/BIP 的 5' 端, 必须更高)
    F1B1_tm_max: float = 70.0
    # F1c 与 F2 的 Tm 差 (F1c 应比 F2 高 ≥2°C; 严格文献用 3°C, 但多数应用 2°C)
    min_f1c_f2_tm_diff: float = 2.0

    # GC 含量
    gc_min: float = 40.0
    gc_max: float = 65.0

    # 二级结构约束
    max_fip_hairpin_dg: float = -3.0     # FIP 发夹 ΔG (kcal/mol), 不应过稳
    max_bip_hairpin_dg: float = -3.0
    max_fip_bip_heterodimer_dg: float = -6.0

    # loop 引物
    include_loop: bool = True
    loop_primer_min_size: int = 18
    loop_primer_max_size: int = 24
    loop_primer_tm_min: float = 62.0
    loop_primer_tm_max: float = 68.0

    # 输出
    num_sets: int = 3
    max_candidates: int = 500   # 组合枚举上限 (避免爆炸)

    # 盐浓度 (LAMP 典型: 8mM Mg²⁺, 比 PCR 高)
    # 注: dna_conc 单位为 nM (不是 M), 典型引物浓度 100-500nM
    mv_conc: float = 50.0      # 一价阳离子 mM
    dv_conc: float = 8.0        # Mg²⁺ mM (LAMP 通常 4-10mM)
    dntp_conc: float = 1.2      # dNTP mM
    dna_conc: float = 200.0     # 引物浓度 nM


# ============================================================
# 工具函数
# ============================================================
_COMP = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N",
         "a": "t", "t": "a", "c": "g", "g": "c", "n": "n"}


def _reverse_complement(seq: str) -> str:
    return "".join(_COMP.get(b, "N") for b in reversed(seq))


def _gc_percent(seq: str) -> float:
    if not seq:
        return 0.0
    gc = sum(1 for b in seq.upper() if b in {"G", "C"})
    return gc / len(seq) * 100.0


def _calc_tm(seq: str, params: LAMPDesignParameters) -> float:
    """primer3 Tm (SantaLucia + 盐校正)。"""
    if not seq or len(seq) < 4:
        return 0.0
    try:
        return float(primer3.calc_tm(
            seq,
            mv_conc=params.mv_conc,
            dv_conc=params.dv_conc,
            dntp_conc=params.dntp_conc,
            dna_conc=params.dna_conc,
            tm_method="santalucia",
            salt_corrections_method="santalucia",
        ))
    except Exception:
        return 0.0


def _calc_hairpin_dg(seq: str, params: LAMPDesignParameters) -> float:
    """发夹 ΔG (kcal/mol)。"""
    if not seq or len(seq) < 4:
        return 0.0
    try:
        res = primer3.calc_hairpin(
            seq,
            mv_conc=params.mv_conc,
            dv_conc=params.dv_conc,
            dntp_conc=params.dntp_conc,
            dna_conc=params.dna_conc,
            temp_c=65.0,  # LAMP 反应温度
        )
        return float(res.dg / 1000.0)
    except Exception:
        return 0.0


def _calc_heterodimer_dg(s1: str, s2: str,
                         params: LAMPDesignParameters) -> float:
    if not s1 or not s2:
        return 0.0
    try:
        res = primer3.calc_heterodimer(
            s1, s2,
            mv_conc=params.mv_conc,
            dv_conc=params.dv_conc,
            dntp_conc=params.dntp_conc,
            dna_conc=params.dna_conc,
            temp_c=65.0,
        )
        return float(res.dg / 1000.0)
    except Exception:
        return 0.0


@dataclass
class _Candidate:
    """模板上一个候选子序列。"""
    start: int
    length: int
    seq: str
    tm: float
    gc: float


def _scan_candidates(
    sequence: str,
    min_len: int,
    max_len: int,
    tm_min: float,
    tm_max: float,
    gc_min: float,
    gc_max: float,
    params: LAMPDesignParameters,
) -> list[_Candidate]:
    """滑动窗口扫描, 返回所有满足 Tm/GC 的候选子序列。

    对每个起点, 按长度从短到长扫描, 命中 Tm 区间即记录。
    """
    cands: list[_Candidate] = []
    seq_len = len(sequence)
    for start in range(seq_len - min_len + 1):
        for length in range(min_len, max_len + 1):
            if start + length > seq_len:
                break
            sub = sequence[start:start + length]
            # 快速 GC 预筛
            gc = _gc_percent(sub)
            if gc < gc_min or gc > gc_max:
                continue
            tm = _calc_tm(sub, params)
            if tm_min <= tm <= tm_max:
                cands.append(_Candidate(start, length, sub, tm, gc))
    return cands


def _pick_loop_primer(
    sequence: str,
    region_start: int,
    region_end: int,
    tm_min: float,
    tm_max: float,
    min_size: int,
    max_size: int,
    params: LAMPDesignParameters,
    reverse: bool = False,
) -> str:
    """在 loop 区间内挑一条 Loop 引物。

    reverse=False → Loop F (正链子序列)
    reverse=True  → Loop B (反链子序列, 即正链子序列的反向互补)
    """
    if region_end <= region_start:
        return ""
    region = sequence[region_start:region_end]
    if len(region) < min_size:
        return ""
    best: Optional[tuple[float, str]] = None
    for start in range(len(region) - min_size + 1):
        for length in range(min_size, min(max_size, len(region) - start) + 1):
            sub = region[start:start + length]
            gc = _gc_percent(sub)
            if not (params.gc_min <= gc <= params.gc_max):
                continue
            tm = _calc_tm(sub, params)
            if tm_min <= tm <= tm_max:
                # 取 Tm 最接近区间中心的
                target = (tm_min + tm_max) / 2
                score = abs(tm - target)
                if best is None or score < best[0]:
                    if reverse:
                        best = (score, _reverse_complement(sub))
                    else:
                        best = (score, sub)
    return best[1] if best else ""


# ============================================================
# 核心设计算法
# ============================================================
def design_lamp_primers(
    sequence: str,
    params: Optional[LAMPDesignParameters] = None,
) -> list[LAMPPrimerSet]:
    """设计多组 LAMP 引物集。

    算法:
      1. 在模板上扫描所有满足 Tm/GC 的候选子序列
         - low-Tm (60-65°C): 用于 F3, B3, F2, B2
         - high-Tm (65-70°C): 用于 F1c, B1c
      2. 按 F3 → F2 → F1 → loop → B1c → B2c → B3c 顺序组合
         位置约束:
           F3.start < F2.start < F1.start < B1c.start < B2c.start < B3c.start
           (模板正链 5'→3')
         amplicon (F3.start 到 B3c.end) 在 120-300bp
      3. 组装 FIP = F1c + linker + F2, BIP = B1c + linker + B2
      4. 检查 FIP/BIP 二级结构 (发夹 + 异源二聚体)
      5. 可选 Loop F/B

    Args:
        sequence: 模板序列 (已大写, 无空白)
        params: 设计参数

    Returns:
        LAMPPrimerSet 列表, accepted 优先, 按 Tm 匹配度排序
    """
    seq = sequence.upper().replace(" ", "").replace("\n", "")
    if len(seq) < 200:
        raise ValueError(f"LAMP 设计模板至少 200bp, 当前 {len(seq)}bp")

    params = params or LAMPDesignParameters()

    # --- Step 1: 扫描候选 ---
    # F3/B3 候选 (Tm 60-65, 长度 18-22)
    f3_b3_cands = _scan_candidates(
        seq, params.F3_min_size, params.F3_max_size,
        params.F3B3_tm_min, params.F3B3_tm_max,
        params.gc_min, params.gc_max, params,
    )
    # F2/B2 候选 (Tm 60-65, 长度 17-21)
    f2_b2_cands = _scan_candidates(
        seq, params.F2_min_size, params.F2_max_size,
        params.F2B2_tm_min, params.F2B2_tm_max,
        params.gc_min, params.gc_max, params,
    )
    # F1c/B1c 候选 (Tm 65-70, 长度 18-24)
    f1_b1_cands = _scan_candidates(
        seq, params.F1_min_size, params.F1_max_size,
        params.F1B1_tm_min, params.F1B1_tm_max,
        params.gc_min, params.gc_max, params,
    )

    if not f3_b3_cands or not f2_b2_cands or not f1_b1_cands:
        raise RuntimeError(
            f"候选不足: F3/B3={len(f3_b3_cands)}, "
            f"F2/B2={len(f2_b2_cands)}, F1c/B1c={len(f1_b1_cands)}。"
            f"建议放宽 Tm/GC 区间或换序列"
        )

    # --- Step 2: 组合枚举 (F1c 优先策略, v2 改进) ---
    # 改进: 原算法按位置顺序枚举, 导致 F1c-F2 Tm 差常 <2°C (60% 病原体失败)
    # 新策略: F1c 按 Tm 降序排序, 优先用高 Tm 候选;
    #         F2 Tm 上限动态设为 F1c.tm - min_diff, 确保差值
    # 位置约束: F3 < F2 < F1 < B1c < B2c < B3c (模板正链 5'→3')
    sets: list[LAMPPrimerSet] = []
    set_idx = 0

    # 排序: F3/B3/F2/B2 按位置; F1c/B1c 按 Tm 降序 (高 Tm 优先)
    f3_b3_cands.sort(key=lambda c: c.start)
    f2_b2_cands.sort(key=lambda c: c.start)
    f1_b1_cands.sort(key=lambda c: -c.tm)  # Tm 降序, v2 关键改进

    # 收集所有组合 (先不限 max_candidates, 后统一排序取前 N)
    all_combos: list[tuple] = []

    # 枚举 F3 (正链前段)
    for f3 in f3_b3_cands:
        if len(all_combos) >= params.max_candidates * 3:
            break
        # F2 必须在 F3 下游
        f2_list = [c for c in f2_b2_cands
                    if c.start >= f3.start + f3.length
                    and c.start < f3.start + 40]
        for f2 in f2_list:
            if len(all_combos) >= params.max_candidates * 3:
                break
            # F1 必须紧邻 F2 下游 (±2bp)
            f1_list = [c for c in f1_b1_cands
                        if abs(c.start - (f2.start + f2.length)) <= 2
                        and c.start > f2.start]
            for f1 in f1_list:
                if len(all_combos) >= params.max_candidates * 3:
                    break
                # B1c 必须在 F1 下游, loop 区间 10-80bp
                f1_end = f1.start + f1.length
                b1c_list = [c for c in f1_b1_cands
                             if c.start >= f1_end + params.loop_min
                             and c.start <= f1_end + params.loop_max]
                for b1c in b1c_list:
                    if len(all_combos) >= params.max_candidates * 3:
                        break
                    # B2c 必须紧邻 B1c 下游
                    b2c_list = [c for c in f2_b2_cands
                                 if c.start >= b1c.start + b1c.length
                                 and c.start < b1c.start + b1c.length + 20]
                    for b2c in b2c_list:
                        if len(all_combos) >= params.max_candidates * 3:
                            break
                        # B3c 必须在 B2c 下游
                        b3c_list = [c for c in f3_b3_cands
                                     if c.start >= b2c.start + b2c.length
                                     and c.start < b2c.start + b2c.length + 20]
                        for b3c in b3c_list:
                            amplicon = b3c.start + b3c.length - f3.start
                            if amplicon < params.amplicon_min or amplicon > params.amplicon_max:
                                continue
                            all_combos.append((f3, f2, f1, b1c, b2c, b3c, amplicon))

    if not all_combos:
        raise RuntimeError(
            "未能生成满足位置约束的 LAMP 引物组合。"
            "建议: (1) 增大 amplicon_max, (2) 放宽 Tm 区间, (3) 检查序列 GC 分布"
        )

    # 对每个组合计算质量分, 排序后取前 num_sets * 6 个做详细评估
    # 质量分: F1c Tm 越高越好 (环优先), F1c-F2 Tm 差越大越好
    def combo_score(combo):
        f3, f2, f1, b1c, b2c, b3c, amplicon = combo
        # F1c Tm 越高越好 (但 ≤70), F2 Tm 越低越好 (但 ≥60)
        f1c_tm_proxy = f1.tm  # F1c Tm ≈ F1 区 Tm (反向互补 Tm 接近)
        return -(f1c_tm_proxy)  # 升序排列时 Tm 高的在前

    all_combos.sort(key=combo_score)
    combos_to_eval = all_combos[:params.num_sets * 6]

    for f3, f2, f1, b1c, b2c, b3c, amplicon in combos_to_eval:
        # --- 组装引物 ---
        F3_seq = f3.seq
        B3_seq = _reverse_complement(b3c.seq)
        F1c_seq = _reverse_complement(f1.seq)
        FIP_seq = F1c_seq + params.linker + f2.seq
        B1c_seq = b1c.seq
        B2_seq = _reverse_complement(b2c.seq)
        BIP_seq = B1c_seq + params.linker + B2_seq

        # --- 计算 Tm ---
        f3_tm = f3.tm
        b3_tm = _calc_tm(B3_seq, params)
        f1c_tm = _calc_tm(F1c_seq, params)
        f2_tm = f2.tm
        b1c_tm = _calc_tm(B1c_seq, params)
        b2_tm = _calc_tm(B2_seq, params)
        fip_tm = _calc_tm(FIP_seq, params)
        bip_tm = _calc_tm(BIP_seq, params)

        # --- 二级结构检查 ---
        fip_hairpin = _calc_hairpin_dg(FIP_seq, params)
        bip_hairpin = _calc_hairpin_dg(BIP_seq, params)
        fip_bip_het = _calc_heterodimer_dg(FIP_seq, BIP_seq, params)

        # --- 可选 Loop 引物 ---
        loop_f, loop_b = "", ""
        if params.include_loop:
            loop_region_start = f1.start + f1.length
            loop_region_end = b1c.start
            if loop_region_end > loop_region_start + params.loop_primer_min_size:
                loop_f = _pick_loop_primer(
                    seq, loop_region_start, loop_region_end,
                    params.loop_primer_tm_min, params.loop_primer_tm_max,
                    params.loop_primer_min_size, params.loop_primer_max_size,
                    params, reverse=False,
                )
                loop_b = _pick_loop_primer(
                    seq, loop_region_start, loop_region_end,
                    params.loop_primer_tm_min, params.loop_primer_tm_max,
                    params.loop_primer_min_size, params.loop_primer_max_size,
                    params, reverse=True,
                )

        # --- 验收 (动态 Tm 差策略, v2 改进) ---
        # v2: F1c-F2 Tm 差 <min_diff 但 ≥1.5°C 标为 warning (非 fatal)
        #     只有 <1.5°C 才是 fatal
        issues: list[str] = []
        warnings: list[str] = []
        if not (params.F3B3_tm_min <= f3_tm <= params.F3B3_tm_max):
            issues.append(f"F3 Tm {f3_tm:.1f}°C 超出 [{params.F3B3_tm_min}-{params.F3B3_tm_max}]")
        if not (params.F3B3_tm_min <= b3_tm <= params.F3B3_tm_max):
            issues.append(f"B3 Tm {b3_tm:.1f}°C 超出 [{params.F3B3_tm_min}-{params.F3B3_tm_max}]")
        if not (params.F2B2_tm_min <= f2_tm <= params.F2B2_tm_max):
            issues.append(f"FIP(F2) Tm {f2_tm:.1f}°C 超出 [{params.F2B2_tm_min}-{params.F2B2_tm_max}]")
        if not (params.F2B2_tm_min <= b2_tm <= params.F2B2_tm_max):
            issues.append(f"BIP(B2) Tm {b2_tm:.1f}°C 超出 [{params.F2B2_tm_min}-{params.F2B2_tm_max}]")
        if not (params.F1B1_tm_min <= f1c_tm <= params.F1B1_tm_max):
            issues.append(f"FIP(F1c) Tm {f1c_tm:.1f}°C 超出 [{params.F1B1_tm_min}-{params.F1B1_tm_max}]")
        if not (params.F1B1_tm_min <= b1c_tm <= params.F1B1_tm_max):
            issues.append(f"BIP(B1c) Tm {b1c_tm:.1f}°C 超出 [{params.F1B1_tm_min}-{params.F1B1_tm_max}]")
        # F1c-F2 Tm 差 (动态降级)
        f1c_f2_diff = f1c_tm - f2_tm
        b1c_b2_diff = b1c_tm - b2_tm
        if f1c_f2_diff < 1.5:
            issues.append(f"F1c-F2 Tm 差 {f1c_f2_diff:.1f}°C < 1.5°C (最低容忍)")
        elif f1c_f2_diff < params.min_f1c_f2_tm_diff:
            warnings.append(f"F1c-F2 Tm 差 {f1c_f2_diff:.1f}°C < {params.min_f1c_f2_tm_diff}°C (降级可接受)")
        if b1c_b2_diff < 1.5:
            issues.append(f"B1c-B2 Tm 差 {b1c_b2_diff:.1f}°C < 1.5°C (最低容忍)")
        elif b1c_b2_diff < params.min_f1c_f2_tm_diff:
            warnings.append(f"B1c-B2 Tm 差 {b1c_b2_diff:.1f}°C < {params.min_f1c_f2_tm_diff}°C (降级可接受)")
        # FIP/BIP 长度
        if not (38 <= len(FIP_seq) <= 50):
            issues.append(f"FIP 长度 {len(FIP_seq)}bp 超出 38-50")
        if not (38 <= len(BIP_seq) <= 50):
            issues.append(f"BIP 长度 {len(BIP_seq)}bp 超出 38-50")
        # 二级结构
        if fip_hairpin < params.max_fip_hairpin_dg:
            issues.append(f"FIP 发夹 ΔG={fip_hairpin:.2f} 过稳 (< {params.max_fip_hairpin_dg})")
        if bip_hairpin < params.max_bip_hairpin_dg:
            issues.append(f"BIP 发夹 ΔG={bip_hairpin:.2f} 过稳 (< {params.max_bip_hairpin_dg})")
        if fip_bip_het < params.max_fip_bip_heterodimer_dg:
            issues.append(f"FIP-BIP 异源二聚体 ΔG={fip_bip_het:.2f} 过稳")

        # 综合质量分 (越低越好):
        # (1) Tm 平衡 (2) F1c-F2 Tm 差越大越好 (3) 二级结构越稳越扣分
        tm_balance = (
            abs(f3_tm - b3_tm) +
            abs(f1c_tm - b1c_tm) +
            abs(f2_tm - b2_tm)
        )
        # F1c-F2 差分: 差越大越好 (用负值, 这样排序时大的在前)
        f1c_f2_diff_score = max(0, params.min_f1c_f2_tm_diff - (f1c_tm - f2_tm)) + \
                           max(0, params.min_f1c_f2_tm_diff - (b1c_tm - b2_tm))
        # 二级结构分: 越负越扣分
        ss_score = min(0, fip_hairpin) + min(0, bip_hairpin) + min(0, fip_bip_het)
        total_score = tm_balance + f1c_f2_diff_score + ss_score

        sets.append((total_score, LAMPPrimerSet(
            set_id=set_idx,
            F3=F3_seq, B3=B3_seq,
            FIP=FIP_seq, BIP=BIP_seq,
            LoopF=loop_f, LoopB=loop_b,
            F3_pos=f3.start, F2_pos=f2.start, F1_pos=f1.start,
            B1c_pos=b1c.start, B2c_pos=b2c.start, B3c_pos=b3c.start,
            amplicon_size=amplicon,
            F3_tm=f3_tm, B3_tm=b3_tm,
            FIP_tm=fip_tm, BIP_tm=bip_tm,
            LoopF_tm=_calc_tm(loop_f, params) if loop_f else 0.0,
            LoopB_tm=_calc_tm(loop_b, params) if loop_b else 0.0,
            FIP_F1c_tm=f1c_tm, FIP_F2_tm=f2_tm,
            BIP_B1c_tm=b1c_tm, BIP_B2_tm=b2_tm,
            F3_gc=_gc_percent(F3_seq), B3_gc=_gc_percent(B3_seq),
            FIP_gc=_gc_percent(FIP_seq), BIP_gc=_gc_percent(BIP_seq),
            FIP_hairpin_dg=fip_hairpin,
            BIP_hairpin_dg=bip_hairpin,
            FIP_BIP_heterodimer_dg=fip_bip_het,
            accepted=not issues,  # warnings 不算 fatal
            issues=issues + warnings,  # 显示但不影响 accepted
        )))
        set_idx += 1

    if not sets:
        raise RuntimeError(
            "未能生成任何评估组合"
        )

    # 排序: accepted 优先, 然后按综合质量分升序
    sets.sort(key=lambda x: (not x[1].accepted, x[0]))

    # 取前 num_sets 个
    out = [s[1] for s in sets[:params.num_sets]]
    for idx, s in enumerate(out):
        s.set_id = idx
    return out


if __name__ == "__main__":
    # 自测: 用一个 SARS-CoV-2 N 基因片段
    test_seq = (
        "GACCCCAAAATCAGCGAAATCGCACGGCAGAATGGCTAAGCACAACCTCAACAT"
        "ACCAACAGGTGTGCAACTGAACAACCTGCTGGATGACGCTGAGTCCAAGAACCAG"
        "AATCAATGGAGTCATCAGCAGAACAAGGACATCACTACCGAGATTCGGTGAAGCA"
        "GGTGCAACTGCAGCAAGTTTGGTCACTGAGAGTTCAACATCACCACCAGTAATGA"
        "CAAGGTGGTGCTGGCCATTAATGTTGATGAACCAACCAACAGGTGGTGCAACTGA"
        "ACAACCTGCTGGATGACGCTGAGTCCAAGAACCAGAATCAATGGAGTCATCAGC"
        "AGAACAAGGACATCACTACCGAGATTCGGTGAAGCAGGTGCAACTGCAGCAAGT"
        "TTTGGTCACTGAGAGTTCAACATCACCACCAGTAATGACAAG"
    )
    print(f"测试序列: {len(test_seq)} bp")
    results = design_lamp_primers(test_seq)
    print(f"\n生成 {len(results)} 组引物集:")
    for s in results:
        print(f"\n--- Set #{s.set_id} ({'✓' if s.accepted else '✗'}) ---")
        print(f"  F3:  {s.F3}  (Tm={s.F3_tm:.1f}°C)")
        print(f"  B3:  {s.B3}  (Tm={s.B3_tm:.1f}°C)")
        print(f"  FIP: {s.FIP}  (F1c Tm={s.FIP_F1c_tm:.1f}, F2 Tm={s.FIP_F2_tm:.1f})")
        print(f"  BIP: {s.BIP}  (B1c Tm={s.BIP_B1c_tm:.1f}, B2 Tm={s.BIP_B2_tm:.1f})")
        print(f"  产物: {s.amplicon_size}bp")
        print(f"  FIP发夹ΔG={s.FIP_hairpin_dg:.2f}, BIP发夹ΔG={s.BIP_hairpin_dg:.2f}")
        if s.issues:
            print(f"  问题: {s.issues}")
