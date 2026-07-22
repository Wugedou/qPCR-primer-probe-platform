"""RPA (Recombinase Polymerase Amplification) 引物设计引擎。

RPA 特点:
- 等温扩增(37~42°C)
- 引物长度 30~35 bp(显著长于 PCR)
- 产物 100~500 bp
- 不需要严格的 Tm 约束
- 避免引物 3' 端发夹与二聚体(对 RPA 极重要)
- 可选 exo-probe (46~52bp) 含 THF/反义 CITE 修饰位点
"""
from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Optional

import primer3


@dataclass
class RPAProbe:
    """RPA exo-probe / nfo-probe 设计结果(可选)。

    Attributes:
        probe_seq: 探针主序列(5'→3',含 dT-FAM / THF / dT-BHQ1 等修饰位点)。
        target_site: 模板上的靶位点 (start, length)。
        reporter_at: 5' 端 dT-FAM(荧光基团)位置。
        quencher_at: 3' 端 dT-BHQ(淬灭基团)位置。
        thf_at: 中间 THF(四氢呋喃)位置 — exonuclease 切割位点。
    """

    probe_seq: str = ""
    target_site: tuple[int, int] = (0, 0)
    reporter_at: int = 0
    quencher_at: int = 0
    thf_at: int = 0
    note: str = ""

    def to_dict(self) -> dict:
        d = asdict(self)
        d["target_site"] = list(self.target_site)
        return d


@dataclass
class RPAPrimerPair:
    """单组 RPA 引物对(或加探针)的设计结果。

    Attributes:
        pair_id: 组号(0-based)。
        forward_seq: 正向引物(5'→3')。
        reverse_seq: 反向引物(5'→3')。
        forward_pos: 正向引物在模板上的起始位置(0-based)。
        reverse_pos: 反向引物(反向互补后)在模板上的起始位置(0-based)。
        product_size: 产物大小(bp)。
        forward_length: 正向引物长度。
        reverse_length: 反向引物长度。
        forward_gc: 正向引物 GC%。
        reverse_gc: 反向引物 GC%。
        forward_tm: 正向引物 Tm(估算)。
        reverse_tm: 反向引物 Tm(估算)。
        forward_3prime_hairpin_dg: 正向 3' 末端发夹 ΔG(kcal/mol)。
        reverse_3prime_hairpin_dg: 反向 3' 末端发夹 ΔG。
        heterodimer_dg: 正反向异源二聚体 ΔG。
        accepted: 是否满足 RPA 基本可用条件。
        issues: 问题列表(空表示通过)。
        probe: 探针对象(若设计则为非空)。
    """

    pair_id: int
    forward_seq: str = ""
    reverse_seq: str = ""
    forward_pos: int = 0
    reverse_pos: int = 0
    product_size: int = 0
    forward_length: int = 0
    reverse_length: int = 0
    forward_gc: float = 0.0
    reverse_gc: float = 0.0
    forward_tm: float = 0.0
    reverse_tm: float = 0.0
    forward_3prime_hairpin_dg: float = 0.0
    reverse_3prime_hairpin_dg: float = 0.0
    heterodimer_dg: float = 0.0
    accepted: bool = True
    issues: list[str] = field(default_factory=list)
    probe: Optional[RPAProbe] = None

    def to_dict(self) -> dict:
        d = asdict(self)
        d["probe"] = self.probe.to_dict() if self.probe else None
        return d


@dataclass
class RPADesignParameters:
    """RPA 引物设计参数。

    Attributes:
        forward_min_size: 正向引物最小长度(默认 30)。
        forward_opt_size: 正向引物最佳长度(默认 32)。
        forward_max_size: 正向引物最大长度(默认 35)。
        reverse_min_size: 反向引物最小长度(默认 30)。
        reverse_opt_size: 反向引物最佳长度(默认 32)。
        reverse_max_size: 反向引物最大长度(默认 35)。
        product_size_min: 产物最小长度(默认 100)。
        product_size_max: 产物最大长度(默认 500)。
        gc_min: 引物 GC% 下限(默认 30)。
        gc_max: 引物 GC% 上限(默认 70)。
        max_3prime_hairpin_dg: 3' 端发夹 ΔG 上限(更负 = 更稳定),默认 -2.0。
        max_heterodimer_dg: 正反向异源二聚体 ΔG 上限,默认 -6.0。
        num_return: 返回候选组数。
        design_probe: 是否同时设计 exo-probe。
        probe_min_size: 探针最小长度(默认 46)。
        probe_max_size: 探针最大长度(默认 52)。
    """

    forward_min_size: int = 30
    forward_opt_size: int = 32
    forward_max_size: int = 35
    reverse_min_size: int = 30
    reverse_opt_size: int = 32
    reverse_max_size: int = 35
    # v2: 产物大小按 Daher 2016 / TwistDx 官方规范 80-400bp (opt 100-200)
    product_size_min: int = 80
    product_size_max: int = 400
    product_size_opt_min: int = 100   # 最优区间下限
    product_size_opt_max: int = 200   # 最优区间上限
    gc_min: float = 30.0
    gc_max: float = 70.0
    # v2: 长串单核苷酸检查 (Daher 2016: 避免 >4-5nt 连续)
    max_homopolymer: int = 5  # >5nt 连续同碱基 = 警告
    max_3prime_hairpin_dg: float = -2.0
    max_heterodimer_dg: float = -6.0
    num_return: int = 5
    design_probe: bool = False
    probe_min_size: int = 46
    probe_max_size: int = 52

    def to_primer3_dict(self, sequence: str) -> dict:
        """转换为 primer3 设计字典。

        注意:primer3-py 对内部 oligo 长度有内置上限(~36bp),
        RPA exo-probe (46-52bp) 超出此限制。因此 primer3 只设计引物对,
        exo-probe 由本模块自定义逻辑生成。
        """
        seq_args = {
            "SEQUENCE_ID": "rpa_target",
            "SEQUENCE_TEMPLATE": sequence,
        }
        global_args = {
            "PRIMER_TASK": "generic",
            "PRIMER_PICK_LEFT_PRIMER": 1,
            "PRIMER_PICK_RIGHT_PRIMER": 1,
            "PRIMER_PICK_INTERNAL_OLIGO": 0,  # 探针走自定义逻辑
            "PRIMER_OPT_SIZE": self.forward_opt_size,
            "PRIMER_MIN_SIZE": self.forward_min_size,
            "PRIMER_MAX_SIZE": self.forward_max_size,
            # RPA 不强制 Tm,但 primer3 要求设置,放宽窗口
            "PRIMER_OPT_TM": 50.0,
            "PRIMER_MIN_TM": 35.0,
            "PRIMER_MAX_TM": 75.0,
            "PRIMER_PAIR_MAX_DIFF_TM": 25.0,
            "PRIMER_MIN_GC": self.gc_min,
            "PRIMER_MAX_GC": self.gc_max,
            "PRIMER_PRODUCT_SIZE_RANGE": [[self.product_size_min, self.product_size_max]],
            "PRIMER_NUM_RETURN": self.num_return,
            "PRIMER_MAX_POLY_X": 5,
            "PRIMER_MAX_NS_ACCEPTED": 0,
            "PRIMER_MIN_THREE_PRIME_DISTANCE": 3,
            "PRIMER_MAX_SELF_ANY": 12,
            "PRIMER_MAX_SELF_END": 6,
            "PRIMER_PAIR_MAX_COMPL_ANY": 12,
            "PRIMER_PAIR_MAX_COMPL_END": 6,
            "PRIMER_MAX_HAIRPIN_TH": 55.0,
        }
        return {"SEQUENCE": seq_args, "GLOBAL": global_args}


def _gc_percent(seq: str) -> float:
    """计算 GC%。"""
    if not seq:
        return 0.0
    gc = sum(1 for b in seq.upper() if b in {"G", "C"})
    return gc / len(seq) * 100.0


def _max_homopolymer_run(seq: str) -> int:
    """返回最长连续同碱基数 (Daher 2016: 避免 >4-5nt 连续)。"""
    if not seq:
        return 0
    seq = seq.upper()
    max_run = 1
    cur = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i-1]:
            cur += 1
            max_run = max(max_run, cur)
        else:
            cur = 1
    return max_run


def _has_long_repeats(seq: str, min_repeat_unit: int = 3, min_repeats: int = 3) -> bool:
    """检查是否有大量小重复 (Daher 2016: 避免 small repeats)。"""
    if not seq or len(seq) < min_repeat_unit * min_repeats:
        return False
    seq = seq.upper()
    for unit_len in range(min_repeat_unit, min(6, len(seq) // min_repeats) + 1):
        for start in range(len(seq) - unit_len * min_repeats + 1):
            unit = seq[start:start + unit_len]
            count = 1
            pos = start + unit_len
            while pos + unit_len <= len(seq) and seq[pos:pos + unit_len] == unit:
                count += 1
                pos += unit_len
            if count >= min_repeats:
                return True
    return False


def _three_prime_hairpin(seq: str, window: int = 10, salt_kw: Optional[dict] = None) -> float:
    """检查 3' 末端 window 区域的发夹 ΔG。"""
    if not seq or len(seq) < window:
        return 0.0
    three_prime = seq[-window:]
    kw = salt_kw or {}
    try:
        res = primer3.calc_hairpin(three_prime, **kw)
        return float(res.dg / 1000.0)
    except Exception:
        return 0.0


def _heterodimer(s1: str, s2: str, salt_kw: Optional[dict] = None) -> float:
    """计算两条引物的异源二聚体 ΔG。"""
    if not s1 or not s2:
        return 0.0
    kw = salt_kw or {}
    try:
        res = primer3.calc_heterodimer(s1, s2, **kw)
        return float(res.dg / 1000.0)
    except Exception:
        return 0.0


def _tm(seq: str, salt_kw: Optional[dict] = None) -> float:
    """估算 Tm。等温 RPA 采用较低温度(38°C)校正。"""
    if not seq:
        return 0.0
    kw = salt_kw or {}
    kw_full = {
        "mv_conc": 50.0, "dv_conc": 3.0, "dntp_conc": 0.8, "dna_conc": 200.0,
        "tm_method": "santalucia", "salt_corrections_method": "santalucia",
    }
    kw_full.update(kw)
    try:
        return float(primer3.calc_tm(seq, **kw_full))
    except Exception:
        try:
            return float(primer3.calc_tm(seq))
        except Exception:
            return 0.0


def _design_exo_probe(
    sequence: str,
    f_seq: str, f_pos: int,
    r_seq: str, r_pos: int,
    min_size: int, max_size: int,
    salt_kw: Optional[dict] = None,
) -> Optional[RPAProbe]:
    """在 F 与 R 引物之间的产物区间内,挑一段作为 exo-probe 区域。

    probe 内含 dT-FAM / THF / dT-BHQ 三个修饰位点;这里仅在序列中用小写标记
    位置,实际实验中替换为修饰碱基。
    """
    f_len = len(f_seq)
    r_len = len(r_seq)
    if r_pos <= f_pos + f_len:
        return None
    inner_start = f_pos + f_len
    inner_end = r_pos
    inner_len = inner_end - inner_start
    if inner_len < max_size:
        return None
    # 取中间一段作为探针区
    probe_len = max(min_size, min(max_size, inner_len // 3))
    probe_start = inner_start + max(0, (inner_len - probe_len) // 2)
    probe_seq = sequence[probe_start:probe_start + probe_len].upper()
    if len(probe_seq) < min_size:
        return None

    # 三段修饰位点:T(thf 前后)各一,中间 THF 1 个
    # 3' 端封闭 (Daher 2016 / TwistDx 规范): 磷酸基团或双脱氧核苷酸
    half = len(probe_seq) // 2
    thf_at = half
    reporter_at = 2  # 第二个碱基位置建议 dT-FAM
    quencher_at = len(probe_seq) - 2  # 倒数第二个位置建议 dT-BHQ
    return RPAProbe(
        probe_seq=probe_seq,
        target_site=(probe_start, len(probe_seq)),
        reporter_at=reporter_at,
        quencher_at=quencher_at,
        thf_at=thf_at,
        note=("实验中: reporter_at 替换为 dT-FAM, quencher_at 替换为 dT-BHQ1, "
              "thf_at 替换为 THF(四氢呋喃), 3' 端加磷酸基团(PO4)或双脱氧核苷酸(ddC)封闭"
              "防止非特异延伸 (Daher 2016/TwistDx 规范)"),
    )


def design_rpa_primers(
    sequence: str,
    params: Optional[RPADesignParameters] = None,
) -> list[RPAPrimerPair]:
    """设计多组 RPA 引物对(可选 exo-probe)。

    Args:
        sequence: 模板序列(已大写,无空白)。
        params: 设计参数;None 使用默认。

    Returns:
        RPAPrimerPair 列表,按产品 acceptance/penalty 排序。

    Raises:
        ValueError: 序列太短。
        RuntimeError: primer3 内部错误或无合格引物。
    """
    if len(sequence) < 200:
        raise ValueError("RPA 设计模板序列至少需要 200 bp")

    params = params or RPADesignParameters()
    args = params.to_primer3_dict(sequence)
    try:
        result = primer3.bindings.design_primers(args["SEQUENCE"], args["GLOBAL"])
    except Exception as exc:
        raise RuntimeError(f"primer3 RPA 设计失败: {exc}") from exc

    count = int(result.get("PRIMER_PAIR_NUM_RETURNED", 0))
    if count == 0:
        raise RuntimeError("未找到满足条件的 RPA 引物对,请放宽参数或检查序列")

    salt_kw = {"mv_conc": 50.0, "dv_conc": 3.0, "dntp_conc": 0.8, "dna_conc": 200.0}
    pairs: list[RPAPrimerPair] = []

    for i in range(count):
        f_seq = result.get(f"PRIMER_LEFT_{i}_SEQUENCE", "")
        r_seq = result.get(f"PRIMER_RIGHT_{i}_SEQUENCE", "")
        f_pos_str = result.get(f"PRIMER_LEFT_{i}", "")
        r_pos_str = result.get(f"PRIMER_RIGHT_{i}", "")
        # primer3 returns "start,length"
        try:
            f_pos = int(str(f_pos_str).split(",")[0])
            f_len = int(str(f_pos_str).split(",")[1])
        except (ValueError, IndexError):
            f_pos, f_len = 0, len(f_seq)
        try:
            r_pos = int(str(r_pos_str).split(",")[0])
            r_len = int(str(r_pos_str).split(",")[1])
        except (ValueError, IndexError):
            r_pos, r_len = 0, len(r_seq)

        # 反向引物位置 primer3 的 r_pos 是反向引物最右端位置;此处保留 0-based 起点
        product_size = int(result.get(f"PRIMER_PAIR_{i}_PRODUCT_SIZE", 0))

        f_gc = _gc_percent(f_seq)
        r_gc = _gc_percent(r_seq)
        f_tm = _tm(f_seq, salt_kw)
        r_tm = _tm(r_seq, salt_kw)

        f_3hp = _three_prime_hairpin(f_seq, window=10, salt_kw=salt_kw)
        r_3hp = _three_prime_hairpin(r_seq, window=10, salt_kw=salt_kw)
        het_dg = _heterodimer(f_seq, r_seq, salt_kw=salt_kw)

        issues: list[str] = []
        warnings: list[str] = []
        if f_gc < params.gc_min or f_gc > params.gc_max:
            issues.append(f"正向引物 GC% {f_gc:.1f} 超出 [{params.gc_min}-{params.gc_max}]")
        if r_gc < params.gc_min or r_gc > params.gc_max:
            issues.append(f"反向引物 GC% {r_gc:.1f} 超出 [{params.gc_min}-{params.gc_max}]")
        if f_3hp < params.max_3prime_hairpin_dg:
            issues.append(f"正向 3' 端发夹 ΔG={f_3hp:.2f} 过稳")
        if r_3hp < params.max_3prime_hairpin_dg:
            issues.append(f"反向 3' 端发夹 ΔG={r_3hp:.2f} 过稳")
        if het_dg < params.max_heterodimer_dg:
            issues.append(f"正反向异源二聚体 ΔG={het_dg:.2f} 过稳")
        # v2: 长串单核苷酸检查 (Daher 2016)
        f_homo = _max_homopolymer_run(f_seq)
        r_homo = _max_homopolymer_run(r_seq)
        if f_homo > params.max_homopolymer:
            warnings.append(f"正向引物最长连续 {f_homo}nt > {params.max_homopolymer} (Daher 2016 警告)")
        if r_homo > params.max_homopolymer:
            warnings.append(f"反向引物最长连续 {r_homo}nt > {params.max_homopolymer} (Daher 2016 警告)")
        # v2: 小重复检查
        if _has_long_repeats(f_seq):
            warnings.append("正向引物含大量小重复序列 (Daher 2016 警告)")
        if _has_long_repeats(r_seq):
            warnings.append("反向引物含大量小重复序列 (Daher 2016 警告)")

        probe: Optional[RPAProbe] = None
        if params.design_probe:
            probe = _design_exo_probe(
                sequence, f_seq, f_pos, r_seq, r_pos,
                params.probe_min_size, params.probe_max_size, salt_kw,
            )

        pairs.append(RPAPrimerPair(
            pair_id=i,
            forward_seq=f_seq,
            reverse_seq=r_seq,
            forward_pos=f_pos,
            reverse_pos=r_pos,
            product_size=product_size,
            forward_length=len(f_seq),
            reverse_length=len(r_seq),
            forward_gc=f_gc,
            reverse_gc=r_gc,
            forward_tm=f_tm,
            reverse_tm=r_tm,
            forward_3prime_hairpin_dg=f_3hp,
            reverse_3prime_hairpin_dg=r_3hp,
            heterodimer_dg=het_dg,
            accepted=not issues,  # warnings 不算 fatal
            issues=issues + warnings,
            probe=probe,
        ))

    # accepted 优先,然后按产品大小从中间向两端排序
    pairs.sort(key=lambda p: (not p.accepted, abs(p.product_size - 200)))
    return pairs
