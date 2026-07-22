"""引物探针设计模块:基于 primer3-py 批量生成候选引物对。

输出每组候选的引物序列、Tm、GC%、探针序列、产物大小、penalty 等参数。
"""
from __future__ import annotations

import json
from dataclasses import dataclass, asdict, field
from pathlib import Path
from typing import Optional

import primer3


@dataclass
class PrimerPairResult:
    """单组引物探针设计结果。

    Attributes:
        pair_id: 组号(0-based)。
        forward_seq: 正向引物序列。
        reverse_seq: 反向引物序列。
        probe_seq: 探针序列(TaqMan/Molecular Beacon)。
        forward_tm: 正向引物 Tm(°C)。
        reverse_tm: 反向引物 Tm(°C)。
        probe_tm: 探针 Tm(°C)。
        forward_gc: 正向引物 GC%。
        reverse_gc: 反向引物 GC%。
        probe_gc: 探针 GC%。
        product_size: 产物大小(bp)。
        penalty: primer3 评估的惩罚值(越低越好)。
        forward_explain: 正向引物各项 penalty 解释。
        reverse_explain: 反向引物各项 penalty 解释。
        pair_penalty: 引物对 penalty。
    """

    pair_id: int
    forward_seq: str = ""
    reverse_seq: str = ""
    probe_seq: str = ""
    forward_tm: float = 0.0
    reverse_tm: float = 0.0
    probe_tm: float = 0.0
    forward_gc: float = 0.0
    reverse_gc: float = 0.0
    probe_gc: float = 0.0
    product_size: int = 0
    penalty: float = 0.0
    forward_explain: str = ""
    reverse_explain: str = ""
    pair_penalty: float = 0.0

    def to_dict(self) -> dict:
        """转换为字典,便于表格化展示。"""
        return asdict(self)


@dataclass
class DesignParameters:
    """引物设计参数。

    Attributes:
        product_size_min: 产物最小长度。
        product_size_max: 产物最大长度。
        primer_tm_min: 引物最小 Tm。
        primer_tm_opt: 引物最佳 Tm。
        primer_tm_max: 引物最大 Tm。
        probe_tm_min: 探针最小 Tm。
        probe_tm_opt: 探针最佳 Tm。
        probe_tm_max: 探针最大 Tm。
        primer_gc_min: 引物 GC% 下限。
        primer_gc_max: 引物 GC% 上限。
        primer_min_size: 引物最小长度。
        primer_opt_size: 引物最佳长度。
        primer_max_size: 引物最大长度。
        num_return: 返回候选组数。
        detection_type: 检测类型(TaqMan / SYBR / Beacon)。
        target_start: 靶标起始(可选,1-based)。
        target_length: 靶标长度(可选)。
    """

    product_size_min: int = 70
    product_size_max: int = 150
    primer_tm_min: float = 58.0
    primer_tm_opt: float = 60.0
    primer_tm_max: float = 62.0
    probe_tm_min: float = 68.0
    probe_tm_opt: float = 70.0
    probe_tm_max: float = 72.0
    primer_gc_min: float = 40.0
    primer_gc_max: float = 60.0
    primer_min_size: int = 18
    primer_opt_size: int = 20
    primer_max_size: int = 25
    num_return: int = 5
    detection_type: str = "TaqMan"
    target_start: Optional[int] = None
    target_length: Optional[int] = None

    def to_primer3_dict(self, sequence: str) -> dict:
        """转换为 primer3.bindings.design_primers 所需的字典。

        Args:
            sequence: 模板序列。

        Returns:
            primer3 参数字典。
        """
        seq_args = {
            "SEQUENCE_ID": "qpcr_target",
            "SEQUENCE_TEMPLATE": sequence,
        }
        # 可选:指定靶标区域
        if self.target_start is not None and self.target_length is not None:
            seq_args["SEQUENCE_TARGET"] = [self.target_start, self.target_length]

        global_args = {
            "PRIMER_TASK": "generic",
            "PRIMER_PICK_LEFT_PRIMER": 1,
            "PRIMER_PICK_INTERNAL_OLIGO": 1 if self.detection_type in {"TaqMan", "Beacon"} else 0,
            "PRIMER_PICK_RIGHT_PRIMER": 1,
            "PRIMER_OPT_SIZE": self.primer_opt_size,
            "PRIMER_MIN_SIZE": self.primer_min_size,
            "PRIMER_MAX_SIZE": self.primer_max_size,
            "PRIMER_OPT_TM": self.primer_tm_opt,
            "PRIMER_MIN_TM": self.primer_tm_min,
            "PRIMER_MAX_TM": self.primer_tm_max,
            "PRIMER_MIN_GC": self.primer_gc_min,
            "PRIMER_MAX_GC": self.primer_gc_max,
            "PRIMER_PRODUCT_SIZE_RANGE": [[self.product_size_min, self.product_size_max]],
            "PRIMER_NUM_RETURN": self.num_return,
            "PRIMER_MAX_POLY_X": 4,
            "PRIMER_MAX_NS_ACCEPTED": 0,
            "PRIMER_MIN_THREE_PRIME_DISTANCE": 3,
            "PRIMER_MAX_SELF_ANY": 8,
            "PRIMER_MAX_SELF_END": 3,
            "PRIMER_PAIR_MAX_COMPL_ANY": 8,
            "PRIMER_PAIR_MAX_COMPL_END": 3,
            "PRIMER_MAX_HAIRPIN_TH": 47.0,
        }
        if self.detection_type in {"TaqMan", "Beacon"}:
            global_args.update({
                "PRIMER_INTERNAL_OPT_SIZE": 22,
                "PRIMER_INTERNAL_MIN_SIZE": 18,
                "PRIMER_INTERNAL_MAX_SIZE": 30,
                "PRIMER_INTERNAL_OPT_TM": self.probe_tm_opt,
                "PRIMER_INTERNAL_MIN_TM": self.probe_tm_min,
                "PRIMER_INTERNAL_MAX_TM": self.probe_tm_max,
                "PRIMER_INTERNAL_MIN_GC": 40.0,
                "PRIMER_INTERNAL_MAX_GC": 60.0,
                "PRIMER_INTERNAL_MAX_POLY_X": 4,
                "PRIMER_INTERNAL_MAX_NS_ACCEPTED": 0,
                "PRIMER_INTERNAL_MAX_SELF_ANY": 8,
                "PRIMER_INTERNAL_MAX_SELF_END": 3,
                "PRIMER_INTERNAL_MAX_HAIRPIN_TH": 47.0,
            })

        return {"SEQUENCE": seq_args, "GLOBAL": global_args}


def _parse_pair(
    pair_id: int,
    primer_dict: dict,
    include_probe: bool,
) -> PrimerPairResult:
    """将 primer3 的单组返回解析为 PrimerPairResult。

    Args:
        pair_id: 候选组号。
        primer_dict: primer3 返回的 dict。
        include_probe: 是否包含探针。

    Returns:
        PrimerPairResult 对象。
    """
    f_seq = primer_dict.get(f"PRIMER_LEFT_{pair_id}_SEQUENCE", "")
    r_seq = primer_dict.get(f"PRIMER_RIGHT_{pair_id}_SEQUENCE", "")
    p_seq = primer_dict.get(f"PRIMER_INTERNAL_{pair_id}_SEQUENCE", "") if include_probe else ""

    f_tm = float(primer_dict.get(f"PRIMER_LEFT_{pair_id}_TM", 0.0))
    r_tm = float(primer_dict.get(f"PRIMER_RIGHT_{pair_id}_TM", 0.0))
    p_tm = float(primer_dict.get(f"PRIMER_INTERNAL_{pair_id}_TM", 0.0)) if include_probe else 0.0

    f_gc = float(primer_dict.get(f"PRIMER_LEFT_{pair_id}_GC_PERCENT", 0.0))
    r_gc = float(primer_dict.get(f"PRIMER_RIGHT_{pair_id}_GC_PERCENT", 0.0))
    p_gc = float(primer_dict.get(f"PRIMER_INTERNAL_{pair_id}_GC_PERCENT", 0.0)) if include_probe else 0.0

    product_size = int(primer_dict.get(f"PRIMER_PAIR_{pair_id}_PRODUCT_SIZE", 0))
    pair_penalty = float(primer_dict.get(f"PRIMER_PAIR_{pair_id}_PENALTY", 0.0))
    f_penalty = float(primer_dict.get(f"PRIMER_LEFT_{pair_id}_PENALTY", 0.0))
    r_penalty = float(primer_dict.get(f"PRIMER_RIGHT_{pair_id}_PENALTY", 0.0))

    return PrimerPairResult(
        pair_id=pair_id,
        forward_seq=f_seq,
        reverse_seq=r_seq,
        probe_seq=p_seq,
        forward_tm=f_tm,
        reverse_tm=r_tm,
        probe_tm=p_tm,
        forward_gc=f_gc,
        reverse_gc=r_gc,
        probe_gc=p_gc,
        product_size=product_size,
        penalty=pair_penalty + f_penalty + r_penalty,
        forward_explain=str(primer_dict.get(f"PRIMER_LEFT_{pair_id}_EXPLAIN", "")),
        reverse_explain=str(primer_dict.get(f"PRIMER_RIGHT_{pair_id}_EXPLAIN", "")),
        pair_penalty=pair_penalty,
    )


def design_primer_pairs(
    sequence: str,
    params: Optional[DesignParameters] = None,
) -> list[PrimerPairResult]:
    """调用 primer3.design_primers 批量设计引物探针。

    Args:
        sequence: 模板序列(已转大写,无空白)。
        params: 设计参数;若为 None 则使用默认值。

    Returns:
        PrimerPairResult 列表,按 penalty 升序排列。

    Raises:
        ValueError: 当序列过短或参数不合理时。
        RuntimeError: 当 primer3 内部异常时。
    """
    if len(sequence) < 100:
        raise ValueError("模板序列至少需要 100 bp 才能设计引物")

    params = params or DesignParameters()
    include_probe = params.detection_type in {"TaqMan", "Beacon"}

    args = params.to_primer3_dict(sequence)
    try:
        result = primer3.bindings.design_primers(args["SEQUENCE"], args["GLOBAL"])
    except Exception as exc:
        raise RuntimeError(f"primer3 设计失败: {exc}") from exc

    count = int(result.get("PRIMER_PAIR_NUM_RETURNED", 0))
    if count == 0:
        raise RuntimeError("未找到满足条件的引物对,请放宽参数或检查序列")

    pairs: list[PrimerPairResult] = []
    for i in range(count):
        pairs.append(_parse_pair(i, result, include_probe))

    pairs.sort(key=lambda x: x.penalty)
    return pairs


def design_single_pair(
    sequence: str,
    forward_primer: str,
    reverse_primer: str,
    probe: str = "",
    detection_type: str = "TaqMan",
) -> PrimerPairResult:
    """对已知的引物对进行基础评估(不调用 design_primers)。

    Args:
        sequence: 模板序列。
        forward_primer: 正向引物。
        reverse_primer: 反向引物。
        probe: 探针(可空)。
        detection_type: 检测类型。

    Returns:
        PrimerPairResult 对象(Tm/GC 等使用 calc_tm/calc_gc 计算)。
    """
    f_tm = primer3.calc_tm(forward_primer)
    r_tm = primer3.calc_tm(reverse_primer)
    f_gc = primer3.calc_gc(forward_primer)
    r_gc = primer3.calc_gc(reverse_primer)

    p_tm = 0.0
    p_gc = 0.0
    if probe:
        p_tm = primer3.calc_tm(probe)
        p_gc = primer3.calc_gc(probe)

    # 简单估算产物大小:在序列中查找引物位置
    f_pos = sequence.find(forward_primer)
    r_comp = str(__import__("Bio.Seq", fromlist=["Seq"]).Seq(reverse_primer).reverse_complement())
    r_pos = sequence.find(r_comp)
    if f_pos >= 0 and r_pos >= 0 and r_pos > f_pos:
        product_size = r_pos + len(r_comp) - f_pos
    else:
        product_size = 0

    return PrimerPairResult(
        pair_id=0,
        forward_seq=forward_primer,
        reverse_seq=reverse_primer,
        probe_seq=probe,
        forward_tm=float(f_tm),
        reverse_tm=float(r_tm),
        probe_tm=float(p_tm),
        forward_gc=float(f_gc),
        reverse_gc=float(r_gc),
        probe_gc=float(p_gc),
        product_size=product_size,
        penalty=abs(f_tm - 60.0) + abs(r_tm - 60.0) + (abs(p_tm - 70.0) if probe else 0),
    )


def load_master_mix_presets(path: Optional[str] = None) -> dict:
    """读取试剂盒预设 JSON。

    Args:
        path: JSON 文件路径;None 时使用默认 data/master_mix_presets.json。

    Returns:
        预设字典。
    """
    if path is None:
        path = str(Path(__file__).resolve().parent.parent / "data" / "master_mix_presets.json")
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)