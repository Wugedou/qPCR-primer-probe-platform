"""序列输入模块：处理手动输入和 GenBank 检索。

支持 FASTA 格式与纯序列文本,并可通过 Biopython Entrez 从 NCBI 下载序列。
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Optional

from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# 常用 IUPAC 碱基字符集合(用于校验)
_VALID_BASES = set("ACGTUNRYSWKMBDHVacgtunryswkmbdhv")


@dataclass
class SequenceInfo:
    """序列信息数据类。

    Attributes:
        sequence: 纯序列字符串(已转大写)。
        description: 序列描述(FASTA 头或自定义)。
        organism: 病原体名称(若可知)。
        gene: 基因名称(若可知)。
        accession: GenBank 检索号(若从 NCBI 下载)。
        length: 序列长度。
    """

    sequence: str
    description: str = ""
    organism: str = ""
    gene: str = ""
    accession: str = ""
    length: int = field(init=False)

    def __post_init__(self) -> None:
        self.sequence = self.sequence.upper().strip()
        self.length = len(self.sequence)


def _clean_sequence(raw: str) -> str:
    """清洗并校验用户输入的序列字符串。

    支持 FASTA(以 > 开头的行)和纯序列两种格式。
    去除空白、数字、空行等非碱基字符后,只保留合法碱基。

    Args:
        raw: 原始输入文本。

    Returns:
        清洗后的大写序列字符串。

    Raises:
        ValueError: 当序列为空或全部字符都不是合法碱基时。
    """
    if not raw:
        raise ValueError("序列不能为空")

    lines = raw.strip().splitlines()
    # FASTA 格式检测:存在 > 开头的行
    if any(line.startswith(">") for line in lines):
        seq_parts: list[str] = []
        for line in lines:
            if line.startswith(">"):
                continue
            seq_parts.append(line)
        cleaned = "".join(seq_parts)
    else:
        # 纯序列模式:去除所有空白和数字(常见 FASTA 行号)
        cleaned = re.sub(r"[\s0-9]+", "", raw)

    # 去除非 IUPAC 字符
    cleaned = "".join(ch for ch in cleaned if ch in _VALID_BASES)

    if not cleaned:
        raise ValueError("清洗后序列为空,请检查输入是否包含合法碱基(A/C/G/T/U/N)")

    if len(cleaned) < 50:
        raise ValueError(f"序列过短({len(cleaned)} bp),引物设计至少需要 50 bp")

    return cleaned


def fetch_sequence_manual(text: str, description: str = "") -> SequenceInfo:
    """从用户粘贴的文本中解析序列。

    Args:
        text: FASTA 文本或纯序列文本。
        description: 序列描述(可由用户指定)。

    Returns:
        SequenceInfo 对象。

    Raises:
        ValueError: 当输入不合法时。
    """
    cleaned = _clean_sequence(text)
    desc = description.strip() if description else "手动输入序列"
    return SequenceInfo(sequence=cleaned, description=desc)


def fetch_sequence_genbank(
    accession: str,
    email: str = "qpcr_tool@example.com",
    organism_hint: str = "",
    gene_hint: str = "",
) -> SequenceInfo:
    """通过 NCBI Entrez 下载指定 accession 的 GenBank 记录并提取序列。

    Args:
        accession: GenBank 检索号(如 NM_001234.1)。
        email: 联系方式邮箱(NCBI 强制要求)。
        organism_hint: 病原体名称提示。
        gene_hint: 基因名称提示。

    Returns:
        SequenceInfo 对象。

    Raises:
        RuntimeError: 当检索或解析失败时。
    """
    if not accession:
        raise ValueError("GenBank accession 不能为空")

    Entrez.email = email
    try:
        # 1) 检索以确定记录存在
        with Entrez.esearch(db="nucleotide", term=accession, retmax=1) as h:
            search = Entrez.read(h)
        if not search.get("IdList"):
            raise RuntimeError(f"未找到 accession: {accession}")

        # 2) 获取 GenBank 格式详细记录
        with Entrez.efetch(
            db="nucleotide",
            id=accession,
            rettype="gb",
            retmode="text",
        ) as h:
            record = SeqIO.read(h, "genbank")

    except Exception as exc:
        raise RuntimeError(f"GenBank 检索失败: {exc}") from exc

    seq_str = str(record.seq).upper()
    if len(seq_str) < 50:
        raise RuntimeError(f"下载到的序列过短({len(seq_str)} bp),请检查 accession")

    # 尝试从注释中解析 organism / gene
    organism = organism_hint or record.annotations.get("organism", "")
    gene = gene_hint
    for feat in record.features:
        if feat.type in {"gene", "CDS", "misc_feature"}:
            if "gene" in feat.qualifiers:
                gene = gene or feat.qualifiers["gene"][0]
                break

    description = record.description or f"GenBank:{accession}"

    return SequenceInfo(
        sequence=seq_str,
        description=description,
        organism=organism,
        gene=gene,
        accession=accession,
    )


def subsequence_by_target(
    full_seq: str, start: int, length: int
) -> str:
    """按位置截取子序列,用于指定靶标区域设计引物。

    Args:
        full_seq: 完整序列。
        start: 起始位置(0-based)。
        length: 子序列长度。

    Returns:
        子序列字符串。

    Raises:
        ValueError: 当位置越界时。
    """
    if start < 0 or start >= len(full_seq):
        raise ValueError(f"起始位置 {start} 超出范围(0-{len(full_seq)-1})")
    if length <= 0:
        raise ValueError("长度必须为正数")
    end = min(start + length, len(full_seq))
    return full_seq[start:end]