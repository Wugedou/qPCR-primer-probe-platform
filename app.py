"""qPCR 引物探针设计 + 仿真验证平台 - Streamlit 主应用。

布局:
- 侧边栏:序列输入(手动 / GenBank)、病原信息、检测类型、参数设置
- 主区:Tabs 1-5 分别对应设计、热力学仿真、扩增仿真、优化循环、报告导出
"""
from __future__ import annotations

import io
import json
import sys
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

# 确保 core 包可导入(无论 streamlit run 的工作目录如何)
_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

from core.amplification import (
    AmplificationParams,
    LAMPIsothermalParams,
    RPAIsothermalParams,
    compute_ct,
    estimate_lod,
    estimate_lamp_tt,
    estimate_rpa_tt,
    fit_standard_curve,
    simulate_amplification,
    simulate_lamp,
    simulate_rpa,
)
from core.lamp_design import (
    LAMPDesignParameters,
    LAMPPrimerSet,
    design_lamp_primers,
)
from core.optimizer import OptimizationEngine, OptimizationResult, OptimizationSuggestion
from core.primer_design import (
    DesignParameters,
    PrimerPairResult,
    design_primer_pairs,
    load_master_mix_presets,
)
from core.report import generate_report, report_from_optimization
from core.rpa_design import (
    RPADesignParameters,
    RPAPrimerPair,
    design_rpa_primers,
)
from core.sequence_input import (
    SequenceInfo,
    fetch_sequence_genbank,
    fetch_sequence_manual,
)
from core.thermo_sim import (
    ReactionSystem,
    ThermoResult,
    compare_master_mixes,
    evaluate_thermodynamics,
    load_preset_systems,
    recommend_best_system,
)


# ============================================================
# 页面配置
# ============================================================
st.set_page_config(
    page_title="qPCR 引物探针设计仿真平台",
    page_icon="🧬",
    layout="wide",
)


# ============================================================
# 工具函数:session_state 操作
# ============================================================
def _init_state() -> None:
    """初始化 session_state 变量。"""
    defaults = {
        "seq_info": None,
        "pairs": None,                # qPCR 候选(PrimerPairResult 列表)
        "rpa_pairs": None,            # RPA 候选(RPAPrimerPair 列表)
        "lamp_sets": None,            # LAMP 候选(LAMPPrimerSet 列表)
        "design_mode": "qpcr",        # 当前设计模式(qpcr / rpa / lamp)
        "thermo_results": None,
        "selected_pair_idx": 0,
        "selected_system_key": "ABI_TaqMan_Fast",
        "amp_curves": None,
        "std_curve": None,
        "lod_info": None,
        "opt_result": None,
        "design_params": None,
        "amp_params": None,
        "copy_numbers": [1e8, 1e7, 1e6, 1e5, 1e4, 1e3, 1e2, 1e1],
        "report_text": None,
        "preset_systems": None,
        "rpa_curves": None,           # RPA 时间-荧光曲线
        "rpa_tt": None,               # RPA Tt 估计
        "lamp_curves": None,          # LAMP 时间-荧光曲线
        "lamp_tt": None,              # LAMP Tt 估计
    }
    for k, v in defaults.items():
        if k not in st.session_state:
            st.session_state[k] = v


_init_state()


# ============================================================
# 设计引擎路由
# ============================================================
def _run_design(
    detection_type: str,
    sequence: str,
    design_params: DesignParameters,
) -> dict:
    """根据 detection_type 调用对应设计引擎。

    Returns:
        字典 {"mode": "qpcr"/"rpa"/"lamp", "qpcr_pairs": ..., "rpa_pairs": ..., "lamp_sets": ...}
    """
    result: dict = {"mode": detection_type.lower()}
    if detection_type == "RPA":
        rpa_params = RPADesignParameters(
            product_size_min=max(100, design_params.product_size_min),
            product_size_max=max(500, design_params.product_size_max),
            num_return=design_params.num_return,
            design_probe=True,
        )
        result["rpa_pairs"] = design_rpa_primers(sequence, rpa_params)
        return result
    if detection_type == "LAMP":
        lamp_params = LAMPDesignParameters(num_sets=3, include_loop=True)
        result["lamp_sets"] = design_lamp_primers(sequence, lamp_params)
        return result
    # qPCR 三种检测类型(TaqMan/SYBR Green/Molecular Beacon)走 primer3
    result["qpcr_pairs"] = design_primer_pairs(sequence, design_params)
    return result


# ============================================================
# 侧边栏
# ============================================================
def render_sidebar() -> tuple[SequenceInfo, DesignParameters, AmplificationParams]:
    """渲染侧边栏输入面板,返回当前 SequenceInfo / 设计 / 扩增参数。"""
    st.sidebar.title("🧪 输入与参数")

    st.sidebar.subheader("1. 序列来源")
    seq_source = st.sidebar.radio(
        "序列输入方式",
        options=["手动输入", "GenBank 检索"],
        index=0,
    )

    seq_info: Optional[SequenceInfo] = st.session_state.seq_info
    if seq_source == "手动输入":
        fasta_text = st.sidebar.text_area(
            "粘贴 FASTA 或纯序列",
            value=(
                st.session_state.seq_info.sequence
                if st.session_state.seq_info is not None
                else ">SARS-CoV-2_N_gene\n"
                "ATGTCTGATAATGGACCCCAAAATCAGCGAAATGCACCCCGCATTACGTTTGGTGGACCCTCAGATTCAACTGGCAGTAACCAGAATGGAGAACGCAGTGGGGCGCGATCAAAACAACGTCGGCCCCAAGGTTTACCCAATAATACTGCGTCTTGGTTCACCGCTCTCACTCAACATGGCAAGGAAGACCTTAAATTCCCTCGAGGACAAGGCGTTCCAATTAACACCAATAGCAGTCCAGATGACCAAATTGGCTACTACCGAAGAGCTACCAGACGAATTCGTGGTGGTGACGGTAAAATGAAAGATCTCAGTCCAAGATGGTATTTCTACTACCTAGGAACTGGGCCAGAAGCTGGACTTCCCTATGGTGCTAACAAAGACGGCATCATATGGGTTGCAACTGAGGGAGCCTTGAATACACCAAAAGATCACATTGGCACCCGCAATCCTGCTAACAATGCTGCAATCGTGCTACAACTTCCTCAAGGAACAACATTGCCAAAAGGCTTCTACGCAGAAGGGAGCAGAGGCGGCAGTCAAGCCTCTTCTCGTTCCTCATCACGTAGTCGCAACAGTTCAAGAAATTCAACTCCAGGCAGCAGTAGGGGAACTTCTCCTGCTAGAATGGCTGGCAATGGCGGTGATGCTGCTCTTGCTTTGCTGCTGCTTGACAGATTGAACCAGCTTGAGAGCAAAATGTCTGGTAAAGGCCAACAACAACAAGGCCAAACTGTCACTAAGAAATCTGCTGCTGAGGCTTCTAAGAAGCCTCGGCAAAAACGTACTGCCACTAAAGCATACAATGTAACACAAGCTTTCGGCAGACGTGGTCCAGAACAAACCCAAGGAAATTTTGGGGACCAGGAACTAATCAGACAAGGAACTGATTACAAACATTGGCCGCAAATTGCACAATTTGCCCCCAGCGCTTCAGCGTTCTTCGGAATGTCGCGCATTGGCATGGAAGTCACACCTTCGGGAACGTGGTTGACCTACACAGGTGCCATCAAATTGGATGACAAAGATCCAAATTTCAAAGATCAAGTCATTTTGCTGAATAAGCATATTGACGCATACAAAACATTCCCACCAACAGAGCCTAAAAAGGACAAAAAGAAGAAGGCTGATGAAACTCAAGCCTTACCGCAGAGACAGAAGAAACAGCAAACTGTGACTCTTCTTCCTGCTGCAGATTTGGATGATTTCTCCAAACAATTGCAACAATCCATGAGCAGTGCTGACTCAACTCAGGCCTAA"
            ),
            height=180,
        )
        if st.sidebar.button("解析序列", type="primary"):
            try:
                seq_info = fetch_sequence_manual(fasta_text)
                st.session_state.seq_info = seq_info
                st.session_state.pairs = None
                st.session_state.thermo_results = None
                st.sidebar.success(f"✓ 解析成功:长度 {seq_info.length} bp")
            except ValueError as exc:
                st.sidebar.error(f"❌ {exc}")
    else:
        accession = st.sidebar.text_input("GenBank accession(如 NM_001234.1)", value="")
        email = st.sidebar.text_input("NCBI 邮箱", value="qpcr_tool@example.com")
        if st.sidebar.button("下载序列", type="primary"):
            if not accession:
                st.sidebar.error("请输入 accession")
            else:
                with st.sidebar.spinner("正在从 NCBI 下载..."):
                    try:
                        seq_info = fetch_sequence_genbank(accession, email=email)
                        st.session_state.seq_info = seq_info
                        st.session_state.pairs = None
                        st.session_state.thermo_results = None
                        st.sidebar.success(
                            f"✓ 下载成功:{seq_info.description[:60]} ({seq_info.length} bp)"
                        )
                    except Exception as exc:
                        st.sidebar.error(f"❌ {exc}")

    if seq_info:
        st.sidebar.markdown("**当前序列**")
        st.sidebar.info(
            f"描述:{seq_info.description[:50]}\n\n"
            f"长度:{seq_info.length} bp\n\n"
            f"GC%:{_calc_gc(seq_info.sequence):.1f}%"
        )

    st.sidebar.divider()
    st.sidebar.subheader("2. 病原信息")
    pathogen_name = st.sidebar.text_input(
        "病原体名称",
        value=st.session_state.get("pathogen_name", ""),
    )
    gene_name = st.sidebar.text_input(
        "靶基因名称",
        value=st.session_state.get("gene_name", ""),
    )
    target_region = st.sidebar.text_input(
        "靶标区域描述",
        value=st.session_state.get("target_region", ""),
    )
    st.session_state.pathogen_name = pathogen_name
    st.session_state.gene_name = gene_name
    st.session_state.target_region = target_region

    st.sidebar.divider()
    st.sidebar.subheader("3. 检测类型与产物")
    detection_type = st.sidebar.selectbox(
        "检测类型",
        options=["TaqMan", "SYBR Green", "Molecular Beacon", "RPA", "LAMP"],
        index=0,
        help="TaqMan/SYBR/Beacon 走 qPCR 循环扩增;RPA 走重组酶等温扩增;LAMP 走环介导等温扩增",
    )
    product_size = st.sidebar.slider(
        "产物大小 (bp)",
        min_value=60, max_value=300,
        value=(70, 150),
        step=5,
    )

    st.sidebar.subheader("4. 引物/探针 Tm 设置")
    primer_tm = st.sidebar.slider(
        "引物 Tm 范围 (°C)",
        min_value=55.0, max_value=65.0,
        value=(58.0, 62.0),
        step=0.5,
    )
    probe_tm = st.sidebar.slider(
        "探针 Tm 范围 (°C)",
        min_value=63.0, max_value=75.0,
        value=(68.0, 72.0),
        step=0.5,
    )

    st.sidebar.subheader("5. 扩增仿真参数")
    cycles = st.sidebar.slider("循环数", min_value=20, max_value=50, value=40, step=1)
    efficiency = st.sidebar.slider(
        "扩增效率 (%)",
        min_value=50, max_value=100, value=95, step=1,
    ) / 100.0
    threshold = st.sidebar.number_input(
        "荧光阈值 (ΔRn)",
        min_value=0.001, max_value=1.0, value=0.1, step=0.01, format="%.3f",
    )
    noise = st.sidebar.slider(
        "噪声水平",
        min_value=0.0, max_value=0.05, value=0.005, step=0.001, format="%.3f",
    )
    num_copies_grad = st.sidebar.multiselect(
        "拷贝数梯度(标准曲线)",
        options=[1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9],
        default=[1e8, 1e7, 1e6, 1e5, 1e4, 1e3, 1e2, 1e1],
    )

    design_params = DesignParameters(
        product_size_min=product_size[0],
        product_size_max=product_size[1],
        primer_tm_min=primer_tm[0],
        primer_tm_opt=(primer_tm[0] + primer_tm[1]) / 2,
        primer_tm_max=primer_tm[1],
        probe_tm_min=probe_tm[0],
        probe_tm_opt=(probe_tm[0] + probe_tm[1]) / 2,
        probe_tm_max=probe_tm[1],
        num_return=5,
        detection_type=detection_type,
    )

    amp_params = AmplificationParams(
        cycles=cycles,
        efficiency=efficiency,
        threshold=threshold,
        noise_std=noise,
    )
    st.session_state.copy_numbers = num_copies_grad or [1e8, 1e6, 1e4, 1e2]
    st.session_state.design_params = design_params
    st.session_state.amp_params = amp_params

    return seq_info, design_params, amp_params


def _calc_gc(seq: str) -> float:
    """计算序列 GC%。"""
    if not seq:
        return 0.0
    gc = sum(1 for b in seq.upper() if b in {"G", "C"})
    return gc / len(seq) * 100.0


# ============================================================
# Tab 1: 引物探针设计
# ============================================================
def tab_primer_design(seq_info: SequenceInfo, design_params: DesignParameters) -> None:
    st.header("🧬 引物探针设计")
    if not seq_info:
        st.warning("请先在左侧输入或下载序列")
        return

    detection = design_params.detection_type
    col_a, col_b = st.columns([3, 1])
    with col_a:
        st.markdown(
            f"**模板长度:** {seq_info.length} bp  |  "
            f"**检测类型:** {detection}  |  "
            f"**产物范围:** {design_params.product_size_min}-{design_params.product_size_max} bp"
        )
    with col_b:
        if st.button("🚀 开始设计引物探针", type="primary", width="stretch"):
            with st.spinner("设计计算中..."):
                try:
                    out = _run_design(detection, seq_info.sequence, design_params)
                    st.session_state.design_mode = out["mode"]
                    if "qpcr_pairs" in out:
                        st.session_state.pairs = out["qpcr_pairs"]
                    if "rpa_pairs" in out:
                        st.session_state.rpa_pairs = out["rpa_pairs"]
                    if "lamp_sets" in out:
                        st.session_state.lamp_sets = out["lamp_sets"]
                    st.session_state.thermo_results = None
                    st.session_state.amp_curves = None
                    st.session_state.rpa_curves = None
                    st.session_state.lamp_curves = None
                    st.session_state.opt_result = None
                    count = (
                        len(out.get("qpcr_pairs", [])) or
                        len(out.get("rpa_pairs", [])) or
                        len(out.get("lamp_sets", []))
                    )
                    st.success(f"✓ 设计完成,共 {count} 组候选")
                except (RuntimeError, ValueError) as exc:
                    st.error(f"❌ 设计失败:{exc}")

    if detection == "RPA":
        _render_rpa_results()
    elif detection == "LAMP":
        _render_lamp_results()
    else:
        _render_qpcr_results()


def _render_qpcr_results() -> None:
    """Tab 1 qPCR 分支(TaqMan / SYBR / Beacon)。"""
    pairs: list[PrimerPairResult] = st.session_state.pairs or []
    if not pairs:
        st.info("点击「开始设计引物探针」生成候选。")
        return

    rows = []
    for p in pairs:
        rows.append({
            "组号": p.pair_id,
            "正向引物": p.forward_seq,
            "反向引物": p.reverse_seq,
            "探针": p.probe_seq or "-",
            "F_Tm(°C)": round(p.forward_tm, 2),
            "R_Tm(°C)": round(p.reverse_tm, 2),
            "P_Tm(°C)": round(p.probe_tm, 2) if p.probe_seq else None,
            "F_GC%": round(p.forward_gc, 1),
            "R_GC%": round(p.reverse_gc, 1),
            "P_GC%": round(p.probe_gc, 1) if p.probe_seq else None,
            "产物(bp)": p.product_size,
            "penalty": round(p.penalty, 4),
        })
    df = pd.DataFrame(rows).sort_values("penalty")
    st.subheader("候选引物探针列表(按 penalty 升序)")
    st.dataframe(df, width="stretch", hide_index=True)

    pair_options = [f"#{p.pair_id} | penalty={p.penalty:.4f} | 产物={p.product_size}bp" for p in pairs]
    current_idx = st.session_state.selected_pair_idx
    if current_idx >= len(pairs):
        current_idx = 0
    sel = st.selectbox(
        "选择要深入分析的候选",
        options=list(range(len(pairs))),
        format_func=lambda i: pair_options[i],
        index=current_idx,
    )
    st.session_state.selected_pair_idx = sel
    selected_pair = pairs[sel]

    with st.expander("📋 选中候选详情", expanded=True):
        c1, c2, c3 = st.columns(3)
        with c1:
            st.markdown("**正向引物**")
            st.code(selected_pair.forward_seq, language="text")
            st.write(
                f"Tm={selected_pair.forward_tm:.2f}°C, "
                f"GC={selected_pair.forward_gc:.1f}%, "
                f"长度={len(selected_pair.forward_seq)} nt"
            )
        with c2:
            st.markdown("**反向引物**")
            st.code(selected_pair.reverse_seq, language="text")
            st.write(
                f"Tm={selected_pair.reverse_tm:.2f}°C, "
                f"GC={selected_pair.reverse_gc:.1f}%, "
                f"长度={len(selected_pair.reverse_seq)} nt"
            )
        with c3:
            if selected_pair.probe_seq:
                st.markdown("**探针**")
                st.code(selected_pair.probe_seq, language="text")
                st.write(
                    f"Tm={selected_pair.probe_tm:.2f}°C, "
                    f"GC={selected_pair.probe_gc:.1f}%, "
                    f"长度={len(selected_pair.probe_seq)} nt"
                )
            else:
                st.markdown("**探针(SYBR 模式无)**")
                st.info("当前为 SYBR Green 检测,无探针序列")


def _render_rpa_results() -> None:
    """Tab 1 RPA 分支。"""
    pairs: list[RPAPrimerPair] = st.session_state.rpa_pairs or []
    if not pairs:
        st.info("点击「开始设计引物探针」生成 RPA 引物。")
        return

    rows = []
    for p in pairs:
        rows.append({
            "组号": p.pair_id,
            "正向引物": p.forward_seq,
            "反向引物": p.reverse_seq,
            "F_len": p.forward_length,
            "R_len": p.reverse_length,
            "F_GC%": round(p.forward_gc, 1),
            "R_GC%": round(p.reverse_gc, 1),
            "F_3'末端发夹dG": round(p.forward_3prime_hairpin_dg, 2),
            "R_3'末端发夹dG": round(p.reverse_3prime_hairpin_dg, 2),
            "异二聚体dG": round(p.heterodimer_dg, 2),
            "产物(bp)": p.product_size,
            "可用": "✓" if p.accepted else "⚠️",
        })
    st.subheader("RPA 引物对候选")
    st.caption("RPA: 等温 37-42°C,引物 30-35 bp,产物 100-500 bp,3' 末端避免发夹/二聚体")
    st.dataframe(pd.DataFrame(rows), width="stretch", hide_index=True)

    pair_options = [
        f"#{p.pair_id} | 产物={p.product_size}bp | F_len={p.forward_length}/R_len={p.reverse_length}"
        for p in pairs
    ]
    current_idx = st.session_state.selected_pair_idx
    if current_idx >= len(pairs):
        current_idx = 0
    sel = st.selectbox(
        "选择候选",
        options=list(range(len(pairs))),
        format_func=lambda i: pair_options[i],
        index=current_idx,
        key="rpa_pair_sel",
    )
    st.session_state.selected_pair_idx = sel
    sp = pairs[sel]

    with st.expander("📋 选中 RPA 引物详情", expanded=True):
        c1, c2 = st.columns(2)
        with c1:
            st.markdown("**正向引物 (5'→3')**")
            st.code(sp.forward_seq, language="text")
            st.write(
                f"长度 {sp.forward_length} nt,GC={sp.forward_gc:.1f}%,Tm={sp.forward_tm:.1f}°C,"
                f"3' 末端发夹 ΔG={sp.forward_3prime_hairpin_dg:.2f} kcal/mol"
            )
        with c2:
            st.markdown("**反向引物 (5'→3')**")
            st.code(sp.reverse_seq, language="text")
            st.write(
                f"长度 {sp.reverse_length} nt,GC={sp.reverse_gc:.1f}%,Tm={sp.reverse_tm:.1f}°C,"
                f"3' 末端发夹 ΔG={sp.reverse_3prime_hairpin_dg:.2f} kcal/mol"
            )
        st.markdown(
            f"**产物大小:** {sp.product_size} bp  |  "
            f"**异二聚体 ΔG:** {sp.heterodimer_dg:.2f} kcal/mol  |  "
            f"**可用:** {'✓' if sp.accepted else '⚠️ ' + '; '.join(sp.issues)}"
        )
        if sp.probe:
            st.markdown("**exo-probe / nfo-probe (含修饰位点建议)**")
            st.code(sp.probe.probe_seq, language="text")
            st.caption(
                f"位点 reporter@{sp.probe.reporter_at} (dT-FAM) · "
                f"THF@{sp.probe.thf_at} · "
                f"quencher@{sp.probe.quencher_at} (dT-BHQ1)"
            )
            st.caption(sp.probe.note)


def _render_lamp_results() -> None:
    """Tab 1 LAMP 分支。"""
    sets: list[LAMPPrimerSet] = st.session_state.lamp_sets or []
    if not sets:
        st.info("点击「开始设计引物探针」生成 LAMP 引物集。")
        return

    rows = []
    for s in sets:
        rows.append({
            "组号": s.set_id,
            "F3": s.F3,
            "B3": s.B3,
            "FIP长度": len(s.FIP),
            "BIP长度": len(s.BIP),
            "F3_Tm": round(s.F3_tm, 1),
            "B3_Tm": round(s.B3_Tm, 1),
            "FIP_F1c_Tm": round(s.FIP_F1c_tm, 1),
            "FIP_F2_Tm": round(s.FIP_F2_tm, 1),
            "BIP_B1c_Tm": round(s.BIP_B1c_tm, 1),
            "BIP_B2_Tm": round(s.BIP_B2_tm, 1),
            "产物(bp)": s.amplicon_size,
            "可用": "✓" if s.accepted else "⚠️",
        })
    st.subheader("LAMP 引物集候选")
    st.caption(
        "LAMP: 等温 60-65°C,4 条核心引物 (F3/B3/FIP/BIP) + 可选 Loop F/B,"
        " 产物 100-1000+ bp 梯状条带"
    )
    st.dataframe(pd.DataFrame(rows), width="stretch", hide_index=True)

    set_options = [f"#{s.set_id} | amplicon={s.amplicon_size}bp" for s in sets]
    current_idx = st.session_state.selected_pair_idx
    if current_idx >= len(sets):
        current_idx = 0
    sel = st.selectbox(
        "选择 LAMP 引物集",
        options=list(range(len(sets))),
        format_func=lambda i: set_options[i],
        index=current_idx,
        key="lamp_set_sel",
    )
    st.session_state.selected_pair_idx = sel
    ss = sets[sel]

    with st.expander("📋 选中 LAMP 引物集详情", expanded=True):
        for label, seq_attr in [
            ("F3", "F3"), ("B3", "B3"), ("FIP (F1c+TTTT+F2)", "FIP"), ("BIP (B1c+TTTT+B2)", "BIP"),
        ]:
            seq = getattr(ss, seq_attr)
            st.markdown(f"**{label}**")
            st.code(seq, language="text")
            st.caption(f"长度 {len(seq)} nt")
        if ss.LoopF:
            st.markdown("**Loop F**")
            st.code(ss.LoopF, language="text")
            st.caption(f"长度 {len(ss.LoopF)} nt")
        if ss.LoopB:
            st.markdown("**Loop B**")
            st.code(ss.LoopB, language="text")
            st.caption(f"长度 {len(ss.LoopB)} nt")
        if ss.issues:
            st.markdown("**问题:** " + "; ".join(ss.issues))
        st.markdown(
            f"**扩增产物大小:** {ss.amplicon_size} bp"
        )


# ============================================================
# Tab 2: 热力学仿真
# ============================================================
def tab_thermo_simulation(seq_info: SequenceInfo) -> None:
    st.header("⚗️ 热力学仿真(多体系对比)")
    # isothermal: thermo tab 大幅降权,但仍展示 3' 末端发夹/异源二聚体(盐参数一致)
    mode = st.session_state.get("design_mode", "qpcr")
    if mode == "rpa":
        _render_rpa_thermo()
        return
    if mode == "lamp":
        _render_lamp_thermo()
        return
    if not st.session_state.pairs:
        st.warning("请先在 Tab 1 设计引物探针")
        return

    pairs = st.session_state.pairs
    presets = st.session_state.preset_systems or load_preset_systems()
    st.session_state.preset_systems = presets

    pair_options = [f"#{p.pair_id} | penalty={p.penalty:.4f}" for p in pairs]
    sel = st.selectbox(
        "选择候选对",
        options=list(range(len(pairs))),
        format_func=lambda i: pair_options[i],
        index=st.session_state.selected_pair_idx,
        key="thermo_pair_sel",
    )
    selected_pair = pairs[sel]
    st.session_state.selected_pair_idx = sel

    if st.button("🔬 评估所有体系", type="primary"):
        with st.spinner("热力学计算中..."):
            results = compare_master_mixes(
                selected_pair.forward_seq,
                selected_pair.reverse_seq,
                selected_pair.probe_seq,
                presets,
            )
            st.session_state.thermo_results = results
        st.success("✓ 完成")

    results: list[ThermoResult] = st.session_state.thermo_results or []
    if not results:
        st.info("点击「评估所有体系」开始计算")
        return

    # 结果表
    rows = []
    for r in results:
        rows.append({
            "体系": r.system_name,
            "F_Tm": round(r.forward_tm, 2),
            "R_Tm": round(r.reverse_tm, 2),
            "P_Tm": round(r.probe_tm, 2) if r.probe_tm else None,
            "Tm差": round(r.tm_difference, 2),
            "F_发夹ΔG": round(r.forward_hairpin_dg, 2),
            "R_发夹ΔG": round(r.reverse_hairpin_dg, 2),
            "P_发夹ΔG": round(r.probe_hairpin_dg, 2) if r.probe_hairpin_dg else None,
            "P_发夹Tm": round(r.probe_hairpin_tm, 2) if r.probe_hairpin_tm else None,
            "异二聚体ΔG": round(r.heterodimer_dg, 2),
            "可用": "✓" if r.acceptable else "⚠️",
            "备注": r.notes,
        })
    df = pd.DataFrame(rows)
    st.subheader("多体系热力学评估")
    st.dataframe(df, width="stretch", hide_index=True)

    # Tm 对比柱状图
    fig = go.Figure()
    systems = [r.system_name for r in results]
    fig.add_trace(go.Bar(name="正向引物 Tm", x=systems, y=[r.forward_tm for r in results]))
    fig.add_trace(go.Bar(name="反向引物 Tm", x=systems, y=[r.reverse_tm for r in results]))
    if any(r.probe_tm for r in results):
        fig.add_trace(go.Bar(name="探针 Tm", x=systems, y=[r.probe_tm for r in results]))
    fig.update_layout(
        title="不同反应体系下的 Tm 对比",
        xaxis_title="反应体系",
        yaxis_title="Tm (°C)",
        barmode="group",
        height=420,
    )
    st.plotly_chart(fig, width="stretch")

    # ΔG 对比
    fig2 = go.Figure()
    fig2.add_trace(go.Bar(
        name="F_发夹ΔG", x=systems, y=[r.forward_hairpin_dg for r in results],
    ))
    fig2.add_trace(go.Bar(
        name="R_发夹ΔG", x=systems, y=[r.reverse_hairpin_dg for r in results],
    ))
    if any(r.probe_hairpin_dg for r in results):
        fig2.add_trace(go.Bar(
            name="P_发夹ΔG", x=systems, y=[r.probe_hairpin_dg for r in results],
        ))
    fig2.add_trace(go.Bar(
        name="异源二聚体ΔG", x=systems, y=[r.heterodimer_dg for r in results],
    ))
    fig2.update_layout(
        title="不同体系下二聚体/发夹 ΔG(越接近 0 越稳定)",
        xaxis_title="反应体系",
        yaxis_title="ΔG (kcal/mol)",
        barmode="group",
        height=420,
    )
    st.plotly_chart(fig2, width="stretch")

    # 推荐体系
    best = recommend_best_system(results)
    if best:
        st.success(f"🎯 推荐体系:**{best.system_name}** —— 引物 Tm 差 {best.tm_difference:.2f}°C")


# ============================================================
# RPA / LAMP 简要热力学评估(无多体系切换,展示发夹/二聚体)
# ============================================================
def _render_rpa_thermo() -> None:
    """RPA 模式的热力学子页。"""
    pairs = st.session_state.rpa_pairs or []
    if not pairs:
        st.info("请先在 Tab 1 设计 RPA 引物")
        return
    sel = st.selectbox(
        "选择候选对",
        options=list(range(len(pairs))),
        format_func=lambda i: f"#{pairs[i].pair_id} | 产物={pairs[i].product_size}bp",
        index=st.session_state.selected_pair_idx,
        key="rpa_thermo_sel",
    )
    st.session_state.selected_pair_idx = sel
    p = pairs[sel]
    rows = [{
        "指标": "正向引物 3' 末端发夹 ΔG",
        "数值": f"{p.forward_3prime_hairpin_dg:.2f} kcal/mol",
    }, {
        "指标": "反向引物 3' 末端发夹 ΔG",
        "数值": f"{p.reverse_3prime_hairpin_dg:.2f} kcal/mol",
    }, {
        "指标": "正反向异源二聚体 ΔG",
        "数值": f"{p.heterodimer_dg:.2f} kcal/mol",
    }, {
        "指标": "RPA 等温温度",
        "数值": "37-42°C (无需严格 Tm)",
    }, {
        "指标": "产物大小",
        "数值": f"{p.product_size} bp",
    }]
    st.caption("RPA 不依赖 master-mix Tm;关键看 3' 末端发夹与异源二聚体")
    st.dataframe(pd.DataFrame(rows), width="stretch", hide_index=True)
    if p.issues:
        st.warning("待优化: " + "; ".join(p.issues))
    else:
        st.success("✓ 关键参数满足 RPA 要求")


def _render_lamp_thermo() -> None:
    """LAMP 模式的热力学子页。"""
    sets = st.session_state.lamp_sets or []
    if not sets:
        st.info("请先在 Tab 1 设计 LAMP 引物集")
        return
    sel = st.selectbox(
        "选择 LAMP 引物集",
        options=list(range(len(sets))),
        format_func=lambda i: f"#{sets[i].set_id} | amplicon={sets[i].amplicon_size}bp",
        index=st.session_state.selected_pair_idx,
        key="lamp_thermo_sel",
    )
    st.session_state.selected_pair_idx = sel
    s = sets[sel]
    rows = [
        {"区段": "F3",  "长度": len(s.F3),  "Tm (°C)": round(s.F3_tm, 1),  "GC (%)": round(s.F3_gc, 1)},
        {"区段": "B3",  "长度": len(s.B3),  "Tm (°C)": round(s.B3_tm, 1),  "GC (%)": round(s.B3_gc, 1)},
        {"区段": "FIP", "长度": len(s.FIP), "F1c_Tm": round(s.FIP_F1c_tm, 1), "F2_Tm": round(s.FIP_F2_tm, 1)},
        {"区段": "BIP", "长度": len(s.BIP), "B1c_Tm": round(s.BIP_B1c_tm, 1), "B2_Tm": round(s.BIP_B2_tm, 1)},
    ]
    st.caption("LAMP 等温扩增 60-65°C, 各段 Tm 是组装建议(预期 FIP/BIP 段 60-65°C, F3/B3 55-60°C)")
    st.dataframe(pd.DataFrame(rows), width="stretch", hide_index=True)
    if s.issues:
        st.warning("待优化: " + "; ".join(s.issues[:5]))
    else:
        st.success("✓ LAMP 引物集组装通过")


# ============================================================
# Tab 3: 扩增曲线仿真
# ============================================================
def tab_amplification(
    seq_info: SequenceInfo,
    amp_params: AmplificationParams,
    copy_numbers: list[float],
    detection_type: str,
) -> None:
    st.header("📈 扩增曲线仿真")
    if not seq_info:
        st.warning("请先在左侧输入或下载序列")
        return
    mode = st.session_state.get("design_mode", "qpcr")
    if detection_type in ("RPA", "LAMP"):
        _render_isothermal(detection_type, copy_numbers)
        return
    # qPCR
    if not st.session_state.pairs:
        st.warning("请先在 Tab 1 设计引物探针")
        return

    pairs = st.session_state.pairs
    sel = st.selectbox(
        "选择候选对",
        options=list(range(len(pairs))),
        format_func=lambda i: f"#{pairs[i].pair_id} | penalty={pairs[i].penalty:.4f}",
        index=st.session_state.selected_pair_idx,
        key="amp_pair_sel",
    )
    selected_pair = pairs[sel]
    st.session_state.selected_pair_idx = sel

    if st.button("▶️ 仿真扩增曲线", type="primary"):
        with st.spinner("仿真中..."):
            curves = simulate_amplification(copy_numbers, amp_params)
            st.session_state.amp_curves = curves
            try:
                std = fit_standard_curve(curves)
                st.session_state.std_curve = std
            except ValueError:
                st.session_state.std_curve = None
            st.session_state.lod_info = estimate_lod(curves)
        st.success("✓ 仿真完成")

    curves = st.session_state.amp_curves or []
    if not curves:
        st.info("点击「仿真扩增曲线」开始模拟")
        return

    # 多条扩增曲线
    fig = go.Figure()
    colors = px.colors.qualitative.Set1
    for i, s in enumerate(curves):
        color = colors[i % len(colors)]
        fig.add_trace(go.Scatter(
            x=s.cycles, y=s.fluorescence,
            mode="lines+markers",
            name=s.label,
            line=dict(color=color, width=2),
            marker=dict(size=4),
        ))
    fig.add_hline(
        y=amp_params.threshold,
        line_dash="dash", line_color="red",
        annotation_text=f"阈值={amp_params.threshold}",
    )
    fig.update_layout(
        title="扩增曲线(各拷贝数梯度)",
        xaxis_title="循环数",
        yaxis_title="荧光强度 (ΔRn)",
        height=480,
        hovermode="x unified",
    )
    st.plotly_chart(fig, width="stretch")

    # Ct 表
    st.subheader("Ct 值表")
    ct_rows = []
    for s in curves:
        ct_rows.append({
            "起始拷贝数": f"{s.start_copies:.0e}",
            "Ct": round(s.ct, 2) if s.ct is not None else "未检出",
        })
    st.dataframe(pd.DataFrame(ct_rows), width="stretch", hide_index=True)

    # 标准曲线
    std = st.session_state.std_curve
    lod = st.session_state.lod_info
    if std is not None:
        st.subheader("标准曲线 (Ct vs log10 拷贝数)")
        fig2 = go.Figure()
        fig2.add_trace(go.Scatter(
            x=std.log_concentrations, y=std.cts,
            mode="markers+lines",
            name="实测点",
            marker=dict(size=10, color="blue"),
            line=dict(color="blue", dash="dot"),
        ))
        # 拟合线
        x_fit = np.linspace(std.log_concentrations.min(), std.log_concentrations.max(), 50)
        y_fit = std.slope * x_fit + std.intercept
        fig2.add_trace(go.Scatter(
            x=x_fit, y=y_fit,
            mode="lines",
            name=f"拟合:y={std.slope:.3f}x+{std.intercept:.2f}",
            line=dict(color="red", width=2),
        ))
        fig2.update_layout(
            xaxis_title="log10(拷贝数)",
            yaxis_title="Ct",
            height=420,
        )
        st.plotly_chart(fig2, width="stretch")

        m1, m2, m3, m4 = st.columns(4)
        m1.metric("斜率", f"{std.slope:.3f}")
        m2.metric("R²", f"{std.r_squared:.4f}")
        m3.metric("扩增效率", f"{std.efficiency_pct:.2f}%")
        m4.metric("线性范围(log10)", f"{std.linear_range[0]:.1f} ~ {std.linear_range[1]:.1f}")

        # 效率判断
        if 90 <= std.efficiency_pct <= 110 and std.r_squared >= 0.98:
            st.success("✅ 标准曲线质量优秀(效率 90~110%, R²≥0.98)")
        elif std.r_squared < 0.95:
            st.error("⚠️ 标准曲线 R² 偏低,建议检查数据点或重新设计")
        else:
            st.warning("⚠️ 标准曲线可优化")

    if lod is not None and lod.get("lod_copies") is not None:
        st.subheader("检测下限 (LOD)")
        st.write(
            f"**估计 LOD:** {lod['lod_copies']:.1f} 拷贝 "
            f"(log10 = {lod['lod_log10']:.2f})"
        )
        st.caption(lod["note"])


# ============================================================
# RPA / LAMP 等温扩增仿真(Tab 3 分支)
# ============================================================
def _render_isothermal(detection_type: str, copy_numbers: list[float]) -> None:
    """Render isothermal amplification model (RPA or LAMP).

    Args:
        detection_type: "RPA" or "LAMP".
        copy_numbers: 拷贝数梯度(从 sidebar 取)。
    """
    if detection_type == "RPA":
        st.caption("RPA: 重组酶聚合酶等温扩增,通常 5-20 分钟出峰")
        # 只支持 RPA 引物
        if not st.session_state.rpa_pairs:
            st.warning("请先在 Tab 1 设计 RPA 引物")
            return
        pairs = st.session_state.rpa_pairs
        sel = st.selectbox(
            "选择 RPA 引物对",
            options=list(range(len(pairs))),
            format_func=lambda i: f"#{pairs[i].pair_id} | 产物={pairs[i].product_size}bp",
            index=st.session_state.selected_pair_idx,
            key="rpa_amp_pair_sel",
        )
        st.session_state.selected_pair_idx = sel

        c1, c2 = st.columns(2)
        with c1:
            duration = st.number_input("反应时长 (分钟)", min_value=5, max_value=60, value=20, step=1)
            efficiency = st.slider("扩增效率 (%)", min_value=50, max_value=100, value=85) / 100.0
        with c2:
            threshold = st.number_input("阈值 (Tt)", min_value=0.001, max_value=1.0, value=0.1, step=0.01, format="%.3f")
            noise = st.slider("噪声水平", min_value=0.0, max_value=0.05, value=0.005, step=0.001, format="%.3f")

        if st.button("▶️ 仿真 RPA 时间-荧光曲线", type="primary"):
            params = RPAIsothermalParams(
                duration_min=float(duration),
                efficiency=efficiency,
                threshold=threshold,
                noise_std=noise,
            )
            curves = simulate_rpa(copy_numbers, params)
            st.session_state.rpa_curves = curves
            st.session_state.rpa_tt = estimate_rpa_tt(curves)
            st.success("✓ RPA 仿真完成")

        curves = st.session_state.rpa_curves or []
        if not curves:
            return
        # 曲线
        fig = go.Figure()
        colors = px.colors.qualitative.Set1
        for i, s in enumerate(curves):
            color = colors[i % len(colors)]
            fig.add_trace(go.Scatter(
                x=s.cycles, y=s.fluorescence,
                mode="lines+markers", name=s.label,
                line=dict(color=color, width=2),
                marker=dict(size=4),
            ))
        fig.add_hline(y=threshold, line_dash="dash", line_color="red",
                      annotation_text=f"阈值={threshold}")
        fig.update_layout(
            title="RPA 时间-荧光曲线",
            xaxis_title="时间 (分钟)", yaxis_title="荧光",
            height=460, hovermode="x unified",
        )
        st.plotly_chart(fig, width="stretch")

        # Tt 表
        rows = []
        for s in curves:
            rows.append({"起始拷贝数": f"{s.start_copies:.0e}",
                         "Tt (分钟)": round(s.ct, 2) if s.ct is not None else "未检出"})
        st.dataframe(pd.DataFrame(rows), width="stretch", hide_index=True)

        tt_info = st.session_state.rpa_tt
        if tt_info and tt_info.get("lod_copies") is not None:
            st.success(
                f"🎯 **最低检出:** {tt_info['lod_copies']:.0e} 拷贝 · "
                f"**Tt:** {tt_info['min_tt']:.2f} 分钟"
            )
            st.caption(tt_info["note"])

    elif detection_type == "LAMP":
        st.caption("LAMP: 环介导等温扩增,通常 15-60 分钟出峰")
        if not st.session_state.lamp_sets:
            st.warning("请先在 Tab 1 设计 LAMP 引物集")
            return
        sets = st.session_state.lamp_sets
        sel = st.selectbox(
            "选择 LAMP 引物集",
            options=list(range(len(sets))),
            format_func=lambda i: f"#{sets[i].set_id} | amplicon={sets[i].amplicon_size}bp",
            index=st.session_state.selected_pair_idx,
            key="lamp_amp_set_sel",
        )
        st.session_state.selected_pair_idx = sel

        c1, c2 = st.columns(2)
        with c1:
            duration = st.number_input("反应时长 (分钟)", min_value=15, max_value=120, value=60, step=5)
            efficiency = st.slider("扩增效率 (%)", min_value=50, max_value=100, value=90) / 100.0
        with c2:
            threshold = st.number_input("阈值 (Tt)", min_value=0.001, max_value=1.0, value=0.1, step=0.01, format="%.3f")
            noise = st.slider("噪声水平", min_value=0.0, max_value=0.05, value=0.005, step=0.001, format="%.3f")

        det_mode = st.radio("检测模式", options=["fluorescence (荧光)", "turbidity (浊度)"],
                            horizontal=True, key="lamp_det_mode")
        det_mode_val = "fluorescence" if det_mode.startswith("fluorescence") else "turbidity"

        if st.button("▶️ 仿真 LAMP 时间-荧光曲线", type="primary"):
            params = LAMPIsothermalParams(
                duration_min=float(duration),
                efficiency=efficiency,
                threshold=threshold,
                noise_std=noise,
                detection_mode=det_mode_val,
            )
            curves = simulate_lamp(copy_numbers, params)
            st.session_state.lamp_curves = curves
            st.session_state.lamp_tt = estimate_lamp_tt(curves)
            st.success("✓ LAMP 仿真完成")

        curves = st.session_state.lamp_curves or []
        if not curves:
            return
        fig = go.Figure()
        colors = px.colors.qualitative.Set1
        for i, s in enumerate(curves):
            color = colors[i % len(colors)]
            fig.add_trace(go.Scatter(
                x=s.cycles, y=s.fluorescence,
                mode="lines+markers", name=s.label,
                line=dict(color=color, width=2),
                marker=dict(size=4),
            ))
        fig.add_hline(y=threshold, line_dash="dash", line_color="red",
                      annotation_text=f"阈值={threshold}")
        fig.update_layout(
            title="LAMP 时间-荧光曲线",
            xaxis_title="时间 (分钟)", yaxis_title="荧光/浊度",
            height=460, hovermode="x unified",
        )
        st.plotly_chart(fig, width="stretch")

        rows = []
        for s in curves:
            rows.append({"起始拷贝数": f"{s.start_copies:.0e}",
                         "Tt (分钟)": round(s.ct, 2) if s.ct is not None else "未检出"})
        st.dataframe(pd.DataFrame(rows), width="stretch", hide_index=True)

        tt_info = st.session_state.lamp_tt
        if tt_info and tt_info.get("lod_copies") is not None:
            st.success(
                f"🎯 **最低检出:** {tt_info['lod_copies']:.0e} 拷贝 · "
                f"**Tt:** {tt_info['min_tt']:.2f} 分钟"
            )
            st.caption(tt_info["note"])


# ============================================================
# Tab 4: 优化循环
# ============================================================
def tab_optimization(
    seq_info: SequenceInfo,
    design_params: DesignParameters,
    amp_params: AmplificationParams,
    copy_numbers: list[float],
    pathogen: str,
    gene: str,
    detection_type: str,
) -> None:
    st.header("🔄 优化循环")
    if detection_type in ("RPA", "LAMP"):
        st.info(
            f"当前检测类型为 {detection_type},优化循环暂仅支持 qPCR。"
            f"等温扩增可在 Tab 3 调整参数重新仿真。"
        )
        return
    if not st.session_state.pairs:
        st.warning("请先在 Tab 1 设计引物探针")
        return

    pairs = st.session_state.pairs
    sel = st.selectbox(
        "选择初始候选",
        options=list(range(len(pairs))),
        format_func=lambda i: f"#{pairs[i].pair_id} | penalty={pairs[i].penalty:.4f}",
        index=st.session_state.selected_pair_idx,
        key="opt_pair_sel",
    )
    initial_pair = pairs[sel]
    st.session_state.selected_pair_idx = sel

    presets = st.session_state.preset_systems or load_preset_systems()
    sys_keys = list(presets.keys())
    sys_labels = [f"{presets[k].name} (Mg={presets[k].mg_conc_mM}mM)" for k in sys_keys]
    sys_sel = st.selectbox(
        "选择反应体系",
        options=sys_keys,
        format_func=lambda k: f"{presets[k].name} (Mg={presets[k].mg_conc_mM}mM, dNTP={presets[k].dntp_conc_mM}mM)",
        index=sys_keys.index(st.session_state.selected_system_key) if st.session_state.selected_system_key in sys_keys else 0,
    )
    st.session_state.selected_system_key = sys_sel
    selected_system = presets[sys_sel]

    if st.button("⚙️ 执行优化循环", type="primary"):
        with st.spinner("优化中..."):
            engine = OptimizationEngine(seq_info.sequence, selected_system)
            opt = engine.optimize(
                initial_pair=initial_pair,
                design_params=design_params,
                copy_numbers=copy_numbers,
                amp_params=amp_params,
            )
            st.session_state.opt_result = opt
            # 重新生成报告
            st.session_state.report_text = report_from_optimization(
                sequence=seq_info.sequence,
                organism=pathogen,
                gene=gene,
                detection_type=detection_type,
                opt_result=opt,
            )
        st.success("✓ 优化完成")

    opt: Optional[OptimizationResult] = st.session_state.opt_result
    if opt is None:
        st.info("点击「执行优化循环」开始")
        return

    # 优化建议
    st.subheader("💡 优化建议")
    for i, s in enumerate(opt.suggestions, 1):
        priority_color = {"high": "🔴", "medium": "🟡", "low": "🟢"}.get(s.priority, "⚪")
        with st.container():
            st.markdown(
                f"{priority_color} **[{s.priority.upper()}] {s.dimension}** — {s.issue}"
            )
            st.markdown(f"&nbsp;&nbsp;&nbsp;&nbsp;👉 建议:{s.action}")

    st.divider()
    st.subheader("📊 优化前后对比")

    # 对比表
    rows = []
    rows.append(_comparison_row("指标", "优化前", "优化后"))
    before_t = opt.before_thermo
    after_t = opt.after_thermo
    rows.append(["F_Tm", f"{opt.before_pair.forward_tm:.2f}",
                 f"{opt.after_pair.forward_tm:.2f}" if opt.after_pair else "-"])
    rows.append(["R_Tm", f"{opt.before_pair.reverse_tm:.2f}",
                 f"{opt.after_pair.reverse_tm:.2f}" if opt.after_pair else "-"])
    rows.append(["P_Tm", f"{opt.before_pair.probe_tm:.2f}" if opt.before_pair.probe_seq else "-",
                 f"{opt.after_pair.probe_tm:.2f}" if opt.after_pair and opt.after_pair.probe_seq else "-"])
    rows.append(["Tm 差", f"{before_t.tm_difference:.2f}",
                 f"{after_t.tm_difference:.2f}" if after_t else "-"])
    rows.append(["产物(bp)", f"{opt.before_pair.product_size}",
                 f"{opt.after_pair.product_size}" if opt.after_pair else "-"])
    if opt.before_std:
        rows.append(["扩增效率(%)",
                     f"{opt.before_std.efficiency_pct:.2f}",
                     f"{opt.after_std.efficiency_pct:.2f}" if opt.after_std else "-"])
        rows.append(["R²",
                     f"{opt.before_std.r_squared:.4f}",
                     f"{opt.after_std.r_squared:.4f}" if opt.after_std else "-"])
    rows.append(["推荐体系", "-", opt.selected_system])

    cmp_df = pd.DataFrame(rows[1:], columns=rows[0])
    st.dataframe(cmp_df, width="stretch", hide_index=True)

    # 对比曲线
    if opt.after_curves:
        st.subheader("📈 优化前后扩增曲线对比")
        fig = go.Figure()
        # 优化前
        colors = px.colors.qualitative.Set1
        for i, s in enumerate(opt.before_curves[:4]):
            fig.add_trace(go.Scatter(
                x=s.cycles, y=s.fluorescence,
                mode="lines",
                name=f"[优化前] {s.label}",
                line=dict(color=colors[i % len(colors)], dash="dot", width=2),
                legendgroup="before",
            ))
        # 优化后
        for i, s in enumerate(opt.after_curves[:4]):
            fig.add_trace(go.Scatter(
                x=s.cycles, y=s.fluorescence,
                mode="lines",
                name=f"[优化后] {s.label}",
                line=dict(color=colors[i % len(colors)], width=3),
                legendgroup="after",
            ))
        fig.add_hline(y=amp_params.threshold, line_dash="dash", line_color="red",
                      annotation_text="阈值")
        fig.update_layout(
            title="优化前后扩增曲线对比",
            xaxis_title="循环数",
            yaxis_title="荧光强度",
            height=480,
            hovermode="x unified",
        )
        st.plotly_chart(fig, width="stretch")


def _comparison_row(a: str, b: str, c: str) -> list:
    return [a, b, c]


# ============================================================
# Tab 5: 报告导出
# ============================================================
def tab_report(
    seq_info: SequenceInfo,
    detection_type: str,
    pathogen: str,
    gene: str,
) -> None:
    st.header("📝 综合报告导出")
    if detection_type == "RPA" and st.session_state.rpa_pairs:
        _render_rpa_report(seq_info, detection_type, pathogen, gene)
        return
    if detection_type == "LAMP" and st.session_state.lamp_sets:
        _render_lamp_report(seq_info, detection_type, pathogen, gene)
        return
    if not st.session_state.pairs:
        st.warning("请先在 Tab 1 设计引物探针")
        return

    if st.session_state.report_text is None:
        # 若尚未执行优化,先基于当前选中直接生成
        pairs = st.session_state.pairs
        sel = st.session_state.selected_pair_idx
        pair = pairs[sel]
        presets = st.session_state.preset_systems or load_preset_systems()
        sys_key = st.session_state.selected_system_key
        if sys_key in presets:
            sys_obj = presets[sys_key]
        else:
            sys_obj = next(iter(presets.values()))
        from core.thermo_sim import evaluate_thermodynamics
        thermo = evaluate_thermodynamics(pair.forward_seq, pair.reverse_seq, pair.probe_seq, sys_obj)
        # 标准曲线
        std = st.session_state.std_curve
        lod = st.session_state.lod_info
        st.session_state.report_text = generate_report(
            sequence=seq_info.sequence if seq_info else "",
            organism=pathogen,
            gene=gene,
            detection_type=detection_type,
            best_pair=pair,
            best_thermo=thermo,
            best_system_name=sys_obj.name,
            std_curve=std,
            lod_info=lod,
        )

    report_text = st.session_state.report_text
    st.text_area("报告内容(可全选复制)", value=report_text, height=520)

    col1, col2 = st.columns(2)
    with col1:
        st.download_button(
            "📥 下载报告 (.txt)",
            data=report_text.encode("utf-8"),
            file_name=f"qpcr_report_{_timestamp()}.txt",
            mime="text/plain",
            width="stretch",
        )
    with col2:
        # 同时导出 CSV:引物探针信息
        pairs = st.session_state.pairs
        if pairs:
            buf = io.StringIO()
            df = pd.DataFrame([p.to_dict() for p in pairs])
            df.to_csv(buf, index=False)
            st.download_button(
                "📊 下载候选表 (.csv)",
                data=buf.getvalue().encode("utf-8"),
                file_name=f"qpcr_candidates_{_timestamp()}.csv",
                mime="text/csv",
                width="stretch",
            )


def _timestamp() -> str:
    from datetime import datetime
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def _render_rpa_report(
    seq_info: SequenceInfo,
    detection_type: str,
    pathogen: str,
    gene: str,
) -> None:
    """RPA 报告导出。"""
    pairs = st.session_state.rpa_pairs
    sel = st.session_state.selected_pair_idx
    pair = pairs[sel] if sel < len(pairs) else pairs[0]
    lines = [
        "=" * 70,
        "        RPA 引物设计方案报告",
        "=" * 70,
        f"病原体: {pathogen or '-'}",
        f"靶基因: {gene or '-'}",
        f"检测类型: RPA",
        f"模板序列长度: {seq_info.length if seq_info else 0} bp",
        "",
        "【最优 RPA 引物对】",
        f"  正向引物 (5'→3'): {pair.forward_seq}",
        f"    长度: {pair.forward_length} nt    GC%: {pair.forward_gc:.1f}    Tm: {pair.forward_tm:.1f}°C",
        f"  反向引物 (5'→3'): {pair.reverse_seq}",
        f"    长度: {pair.reverse_length} nt    GC%: {pair.reverse_gc:.1f}    Tm: {pair.reverse_tm:.1f}°C",
        f"  产物大小: {pair.product_size} bp",
        f"  3' 末端发夹 ΔG: F={pair.forward_3prime_hairpin_dg:.2f}  R={pair.reverse_3prime_hairpin_dg:.2f}",
        f"  异源二聚体 ΔG: {pair.heterodimer_dg:.2f}",
    ]
    if pair.probe:
        lines += [
            "",
            "【exo-probe (含修饰位点建议)】",
            f"  序列 (5'→3'): {pair.probe.probe_seq}",
            f"  长度: {len(pair.probe.probe_seq)} nt",
            f"  修饰位点: dT-FAM@{pair.probe.reporter_at}, THF@{pair.probe.thf_at}, dT-BHQ1@{pair.probe.quencher_at}",
        ]
    if pair.issues:
        lines += ["", "【待优化项】 " + "; ".join(pair.issues)]
    lines.append("=" * 70)
    text = "\n".join(lines)
    st.text_area("RPA 报告内容", value=text, height=420)
    st.download_button(
        "📥 下载 RPA 报告 (.txt)",
        data=text.encode("utf-8"),
        file_name=f"rpa_report_{_timestamp()}.txt",
        mime="text/plain",
        width="stretch",
    )


def _render_lamp_report(
    seq_info: SequenceInfo,
    detection_type: str,
    pathogen: str,
    gene: str,
) -> None:
    """LAMP 报告导出。"""
    sets = st.session_state.lamp_sets
    sel = st.session_state.selected_pair_idx
    s = sets[sel] if sel < len(sets) else sets[0]
    lines = [
        "=" * 70,
        "        LAMP 引物集方案报告",
        "=" * 70,
        f"病原体: {pathogen or '-'}",
        f"靶基因: {gene or '-'}",
        f"检测类型: LAMP",
        f"模板序列长度: {seq_info.length if seq_info else 0} bp",
        "",
        "【LAMP 引物集】",
        f"  F3 (5'→3',{len(s.F3)} nt, Tm={s.F3_tm:.1f}°C): {s.F3}",
        f"  B3 (5'→3',{len(s.B3)} nt, Tm={s.B3_tm:.1f}°C): {s.B3}",
        f"  FIP (5'→3',{len(s.FIP)} nt): {s.FIP}",
        f"    F1c 段 Tm={s.FIP_F1c_tm:.1f}°C, F2 段 Tm={s.FIP_F2_tm:.1f}°C",
        f"  BIP (5'→3',{len(s.BIP)} nt): {s.BIP}",
        f"    B1c 段 Tm={s.BIP_B1c_tm:.1f}°C, B2 段 Tm={s.BIP_B2_tm:.1f}°C",
    ]
    if s.LoopF:
        lines.append(f"  Loop F (5'→3',{len(s.LoopF)} nt): {s.LoopF}")
    if s.LoopB:
        lines.append(f"  Loop B (5'→3',{len(s.LoopB)} nt): {s.LoopB}")
    lines += [
        f"  扩增产物大小: {s.amplicon_size} bp",
        f"  可用: {'✓' if s.accepted else '⚠️ ' + '; '.join(s.issues)}",
        "=" * 70,
    ]
    text = "\n".join(lines)
    st.text_area("LAMP 报告内容", value=text, height=440)
    st.download_button(
        "📥 下载 LAMP 报告 (.txt)",
        data=text.encode("utf-8"),
        file_name=f"lamp_report_{_timestamp()}.txt",
        mime="text/plain",
        width="stretch",
    )


# ============================================================
# 主流程
# ============================================================
def main() -> None:
    st.title("🧬 qPCR 引物探针设计 + 仿真验证平台")
    st.caption(
        "给定核酸序列 → 自动设计引物探针 → 多体系热力学仿真 → 扩增曲线仿真 → "
        "自动优化循环 → 输出综合报告"
    )

    seq_info, design_params, amp_params = render_sidebar()

    tab1, tab2, tab3, tab4, tab5 = st.tabs([
        "🧬 1. 引物探针设计",
        "⚗️ 2. 热力学仿真",
        "📈 3. 扩增曲线仿真",
        "🔄 4. 优化循环",
        "📝 5. 综合报告",
    ])

    pathogen_name = st.session_state.get("pathogen_name", "")
    gene_name = st.session_state.get("gene_name", "")

    with tab1:
        tab_primer_design(seq_info, design_params)

    with tab2:
        tab_thermo_simulation(seq_info)

    with tab3:
        tab_amplification(
            seq_info,
            amp_params,
            st.session_state.copy_numbers,
            design_params.detection_type,
        )

    with tab4:
        # 读取 sidebar 文本输入(从 state 中拿,避免重复访问)
        tab_optimization(
            seq_info,
            design_params,
            amp_params,
            st.session_state.copy_numbers,
            pathogen=pathogen_name,
            gene=gene_name,
            detection_type=design_params.detection_type,
        )

    with tab5:
        tab_report(seq_info, design_params.detection_type, pathogen_name, gene_name)


if __name__ == "__main__":
    main()