"""多病原体 LAMP/RPA 设计验证 — 改进前基线 + 改进后对比。

用 5 个真实病原体靶序列测试设计成功率, 与文献基准对照 Tt 值。
"""
from __future__ import annotations
import sys, json, time
from pathlib import Path

HERE = Path(__file__).resolve().parent
PROJ_ROOT = HERE.parent
sys.path.insert(0, str(PROJ_ROOT))

from core.lamp_design import design_lamp_primers, LAMPDesignParameters
from core.rpa_design import design_rpa_primers, RPADesignParameters
from core.amplification import (
    simulate_lamp, simulate_rpa,
    LAMPIsothermalParams, RPAIsothermalParams,
)


def load_targets() -> list[dict]:
    with open(HERE / "real_pathogen_sequences.json") as f:
        return json.load(f)


def test_lamp(seq: str, params: LAMPDesignParameters | None = None):
    params = params or LAMPDesignParameters(num_sets=3, max_candidates=500)
    t0 = time.time()
    try:
        sets = design_lamp_primers(seq, params)
        elapsed = time.time() - t0
        accepted = [s for s in sets if s.accepted]
        return {
            "total": len(sets), "accepted": len(accepted), "time": elapsed,
            "best": accepted[0] if accepted else (sets[0] if sets else None),
        }
    except Exception as e:
        return {"total": 0, "accepted": 0, "time": 0, "error": str(e)}


def test_rpa(seq: str, params: RPADesignParameters | None = None):
    params = params or RPADesignParameters(num_return=5, design_probe=True)
    t0 = time.time()
    try:
        pairs = design_rpa_primers(seq, params)
        elapsed = time.time() - t0
        accepted = [p for p in pairs if p.accepted]
        return {
            "total": len(pairs), "accepted": len(accepted), "time": elapsed,
            "best": accepted[0] if accepted else (pairs[0] if pairs else None),
        }
    except Exception as e:
        return {"total": 0, "accepted": 0, "time": 0, "error": str(e)}


def test_isothermal(copies: list[float] | None = None):
    """等温扩增仿真, 返回 LAMP/RPA 的 Tt 值。"""
    copies = copies or [10, 100, 1000, 10000, 100000]
    lamp_params = LAMPIsothermalParams()
    rpa_params = RPAIsothermalParams()
    lamp_samples = simulate_lamp(copies, lamp_params)
    rpa_samples = simulate_rpa(copies, rpa_params)
    return {
        "lamp": {s.start_copies: s.ct for s in lamp_samples},
        "rpa": {s.start_copies: s.ct for s in rpa_samples},
    }


def main():
    print(f"# 多病原体 LAMP/RPA 设计验证")
    print(f"# 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    targets = load_targets()
    print(f"# 病原体数: {len(targets)}")

    print("\n" + "━" * 72)
    print(f"{'病原体':20s} | {'GC%':>5s} | {'LAMP':>12s} | {'RPA':>12s}")
    print("-" * 72)

    lamp_success = 0
    rpa_success = 0
    lamp_results = []
    rpa_results = []

    for t in targets:
        name = t["name"]
        seq = t["sequence"]
        gc = t["gc"]

        l = test_lamp(seq)
        r = test_rpa(seq)
        lamp_results.append((name, l))
        rpa_results.append((name, r))

        l_str = f"{l['accepted']}/{l['total']}" if l["total"] > 0 else "ERR"
        r_str = f"{r['accepted']}/{r['total']}" if r["total"] > 0 else "ERR"
        print(f"{name:20s} | {gc:>5.1f} | {l_str:>12s} | {r_str:>12s}")

        if l.get("accepted", 0) >= 1:
            lamp_success += 1
        if r.get("accepted", 0) >= 1:
            rpa_success += 1

    print("\n" + "━" * 72)
    print(f"\nLAMP 设计成功率: {lamp_success}/{len(targets)} = {lamp_success/len(targets)*100:.0f}%")
    print(f"RPA 设计成功率:  {rpa_success}/{len(targets)} = {rpa_success/len(targets)*100:.0f}%")
    print(f"验收标准: ≥80%")

    # 等温扩增仿真
    print(f"\n{'━' * 72}")
    print("等温扩增仿真 Tt 值:")
    iso = test_isothermal()
    print(f"{'拷贝数':>12s} | {'LAMP Tt':>10s} | {'RPA Tt':>10s}")
    for c in sorted(iso["lamp"].keys()):
        l_t = iso["lamp"][c]
        r_t = iso["rpa"].get(c)
        l_str = f"{l_t:.1f} min" if l_t else "—"
        r_str = f"{r_t:.1f} min" if r_t else "—"
        print(f"{c:>12,.0f} | {l_str:>10s} | {r_str:>10s}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
