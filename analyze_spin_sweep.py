#!/usr/bin/env python3
"""Analyze PHANTOM spin-period Sobol sweep results for numerical stability tests.

Usage:
    python3 analyze_spin_sweep.py sobol_mass_runs/<batch_dir>
    python3 analyze_spin_sweep.py sobol_mass_runs/<batch_a> --compare sobol_mass_runs/<batch_b>

Writes figures and summary tables into <batch_dir>/analysis/
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

P_CRIT_DEFAULT_S = 7233.0  # fluid spin limit from setup_solarsystem (override via setup.log)
STABLE_DISPERSION_MAX = 2.0
STABLE_UNBOUND_MAX = 0.01  # 1%


def _read_p_crit(run_dir: Path) -> float:
    log = run_dir / "setup.log"
    if not log.is_file():
        return P_CRIT_DEFAULT_S
    text = log.read_bytes().decode("utf-8", errors="replace")
    m = re.search(r"fluid spin-limit period\s*=\s*([\d.E+-]+)\s*s", text, re.I)
    if m:
        return float(m.group(1))
    return P_CRIT_DEFAULT_S


def _read_np_apophis(run_dir: Path) -> Optional[int]:
    setup = run_dir / "apophis.setup"
    if not setup.is_file():
        return None
    m = re.search(r"np_apophis\s*=\s*(\d+)", setup.read_text(encoding="utf-8"))
    return int(m.group(1)) if m else None


def load_batch(batch_dir: Path) -> Tuple[pd.DataFrame, Dict[str, object]]:
    csv_path = batch_dir / "sobol_mass_outputs.csv"
    if not csv_path.is_file():
        raise FileNotFoundError(f"Missing {csv_path}")

    df = pd.read_csv(csv_path)
    ok = df[df["status"] == "ok"].copy()
    if ok.empty:
        raise ValueError(f"No successful runs in {csv_path}")

    run1_csv = Path(ok.iloc[0]["run_dir"])
    run1 = run1_csv if run1_csv.is_dir() else (batch_dir / "run_0001")
    p_crit = _read_p_crit(run1)
    np_apophis = _read_np_apophis(run1)
    if np_apophis is None:
        for candidate in sorted(batch_dir.glob("run_*/apophis.setup")):
            np_apophis = _read_np_apophis(candidate.parent)
            if np_apophis is not None:
                break

    if "apophis_spin_period" not in ok.columns:
        cols = ", ".join(ok.columns)
        raise ValueError(
            f"{csv_path} has no apophis_spin_period column (found: {cols}). "
            "Use a spin-period Sobol batch (--spin-period-min/max), not a cohesion or Roche sweep."
        )
    ok["P_in_s"] = ok["apophis_spin_period"].astype(float)
    ok["P_over_Pcrit"] = ok["P_in_s"] / p_crit
    ok["stable"] = (
        (ok["dispersion_ratio"] <= STABLE_DISPERSION_MAX)
        & (ok["unbound_fraction"] <= STABLE_UNBOUND_MAX)
    )
    ok["broken_up"] = ~ok["stable"]

    meta = {
        "batch_dir": str(batch_dir.resolve()),
        "n_runs": int(len(ok)),
        "n_failed": int((df["status"] != "ok").sum()),
        "np_apophis": np_apophis,
        "P_crit_s": p_crit,
        "P_in_min_s": float(ok["P_in_s"].min()),
        "P_in_max_s": float(ok["P_in_s"].max()),
    }
    return ok, meta


def estimate_p_crit_transition(df: pd.DataFrame, p_crit: float) -> Dict[str, float]:
    """Bracket numerical spin limit from stable/unstable labels."""
    d = df.sort_values("P_in_s")
    stable = d["stable"].to_numpy()
    periods = d["P_in_s"].to_numpy()

    # Largest stable P and smallest unstable P
    stable_p = periods[stable] if stable.any() else np.array([np.nan])
    unstable_p = periods[~stable] if (~stable).any() else np.array([np.nan])

    p_low = float(np.nanmax(stable_p)) if stable.any() else float("nan")
    p_high = float(np.nanmin(unstable_p)) if (~stable).any() else float("nan")
    p_mid = (
        0.5 * (p_low + p_high)
        if np.isfinite(p_low) and np.isfinite(p_high)
        else float("nan")
    )

    # Fraction unstable in bins across P/P_crit
    ratio = d["P_over_Pcrit"].to_numpy()
    frac_unstable: List[Dict[str, float]] = []
    edges = np.linspace(ratio.min(), ratio.max(), 9)
    for lo, hi in zip(edges[:-1], edges[1:]):
        mask = (ratio >= lo) & (ratio < hi)
        if mask.sum() == 0:
            continue
        frac_unstable.append(
            {
                "P_over_Pcrit_lo": float(lo),
                "P_over_Pcrit_hi": float(hi),
                "n": int(mask.sum()),
                "frac_unstable": float((~stable[mask]).mean()),
            }
        )

    return {
        "P_stable_max_s": p_low,
        "P_unstable_min_s": p_high,
        "P_numerical_bracket_mid_s": p_mid,
        "P_numerical_over_Pcrit": p_mid / p_crit if np.isfinite(p_mid) else float("nan"),
        "frac_stable": float(stable.mean()),
        "frac_unstable": float((~stable).mean()),
        "binned_frac_unstable": frac_unstable,
    }


def summary_table(df: pd.DataFrame, meta: Dict[str, object]) -> pd.DataFrame:
    p_crit = float(meta["P_crit_s"])
    below = df[df["P_over_Pcrit"] < 0.98]
    above = df[df["P_over_Pcrit"] > 1.02]
    near = df[(df["P_over_Pcrit"] >= 0.98) & (df["P_over_Pcrit"] <= 1.02)]

    def _block(label: str, sub: pd.DataFrame) -> Dict[str, object]:
        if sub.empty:
            return {"region": label, "n": 0}
        return {
            "region": label,
            "n": len(sub),
            "P_over_Pcrit_min": float(sub["P_over_Pcrit"].min()),
            "P_over_Pcrit_max": float(sub["P_over_Pcrit"].max()),
            "dispersion_median": float(sub["dispersion_ratio"].median()),
            "dispersion_max": float(sub["dispersion_ratio"].max()),
            "unbound_median": float(sub["unbound_fraction"].median()),
            "unbound_max": float(sub["unbound_fraction"].max()),
            "frac_stable": float(sub["stable"].mean()),
        }

    rows = [
        _block("below_0.98_Pcrit", below),
        _block("near_critical_0.98-1.02", near),
        _block("above_1.02_Pcrit", above),
        {
            "region": "all",
            "n": len(df),
            "P_over_Pcrit_min": float(df["P_over_Pcrit"].min()),
            "P_over_Pcrit_max": float(df["P_over_Pcrit"].max()),
            "dispersion_median": float(df["dispersion_ratio"].median()),
            "dispersion_max": float(df["dispersion_ratio"].max()),
            "unbound_median": float(df["unbound_fraction"].median()),
            "unbound_max": float(df["unbound_fraction"].max()),
            "frac_stable": float(df["stable"].mean()),
        },
    ]
    out = pd.DataFrame(rows)
    out["np_apophis"] = meta.get("np_apophis")
    out["P_crit_s"] = p_crit
    return out


def plot_sweep(
    df: pd.DataFrame,
    meta: Dict[str, object],
    out_dir: Path,
    transition: Optional[Dict[str, object]] = None,
) -> None:
    p_crit = float(meta["P_crit_s"])
    np_n = meta.get("np_apophis") or "?"
    label = f"n={np_n}"

    fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True)

    ax = axes[0]
    colors = np.where(df["stable"], "#2a9d8f", "#e76f51")
    ax.scatter(df["P_over_Pcrit"], df["dispersion_ratio"], c=colors, s=55, edgecolors="k", linewidths=0.4)
    ax.axhline(STABLE_DISPERSION_MAX, color="gray", ls="--", lw=1, label=f"dispersion = {STABLE_DISPERSION_MAX}")
    ax.axvline(1.0, color="black", ls=":", lw=1.2, label=r"$P_{\mathrm{in}}/P_{\mathrm{crit}}=1$")
    ax.set_ylabel("Dispersion ratio")
    ax.set_title(f"Spin stability sweep ({label}, {meta['n_runs']} runs)")
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    ax.scatter(df["P_over_Pcrit"], df["unbound_fraction"] * 100, c=colors, s=55, edgecolors="k", linewidths=0.4)
    ax.axhline(STABLE_UNBOUND_MAX * 100, color="gray", ls="--", lw=1, label=f"unbound = {STABLE_UNBOUND_MAX*100:.0f}%")
    ax.axvline(1.0, color="black", ls=":", lw=1.2)
    ax.set_xlabel(r"$P_{\mathrm{spin\,in}}\,/\,P_{\mathrm{crit}}$" + f"  ($P_{{\\mathrm{{crit}}}}$={p_crit:.0f} s)")
    ax.set_ylabel("Unbound fraction (%)")
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_dir / "spin_stability_scatter.png", dpi=160)
    plt.close(fig)

    # Sorted strip chart
    d = df.sort_values("P_in_s")
    fig, ax = plt.subplots(figsize=(9, 4))
    x = np.arange(len(d))
    ax.bar(x - 0.2, d["dispersion_ratio"], width=0.4, label="Dispersion", color="#457b9d", alpha=0.85)
    ax2 = ax.twinx()
    ax2.bar(x + 0.2, d["unbound_fraction"] * 100, width=0.4, label="Unbound (%)", color="#e9c46a", alpha=0.9)
    ax.axhline(STABLE_DISPERSION_MAX, color="#457b9d", ls="--", lw=1)
    ax2.axhline(STABLE_UNBOUND_MAX * 100, color="#e9c46a", ls="--", lw=1)
    tick_step = max(1, len(d) // 16)
    ticks = x[::tick_step]
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{p:.0f}" for p in d["P_in_s"].iloc[::tick_step]], rotation=45, ha="right")
    ax.set_xlabel("Input spin period (s)")
    ax.set_ylabel("Dispersion ratio")
    ax2.set_ylabel("Unbound fraction (%)")
    ax.set_title(f"Runs sorted by spin period ({label})")
    fig.tight_layout()
    fig.savefig(out_dir / "spin_metrics_sorted.png", dpi=160)
    plt.close(fig)

    # Binned breakup probability (useful for numerical P_crit estimate)
    bins = transition.get("binned_frac_unstable") or []
    if bins:
        mids = [(b["P_over_Pcrit_lo"] + b["P_over_Pcrit_hi"]) / 2 for b in bins]
        frac = [b["frac_unstable"] for b in bins]
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.plot(mids, frac, "o-", color="#e76f51", lw=2, markersize=8)
        ax.axvline(1.0, color="k", ls=":", lw=1.2)
        ax.set_xlabel(r"$P_{\mathrm{in}}/P_{\mathrm{crit}}$")
        ax.set_ylabel("Fraction of runs unstable")
        ax.set_ylim(-0.05, 1.05)
        ax.set_title(f"Breakup probability vs spin ({label})")
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(out_dir / "breakup_probability.png", dpi=160)
        plt.close(fig)


def plot_resolution_compare(
    batches: List[Tuple[str, pd.DataFrame, Dict[str, object]]], out_dir: Path
) -> None:
    if len(batches) < 2:
        return
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    for name, df, meta in batches:
        axes[0].scatter(
            df["P_over_Pcrit"],
            df["dispersion_ratio"],
            s=45,
            alpha=0.75,
            label=name,
        )
        axes[1].scatter(
            df["P_over_Pcrit"],
            df["unbound_fraction"] * 100,
            s=45,
            alpha=0.75,
            label=name,
        )
    for ax in axes:
        ax.axvline(1.0, color="k", ls=":", lw=1)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=9)
    axes[0].set_ylabel("Dispersion ratio")
    axes[0].set_xlabel(r"$P_{\mathrm{in}}/P_{\mathrm{crit}}$")
    axes[1].set_ylabel("Unbound fraction (%)")
    axes[1].set_xlabel(r"$P_{\mathrm{in}}/P_{\mathrm{crit}}$")
    axes[0].set_title("Resolution comparison")
    fig.tight_layout()
    fig.savefig(out_dir / "resolution_comparison.png", dpi=160)
    plt.close(fig)


def write_report(
    df: pd.DataFrame,
    meta: Dict[str, object],
    transition: Dict[str, object],
    table: pd.DataFrame,
    out_dir: Path,
) -> None:
    table.to_csv(out_dir / "region_summary.csv", index=False)
    df.sort_values("P_in_s").to_csv(out_dir / "runs_enriched.csv", index=False)

    report = {
        "meta": meta,
        "stability_thresholds": {
            "dispersion_ratio_max": STABLE_DISPERSION_MAX,
            "unbound_fraction_max": STABLE_UNBOUND_MAX,
        },
        "transition_estimate": transition,
        "region_summary": table.to_dict(orient="records"),
    }
    (out_dir / "analysis_summary.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8"
    )

    lines = [
        "PHANTOM spin stability analysis",
        "=" * 40,
        f"Batch: {meta['batch_dir']}",
        f"np_apophis: {meta.get('np_apophis', '?')}",
        f"Successful runs: {meta['n_runs']} (failed: {meta.get('n_failed', 0)})",
        f"Analytic P_crit (fluid): {meta['P_crit_s']:.1f} s",
        f"Input period range: {meta['P_in_min_s']:.1f} – {meta['P_in_max_s']:.1f} s",
        "",
        "Stability criteria:",
        f"  stable if dispersion <= {STABLE_DISPERSION_MAX} AND unbound <= {STABLE_UNBOUND_MAX*100:.0f}%",
        "",
        f"Fraction stable: {transition['frac_stable']*100:.1f}%",
        f"Fraction unstable: {transition['frac_unstable']*100:.1f}%",
    ]
    if np.isfinite(transition.get("P_numerical_bracket_mid_s", float("nan"))):
        lines += [
            "",
            "Numerical transition bracket:",
            f"  max stable P_in: {transition['P_stable_max_s']:.1f} s",
            f"  min unstable P_in: {transition['P_unstable_min_s']:.1f} s",
            f"  midpoint: {transition['P_numerical_bracket_mid_s']:.1f} s "
            f"({transition['P_numerical_over_Pcrit']:.3f} P_crit)",
        ]
    lines += ["", "Region summary:", table.to_string(index=False)]
    (out_dir / "analysis_report.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def analyze_batch(batch_dir: Path, out_subdir: str = "analysis") -> Tuple[pd.DataFrame, Dict[str, object]]:
    df, meta = load_batch(batch_dir)
    transition = estimate_p_crit_transition(df, float(meta["P_crit_s"]))
    table = summary_table(df, meta)
    out_dir = batch_dir / out_subdir
    out_dir.mkdir(parents=True, exist_ok=True)
    plot_sweep(df, meta, out_dir, transition=transition)
    write_report(df, meta, transition, table, out_dir)
    print(f"[OK] {batch_dir.name} -> {out_dir}")
    return df, meta


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("batch_dir", type=Path, help="Sweep output folder with sobol_mass_outputs.csv")
    parser.add_argument(
        "--compare",
        type=Path,
        action="append",
        default=[],
        help="Additional batch(es) for resolution comparison plots",
    )
    args = parser.parse_args()

    batch_dir = args.batch_dir.resolve()
    if not batch_dir.is_dir():
        print(f"Not a directory: {batch_dir}", file=sys.stderr)
        return 1

    df, meta = analyze_batch(batch_dir)
    compare_batches: List[Tuple[str, pd.DataFrame, Dict[str, object]]] = [
        (f"n={meta.get('np_apophis', '?')}", df, meta)
    ]
    for other in args.compare:
        try:
            odf, ometa = load_batch(other.resolve())
            compare_batches.append((f"n={ometa.get('np_apophis', '?')}", odf, ometa))
        except (FileNotFoundError, ValueError) as exc:
            print(f"[WARN] skip compare {other}: {exc}", file=sys.stderr)

    out_dir = batch_dir / "analysis"
    if len(compare_batches) > 1:
        plot_resolution_compare(compare_batches, out_dir)

    print(f"Figures and tables in: {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
