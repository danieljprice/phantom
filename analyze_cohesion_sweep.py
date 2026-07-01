#!/usr/bin/env python3
"""Analyze PHANTOM cohesion (kt_cgs) calibration sweep results.

Typical use: apophis_only=T, fixed fast spin (e.g. P=6500 s), Sobol over kt_cgs.

Usage:
    python3 analyze_cohesion_sweep.py sobol_mass_runs/<batch_dir>
    python3 analyze_cohesion_sweep.py <batch> --kt-literature 3e7

Writes figures and tables into <batch_dir>/analysis/
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
from matplotlib.ticker import LogLocator, LogFormatterMathtext

STABLE_DISPERSION_MAX = 2.0
STABLE_UNBOUND_MAX = 0.01
KT_LITERATURE_25PA_DEFAULT = 3.0e7  # σ≈25 Pa → kt for ~14 m DEM grains (thesis estimate)


def _format_kt_log_axis(ax, kt_min: float, kt_max: float) -> None:
    """Readable log ticks for kt_cgs (matplotlib often shows only one label by default)."""
    lo = max(float(kt_min), 1.0)
    hi = max(float(kt_max), lo * 10)
    exp_lo = int(np.floor(np.log10(lo)))
    exp_hi = int(np.ceil(np.log10(hi)))
    ticks = [10.0**e for e in range(exp_lo, exp_hi + 1)]
    ax.set_xscale("log")
    ax.set_xlim(max(lo * 0.8, 10.0**exp_lo), min(hi * 1.2, 10.0**exp_hi))
    ax.set_xticks(ticks)
    ax.xaxis.set_major_formatter(LogFormatterMathtext(base=10))
    ax.xaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1))


def _read_spin_period_s(run_dir: Path) -> Optional[float]:
    setup = run_dir / "apophis.setup"
    if not setup.is_file():
        return None
    m = re.search(r"apophis_spin_period\s*=\s*([\d.E+-]+)", setup.read_text(encoding="utf-8"))
    return float(m.group(1)) if m else None


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

    if "kt_cgs" not in ok.columns:
        raise ValueError("CSV has no kt_cgs column — is this a cohesion sweep batch?")

    run1 = Path(ok.iloc[0]["run_dir"])
    if not run1.is_dir():
        run1 = batch_dir / "run_0001"

    ok["kt_cgs"] = pd.to_numeric(ok["kt_cgs"], errors="coerce")
    ok["stable"] = (
        (ok["dispersion_ratio"] <= STABLE_DISPERSION_MAX)
        & (ok["unbound_fraction"] <= STABLE_UNBOUND_MAX)
    )
    ok["broken_up"] = ~ok["stable"]

    p_spin = _read_spin_period_s(run1)
    np_apophis = _read_np_apophis(run1)

    meta = {
        "batch_dir": str(batch_dir.resolve()),
        "n_runs": int(len(ok)),
        "n_failed": int((df["status"] != "ok").sum()),
        "np_apophis": np_apophis,
        "P_spin_s": p_spin,
        "kt_min": float(ok["kt_cgs"].min()),
        "kt_max": float(ok["kt_cgs"].max()),
    }
    return ok, meta


def estimate_kt_transition(df: pd.DataFrame) -> Dict[str, object]:
    """Bracket minimum kt_cgs for stability (sorted ascending by cohesion)."""
    d = df.sort_values("kt_cgs")
    stable = d["stable"].to_numpy()
    kt = d["kt_cgs"].to_numpy()

    stable_kt = kt[stable] if stable.any() else np.array([np.nan])
    unstable_kt = kt[~stable] if (~stable).any() else np.array([np.nan])

    kt_stable_min = float(np.nanmin(stable_kt)) if stable.any() else float("nan")
    kt_unstable_max = float(np.nanmax(unstable_kt)) if (~stable).any() else float("nan")
    kt_mid = (
        0.5 * (kt_stable_min + kt_unstable_max)
        if np.isfinite(kt_stable_min) and np.isfinite(kt_unstable_max)
        else float("nan")
    )

    return {
        "kt_stable_min": kt_stable_min,
        "kt_unstable_max": kt_unstable_max,
        "kt_bracket_mid": kt_mid,
        "frac_stable": float(stable.mean()),
        "frac_unstable": float((~stable).mean()),
    }


def plot_sweep(
    df: pd.DataFrame,
    meta: Dict[str, object],
    out_dir: Path,
    kt_lit: float,
    transition: Dict[str, object],
) -> None:
    np_n = meta.get("np_apophis") or "?"
    p_spin = meta.get("P_spin_s")
    spin_note = f"P={p_spin:.0f} s" if p_spin else "fixed spin"
    label = f"n={np_n}, {spin_note}"

    colors = np.where(df["stable"], "#2a9d8f", "#e76f51")
    kt_pos = df["kt_cgs"].clip(lower=1.0)  # log axis; treat 0 as 1 for display

    fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True)

    axes[0].scatter(kt_pos, df["dispersion_ratio"], c=colors, s=55, edgecolors="k", linewidths=0.4)
    axes[0].axhline(STABLE_DISPERSION_MAX, color="gray", ls="--", lw=1)
    axes[0].axvline(kt_lit, color="k", ls=":", lw=1.2, label=rf"$\sigma$=25 Pa $\rightarrow$ $k_t={kt_lit:.2g}$")
    axes[0].set_ylabel("Dispersion ratio")
    axes[0].set_title(f"Cohesion calibration — {label} ({meta['n_runs']} runs)")
    _format_kt_log_axis(axes[0], float(df["kt_cgs"].min()), float(df["kt_cgs"].max()))
    axes[0].legend(fontsize=9)
    axes[0].grid(True, alpha=0.3, which="both")

    axes[1].scatter(kt_pos, df["unbound_fraction"] * 100, c=colors, s=55, edgecolors="k", linewidths=0.4)
    axes[1].axhline(STABLE_UNBOUND_MAX * 100, color="gray", ls="--", lw=1)
    axes[1].axvline(kt_lit, color="k", ls=":", lw=1.2)
    axes[1].set_xlabel(r"$k_{t,\mathrm{cgs}}$ (g/s$^2$ per cm gap)")
    axes[1].set_ylabel("Unbound fraction (%)")
    _format_kt_log_axis(axes[1], float(df["kt_cgs"].min()), float(df["kt_cgs"].max()))
    axes[1].grid(True, alpha=0.3, which="both")

    fig.tight_layout()
    fig.savefig(out_dir / "cohesion_stability_scatter.png", dpi=160)
    plt.close(fig)

    d = df.sort_values("kt_cgs")
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
    ax.set_xticklabels([f"{v:.2g}" for v in d["kt_cgs"].iloc[::tick_step]], rotation=45, ha="right")
    ax.set_xlabel(r"$k_{t,\mathrm{cgs}}$ (sorted)")
    ax.set_ylabel("Dispersion ratio")
    ax2.set_ylabel("Unbound fraction (%)")
    ax.set_title(f"Runs sorted by cohesion ({label})")
    fig.tight_layout()
    fig.savefig(out_dir / "cohesion_metrics_sorted.png", dpi=160)
    plt.close(fig)

    # Binned stability vs log(kt)
    d2 = df.dropna(subset=["kt_cgs"]).sort_values("kt_cgs")
    if len(d2) >= 4:
        log_kt = np.log10(d2["kt_cgs"].clip(lower=1.0))
        edges = np.linspace(log_kt.min(), log_kt.max(), min(8, len(d2)))
        mids, frac = [], []
        for lo, hi in zip(edges[:-1], edges[1:]):
            sub = d2[(log_kt >= lo) & (log_kt < hi)]
            if len(sub) == 0:
                continue
            mids.append(10 ** (0.5 * (lo + hi)))
            frac.append(float((~sub["stable"]).mean()))
        if mids:
            fig, ax = plt.subplots(figsize=(7, 4))
            ax.plot(mids, frac, "o-", color="#e76f51", lw=2, markersize=8)
            ax.axvline(kt_lit, color="k", ls=":", lw=1.2, label=rf"$\sigma$=25 Pa ($k_t={kt_lit:.2g}$)")
            _format_kt_log_axis(ax, float(d2["kt_cgs"].min()), float(d2["kt_cgs"].max()))
            ax.set_xlabel(r"$k_{t,\mathrm{cgs}}$ (g/s$^2$ per cm gap)")
            ax.set_ylabel("Fraction unstable")
            ax.set_ylim(-0.05, 1.05)
            ax.set_title(f"Instability vs cohesion ({label})")
            ax.legend(fontsize=9)
            ax.grid(True, alpha=0.3, which="both")
            fig.tight_layout()
            fig.savefig(out_dir / "cohesion_breakup_probability.png", dpi=160)
            plt.close(fig)


def write_report(
    df: pd.DataFrame,
    meta: Dict[str, object],
    transition: Dict[str, object],
    kt_lit: float,
    out_dir: Path,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    df.sort_values("kt_cgs").to_csv(out_dir / "runs_enriched.csv", index=False)

    report = {
        "meta": meta,
        "stability_thresholds": {
            "dispersion_max": STABLE_DISPERSION_MAX,
            "unbound_max": STABLE_UNBOUND_MAX,
        },
        "kt_literature_25Pa": kt_lit,
        "transition_estimate": transition,
    }
    (out_dir / "analysis_summary.json").write_text(json.dumps(report, indent=2), encoding="utf-8")

    lines = [
        "PHANTOM cohesion calibration analysis",
        "=" * 40,
        f"Batch: {meta['batch_dir']}",
        f"np_apophis: {meta.get('np_apophis', '?')}",
        f"Fixed spin period: {meta.get('P_spin_s', '?')} s",
        f"Successful runs: {meta['n_runs']} (failed: {meta.get('n_failed', 0)})",
        f"kt_cgs range: {meta['kt_min']:.6g} – {meta['kt_max']:.6g}",
        "",
        "Stability criteria:",
        f"  stable if dispersion <= {STABLE_DISPERSION_MAX} AND unbound <= {STABLE_UNBOUND_MAX*100:.0f}%",
        "",
        f"Literature reference (σ=25 Pa): kt ≈ {kt_lit:.6g}",
        f"Fraction stable: {transition['frac_stable']*100:.1f}%",
    ]
    if np.isfinite(transition.get("kt_stable_min", float("nan"))):
        lines += [
            "",
            "Numerical kt bracket (minimum cohesion for stability):",
            f"  min stable kt_cgs: {transition['kt_stable_min']:.6g}",
            f"  max unstable kt_cgs: {transition['kt_unstable_max']:.6g}",
            f"  midpoint: {transition['kt_bracket_mid']:.6g}",
        ]
    (out_dir / "analysis_report.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("batch_dir", type=Path)
    parser.add_argument(
        "--kt-literature",
        type=float,
        default=KT_LITERATURE_25PA_DEFAULT,
        help="Reference kt_cgs from literature σ≈25 Pa (vertical line on plots)",
    )
    args = parser.parse_args()

    batch_dir = args.batch_dir.resolve()
    df, meta = load_batch(batch_dir)
    transition = estimate_kt_transition(df)
    out_dir = batch_dir / "analysis"
    out_dir.mkdir(parents=True, exist_ok=True)
    plot_sweep(df, meta, out_dir, args.kt_literature, transition)
    write_report(df, meta, transition, args.kt_literature, out_dir)
    print(f"[OK] Figures and tables in: {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
