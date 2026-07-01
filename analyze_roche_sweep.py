#!/usr/bin/env python3
"""Analyze Roche-limit (Earth distance) sweep results.

Usage:
    python3 analyze_roche_sweep.py sobol_mass_runs/<batch_dir>
    python3 analyze_roche_sweep.py <batch_a> --compare <batch_b> [--label "kc=0"] [--compare "kt=1e5:path/to/batch_b"]

Writes figures and tables into <batch_dir>/analysis/
"""
from __future__ import annotations

import argparse
import json
import math
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

STABLE_DISPERSION_MAX = 2.0
STABLE_UNBOUND_MAX = 0.01

# Distinct colours for overlaying several batches on one figure.
COMPARE_COLORS = ("#264653", "#e76f51", "#2a9d8f", "#e9c46a", "#9b5de5", "#f15bb5")

_EARTH_MSUN = 3.0035e-6
_BODY_BLOCK_RE = re.compile(
    r">\s*reading ephemeris from\s+(\w+)\.txt\s*\n(?:.*\n)*?"
    r"\s*x\s*=\s*([\d.E+-]+)\s+([\d.E+-]+)\s+([\d.E+-]+)\n"
    r"\s*v\s*=\s*([\d.E+-]+)\s+([\d.E+-]+)\s+([\d.E+-]+)",
    re.I,
)


def _pericentre_km_from_state(dr: Tuple[float, float, float], dv: Tuple[float, float, float], mu: float) -> float:
    """Geocentric two-body pericentre (km) from relative state in Phantom code units."""
    hx = dr[1] * dv[2] - dr[2] * dv[1]
    hy = dr[2] * dv[0] - dr[0] * dv[2]
    hz = dr[0] * dv[1] - dr[1] * dv[0]
    h2 = hx * hx + hy * hy + hz * hz
    r = math.sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2])
    if r <= 0.0 or mu <= 0.0:
        return float("nan")
    ex = (dv[1] * hz - dv[2] * hy) / mu - dr[0] / r
    ey = (dv[2] * hx - dv[0] * hz) / mu - dr[1] / r
    ez = (dv[0] * hy - dv[1] * hx) / mu - dr[2] / r
    e = math.sqrt(ex * ex + ey * ey + ez * ez)
    return h2 / (mu * (1.0 + e))


def _parse_setup_roche_metrics(setup_log: Path) -> Dict[str, str]:
    """Read Roche metrics from phantomsetup log (matches run_mass_sobol_phantom parser)."""
    out: Dict[str, str] = {}
    try:
        text = setup_log.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return out
    m_peri = re.search(
        r"pericentre distance\s*=\s*([\d.E+-]+)\s*km\s*\(\s*([\d.E+-]+)\s*R_Earth\s*\)",
        text,
        re.I,
    )
    if m_peri:
        out["pericentre_distance_km"] = f"{float(m_peri.group(1)):.12g}"
        out["pericentre_distance_re"] = f"{float(m_peri.group(2)):.12g}"
    m_geo = re.search(
        r"geocentric separation(?: \(initial\))?\s*=\s*([\d.E+-]+)\s*km\s*\(\s*([\d.E+-]+)\s*R_Earth\s*\)",
        text,
        re.I,
    )
    if m_geo:
        out["geocentric_separation_km"] = f"{float(m_geo.group(1)):.12g}"
        out["geocentric_separation_re"] = f"{float(m_geo.group(2)):.12g}"
    m2 = re.search(
        r"r_tidal is\s+[\d.E+-]+\s+au,\s*([\d.E+-]+)\s*km",
        text,
        re.I,
    )
    if m2:
        r_tidal_km = float(m2.group(1))
        out["r_tidal_km"] = f"{r_tidal_km:.12g}"
        if r_tidal_km > 0 and "pericentre_distance_km" in out:
            out["pericentre_over_r_tidal"] = f"{float(out['pericentre_distance_km']) / r_tidal_km:.12g}"
    return out


def _compute_pericentre_from_ephemeris(setup_log: Path, scale_earth_sep: float) -> Dict[str, float]:
    """Backfill pericentre for older setup logs without formatted pericentre output."""
    try:
        text = setup_log.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return {}
    states: Dict[str, Tuple[Tuple[float, float, float], Tuple[float, float, float]]] = {}
    for m in _BODY_BLOCK_RE.finditer(text):
        name = m.group(1).lower()
        pos = (float(m.group(2)), float(m.group(3)), float(m.group(4)))
        vel = (float(m.group(5)), float(m.group(6)), float(m.group(7)))
        states[name] = (pos, vel)
    if "earth" not in states or "apophis" not in states:
        return {}
    (r_earth, v_earth) = states["earth"]
    (r_apophis, v_apophis) = states["apophis"]
    dr = tuple(r_apophis[i] - r_earth[i] for i in range(3))
    dv = tuple(v_apophis[i] - v_earth[i] for i in range(3))
    if scale_earth_sep != 1.0:
        dr = tuple(scale_earth_sep * dr[i] for i in range(3))
    rp_km = _pericentre_km_from_state(dr, dv, _EARTH_MSUN)
    if not math.isfinite(rp_km):
        return {}
    out = {"pericentre_distance_km": rp_km, "pericentre_distance_re": rp_km / 6371.0}
    m2 = re.search(r"r_tidal is\s+[\d.E+-]+\s+au,\s*([\d.E+-]+)\s*km", text, re.I)
    if m2:
        r_tidal_km = float(m2.group(1))
        if r_tidal_km > 0:
            out["pericentre_over_r_tidal"] = rp_km / r_tidal_km
            out["r_tidal_km"] = r_tidal_km
    return out


def _backfill_pericentre_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Fill missing pericentre metrics from each run's setup.log when possible."""
    if df["peri_km"].notna().all():
        return df
    for idx, row in df.iterrows():
        if pd.notna(row.get("peri_km")):
            continue
        run_dir = row.get("run_dir")
        if not isinstance(run_dir, str) or not run_dir:
            continue
        setup_log = Path(run_dir) / "setup.log"
        if not setup_log.is_file():
            continue
        parsed = _parse_setup_roche_metrics(setup_log)
        scale = float(row["scale_earth_sep"]) if pd.notna(row.get("scale_earth_sep")) else 1.0
        if "pericentre_distance_km" not in parsed:
            computed = _compute_pericentre_from_ephemeris(setup_log, scale)
            parsed = {**computed, **parsed}
        if "pericentre_distance_km" in parsed:
            df.at[idx, "pericentre_distance_km"] = parsed["pericentre_distance_km"]
            df.at[idx, "peri_km"] = float(parsed["pericentre_distance_km"])
        if "pericentre_distance_re" in parsed:
            df.at[idx, "pericentre_distance_re"] = parsed["pericentre_distance_re"]
        if "pericentre_over_r_tidal" in parsed:
            df.at[idx, "pericentre_over_r_tidal"] = parsed["pericentre_over_r_tidal"]
            df.at[idx, "peri_over_rtidal"] = float(parsed["pericentre_over_r_tidal"])
        if pd.isna(row.get("r_tidal_km")) and "r_tidal_km" in parsed:
            df.at[idx, "r_tidal_km"] = parsed["r_tidal_km"]
    return df


def _read_kt_cgs(batch_dir: Path) -> Optional[float]:
    """Read kt_cgs from the first run's apophis.in, if present."""
    for in_path in sorted(batch_dir.glob("run_*/apophis.in")):
        m = re.search(r"^\s*kt_cgs\s*=\s*([\d.E+-]+)", in_path.read_text(encoding="utf-8"), re.M)
        if m:
            return float(m.group(1))
    return None


def default_batch_label(batch_dir: Path) -> str:
    kt = _read_kt_cgs(batch_dir)
    if kt is None:
        return batch_dir.name
    if kt <= 0.0:
        return "kc=0"
    return f"kt={kt:.6g}"


def parse_compare_entry(entry: str) -> Tuple[str, Path]:
    """Parse 'label=path' or plain 'path' (auto-label from kt_cgs or folder name)."""
    if "=" in entry:
        label, raw = entry.split("=", 1)
        path = Path(raw).resolve()
        return label.strip(), path
    path = Path(entry).resolve()
    return default_batch_label(path), path


def x_column(df: pd.DataFrame) -> str:
    if df["peri_over_rtidal"].notna().any():
        return "peri_over_rtidal"
    if df["sep_over_rtidal"].notna().any():
        return "sep_over_rtidal"
    return "scale_earth_sep"


def x_label(xcol: str) -> str:
    if xcol == "peri_over_rtidal":
        return r"$r_p / r_{\mathrm{tidal}}$"
    if xcol == "sep_over_rtidal":
        return r"$d / r_{\mathrm{tidal}}$ (initial separation)"
    return "scale_earth_sep"


def breakup_curve(df: pd.DataFrame, xcol: str) -> Tuple[List[float], List[float]]:
    """Binned fraction disrupted vs separation."""
    d = df.dropna(subset=[xcol]).sort_values(xcol)
    if len(d) < 3:
        return [], []
    edges = np.linspace(d[xcol].min(), d[xcol].max(), min(9, len(d) + 1))
    mids, frac = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        sub = d[(d[xcol] >= lo) & (d[xcol] < hi)]
        if len(sub) == 0:
            continue
        mids.append(float(0.5 * (lo + hi)))
        frac.append(float((~sub["stable"]).mean()))
    return mids, frac


def load_batch(batch_dir: Path) -> Tuple[pd.DataFrame, Dict[str, object]]:
    csv_path = batch_dir / "sobol_mass_outputs.csv"
    if not csv_path.is_file():
        raise FileNotFoundError(f"Missing {csv_path}")
    df = pd.read_csv(csv_path)
    ok = df[df["status"] == "ok"].copy()
    if ok.empty:
        raise ValueError(f"No successful runs in {csv_path}")

    for col in (
        "pericentre_distance_km",
        "pericentre_over_r_tidal",
        "geocentric_separation_km",
        "separation_over_r_tidal",
        "scale_earth_sep",
    ):
        if col not in ok.columns:
            ok[col] = float("nan")

    ok["peri_km"] = pd.to_numeric(ok["pericentre_distance_km"], errors="coerce")
    ok["peri_over_rtidal"] = pd.to_numeric(ok["pericentre_over_r_tidal"], errors="coerce")
    ok["sep_km"] = pd.to_numeric(ok["geocentric_separation_km"], errors="coerce")
    ok["sep_over_rtidal"] = pd.to_numeric(ok["separation_over_r_tidal"], errors="coerce")
    ok["scale_earth_sep"] = pd.to_numeric(ok["scale_earth_sep"], errors="coerce")
    ok["stable"] = (
        (ok["dispersion_ratio"] <= STABLE_DISPERSION_MAX)
        & (ok["unbound_fraction"] <= STABLE_UNBOUND_MAX)
    )
    ok = _backfill_pericentre_columns(ok)

    meta = {
        "batch_dir": str(batch_dir.resolve()),
        "n_runs": int(len(ok)),
        "n_failed": int((df["status"] != "ok").sum()),
        "peri_km_min": float(ok["peri_km"].min()) if ok["peri_km"].notna().any() else float("nan"),
        "peri_km_max": float(ok["peri_km"].max()) if ok["peri_km"].notna().any() else float("nan"),
        "sep_km_min": float(ok["sep_km"].min()) if ok["sep_km"].notna().any() else float("nan"),
        "sep_km_max": float(ok["sep_km"].max()) if ok["sep_km"].notna().any() else float("nan"),
        "r_tidal_km": float(pd.to_numeric(ok.get("r_tidal_km", pd.Series(dtype=float)), errors="coerce").iloc[0])
        if "r_tidal_km" in ok.columns and ok["r_tidal_km"].notna().any()
        else float("nan"),
    }
    return ok, meta


def estimate_roche_bracket(df: pd.DataFrame) -> Dict[str, object]:
    xcol = x_column(df)
    d = df.sort_values(xcol)
    stable = d["stable"].to_numpy()
    x = d[xcol].to_numpy()
    stable_x = x[stable] if stable.any() else np.array([np.nan])
    unstable_x = x[~stable] if (~stable).any() else np.array([np.nan])
    out = {
        "x_column": xcol,
        "frac_stable": float(stable.mean()),
    }
    if xcol == "peri_over_rtidal":
        out["max_stable_peri_over_rtidal"] = float(np.nanmax(stable_x)) if stable.any() else float("nan")
        out["min_unstable_peri_over_rtidal"] = float(np.nanmin(unstable_x)) if (~stable).any() else float("nan")
    else:
        out["max_stable_sep_over_rtidal"] = float(np.nanmax(stable_x)) if stable.any() else float("nan")
        out["min_unstable_sep_over_rtidal"] = float(np.nanmin(unstable_x)) if (~stable).any() else float("nan")
    return out


def plot_roche(df: pd.DataFrame, meta: Dict[str, object], out_dir: Path, label: str) -> None:
    xcol = x_column(df)
    xlabel = x_label(xcol)
    colors = np.where(df["stable"], "#2a9d8f", "#e76f51")

    fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
    axes[0].scatter(df[xcol], df["dispersion_ratio"], c=colors, s=55, edgecolors="k", linewidths=0.4)
    axes[0].axhline(STABLE_DISPERSION_MAX, color="gray", ls="--", lw=1)
    axes[0].axvline(1.0, color="k", ls=":", lw=1.2, label=r"$r_p=r_{\mathrm{tidal}}$ (fluid)")
    axes[0].set_ylabel("Dispersion ratio")
    axes[0].set_title(f"Roche sweep — {label} ({meta['n_runs']} runs)")
    axes[0].legend(fontsize=9)
    axes[0].grid(True, alpha=0.3)

    axes[1].scatter(df[xcol], df["unbound_fraction"] * 100, c=colors, s=55, edgecolors="k", linewidths=0.4)
    axes[1].axhline(STABLE_UNBOUND_MAX * 100, color="gray", ls="--", lw=1)
    axes[1].axvline(1.0, color="k", ls=":", lw=1.2)
    axes[1].set_xlabel(xlabel)
    axes[1].set_ylabel("Unbound fraction (%)")
    axes[1].grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_dir / "roche_stability_scatter.png", dpi=160)
    plt.close(fig)

    # Breakup probability vs separation bins
    if df[xcol].notna().sum() >= 3:
        d = df.dropna(subset=[xcol]).sort_values(xcol)
        edges = np.linspace(d[xcol].min(), d[xcol].max(), min(9, len(d) + 1))
        mids, frac = [], []
        for lo, hi in zip(edges[:-1], edges[1:]):
            sub = d[(d[xcol] >= lo) & (d[xcol] < hi)]
            if len(sub) == 0:
                continue
            mids.append(0.5 * (lo + hi))
            frac.append(float((~sub["stable"]).mean()))
        if mids:
            fig, ax = plt.subplots(figsize=(7, 4))
            ax.plot(mids, frac, "o-", color="#e76f51", lw=2, markersize=8)
            ax.axvline(1.0, color="k", ls=":", lw=1.2)
            ax.set_xlabel(xlabel)
            ax.set_ylabel("Fraction disrupted")
            ax.set_ylim(-0.05, 1.05)
            ax.set_title(f"Disruption probability — {label}")
            ax.grid(True, alpha=0.3)
            fig.tight_layout()
            fig.savefig(out_dir / "roche_breakup_probability.png", dpi=160)
            plt.close(fig)


def plot_batch_compare(
    batches: List[Tuple[str, pd.DataFrame, Dict[str, object]]], out_dir: Path
) -> None:
    """Overlay dispersion / unbound / breakup probability for multiple Roche sweeps."""
    if len(batches) < 2:
        return

    xcol = x_column(batches[0][1])
    xlabel = x_label(xcol)

    fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
    for i, (name, df, _meta) in enumerate(batches):
        color = COMPARE_COLORS[i % len(COMPARE_COLORS)]
        axes[0].scatter(
            df[xcol],
            df["dispersion_ratio"],
            s=55,
            alpha=0.85,
            color=color,
            edgecolors="k",
            linewidths=0.35,
            label=name,
        )
        axes[1].scatter(
            df[xcol],
            df["unbound_fraction"] * 100,
            s=55,
            alpha=0.85,
            color=color,
            edgecolors="k",
            linewidths=0.35,
            label=name,
        )

    axes[0].axhline(STABLE_DISPERSION_MAX, color="gray", ls="--", lw=1)
    axes[1].axhline(STABLE_UNBOUND_MAX * 100, color="gray", ls="--", lw=1)
    for ax in axes:
        ax.axvline(1.0, color="k", ls=":", lw=1.2)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=9)

    axes[0].set_ylabel("Dispersion ratio")
    axes[0].set_title("Roche sweep comparison")
    axes[1].set_xlabel(xlabel)
    axes[1].set_ylabel("Unbound fraction (%)")
    fig.tight_layout()
    fig.savefig(out_dir / "roche_cohesion_comparison.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4))
    for i, (name, df, _meta) in enumerate(batches):
        mids, frac = breakup_curve(df, xcol)
        if not mids:
            continue
        color = COMPARE_COLORS[i % len(COMPARE_COLORS)]
        ax.plot(mids, frac, "o-", lw=2, markersize=7, color=color, label=name)
    ax.axvline(1.0, color="k", ls=":", lw=1.2)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Fraction disrupted")
    ax.set_ylim(-0.05, 1.05)
    ax.set_title("Disruption probability comparison")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_dir / "roche_cohesion_breakup_comparison.png", dpi=160)
    plt.close(fig)


def write_outputs(df: pd.DataFrame, meta: Dict[str, object], bracket: Dict[str, object], out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    sort_col = "peri_km" if df["peri_km"].notna().any() else "sep_km"
    df.sort_values(sort_col).to_csv(out_dir / "runs_enriched.csv", index=False)
    report = {"meta": meta, "bracket": bracket, "stability_thresholds": {
        "dispersion_max": STABLE_DISPERSION_MAX, "unbound_max": STABLE_UNBOUND_MAX}}
    (out_dir / "analysis_summary.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    lines = [
        "Roche limit sweep analysis",
        "=" * 40,
        f"Runs OK: {meta['n_runs']}  failed: {meta.get('n_failed', 0)}",
        f"Fluid r_tidal (setup): {meta.get('r_tidal_km', float('nan')):.1f} km",
        f"Fraction stable: {bracket['frac_stable']*100:.1f}%",
    ]
    if np.isfinite(meta.get("peri_km_min", float("nan"))):
        lines.insert(
            4,
            f"Pericentre range: {meta['peri_km_min']:.1f} – {meta['peri_km_max']:.1f} km",
        )
    if np.isfinite(meta.get("sep_km_min", float("nan"))):
        lines.append(
            f"Initial separation range: {meta['sep_km_min']:.1f} – {meta['sep_km_max']:.1f} km"
        )
    xcol = bracket.get("x_column", "sep_over_rtidal")
    if xcol == "peri_over_rtidal":
        if np.isfinite(bracket.get("max_stable_peri_over_rtidal", float("nan"))):
            lines.append(
                f"Max stable r_p/r_tidal: {bracket['max_stable_peri_over_rtidal']:.3f}; "
                f"Min unstable r_p/r_tidal: {bracket['min_unstable_peri_over_rtidal']:.3f}"
            )
    elif np.isfinite(bracket.get("max_stable_sep_over_rtidal", float("nan"))):
        lines.append(
            f"Max stable d/r_tidal: {bracket['max_stable_sep_over_rtidal']:.3f}; "
            f"Min unstable d/r_tidal: {bracket['min_unstable_sep_over_rtidal']:.3f}"
        )
    (out_dir / "analysis_report.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("batch_dir", type=Path)
    parser.add_argument(
        "--compare",
        action="append",
        default=[],
        metavar="[LABEL=]BATCH_DIR",
        help="Additional batch to overlay (optional 'label=path'; auto-labels kt_cgs from apophis.in)",
    )
    parser.add_argument("--label", default=None, help="Plot label for primary batch")
    args = parser.parse_args()

    batch_dir = args.batch_dir.resolve()
    df, meta = load_batch(batch_dir)
    bracket = estimate_roche_bracket(df)
    label = args.label or batch_dir.name
    out_dir = batch_dir / "analysis"
    out_dir.mkdir(parents=True, exist_ok=True)
    plot_roche(df, meta, out_dir, label)
    write_outputs(df, meta, bracket, out_dir)
    print(f"[OK] {out_dir}")

    if args.compare:
        compare_batches: List[Tuple[str, pd.DataFrame, Dict[str, object]]] = [
            (label, df, meta),
        ]
        for entry in args.compare:
            try:
                clabel, cpath = parse_compare_entry(entry)
                cdf, cmeta = load_batch(cpath)
                compare_batches.append((clabel, cdf, cmeta))
            except (FileNotFoundError, ValueError) as exc:
                print(f"[WARN] skip compare {entry!r}: {exc}", file=sys.stderr)
        if len(compare_batches) > 1:
            plot_batch_compare(compare_batches, out_dir)
            print(f"[OK] comparison figures in {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
