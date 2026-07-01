#!/usr/bin/env python3
"""Recompute dispersion_ratio and unbound_fraction from existing sink dumps.

Use after updating extract_breakup_metrics (e.g. main-body unbound fix)
without rerunning PHANTOM.

Usage:
    python3 recompute_sweep_metrics.py sobol_mass_runs/<batch_dir>
    python3 recompute_sweep_metrics.py sobol_mass_runs/<batch_dir> --sink-apophis-id 1
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_mass_sobol_phantom import APOPHIS_SINK_ID_DEFAULT, extract_breakup_metrics  # noqa: E402


def _apophis_sink_id(run_dir: Path, prefix: str, configured: int) -> int:
    setup = run_dir / f"{prefix}.setup"
    if setup.is_file() and re.search(r"apophis_only\s*=\s*T", setup.read_text(encoding="utf-8"), re.I):
        if configured == APOPHIS_SINK_ID_DEFAULT:
            return 1
    return configured


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("batch_dir", type=Path)
    parser.add_argument("--prefix", default="apophis")
    parser.add_argument("--sink-apophis-id", type=int, default=11)
    args = parser.parse_args()

    batch_dir = args.batch_dir.resolve()
    csv_path = batch_dir / "sobol_mass_outputs.csv"
    if not csv_path.is_file():
        print(f"Missing {csv_path}", file=sys.stderr)
        return 1

    with csv_path.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        print("Empty CSV", file=sys.stderr)
        return 1

    fieldnames = list(rows[0].keys())
    updated = 0
    for row in rows:
        if row.get("status") != "ok":
            continue
        run_dir = Path(row["run_dir"])
        if not run_dir.is_dir():
            run_dir = batch_dir / f"run_{int(row['run_id']):04d}"
        sink_id = _apophis_sink_id(run_dir, args.prefix, args.sink_apophis_id)
        disp, unb = extract_breakup_metrics(run_dir, args.prefix, sink_id)
        row["dispersion_ratio"] = f"{disp:.12g}" if disp == disp else ""
        row["unbound_fraction"] = f"{unb:.12g}" if unb == unb else ""
        updated += 1
        print(f"run {row['run_id']}: dispersion={disp:.4g} unbound={unb:.4g}")

    with csv_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)

    print(f"[OK] Updated {updated} rows in {csv_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
