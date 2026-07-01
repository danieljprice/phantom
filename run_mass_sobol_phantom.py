#!/usr/bin/env python3
"""Sobol parameter sweep runner for PHANTOM solarsystem-style setups (mass + optional scales).

Runs PHANTOM once per sample; use ``--jobs N`` to run independent cases in parallel processes."""

from __future__ import annotations

import argparse
import concurrent.futures
from concurrent.futures import as_completed
import csv
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, NamedTuple, Optional, Sequence, Tuple

import numpy as np
from scipy.spatial import cKDTree


AU_IN_KM = 149_597_870.7
SPIN_MEAN_MAX_DUMPS = 50
MAIN_BODY_INTACT_DISP_RATIO = 1.15
MAIN_BODY_MIN_CLUSTER_FRAC = 0.50
MAIN_BODY_LINK_FACTOR_INIT = 2.5
MAIN_BODY_LINK_FACTOR_MAX = 40.0
SPIN_INTRINSIC_MAX_HOURS = 24.0
SPIN_INTRINSIC_MAX_CA_FRACTION = 0.3
SPIN_APPROACH_HOURS_BEFORE_CA = 24.0
EARTH_SINK_ID_DEFAULT = 4
APOPHIS_SINK_ID_DEFAULT = 11

# Literature Apophis shape assets (Shapes/apophis_v233s7.obj longest axis ~0.409741 km).
REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_APOPHIS_SHAPE_CONFIG = REPO_ROOT / "Shapes" / "apophis.shape"
DEFAULT_APOPHIS_OBJ = REPO_ROOT / "Shapes" / "apophis_v233s7.obj"
LITERATURE_MESH_SCALE_KM = 0.205
HORIZONS_APOPHIS_RADIUS_KM = 0.170
LITERATURE_SCALE_R_APOPHIS = 0.409741 / 2.0 / HORIZONS_APOPHIS_RADIUS_KM

# MAXPTMASS is a compile-time array bound on the number of sink particles (config.F90).
# The default upstream value is 1000. The close-packed lattice placement routine overshoots
# the requested np_apophis by ~2-3% (observed: np_apophis=1000 -> 1012 actual sinks), so
# any run with np_apophis >= ~975 will exceed MAXPTMASS=1000 and phantomsetup will abort:
#   "ERROR: nptmass=1012 exceeds ptmass array dimensions of 1000"
# Fix: compile with MAXPTMASS=2000 (set in the Makefile 'ifndef MAXPTMASS' block).
_MAXPTMASS_DEFAULT = 1000
_LATTICE_OVERSHOOT_FRACTION = 0.025
_MAXPTMASS_WARN_THRESHOLD = int(_MAXPTMASS_DEFAULT * (1.0 - _LATTICE_OVERSHOOT_FRACTION))  # 975

# Matches the PHANTOM error line emitted when nptmass exceeds the compiled array bound.
_MAXPTMASS_ERR_RE = re.compile(
    r"nptmass\s*=\s*(\d+)\s+exceeds\s+ptmass\s+array\s+dimensions\s+of\s+(\d+)",
    re.IGNORECASE,
)
_MAXP_OVERFLOW_RE = re.compile(
    r"STOP ipart > maxp|ipart > maxp",
    re.IGNORECASE,
)

# Auto-generated batch folder suffix (after "{prefix}_{timestamp}_"); keep paths portable.
_BATCH_SLUG_DEFAULT_MAX_LEN = 120
_SLUG_SAFE_RE = re.compile(r"[^A-Za-z0-9_.\-]+")

# Default install root: this directory. Override with PHANTOM_DIR.
_DEFAULT_PHANTOM_DIR = Path(__file__).resolve().parent
_DEFAULT_PREFIX = "apophis"
# Template apophis.setup / apophis.in live in the sibling apophis/ project (Honours2026/apophis).
_DEFAULT_BASE_DIR = _DEFAULT_PHANTOM_DIR.parent / "apophis"

_NP_APOPHIS_RE = re.compile(r"^(\s*np_apophis\s*=\s*)\d+", re.MULTILINE | re.IGNORECASE)


def resolve_phantom_executable(phantom_dir: Path, name: str, *, must_exist: bool) -> Path:
    """Prefer ``phantom_dir/bin/<name>``, then ``phantom_dir/<name>`` (flat install)."""
    bin_path = phantom_dir / "bin" / name
    root_path = phantom_dir / name
    for p in (bin_path, root_path):
        if p.is_file():
            return p
    if not must_exist:
        return bin_path
    raise FileNotFoundError(f"Missing PHANTOM binary {name!r}: tried {bin_path}, {root_path}")


def _read_np_apophis_from_setup(setup_path: Path) -> Optional[int]:
    m = re.search(r"^\s*np_apophis\s*=\s*(\d+)", setup_path.read_text(encoding="utf-8"), re.MULTILINE)
    return int(m.group(1)) if m else None


def _shape_crop_enabled_in_setup(setup_path: Path) -> bool:
    """True when the setup file references a non-empty apophis_shape_file."""
    match = re.search(
        r"^\s*apophis_shape_file\s*=\s*(\S+)",
        setup_path.read_text(encoding="utf-8"),
        re.MULTILINE,
    )
    if not match:
        return False
    ref = match.group(1).strip()
    return bool(ref) and ref.lower() not in ("none", "''", '""', "0", "f", "false")


def _resolve_use_shape_crop(sample: "RunSample", setup_path: Path) -> bool:
    """Use explicit sample flag when set; otherwise infer from template apophis_shape_file."""
    if sample.use_shape_crop is not None:
        return sample.use_shape_crop
    return _shape_crop_enabled_in_setup(setup_path)


def _dem_maxp_for_np(np_apophis: int, use_shape_crop: bool) -> int:
    """Minimum --maxp for phantomsetup's temporary SPH lattice (before DEM conversion).

    The close-packed grid is allocated at full lattice size before shape cropping.
    With apophis_shape_file set, phantomsetup runs a second, slightly larger lattice
    (np=8000 -> 20×24×25=12000 then 21×24×26=13104 sites).
    """
    factor = 1.70 if use_shape_crop else 1.55
    floor = 4000 if use_shape_crop else 2000
    return max(floor, int(np_apophis * factor + 0.5))


def _diagnose_setup_log(log_path: Path) -> Optional[str]:
    """Scan a phantomsetup log and return an actionable message for known failure modes."""
    try:
        text = log_path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None
    m = _MAXPTMASS_ERR_RE.search(text)
    if m:
        actual = int(m.group(1))
        compiled = int(m.group(2))
        needed = int(actual * 1.15) + 1
        return (
            f"phantomsetup: nptmass={actual} exceeds compiled MAXPTMASS={compiled}. "
            f"Rebuild phantom and phantomsetup with MAXPTMASS>={needed}: "
            f"add 'MAXPTMASS={needed}' to the Makefile 'ifndef MAXPTMASS' block, "
            "then run 'make setup && make'."
        )
    if _MAXP_OVERFLOW_RE.search(text):
        m = re.search(r"Total memory allocated to arrays:.*?n\s*=\s*(\d+)", text, re.DOTALL)
        alloc = f" (log shows --maxp={m.group(1)})" if m else ""
        return (
            "phantomsetup: temporary SPH lattice exceeds --maxp (ipart > maxp)"
            f"{alloc}. The runner scales --maxp from np_apophis and apophis_shape_file; "
            "if this persists, raise np_apophis or rebuild with a newer run_mass_sobol_phantom.py."
        )
    if "FATAL ERROR" in text:
        return f"phantomsetup failed; see {log_path}"
    missing = sorted(set(re.findall(r"missing A in (\S+)", text)))
    if missing or "ERRORS during setup" in text:
        if missing:
            preview = ", ".join(missing[:6])
            if len(missing) > 6:
                preview += ", ..."
            return (
                "phantomsetup: missing or invalid JPL Horizons ephemeris .txt files "
                f"({preview}). setup_solarsystem still reads all planet files even when "
                "apophis_only=T. Copy *.txt into each run via --ephemeris-cache-dir "
                "(defaults to --base-dir when it contains *.txt)."
            )
    return None


def _verify_phantomsetup_outputs(
    run_dir: Path, prefix: str, setup_log: Path
) -> None:
    """phantomsetup can return exit code 0 on Fortran STOP; validate log and dump file."""
    diag = _diagnose_setup_log(setup_log)
    if diag:
        raise RuntimeError(diag)
    if not any(run_dir.glob(f"{prefix}_00000*")):
        raise RuntimeError(
            f"phantomsetup produced no initial dump ({prefix}_00000*); see {setup_log}"
        )


def _diagnose_maxptmass_error(log_path: Path) -> Optional[str]:
    return _diagnose_setup_log(log_path)


def _np_apophis_maxptmass_warning(sample: "RunSample") -> Optional[str]:
    """Warn before running when np_apophis is close enough to MAXPTMASS that lattice overshoot
    will push the actual sink count over the compiled limit, causing phantomsetup to abort.
    """
    if sample.np_apophis is None or sample.np_apophis < _MAXPTMASS_WARN_THRESHOLD:
        return None
    return (
        f"np_apophis={sample.np_apophis} is at or above the safe limit ({_MAXPTMASS_WARN_THRESHOLD}) "
        f"for default MAXPTMASS={_MAXPTMASS_DEFAULT}. The close-packed lattice overshoots "
        "np_apophis by ~2-3%, which will cause phantomsetup to abort. "
        "Ensure PHANTOM is compiled with MAXPTMASS>np_apophis "
        "(e.g. set MAXPTMASS=2000 in the Makefile 'ifndef MAXPTMASS' block, then run 'make setup && make')."
    )


# Scale parameters varied via CLI min/max: (RunSample/.setup attribute, argparse lo/hi attrs, batch slug token).
# Order MUST match Sobol dimension ordering used in build_run_samples and CSV columns.
# Spin parameters use the full setup-file key as the RunSample attribute so apply_run_sample_to_setup
# can patch them with a single generic loop (same mechanism as scale_vel etc.).
_SCALE_VARIATION_SPEC: Tuple[Tuple[str, str, str, str], ...] = (
    ("scale_vel", "scale_vel_min", "scale_vel_max", "sv"),
    ("scale_pos", "scale_pos_min", "scale_pos_max", "sp"),
    ("scale_earth_sep", "scale_earth_sep_min", "scale_earth_sep_max", "ses"),
    ("scale_r_apophis", "scale_r_apophis_min", "scale_r_apophis_max", "sra"),
    ("scale_rho", "scale_rho_min", "scale_rho_max", "srho"),
    # Spin parameters (setup_solarsystem.f90): period in seconds; axis as unit vector components.
    ("apophis_spin_period",   "spin_period_min",   "spin_period_max",   "sspd"),
    ("apophis_spin_axis_x",   "spin_axis_x_min",   "spin_axis_x_max",   "sax"),
    ("apophis_spin_axis_y",   "spin_axis_y_min",   "spin_axis_y_max",   "say"),
    ("apophis_spin_axis_z",   "spin_axis_z_min",   "spin_axis_z_max",   "saz"),
)

# DEM contact parameters varied via CLI min/max: patched into the run's .in file AFTER
# phantomsetup generates it (phantomsetup overwrites .in, so setup-time patching is lost).
# (RunSample attribute == .in file key, argparse lo/hi dest stems, batch slug token)
_IN_VARIATION_SPEC: Tuple[Tuple[str, str, str, str], ...] = (
    ("kt_cgs",        "kc_min",    "kc_max",    "kc"),
    ("ct_dem",        "ct_min",    "ct_max",    "ct"),
    ("epsilon_n_dem", "eps_n_min", "eps_n_max", "eps"),
    ("kn_cgs",        "kn_min",    "kn_max",    "kn"),
)

# Simulation timeframe parameters: stored in hours, written as "X hr" strings to the .setup file.
# Fixed override (--tmax-hours / --dtmax-hours) applies to every run without consuming a Sobol dimension.
# Sweep bounds (--tmax-hours-min/max / --dtmax-hours-min/max) make the parameter a Sobol dimension.
# The two modes are mutually exclusive per parameter.
_TIME_VARIATION_SPEC: Tuple[Tuple[str, str, str, str], ...] = (
    ("tmax_hours",  "tmax_hours_min",  "tmax_hours_max",  "tmx"),
    ("dtmax_hours", "dtmax_hours_min", "dtmax_hours_max", "dtmx"),
)


def _active_scale_variations(args: argparse.Namespace) -> List[Tuple[str, float, float, str]]:
    """Scale dimensions that are active (both bounds set), in canonical order."""
    out: List[Tuple[str, float, float, str]] = []
    for param, lo_attr, hi_attr, slug_tok in _SCALE_VARIATION_SPEC:
        lo = getattr(args, lo_attr)
        hi = getattr(args, hi_attr)
        if lo is not None and hi is not None:
            out.append((param, lo, hi, slug_tok))
    return out


def _active_in_variations(args: argparse.Namespace) -> List[Tuple[str, float, float, str]]:
    """In-file DEM contact dimensions that are active (both bounds set), in canonical order."""
    out: List[Tuple[str, float, float, str]] = []
    for param, lo_attr, hi_attr, slug_tok in _IN_VARIATION_SPEC:
        lo = getattr(args, lo_attr, None)
        hi = getattr(args, hi_attr, None)
        if lo is not None and hi is not None:
            out.append((param, lo, hi, slug_tok))
    return out


def _active_time_variations(args: argparse.Namespace) -> List[Tuple[str, float, float, str]]:
    """Timeframe dimensions that are being swept (both min and max set), in canonical order."""
    out: List[Tuple[str, float, float, str]] = []
    for param, lo_attr, hi_attr, slug_tok in _TIME_VARIATION_SPEC:
        lo = getattr(args, lo_attr, None)
        hi = getattr(args, hi_attr, None)
        if lo is not None and hi is not None:
            out.append((param, lo, hi, slug_tok))
    return out


def _mass_bounds_active(args: argparse.Namespace) -> bool:
    """True when mass is a Sobol dimension: both kg bounds set (same contract as optional scale min/max)."""
    return args.mass_min_kg is not None and args.mass_max_kg is not None


# ---------------------------------------------------------------------------
# SciPy multi-dimensional Sobol (preferred) vs independent 1D fallback
# ---------------------------------------------------------------------------
def _sobol_nd_scipy(n: int, dim: int, seed: int) -> List[List[float]]:
    from scipy.stats import qmc

    engine = qmc.Sobol(d=dim, scramble=True, seed=seed)
    raw = engine.random(n=n)
    # scipy returns [0,1); keep as list of rows
    return [list(map(float, row)) for row in raw]


def _sobol_nd_fallback(n: int, dim: int, seed: int) -> List[List[float]]:
    """
    Independent 1D quasi-random sequences per dimension (different seeds).
    Joint space coverage is weaker than true multi-d Sobol; install scipy for better designs.
    """
    cols: List[List[float]] = []
    for d in range(dim):
        cols.append(sobol_1d_samples(n, seed + 1_000_003 * (d + 1)))
    return [[cols[d][i] for d in range(dim)] for i in range(n)]


def sobol_nd_unit_samples(n: int, dim: int, seed: int) -> List[List[float]]:
    if n <= 0:
        return []
    if dim <= 0:
        return [[] for _ in range(n)]
    try:
        return _sobol_nd_scipy(n, dim, seed)
    except (ImportError, ValueError) as exc:
        # SciPy missing or Sobol engine rejected parameters — use weaker independent 1D fallback.
        label = "SciPy not installed" if isinstance(exc, ImportError) else f"Sobol engine: {exc}"
        print(f"[WARN] {label}; using independent 1D quasi-random fallback.", file=sys.stderr)
        return _sobol_nd_fallback(n, dim, seed)


@dataclass
class RunRecord:
    run_id: int
    mass_input_kg: float
    run_dir: str
    status: str
    closest_approach_km: float
    closest_approach_au: float
    error: str
    dispersion_ratio: float = float("nan")
    unbound_fraction: float = float("nan")
    intrinsic_spin_period_hr: float = float("nan")
    approach_spin_period_hr: float = float("nan")
    post_flyby_spin_period_hr: float = float("nan")
    # Values written for CSV secondary columns (setup-related names e.g. scale_vel; mass uses mass_input_kg).
    param_columns: Dict[str, str] = field(default_factory=dict)


@dataclass
class RunSample:
    """Per-run values: None means leave the copied template unchanged for that key."""

    mass_kg: Optional[float] = None
    scale_vel: Optional[float] = None
    scale_pos: Optional[float] = None
    scale_earth_sep: Optional[float] = None
    scale_r_apophis: Optional[float] = None
    scale_rho: Optional[float] = None
    use_dem: Optional[bool] = None
    use_shape_crop: Optional[bool] = None
    apophis_only: Optional[bool] = None
    np_apophis: Optional[int] = None
    # Spin: patched into the setup file; Fortran applies them only when use_dem=T, np_apophis>1.
    apophis_spin_period:  Optional[float] = None  # seconds (0 = no spin)
    apophis_spin_axis_x: Optional[float] = None
    apophis_spin_axis_y: Optional[float] = None
    apophis_spin_axis_z: Optional[float] = None
    # DEM contact params: patched into the .in file after phantomsetup; only active when isink_potential=2.
    kt_cgs:        Optional[float] = None  # DEM tensile spring (g/s^2 per cm gap); 0 = no cohesion
    ct_dem:        Optional[float] = None  # tangential damping coefficient
    epsilon_n_dem: Optional[float] = None  # normal restitution coefficient [0,1]
    kn_cgs:        Optional[float] = None  # normal spring constant (dyne/cm)
    # Timeframe: patched into the .setup file before phantomsetup as "X hr" strings.
    tmax_hours:  Optional[float] = None   # simulation end time in hours
    dtmax_hours: Optional[float] = None   # dump interval in hours


class RunWorkerPayload(NamedTuple):
    """Picklable bundle for ProcessPoolExecutor workers (paths as str)."""

    run_id: int
    sample: RunSample
    base_setup: str
    base_input: str
    output_root: str
    prefix: str
    phantomsetup_bin: str
    phantom_bin: str
    ref_mass_kg: Optional[float]
    dry_run: bool
    earth_sink_id: int
    apophis_sink_id: int
    ephemeris_cache_dir: Optional[str]
    shape_file: Optional[str]
    literature_scale_r_allowed: bool


def _shape_crop_may_be_enabled(args: argparse.Namespace) -> bool:
    return bool(args.vary_use_shape_crop or args.use_shape_crop_fixed == "true")


def resolve_default_shape_file() -> Path:
    """Default shape config for literature-scale Apophis mesh cropping."""
    if DEFAULT_APOPHIS_SHAPE_CONFIG.is_file():
        return DEFAULT_APOPHIS_SHAPE_CONFIG
    if DEFAULT_APOPHIS_OBJ.is_file():
        print(
            "[WARN] Shapes/apophis.shape missing; using bare .obj — PHANTOM forces mesh scale=1.0. "
            "Add apophis.shape with 'mesh apophis_v233s7.obj 0.205' for literature sizing.",
            file=sys.stderr,
        )
        return DEFAULT_APOPHIS_OBJ
    raise FileNotFoundError(
        f"Default Apophis shape assets not found under {REPO_ROOT / 'Shapes'} "
        "(expected apophis.shape and/or apophis_v233s7.obj)"
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Generate Sobol samples over Apophis mass and optional setup_solarsystem.f90 "
            "parameters, then run PHANTOM (one subprocess per sample; use --jobs for parallelism). "
            "Multi-d Sobol uses scipy when installed; otherwise independent 1D sequences per dimension "
            "are used (weaker joint coverage — see pip install -r requirements.txt)."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--prefix", default=_DEFAULT_PREFIX, help="PHANTOM file prefix.")
    parser.add_argument(
        "--base-dir",
        default=str(_DEFAULT_BASE_DIR),
        help="Directory containing <prefix>.in and <prefix>.setup.",
    )
    parser.add_argument(
        "--ephemeris-cache-dir",
        default=None,
        metavar="DIR",
        help=(
            "Copy every *.txt from this directory into each run folder before phantomsetup. "
            "Defaults to --base-dir when that folder contains *.txt. PHANTOM then skips JPL "
            "Horizons download when <object>.txt already exists. Files must match the epoch "
            "in your .setup."
        ),
    )
    parser.add_argument(
        "--num-samples",
        type=int,
        default=8,
        help="Number of Sobol samples (runs) to generate.",
    )
    parser.add_argument(
        "--mass-min-kg",
        type=float,
        default=None,
        help="Lower bound for Apophis mass in kg; omit with --mass-max-kg to leave mass unvaried.",
    )
    parser.add_argument(
        "--mass-max-kg",
        type=float,
        default=None,
        help="Upper bound for Apophis mass in kg; omit with --mass-min-kg to leave mass unvaried.",
    )
    parser.add_argument(
        "--apophis-ref-mass-kg",
        type=float,
        default=None,
        metavar="KG",
        help=(
            "Baseline Apophis mass in kg at scale_rho=1 (i.e. using the template density). "
            "Required when --mass-min-kg and --mass-max-kg are set. "
            "The sampled mass is converted to scale_rho = mass_kg / ref_mass_kg and patched "
            "into the setup file."
        ),
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Scramble seed for Sobol sequence.",
    )
    # Real scalings and spin parameters — optional Sobol dimensions when both min and max are set.
    # Spin only takes effect in DEM runs (use_dem=T, np_apophis>1); setting bounds without DEM
    # is allowed but a warning is issued at runtime.
    for key, helpt in (
        ("scale-vel", "scale_vel (Apophis velocity scale)"),
        ("scale-pos", "scale_pos (Apophis heliocentric position scale; not Earth distance)"),
        ("scale-earth-sep", "scale_earth_sep (Earth–Apophis geocentric separation; 1=ephemeris; needs apophis_only=F)"),
        ("scale-r-apophis", "scale_r_apophis (Apophis radius scale)"),
        ("scale-rho", "scale_rho (bulk density scale when mass is density-derived)"),
        ("spin-period",  "apophis_spin_period in seconds (0 = no spin; DEM runs only)"),
        ("spin-axis-x",  "apophis_spin_axis_x (DEM runs only; PHANTOM normalizes the axis vector)"),
        ("spin-axis-y",  "apophis_spin_axis_y (DEM runs only)"),
        ("spin-axis-z",  "apophis_spin_axis_z (DEM runs only)"),
    ):
        parser.add_argument(
            f"--{key}-min",
            type=float,
            default=None,
            help=f"Lower bound for {helpt}; omit with max to leave unvaried.",
        )
        parser.add_argument(
            f"--{key}-max",
            type=float,
            default=None,
            help=f"Upper bound for {helpt}; omit with min to leave unvaried.",
        )
    # DEM contact parameters — optional Sobol dimensions patched into .in after phantomsetup.
    # Only active when both min and max are set and the run uses DEM (isink_potential=2).
    for key, helpt in (
        ("kc",    "kt_cgs (DEM tensile spring constant, g/s^2 per cm gap; 0=off)"),
        ("ct",    "ct_dem (tangential damping coefficient; DEM only)"),
        ("eps-n", "epsilon_n_dem (normal restitution coefficient [0,1]; DEM only)"),
        ("kn",    "kn_cgs (normal spring constant in dyne/cm; DEM only)"),
    ):
        parser.add_argument(
            f"--{key}-min",
            type=float,
            default=None,
            help=f"Lower bound for {helpt}; omit with max to leave unvaried.",
        )
        parser.add_argument(
            f"--{key}-max",
            type=float,
            default=None,
            help=f"Upper bound for {helpt}; omit with min to leave unvaried.",
        )
    parser.add_argument(
        "--kc-fixed",
        type=float,
        default=None,
        metavar="DYNE_CM",
        help=(
            "Fix kt_cgs (DEM tensile spring constant, g/s^2 per cm gap) for every run. "
            "Patched into each run's .in after phantomsetup. "
            "Mutually exclusive with --kc-min/--kc-max."
        ),
    )
    parser.add_argument(
        "--spin-period-fixed",
        type=float,
        default=None,
        metavar="SECONDS",
        help=(
            "Fix apophis_spin_period (seconds) for every run. "
            "Mutually exclusive with --spin-period-min/--spin-period-max."
        ),
    )
    # Timeframe parameters — fix for all runs or sweep as Sobol dimensions.
    # Fixed (--tmax-hours / --dtmax-hours) and min/max bounds are mutually exclusive per parameter.
    parser.add_argument(
        "--tmax-hours",
        type=float,
        default=None,
        metavar="H",
        help=(
            "Fix tmax_in for every run in this batch (hours). "
            "E.g. --tmax-hours 12 sets tmax_in = '12 hr' in each run's .setup. "
            "Mutually exclusive with --tmax-hours-min/max."
        ),
    )
    parser.add_argument(
        "--tmax-hours-min",
        type=float,
        default=None,
        metavar="H",
        help="Lower bound for tmax_in sweep (hours). Requires --tmax-hours-max. Mutually exclusive with --tmax-hours.",
    )
    parser.add_argument(
        "--tmax-hours-max",
        type=float,
        default=None,
        metavar="H",
        help="Upper bound for tmax_in sweep (hours). Requires --tmax-hours-min.",
    )
    parser.add_argument(
        "--dtmax-hours",
        type=float,
        default=None,
        metavar="H",
        help=(
            "Fix dtmax_in (dump interval) for every run in this batch (hours). "
            "Mutually exclusive with --dtmax-hours-min/max."
        ),
    )
    parser.add_argument(
        "--dtmax-hours-min",
        type=float,
        default=None,
        metavar="H",
        help="Lower bound for dtmax_in sweep (hours). Requires --dtmax-hours-max. Mutually exclusive with --dtmax-hours.",
    )
    parser.add_argument(
        "--dtmax-hours-max",
        type=float,
        default=None,
        metavar="H",
        help="Upper bound for dtmax_in sweep (hours). Requires --dtmax-hours-min.",
    )
    parser.add_argument(
        "--scale-earth-sep-list",
        nargs="+",
        type=float,
        default=None,
        metavar="S",
        help=(
            "Run one simulation per scale_earth_sep value "
            "(e.g. 0.3 0.5 0.7 1.0 1.2 → Roche-limit distance sweep). "
            "Mutually exclusive with --scale-earth-sep-min/max. "
            "Requires apophis_only=F in template (Earth present). "
            "Run count = len(list); use --spin-period-fixed 0 for no spin."
        ),
    )
    parser.add_argument(
        "--scale-earth-sep-fixed",
        type=float,
        default=None,
        metavar="S",
        help="Fix scale_earth_sep for every run. Mutually exclusive with --scale-earth-sep-min/max and --scale-earth-sep-list.",
    )
    parser.add_argument(
        "--vary-use-dem",
        action="store_true",
        help="Sample use_dem (logical) with one Sobol dimension (u>=0.5 -> T).",
    )
    parser.add_argument(
        "--use-dem-fixed",
        choices=("true", "false"),
        default=None,
        metavar="{true,false}",
        help=(
            "Force use_dem to a fixed value in every run. "
            "Mutually exclusive with --vary-use-dem; omit to leave the template value unchanged."
        ),
    )
    parser.add_argument(
        "--np-apophis",
        type=int,
        default=None,
        metavar="N",
        help=(
            "Fix np_apophis to a single value in every run's setup (0=none, 1=sink, N>1=gas/DEM). "
            "Mutually exclusive with --np-apophis-list. Omit to leave the template value unchanged."
        ),
    )
    parser.add_argument(
        "--np-apophis-list",
        nargs="+",
        type=int,
        default=None,
        metavar="N",
        help=(
            "Run one simulation per np_apophis value within a single batch "
            "(e.g. --np-apophis-list 250 500 1000 → run_0001=250, run_0002=500, run_0003=1000). "
            "Mutually exclusive with --np-apophis and --num-samples (run count = len of list). "
            "Other parameters stay at template or fixed values (--use-dem-fixed etc. still apply)."
        ),
    )
    parser.add_argument(
        "--vary-use-shape-crop",
        action="store_true",
        help=(
            "Sample shape cropping with one Sobol dimension (u>=0.5 -> T, writes --shape-file "
            "to apophis_shape_file in setup; else blanks apophis_shape_file to disable cropping)."
        ),
    )
    parser.add_argument(
        "--use-shape-crop-fixed",
        choices=("true", "false"),
        default=None,
        metavar="{true,false}",
        help=(
            "Force shape cropping to a fixed value in every run: 'true' writes --shape-file to "
            "apophis_shape_file; 'false' blanks apophis_shape_file to disable cropping. "
            "Mutually exclusive with --vary-use-shape-crop; omit to leave the template apophis_shape_file unchanged."
        ),
    )
    parser.add_argument(
        "--shape-file",
        default=None,
        metavar="PATH",
        help=(
            "Path to the shape config file (or .obj) staged into each run when shape cropping is "
            "enabled. Default when cropping is on: Shapes/apophis.shape (mesh apophis_v233s7.obj "
            f"at literature scale {LITERATURE_MESH_SCALE_KM} km)."
        ),
    )
    parser.add_argument(
        "--vary-apophis-only",
        action="store_true",
        help=(
            "Sample apophis_only with one Sobol dimension. When True for a run, Earth is not in "
            "the setup and closest-approach extraction is skipped (NaN). See also --sink-*-id."
        ),
    )
    parser.add_argument(
        "--apophis-only-fixed",
        choices=("true", "false"),
        default=None,
        metavar="{true,false}",
        help=(
            "Force apophis_only to a fixed value in every run. "
            "Mutually exclusive with --vary-apophis-only. "
            "When true, auto-sets --sink-apophis-id 1 (override explicitly if needed)."
        ),
    )
    parser.add_argument(
        "--sink-earth-id",
        type=int,
        default=EARTH_SINK_ID_DEFAULT,
        help="Sink index for Earth in .ev filenames (default matches full solar system).",
    )
    parser.add_argument(
        "--sink-apophis-id",
        type=int,
        default=APOPHIS_SINK_ID_DEFAULT,
        help="Sink index for Apophis in .ev filenames (default matches full solar system).",
    )
    parser.add_argument(
        "--output-root",
        default="sobol_mass_runs",
        help="Directory for run folders and summary outputs.",
    )
    parser.add_argument(
        "--batch-label",
        default=None,
        metavar="TEXT",
        help=(
            "If set, use this string as the sweep suffix (sanitized) instead of the auto slug "
            "built from varied parameters and bounds. Batch folder: <prefix>_<timestamp>_<suffix>."
        ),
    )
    parser.add_argument(
        "--batch-slug-max-len",
        type=int,
        default=_BATCH_SLUG_DEFAULT_MAX_LEN,
        metavar="N",
        help=(
            "Maximum length of the sweep suffix for auto slugs (and for --batch-label after "
            "sanitization); longer values are truncated with a short hash appended."
        ),
    )
    parser.add_argument(
        "--phantom-dir",
        default=os.environ.get("PHANTOM_DIR", str(_DEFAULT_PHANTOM_DIR)),
        help=(
            "PHANTOM install root: resolves each binary as <root>/bin/<name> first, "
            "then <root>/<name> (upstream layout vs flat copy)."
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Prepare run directories and setup files without executing PHANTOM.",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=1,
        metavar="N",
        help=(
            "Number of parallel worker processes for independent PHANTOM runs. "
            "If PHANTOM is built with OpenMP, set OMP_NUM_THREADS=1 (or ensure jobs×threads "
            "does not oversubscribe CPU cores)."
        ),
    )
    parser.add_argument(
        "--no-cleanup",
        action="store_true",
        dest="no_cleanup",
        help=(
            "Keep heavy output files (binary dumps, .ev files, phantom.log) after each run "
            "completes. By default these are deleted once metrics are extracted to save disk space. "
            "Use this flag when you need the files for post-processing or visualisation."
        ),
    )
    parser.add_argument(
        "--saltelli-n",
        type=int,
        default=None,
        metavar="N",
        help=(
            "Use SALib Saltelli/Sobol sample design with base size N (ignores --num-samples for layout). "
            "Approximate PHANTOM runs: N*(D+2), or N*(2*D+2) with --saltelli-calc-second-order (D=varying dims)."
        ),
    )
    parser.add_argument(
        "--saltelli-calc-second-order",
        action="store_true",
        help="Also estimate second-order Sobol indices (requires more runs). For use with Analysis.py.",
    )
    parser.add_argument(
        "--interactive",
        "-i",
        action="store_true",
        help="Prompt for each parameter (then run); CLI flags set defaults shown at each prompt.",
    )
    return parser


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    return build_parser().parse_args(argv)


def _strip_interactive_flags(argv: Sequence[str]) -> Tuple[List[str], bool]:
    out: List[str] = []
    interactive = False
    for tok in argv:
        if tok in ("--interactive", "-i"):
            interactive = True
            continue
        out.append(tok)
    return out, interactive


def validate_args(args: argparse.Namespace) -> None:
    if args.jobs < 1:
        raise ValueError("jobs must be >= 1")

    np_list = getattr(args, "np_apophis_list", None)
    sep_list = getattr(args, "scale_earth_sep_list", None)

    # np_apophis_list / scale_earth_sep_list drive run count directly.
    if np_list is None and sep_list is None:
        if args.saltelli_n is not None:
            if args.saltelli_n < 2:
                raise ValueError("saltelli-n must be >= 2")
        elif args.num_samples < 1:
            raise ValueError("num-samples must be >= 1")

    if args.batch_slug_max_len < 17:
        raise ValueError("batch-slug-max-len must be >= 17 (timestamp prefix is 15 chars + separators)")
    m_lo, m_hi = args.mass_min_kg, args.mass_max_kg
    if (m_lo is None) ^ (m_hi is None):
        raise ValueError("mass bounds: set both --mass-min-kg and --mass-max-kg, or neither")
    if _mass_bounds_active(args):
        if m_lo <= 0 or m_hi <= 0 or m_hi <= m_lo:
            raise ValueError("mass bounds must satisfy 0 < mass-min-kg < mass-max-kg")
        if args.apophis_ref_mass_kg is None:
            raise ValueError("--apophis-ref-mass-kg is required when --mass-min-kg and --mass-max-kg are set")
        if args.apophis_ref_mass_kg <= 0:
            raise ValueError("--apophis-ref-mass-kg must be > 0")
        active_scale = _active_scale_variations(args)
        if any(p == "scale_rho" for p, *_ in active_scale):
            raise ValueError(
                "--mass-min-kg/max-kg and --scale-rho-min/max cannot both be active: "
                "both patch scale_rho in the setup file"
            )
    for param, lo_attr, hi_attr, _ in _SCALE_VARIATION_SPEC:
        lo = getattr(args, lo_attr)
        hi = getattr(args, hi_attr)
        if (lo is None) ^ (hi is None):
            raise ValueError(f"{param}: set both min and max, or neither")
        if lo is not None and hi is not None and hi <= lo:
            raise ValueError(f"{param}: require min < max")
    for param, lo_attr, hi_attr, _ in _IN_VARIATION_SPEC:
        lo = getattr(args, lo_attr, None)
        hi = getattr(args, hi_attr, None)
        if (lo is None) ^ (hi is None):
            raise ValueError(f"{param}: set both min and max, or neither")
        if lo is not None and hi is not None and hi <= lo:
            raise ValueError(f"{param}: require min < max")
    if getattr(args, "kc_fixed", None) is not None:
        if args.kc_min is not None or args.kc_max is not None:
            raise ValueError("--kc-fixed is mutually exclusive with --kc-min and --kc-max")
        if args.kc_fixed < 0:
            raise ValueError("--kc-fixed must be >= 0")
    if getattr(args, "spin_period_fixed", None) is not None:
        if args.spin_period_min is not None or args.spin_period_max is not None:
            raise ValueError(
                "--spin-period-fixed is mutually exclusive with --spin-period-min and --spin-period-max"
            )
        if args.spin_period_fixed < 0:
            raise ValueError("--spin-period-fixed must be >= 0 (0 = no spin)")
    for param, lo_attr, hi_attr, _ in _TIME_VARIATION_SPEC:
        fixed_val = getattr(args, param, None)
        lo = getattr(args, lo_attr, None)
        hi = getattr(args, hi_attr, None)
        if fixed_val is not None and (lo is not None or hi is not None):
            flag = param.replace("_", "-")
            raise ValueError(
                f"{param}: use either fixed (--{flag}) or min/max sweep bounds, not both."
            )
        if (lo is None) ^ (hi is None):
            raise ValueError(f"{param}: set both --{lo_attr.replace('_', '-')} and "
                             f"--{hi_attr.replace('_', '-')}, or neither.")
        if lo is not None and hi is not None and hi <= lo:
            raise ValueError(f"{param}: require min < max ({lo} >= {hi}).")
        if fixed_val is not None and fixed_val <= 0:
            raise ValueError(f"{param}: must be positive (got {fixed_val}).")
        if lo is not None and lo <= 0:
            raise ValueError(f"{param}: min must be positive (got {lo}).")

    if args.vary_use_dem and args.use_dem_fixed is not None:
        raise ValueError("--vary-use-dem and --use-dem-fixed are mutually exclusive")

    if args.vary_apophis_only and getattr(args, "apophis_only_fixed", None) is not None:
        raise ValueError("--vary-apophis-only and --apophis-only-fixed are mutually exclusive")
    if getattr(args, "apophis_only_fixed", None) == "true" and args.sink_apophis_id == APOPHIS_SINK_ID_DEFAULT:
        args.sink_apophis_id = 1
        print("[INFO] --apophis-only-fixed true: auto-set --sink-apophis-id 1 "
              "(no solar-system bodies; Apophis sinks start at 1). "
              "Override with explicit --sink-apophis-id if needed.", file=sys.stderr)
    else:
        template_setup = Path(args.base_dir) / f"{args.prefix}.setup"
        if (
            template_setup.is_file()
            and read_setup_logical(template_setup, "apophis_only")
            and args.sink_apophis_id == APOPHIS_SINK_ID_DEFAULT
            and not args.vary_apophis_only
            and getattr(args, "apophis_only_fixed", None) != "false"
        ):
            args.sink_apophis_id = 1
            print(
                "[INFO] Template apophis_only=T: auto-set --sink-apophis-id 1 "
                "(Apophis DEM sinks are numbered from 1, not 11). "
                "Override with explicit --sink-apophis-id if needed.",
                file=sys.stderr,
            )

    if args.np_apophis is not None and args.np_apophis < 0:
        raise ValueError("--np-apophis must be >= 0")

    # --np-apophis-list and --np-apophis are mutually exclusive.
    if np_list is not None and args.np_apophis is not None:
        raise ValueError("--np-apophis-list and --np-apophis are mutually exclusive")
    if np_list is not None:
        if len(np_list) < 1:
            raise ValueError("--np-apophis-list requires at least one value")
        if any(n < 0 for n in np_list):
            raise ValueError("all --np-apophis-list values must be >= 0")
        # Warn early for values that risk MAXPTMASS overflow after lattice overshoot.
        high = [n for n in np_list if n >= _MAXPTMASS_WARN_THRESHOLD]
        if high:
            print(
                f"[WARN] --np-apophis-list values {high} >= {_MAXPTMASS_WARN_THRESHOLD}: "
                f"lattice overshoot (~2-3%) may exceed default MAXPTMASS={_MAXPTMASS_DEFAULT}. "
                "Ensure PHANTOM is compiled with MAXPTMASS>max(np_apophis_list) "
                "(e.g. MAXPTMASS=2000 in the Makefile).",
                file=sys.stderr,
            )

    if sep_list is not None:
        if args.scale_earth_sep_min is not None or args.scale_earth_sep_max is not None:
            raise ValueError("--scale-earth-sep-list is mutually exclusive with --scale-earth-sep-min/max")
        if getattr(args, "scale_earth_sep_fixed", None) is not None:
            raise ValueError("--scale-earth-sep-list is mutually exclusive with --scale-earth-sep-fixed")
        if len(sep_list) < 1:
            raise ValueError("--scale-earth-sep-list requires at least one value")
        if any(s <= 0 for s in sep_list):
            raise ValueError("all --scale-earth-sep-list values must be > 0")
        template_setup = Path(args.base_dir) / f"{args.prefix}.setup"
        if template_setup.is_file() and read_setup_logical(template_setup, "apophis_only"):
            print(
                "[WARN] Roche / scale_earth_sep sweeps require Earth in the simulation. "
                "Set apophis_only=F in the template or pass --apophis-only-fixed false.",
                file=sys.stderr,
            )

    if getattr(args, "scale_earth_sep_fixed", None) is not None:
        if args.scale_earth_sep_min is not None or args.scale_earth_sep_max is not None:
            raise ValueError("--scale-earth-sep-fixed is mutually exclusive with --scale-earth-sep-min/max")
        if args.scale_earth_sep_fixed <= 0:
            raise ValueError("--scale-earth-sep-fixed must be > 0")

    if args.vary_use_shape_crop and args.use_shape_crop_fixed is not None:
        raise ValueError("--vary-use-shape-crop and --use-shape-crop-fixed are mutually exclusive")
    if _shape_crop_may_be_enabled(args) and not args.shape_file:
        args.shape_file = str(resolve_default_shape_file())
    if _shape_crop_may_be_enabled(args):
        if not Path(args.shape_file).is_file():
            raise ValueError(f"--shape-file not found: {args.shape_file}")

    # Warn if spin bounds are active but DEM is not enabled; spin only affects DEM particles.
    _spin_lo_dests = frozenset({
        "spin_period_min",
        "spin_axis_x_min",
        "spin_axis_y_min",
        "spin_axis_z_min",
    })
    if any(getattr(args, d, None) is not None for d in _spin_lo_dests):
        dem_on = args.vary_use_dem or getattr(args, "use_dem_fixed", None) == "true"
        if not dem_on:
            print(
                "[WARN] Spin bounds are set but DEM is not enabled "
                "(--vary-use-dem or --use-dem-fixed true). "
                "Spin velocity kicks only apply when use_dem=T and np_apophis>1 in the setup file.",
                file=sys.stderr,
            )

    # np_apophis_list itself constitutes variation; allow dim=0 in that mode.
    dim = count_dimensions(args)
    if dim == 0 and not np_list and not sep_list:
        raise ValueError(
            "No varying dimensions: set both --mass-min-kg and --mass-max-kg and/or pass scale */ "
            "vary-use-dem / vary-apophis-only / --np-apophis-list / --scale-earth-sep-list."
        )


def count_dimensions(args: argparse.Namespace) -> int:
    n = 0
    if _mass_bounds_active(args):
        n += 1
    n += len(_active_scale_variations(args))
    n += len(_active_in_variations(args))
    n += len(_active_time_variations(args))
    if args.vary_use_dem:
        n += 1
    if args.vary_use_shape_crop:
        n += 1
    if args.vary_apophis_only:
        n += 1
    return n


def sobol_1d_samples(n: int, seed: int) -> List[float]:
    """Generate 1D scrambled Sobol-like samples in [0, 1)."""
    if n <= 0:
        return []
    max_bits = max(1, math.ceil(math.log2(n + 1)))
    direction = [1 << (32 - i) for i in range(1, max_bits + 1)]

    scramble = (seed * 2654435761) & 0xFFFFFFFF
    x = scramble
    out: List[float] = []
    for i in range(1, n + 1):
        c = (i & -i).bit_length() - 1
        x ^= direction[c]
        out.append((x & 0xFFFFFFFF) / 2**32)
    return out



def format_logical_token(val: bool) -> str:
    """Right-padded T/F for list-directed Fortran logical reads in .setup files."""
    ch = "T" if val else "F"
    return f"{ch:>10}"


def format_real_token(val: float) -> str:
    return f"{val:.10g}"


def hours_to_phantom_time_string(h: float) -> str:
    """Format hours as a PHANTOM-readable time token (e.g. 12.5 -> '12.5 hr')."""
    return f"{h:.6g} hr"


def replace_setup_assignment(setup_text: str, key: str, value_str: str) -> str:
    pattern = re.compile(rf"^(\s*{re.escape(key)}\s*=\s*)([^!]*)(!.*)?$", re.MULTILINE)
    matches = list(pattern.finditer(setup_text))
    if not matches:
        raise RuntimeError(f"Setup key {key!r} not found in .setup file")
    if len(matches) > 1:
        raise RuntimeError(
            f"Setup key {key!r} appears on {len(matches)} lines; expected exactly one assignment."
        )
    match = matches[0]
    comment = match.group(3) if match.group(3) is not None else ""
    updated_line = f"{match.group(1)}{value_str}"
    if comment:
        updated_line += f" {comment.strip()}"
    return pattern.sub(updated_line, setup_text, count=1)


def validate_assignment(setup_text: str, key: str, expected: str) -> None:
    pattern = re.compile(rf"^\s*{re.escape(key)}\s*=\s*([^!]*?)\s*(?:!.*)?$", re.MULTILINE)
    match = pattern.search(setup_text)
    if not match:
        raise RuntimeError(f"{key} missing after setup update")
    assigned = match.group(1).strip()
    if assigned != expected.strip():
        raise RuntimeError(f"{key} mismatch after setup update (expected {expected!r}, got {assigned!r})")


def apply_run_sample_to_setup(
    setup_path: Path,
    sample: RunSample,
    ref_mass_kg: Optional[float] = None,
    shape_setup_value: Optional[str] = None,
) -> Dict[str, str]:
    """Patch setup file; returns map of column_name -> value string for CSV."""
    text = setup_path.read_text(encoding="utf-8")
    columns: Dict[str, str] = {}

    if sample.mass_kg is not None:
        if ref_mass_kg is None:
            raise RuntimeError("ref_mass_kg required to convert mass sample to scale_rho")
        scale_rho_val = sample.mass_kg / ref_mass_kg
        tok = format_real_token(scale_rho_val)
        text = replace_setup_assignment(text, "scale_rho", tok)
        validate_assignment(text, "scale_rho", tok)
    for param, _, _, _ in _SCALE_VARIATION_SPEC:
        v = getattr(sample, param)
        if v is not None:
            tok = format_real_token(v)
            text = replace_setup_assignment(text, param, tok)
            validate_assignment(text, param, tok)
            columns[param] = f"{v:.12g}"

    if sample.use_dem is not None:
        tok = format_logical_token(sample.use_dem)
        text = replace_setup_assignment(text, "use_dem", tok)
        validate_assignment(text, "use_dem", tok)
        columns["use_dem"] = "T" if sample.use_dem else "F"

    if sample.use_shape_crop is not None:
        path_tok = (shape_setup_value or "") if sample.use_shape_crop else ""
        text = replace_setup_assignment(text, "apophis_shape_file", path_tok)
        validate_assignment(text, "apophis_shape_file", path_tok)
        columns["use_shape_crop"] = "T" if sample.use_shape_crop else "F"

    if sample.apophis_only is not None:
        tok = format_logical_token(sample.apophis_only)
        text = replace_setup_assignment(text, "apophis_only", tok)
        validate_assignment(text, "apophis_only", tok)
        columns["apophis_only"] = "T" if sample.apophis_only else "F"

    if sample.np_apophis is not None:
        if not _NP_APOPHIS_RE.search(text):
            raise RuntimeError("np_apophis assignment not found in .setup file")
        text = _NP_APOPHIS_RE.sub(lambda m: m.group(1) + str(sample.np_apophis), text, count=1)
        columns["np_apophis"] = str(sample.np_apophis)

    # Timeframe parameters: stored as hours, written as "X hr" strings to the .setup file.
    _SETUP_TIME_KEYS = {"tmax_hours": "tmax_in", "dtmax_hours": "dtmax_in"}
    for param, setup_key in _SETUP_TIME_KEYS.items():
        h = getattr(sample, param)
        if h is not None:
            tok = hours_to_phantom_time_string(h)
            text = replace_setup_assignment(text, setup_key, tok)
            validate_assignment(text, setup_key, tok)
            columns[param] = f"{h:.12g}"

    setup_path.write_text(text, encoding="utf-8")
    return columns


def apply_run_sample_to_in(in_path: Path, sample: RunSample) -> Dict[str, str]:
    """Patch DEM contact parameters into a run's .in file; returns column map for CSV.

    Only patches when the .in file contains isink_potential = 2 (DEM mode active) AND at
    least one in-file parameter is set on the sample.  No-op otherwise — safe to call
    unconditionally after phantomsetup.
    """
    columns: Dict[str, str] = {}
    active = [(p, getattr(sample, p)) for p, *_ in _IN_VARIATION_SPEC if getattr(sample, p) is not None]
    if not active:
        return columns

    text = in_path.read_text(encoding="utf-8")
    if not re.search(r"^\s*isink_potential\s*=\s*2\b", text, re.MULTILINE):
        return columns  # not DEM mode; DEM keys absent from .in

    for param, v in active:
        tok = format_real_token(v)
        try:
            text = replace_setup_assignment(text, param, tok)
            validate_assignment(text, param, tok)
        except RuntimeError as exc:
            raise RuntimeError(
                f"apply_run_sample_to_in: {exc}. "
                "Ensure the binary was compiled from current source (make setup && make) "
                "so that DEM contact keys appear in the generated .in file."
            ) from None
        columns[param] = f"{v:.12g}"

    in_path.write_text(text, encoding="utf-8")
    return columns


def _parse_setup_roche_metrics(setup_log: Path) -> Dict[str, str]:
    """Read pericentre, geocentric separation, and fluid Roche radius from phantomsetup log."""
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
    m_vrel = re.search(
        r"Apophis-Earth relative velocity\s*=\s*([\d.E+-]+)\s*km/s",
        text,
        re.I,
    )
    if m_vrel:
        out["relative_velocity_kms"] = f"{float(m_vrel.group(1)):.12g}"
    m2 = re.search(
        r"r_tidal is\s+[\d.E+-]+\s+au,\s*([\d.E+-]+)\s*km",
        text,
        re.I,
    )
    if m2:
        r_tidal_km = float(m2.group(1))
        out["r_tidal_km"] = f"{r_tidal_km:.12g}"
        if r_tidal_km > 0:
            if "pericentre_distance_km" in out:
                peri = float(out["pericentre_distance_km"])
                out["pericentre_over_r_tidal"] = f"{peri / r_tidal_km:.12g}"
            if "geocentric_separation_km" in out:
                sep = float(out["geocentric_separation_km"])
                out["separation_over_r_tidal"] = f"{sep / r_tidal_km:.12g}"
    return out


def _apply_list_mode_fixed_flags(s: RunSample, args: argparse.Namespace) -> None:
    if args.use_dem_fixed is not None:
        s.use_dem = args.use_dem_fixed == "true"
    if args.use_shape_crop_fixed is not None:
        s.use_shape_crop = args.use_shape_crop_fixed == "true"
    if getattr(args, "apophis_only_fixed", None) is not None:
        s.apophis_only = args.apophis_only_fixed == "true"
    if args.np_apophis is not None:
        s.np_apophis = args.np_apophis
    if getattr(args, "scale_earth_sep_fixed", None) is not None:
        s.scale_earth_sep = args.scale_earth_sep_fixed
    _apply_fixed_run_sample_overrides(s, args)


def _apply_fixed_run_sample_overrides(s: RunSample, args: argparse.Namespace) -> None:
    """Apply CLI fixed overrides that are not Sobol dimensions."""
    if getattr(args, "tmax_hours", None) is not None:
        s.tmax_hours = args.tmax_hours
    if getattr(args, "dtmax_hours", None) is not None:
        s.dtmax_hours = args.dtmax_hours
    if getattr(args, "kc_fixed", None) is not None:
        s.kt_cgs = args.kc_fixed
    if getattr(args, "spin_period_fixed", None) is not None:
        s.apophis_spin_period = args.spin_period_fixed
    if getattr(args, "scale_earth_sep_fixed", None) is not None:
        s.scale_earth_sep = args.scale_earth_sep_fixed


def build_scale_earth_sep_list_samples(args: argparse.Namespace) -> List[RunSample]:
    """One RunSample per --scale-earth-sep-list value (Roche distance sweep)."""
    out: List[RunSample] = []
    for sep in args.scale_earth_sep_list:
        s = RunSample(scale_earth_sep=sep)
        _apply_list_mode_fixed_flags(s, args)
        out.append(s)
    return out


def build_np_list_samples(args: argparse.Namespace) -> List[RunSample]:
    """One RunSample per value in --np-apophis-list; all other params at template or fixed values.

    This lets a single batch produce run_0001=np1, run_0002=np2, ... without the Sobol dimension
    requirement, replacing the previous approach of spawning one subprocess per count.
    """
    out: List[RunSample] = []
    for n in args.np_apophis_list:
        s = RunSample(np_apophis=n)
        _apply_list_mode_fixed_flags(s, args)
        out.append(s)
    return out


def build_run_samples(num_samples: int, args: argparse.Namespace) -> List[RunSample]:
    dim = count_dimensions(args)
    unit_rows = sobol_nd_unit_samples(num_samples, dim, args.seed)
    out: List[RunSample] = []
    for row in unit_rows:
        di = 0
        s = RunSample()
        if _mass_bounds_active(args):
            s.mass_kg = args.mass_min_kg + row[di] * (args.mass_max_kg - args.mass_min_kg)
            di += 1
        for param, lo, hi, _ in _active_scale_variations(args):
            setattr(s, param, lo + row[di] * (hi - lo))
            di += 1
        for param, lo, hi, _ in _active_in_variations(args):
            setattr(s, param, lo + row[di] * (hi - lo))
            di += 1
        for param, lo, hi, _ in _active_time_variations(args):
            setattr(s, param, lo + row[di] * (hi - lo))
            di += 1
        _apply_fixed_run_sample_overrides(s, args)
        if args.vary_use_dem:
            s.use_dem = row[di] >= 0.5
            di += 1
        elif args.use_dem_fixed is not None:
            s.use_dem = args.use_dem_fixed == "true"
        if args.vary_use_shape_crop:
            s.use_shape_crop = row[di] >= 0.5
            di += 1
        elif args.use_shape_crop_fixed is not None:
            s.use_shape_crop = args.use_shape_crop_fixed == "true"
        if args.vary_apophis_only:
            s.apophis_only = row[di] >= 0.5
            di += 1
        elif getattr(args, "apophis_only_fixed", None) is not None:
            s.apophis_only = args.apophis_only_fixed == "true"
        if args.np_apophis is not None:
            s.np_apophis = args.np_apophis
        if di != dim:
            raise RuntimeError("internal error: dimension index mismatch")
        out.append(s)
    return out


def sample_column_order(args: argparse.Namespace) -> List[str]:
    """Stable CSV column order for varied parameters."""
    order: List[str] = []
    if _mass_bounds_active(args):
        order.append("mass_input_kg")
    for param, _, _, _ in _active_scale_variations(args):
        order.append(param)
    for param, _, _, _ in _active_in_variations(args):
        order.append(param)
    for param, _, _, _ in _active_time_variations(args):
        order.append(param)
    if args.vary_use_dem:
        order.append("use_dem")
    if args.vary_use_shape_crop:
        order.append("use_shape_crop")
    if args.vary_apophis_only:
        order.append("apophis_only")
    # np_apophis_list mode: particle count is the varying quantity, record it as a column.
    if getattr(args, "np_apophis_list", None):
        order.append("np_apophis")
    if getattr(args, "scale_earth_sep_list", None):
        order.append("scale_earth_sep")
    if (
        getattr(args, "scale_earth_sep_list", None)
        or args.scale_earth_sep_min is not None
        or getattr(args, "scale_earth_sep_fixed", None) is not None
    ):
        for key in (
            "pericentre_distance_km",
            "pericentre_distance_re",
            "pericentre_over_r_tidal",
            "relative_velocity_kms",
            "geocentric_separation_km",
            "geocentric_separation_re",
            "r_tidal_km",
            "separation_over_r_tidal",
        ):
            order.append(key)
    return order


def build_salib_problem(args: argparse.Namespace) -> Dict[str, object]:
    """SALib problem dict; parameter order matches Saltelli rows and RunSample mapping."""
    names: List[str] = []
    bounds: List[List[float]] = []
    if _mass_bounds_active(args):
        names.append("mass_input_kg")
        bounds.append([float(args.mass_min_kg), float(args.mass_max_kg)])
    for param, lo, hi, _ in _active_scale_variations(args):
        names.append(param)
        bounds.append([float(lo), float(hi)])
    for param, lo, hi, _ in _active_in_variations(args):
        names.append(param)
        bounds.append([float(lo), float(hi)])
    for param, lo, hi, _ in _active_time_variations(args):
        names.append(param)
        bounds.append([float(lo), float(hi)])
    if args.vary_use_dem:
        names.append("use_dem")
        bounds.append([0.0, 1.0])
    if args.vary_use_shape_crop:
        names.append("use_shape_crop")
        bounds.append([0.0, 1.0])
    if args.vary_apophis_only:
        names.append("apophis_only")
        bounds.append([0.0, 1.0])
    dim = len(names)
    if dim == 0:
        raise ValueError("Saltelli mode needs at least one varying dimension")
    return {"num_vars": dim, "names": names, "bounds": bounds}


def run_sample_from_salib_row(row: Sequence[float], args: argparse.Namespace) -> RunSample:
    """Build RunSample from one SALib Sobol row (physical bounds)."""
    s = RunSample()
    i = 0
    if _mass_bounds_active(args):
        s.mass_kg = float(row[i])
        i += 1
    for param, _, _, _ in _active_scale_variations(args):
        setattr(s, param, float(row[i]))
        i += 1
    for param, _, _, _ in _active_in_variations(args):
        setattr(s, param, float(row[i]))
        i += 1
    for param, _, _, _ in _active_time_variations(args):
        setattr(s, param, float(row[i]))
        i += 1
    _apply_fixed_run_sample_overrides(s, args)
    if args.vary_use_dem:
        s.use_dem = float(row[i]) >= 0.5
        i += 1
    elif args.use_dem_fixed is not None:
        s.use_dem = args.use_dem_fixed == "true"
    if args.vary_use_shape_crop:
        s.use_shape_crop = float(row[i]) >= 0.5
        i += 1
    elif args.use_shape_crop_fixed is not None:
        s.use_shape_crop = args.use_shape_crop_fixed == "true"
    if args.vary_apophis_only:
        s.apophis_only = float(row[i]) >= 0.5
        i += 1
    elif getattr(args, "apophis_only_fixed", None) is not None:
        s.apophis_only = args.apophis_only_fixed == "true"
    if args.np_apophis is not None:
        s.np_apophis = args.np_apophis
    if i != len(row):
        raise RuntimeError(f"internal error: saltelli row consumed {i} values but width is {len(row)}")
    return s


def expected_saltelli_num_evals(num_vars: int, base_n: int, calc_second_order: bool) -> int:
    """Match SALib.sample.sobol.sample row count."""
    if calc_second_order:
        return base_n * (2 * num_vars + 2)
    return base_n * (num_vars + 2)


def _fmt_slug_float(x: float) -> str:
    """Format a float for use in batch directory names (letters, digits, ., -, _)."""
    return _SLUG_SAFE_RE.sub("", f"{x:.10g}")


def canonical_sweep_descriptor(args: argparse.Namespace) -> str:
    """Stable string for hashing when the auto slug must be truncated."""
    parts = [
        f"num_samples={args.num_samples}",
        f"seed={args.seed}",
        f"saltelli_n={getattr(args, 'saltelli_n', None)}",
        f"saltelli_second={getattr(args, 'saltelli_calc_second_order', False)}",
        f"mass_active={_mass_bounds_active(args)}",
    ]
    if _mass_bounds_active(args):
        parts.append(f"mass_min_kg={args.mass_min_kg}")
        parts.append(f"mass_max_kg={args.mass_max_kg}")
    for param, lo, hi, _ in _active_scale_variations(args):
        parts.append(f"{param}={lo}:{hi}")
    for param, lo, hi, _ in _active_in_variations(args):
        parts.append(f"{param}={lo}:{hi}")
    for param, lo, hi, _ in _active_time_variations(args):
        parts.append(f"{param}={lo}:{hi}")
    for param, _, _, _ in _TIME_VARIATION_SPEC:
        fixed_val = getattr(args, param, None)
        if fixed_val is not None:
            parts.append(f"{param}_fixed={fixed_val}")
    if getattr(args, "kc_fixed", None) is not None:
        parts.append(f"kc_fixed={args.kc_fixed}")
    if getattr(args, "spin_period_fixed", None) is not None:
        parts.append(f"spin_period_fixed={args.spin_period_fixed}")
    parts.append(f"vary_use_dem={args.vary_use_dem}")
    parts.append(f"vary_use_shape_crop={args.vary_use_shape_crop}")
    parts.append(f"vary_apophis_only={args.vary_apophis_only}")
    parts.append(f"apophis_only_fixed={getattr(args, 'apophis_only_fixed', None)}")
    # np_apophis_list must be included so sweeps that differ only in particle counts
    # produce distinct hashes when the slug is truncated.
    np_list = getattr(args, "np_apophis_list", None)
    parts.append(f"np_apophis_list={np_list}")
    return "|".join(parts)


def build_auto_batch_sweep_slug(args: argparse.Namespace, max_len: int) -> str:
    """Sweep suffix from CLI: sample count, seed, each varied dimension and its bounds."""
    np_list = getattr(args, "np_apophis_list", None)
    sep_list = getattr(args, "scale_earth_sep_list", None)
    if np_list is not None:
        if len(np_list) <= 5:
            tok_np = "np" + "_".join(str(n) for n in np_list)
        else:
            tok_np = f"np{len(np_list)}vals_{min(np_list)}-{max(np_list)}"
        tokens: List[str] = [tok_np]
    elif sep_list is not None:
        if len(sep_list) <= 5:
            tok_sep = "ses" + "_".join(_fmt_slug_float(s) for s in sep_list)
        else:
            tok_sep = f"ses{len(sep_list)}vals_{_fmt_slug_float(min(sep_list))}-{_fmt_slug_float(max(sep_list))}"
        tokens = [tok_sep]
    elif args.saltelli_n is not None:
        tok_n = f"salt{args.saltelli_n}"
        if args.saltelli_calc_second_order:
            tok_n += "s2"
        tokens = [tok_n, f"s{args.seed}"]
    else:
        tokens = [f"n{args.num_samples}", f"s{args.seed}"]
    if _mass_bounds_active(args):
        tokens.append(
            f"m{_fmt_slug_float(args.mass_min_kg)}-{_fmt_slug_float(args.mass_max_kg)}"
        )
    for _, lo, hi, slug_tok in _active_scale_variations(args):
        tokens.append(f"{slug_tok}{_fmt_slug_float(lo)}-{_fmt_slug_float(hi)}")
    for _, lo, hi, slug_tok in _active_in_variations(args):
        tokens.append(f"{slug_tok}{_fmt_slug_float(lo)}-{_fmt_slug_float(hi)}")
    for _, lo, hi, slug_tok in _active_time_variations(args):
        tokens.append(f"{slug_tok}{_fmt_slug_float(lo)}-{_fmt_slug_float(hi)}")
    for param, _, _, slug_tok in _TIME_VARIATION_SPEC:
        fixed_val = getattr(args, param, None)
        if fixed_val is not None:
            tokens.append(f"{slug_tok}{_fmt_slug_float(fixed_val)}")
    if args.vary_use_dem:
        tokens.append("dem")
    if args.vary_use_shape_crop:
        tokens.append("shapecrop")
    if args.vary_apophis_only:
        tokens.append("ao")
    if getattr(args, "apophis_only_fixed", None) == "true":
        tokens.append("aoT")
    slug = "_".join(tokens)
    if len(slug) <= max_len:
        return slug
    h = hashlib.sha256(canonical_sweep_descriptor(args).encode()).hexdigest()[:8]
    keep = max_len - (len(h) + 1)
    if keep < 1:
        return h[:max_len]
    return f"{slug[:keep]}_{h}"


def sanitize_batch_label(label: str) -> str:
    label = label.strip()
    if not label:
        raise ValueError("--batch-label cannot be empty")
    s = _SLUG_SAFE_RE.sub("_", label)
    s = re.sub(r"_+", "_", s).strip("_")
    if not s:
        raise ValueError("--batch-label has no valid characters after sanitization")
    return s


def build_batch_directory_basename(args: argparse.Namespace, timestamp: str) -> str:
    """
    Batch folder basename: <prefix>_<timestamp>_<sweep_suffix>.
    Per-run folders remain run_XXXX inside this directory.
    """
    max_len = int(args.batch_slug_max_len)
    if args.batch_label is not None:
        slug = sanitize_batch_label(args.batch_label)
        if len(slug) > max_len:
            h = hashlib.sha256(slug.encode()).hexdigest()[:8]
            keep = max(1, max_len - (len(h) + 1))
            slug = f"{slug[:keep]}_{h}"
    else:
        slug = build_auto_batch_sweep_slug(args, max_len=max_len)
    return f"{args.prefix}_{timestamp}_{slug}"


def run_command(cmd: Sequence[str], cwd: Path, log_path: Path) -> None:
    with log_path.open("w", encoding="utf-8") as log_file:
        proc = subprocess.run(cmd, cwd=str(cwd), stdout=log_file, stderr=subprocess.STDOUT)
    if proc.returncode != 0:
        raise RuntimeError(f"Command failed ({proc.returncode}): {' '.join(cmd)}")


def resolve_ephemeris_cache_dir(
    base_dir: Path, explicit: Optional[str]
) -> Optional[Path]:
    """Return ephemeris cache directory: explicit path, else base_dir if it has *.txt."""
    if explicit:
        return Path(explicit).expanduser().resolve()
    if any(base_dir.glob("*.txt")):
        return base_dir
    return None


def copy_ephemeris_txt_cache(cache_dir: Path, run_dir: Path) -> int:
    """Copy *.txt ephemeris snippets from cache_dir into run_dir (non-recursive). Returns file count."""
    n = 0
    for src in sorted(cache_dir.glob("*.txt")):
        if src.is_file():
            shutil.copy2(src, run_dir / src.name)
            n += 1
    return n


def stage_shape_from_setup_template(setup_path: Path, base_dir: Path, run_dir: Path) -> None:
    """Copy apophis_shape_file from base_dir when the setup references a local relative path."""
    match = re.search(
        r"^\s*apophis_shape_file\s*=\s*(\S+)",
        setup_path.read_text(encoding="utf-8"),
        re.MULTILINE,
    )
    if not match:
        return
    ref = match.group(1).strip()
    if not ref:
        return
    shape_path = Path(ref)
    if shape_path.is_absolute():
        candidate = shape_path
    else:
        candidate = (base_dir / shape_path).resolve()
    if candidate.is_file() and not (run_dir / candidate.name).is_file():
        stage_shape_assets(candidate, run_dir)


def stage_shape_assets(shape_file: Path, run_dir: Path) -> str:
    """Copy shape config and referenced mesh into run_dir; return basename for apophis_shape_file."""
    src = shape_file.resolve()
    if not src.is_file():
        raise FileNotFoundError(f"shape file not found: {src}")

    if src.suffix.lower() == ".obj":
        shutil.copy2(src, run_dir / src.name)
        print(
            "[WARN] Shape crop uses bare .obj; PHANTOM applies mesh scale=1.0 (not literature sizing). "
            "Prefer Shapes/apophis.shape.",
            file=sys.stderr,
        )
        return src.name

    dest_cfg = run_dir / src.name
    shutil.copy2(src, dest_cfg)
    for line in dest_cfg.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        parts = stripped.split()
        if len(parts) >= 2 and parts[0].lower() == "mesh":
            obj_ref = parts[1]
            obj_path = Path(obj_ref)
            if not obj_path.is_absolute():
                obj_path = (src.parent / obj_path).resolve()
            if not obj_path.is_file():
                raise FileNotFoundError(
                    f"mesh OBJ '{obj_ref}' not found (resolved {obj_path})"
                )
            shutil.copy2(obj_path, run_dir / obj_path.name)
            break
    return src.name


def _maybe_apply_literature_scale_r(
    sample: RunSample, *, literature_scale_r_allowed: bool
) -> None:
    """Match lattice r_apophis to literature mesh half-extent when shape crop is on."""
    if not sample.use_shape_crop or not literature_scale_r_allowed:
        return
    if sample.scale_r_apophis is not None:
        return
    sample.scale_r_apophis = LITERATURE_SCALE_R_APOPHIS


def parse_sink_rows(path: Path) -> Tuple[List[float], List[Tuple[float, float, float]]]:
    rows_t: List[float] = []
    rows_xyz: List[Tuple[float, float, float]] = []
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            parts = stripped.split()
            if len(parts) < 4:
                continue
            try:
                t = float(parts[0])
                x = float(parts[1])
                y = float(parts[2])
                z = float(parts[3])
            except ValueError:
                continue
            rows_t.append(t)
            rows_xyz.append((x, y, z))
    if not rows_t:
        raise RuntimeError(f"No numeric sink rows found in {path}")
    return rows_t, rows_xyz


def _sink_times_close(t_a: float, t_b: float) -> bool:
    """Whether two sink row times refer to the same output step (robust float compare)."""
    return math.isclose(t_a, t_b, rel_tol=1e-9, abs_tol=1e-12)


def _pair_sink_rows_by_time(
    t_earth: Sequence[float],
    xyz_earth: Sequence[Tuple[float, float, float]],
    t_apophis: Sequence[float],
    xyz_apophis: Sequence[Tuple[float, float, float]],
) -> List[Tuple[Tuple[float, float, float], Tuple[float, float, float]]]:
    """Pair Earth and Apophis positions at matching simulation times (sorted merge).

    Rows are matched only when ``t`` agrees within ``_sink_times_close``; index order
    is not assumed, so mis-ordered or unequal-length files still pair correctly when
    they share the same time stamps.
    """
    earth = sorted(zip(t_earth, xyz_earth), key=lambda r: r[0])
    apo = sorted(zip(t_apophis, xyz_apophis), key=lambda r: r[0])
    i, j = 0, 0
    out: List[Tuple[Tuple[float, float, float], Tuple[float, float, float]]] = []
    while i < len(earth) and j < len(apo):
        te, ve = earth[i]
        ta, va = apo[j]
        if _sink_times_close(te, ta):
            out.append((ve, va))
            i += 1
            j += 1
        elif te < ta:
            i += 1
        else:
            j += 1
    return out


def _sink_ev_dump_sort_key(path: Path) -> Tuple[int, str]:
    """Sort sink `.ev` dumps by trailing `N<number>` before extension (latest dump last)."""
    m = re.search(r"N(\d+)\.ev$", path.name, re.IGNORECASE)
    if m:
        return (int(m.group(1)), path.name)
    return (-1, path.name)


def _combined_sink_sort_key(path: Path) -> Tuple[int, int]:
    m = re.search(r"(\d+)\.sink$", path.name, re.IGNORECASE)
    return (int(m.group(1)) if m else 0, path.stat().st_size)


# Cap combined .sink scans for closest-approach (12 hr @ 10 min dumps ≈ 72 steps).
_COMBINED_SINK_MAX_TIMESTEPS = 250


def _closest_approach_from_combined_sink(
    run_dir: Path, prefix: str, earth_sink_id: int, apophis_sink_id: int
) -> Tuple[float, float]:
    """Minimum Earth–Apophis separation from combined ``prefixNN.sink`` (km, au)."""
    sink_files = sorted(run_dir.glob(f"{prefix}*.sink"), key=_combined_sink_sort_key)
    if not sink_files:
        return float("nan"), float("nan")

    closest_km = float("inf")
    cur_t_key: Optional[str] = None
    earth_pos: Optional[Tuple[float, float, float]] = None
    apo_particles: List[Tuple[float, float, float, float]] = []
    n_timesteps = 0

    def _flush_timestep() -> None:
        nonlocal closest_km, earth_pos, apo_particles, n_timesteps
        if earth_pos is None or not apo_particles:
            return
        n_timesteps += 1
        masses = [p[3] for p in apo_particles]
        m_tot = float(sum(masses))
        if m_tot <= 0.0:
            return
        com = (
            sum(p[0] * m for p, m in zip(apo_particles, masses)) / m_tot,
            sum(p[1] * m for p, m in zip(apo_particles, masses)) / m_tot,
            sum(p[2] * m for p, m in zip(apo_particles, masses)) / m_tot,
        )
        dx = com[0] - earth_pos[0]
        dy = com[1] - earth_pos[1]
        dz = com[2] - earth_pos[2]
        d = math.sqrt(dx * dx + dy * dy + dz * dz)
        if d < closest_km:
            closest_km = d

    for path in sink_files:
        with path.open("r", encoding="utf-8", errors="replace") as f:
            for line in f:
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    continue
                parts = stripped.split()
                if len(parts) < 19:
                    continue
                try:
                    vals = [float(parts[i]) for i in range(19)]
                except ValueError:
                    continue
                sid = int(vals[18])
                key = f"{vals[0]:.9g}"
                if cur_t_key is not None and key != cur_t_key:
                    _flush_timestep()
                    earth_pos = None
                    apo_particles = []
                cur_t_key = key
                if sid == earth_sink_id:
                    earth_pos = (vals[1], vals[2], vals[3])
                elif sid >= apophis_sink_id:
                    apo_particles.append((vals[1], vals[2], vals[3], vals[4]))
                if n_timesteps >= _COMBINED_SINK_MAX_TIMESTEPS:
                    break
            if n_timesteps >= _COMBINED_SINK_MAX_TIMESTEPS:
                break
    _flush_timestep()

    if not math.isfinite(closest_km):
        return float("nan"), float("nan")
    return closest_km, closest_km / AU_IN_KM


def _closest_approach_time_from_combined_sink(
    run_dir: Path, prefix: str, earth_sink_id: int, apophis_sink_id: int
) -> float:
    """Simulation time at minimum Earth–Apophis separation from combined ``.sink`` file."""
    sink_files = sorted(run_dir.glob(f"{prefix}*.sink"), key=_combined_sink_sort_key)
    if not sink_files:
        return float("nan")

    closest_km = float("inf")
    t_ca = float("nan")
    cur_t = float("nan")
    cur_t_key: Optional[str] = None
    earth_pos: Optional[Tuple[float, float, float]] = None
    apo_particles: List[Tuple[float, float, float, float]] = []
    n_timesteps = 0

    def _flush_timestep() -> None:
        nonlocal closest_km, t_ca, earth_pos, apo_particles, cur_t, n_timesteps
        if earth_pos is None or not apo_particles:
            return
        n_timesteps += 1
        masses = [p[3] for p in apo_particles]
        m_tot = float(sum(masses))
        if m_tot <= 0.0:
            return
        com = (
            sum(p[0] * m for p, m in zip(apo_particles, masses)) / m_tot,
            sum(p[1] * m for p, m in zip(apo_particles, masses)) / m_tot,
            sum(p[2] * m for p, m in zip(apo_particles, masses)) / m_tot,
        )
        dx = com[0] - earth_pos[0]
        dy = com[1] - earth_pos[1]
        dz = com[2] - earth_pos[2]
        d = math.sqrt(dx * dx + dy * dy + dz * dz)
        if d < closest_km:
            closest_km = d
            t_ca = cur_t

    for path in sink_files:
        with path.open("r", encoding="utf-8", errors="replace") as f:
            for line in f:
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    continue
                parts = stripped.split()
                if len(parts) < 19:
                    continue
                try:
                    vals = [float(parts[i]) for i in range(19)]
                except ValueError:
                    continue
                sid = int(vals[18])
                key = f"{vals[0]:.9g}"
                if cur_t_key is not None and key != cur_t_key:
                    _flush_timestep()
                    earth_pos = None
                    apo_particles = []
                cur_t_key = key
                cur_t = vals[0]
                if sid == earth_sink_id:
                    earth_pos = (vals[1], vals[2], vals[3])
                elif sid >= apophis_sink_id:
                    apo_particles.append((vals[1], vals[2], vals[3], vals[4]))
                if n_timesteps >= _COMBINED_SINK_MAX_TIMESTEPS:
                    break
            if n_timesteps >= _COMBINED_SINK_MAX_TIMESTEPS:
                break
    _flush_timestep()
    return t_ca


def _resolve_closest_approach(
    run_dir: Path,
    prefix: str,
    earth_sink_id: int,
    apophis_sink_id: int,
    param_columns: Dict[str, str],
) -> Tuple[float, float]:
    """Closest approach from sink trajectories, else pericentre from ``setup.log``."""
    try:
        closest_km, closest_au = extract_closest_approach(
            run_dir, prefix, earth_sink_id, apophis_sink_id
        )
        if math.isfinite(closest_km):
            return closest_km, closest_au
    except RuntimeError:
        pass
    for key in ("pericentre_distance_km", "geocentric_separation_km"):
        fallback = param_columns.get(key)
        if fallback:
            closest_km = float(fallback)
            return closest_km, closest_km / AU_IN_KM
    raise RuntimeError(
        "Could not determine Earth–Apophis separation from sink output or setup.log "
        f"(sink IDs {earth_sink_id} and {apophis_sink_id})."
    )


def _min_separation_from_pairs(
    pairs: Sequence[Tuple[Tuple[float, float, float], Tuple[float, float, float]]]
) -> float:
    closest_km = float("inf")
    for ve, va in pairs:
        dx = va[0] - ve[0]
        dy = va[1] - ve[1]
        dz = va[2] - ve[2]
        d = math.sqrt(dx * dx + dy * dy + dz * dz)
        if d < closest_km:
            closest_km = d
    return closest_km


def extract_closest_approach(
    run_dir: Path, prefix: str, earth_sink_id: int, apophis_sink_id: int
) -> Tuple[float, float]:
    earth_candidates = sorted(
        run_dir.glob(f"{prefix}Sink{earth_sink_id:04d}N*.ev"),
        key=_sink_ev_dump_sort_key,
    )
    apophis_candidates = sorted(
        run_dir.glob(f"{prefix}Sink{apophis_sink_id:04d}N*.ev"),
        key=_sink_ev_dump_sort_key,
    )

    if earth_candidates and apophis_candidates:
        t_earth, xyz_earth = parse_sink_rows(earth_candidates[-1])
        t_apophis, xyz_apophis = parse_sink_rows(apophis_candidates[-1])
        pairs = _pair_sink_rows_by_time(t_earth, xyz_earth, t_apophis, xyz_apophis)
        if pairs:
            closest_km = _min_separation_from_pairs(pairs)
            if math.isfinite(closest_km):
                return closest_km, closest_km / AU_IN_KM

    combined_km, combined_au = _closest_approach_from_combined_sink(
        run_dir, prefix, earth_sink_id, apophis_sink_id
    )
    if math.isfinite(combined_km):
        return combined_km, combined_au

    raise RuntimeError(
        "Could not find Earth/Apophis sink trajectories "
        f"(expected per-sink {prefix}Sink{earth_sink_id:04d}N*.ev / "
        f"{prefix}Sink{apophis_sink_id:04d}N*.ev or combined {prefix}*.sink)."
    )


def _closest_approach_time(
    run_dir: Path, prefix: str, earth_sink_id: int, apophis_sink_id: int
) -> float:
    """Simulation time (code units) at minimum Earth–Apophis separation."""
    earth_candidates = sorted(
        run_dir.glob(f"{prefix}Sink{earth_sink_id:04d}N*.ev"),
        key=_sink_ev_dump_sort_key,
    )
    apophis_candidates = sorted(
        run_dir.glob(f"{prefix}Sink{apophis_sink_id:04d}N*.ev"),
        key=_sink_ev_dump_sort_key,
    )
    if earth_candidates and apophis_candidates:
        t_earth, xyz_earth = parse_sink_rows(earth_candidates[-1])
        t_apophis, xyz_apophis = parse_sink_rows(apophis_candidates[-1])
        earth = sorted(zip(t_earth, xyz_earth), key=lambda r: r[0])
        apo = sorted(zip(t_apophis, xyz_apophis), key=lambda r: r[0])
        i, j = 0, 0
        closest_km = float("inf")
        t_ca = float("nan")
        while i < len(earth) and j < len(apo):
            te, ve = earth[i]
            ta, va = apo[j]
            if _sink_times_close(te, ta):
                dx = va[0] - ve[0]
                dy = va[1] - ve[1]
                dz = va[2] - ve[2]
                d = math.sqrt(dx * dx + dy * dy + dz * dz)
                if d < closest_km:
                    closest_km = d
                    t_ca = te
                i += 1
                j += 1
            elif te < ta:
                i += 1
            else:
                j += 1
        if math.isfinite(closest_km) and math.isfinite(t_ca):
            return t_ca

    t_ca = _closest_approach_time_from_combined_sink(
        run_dir, prefix, earth_sink_id, apophis_sink_id
    )
    if math.isfinite(t_ca):
        return t_ca

    raise RuntimeError(
        "Could not find Earth/Apophis sink trajectories "
        f"(expected per-sink {prefix}Sink{earth_sink_id:04d}N*.ev / "
        f"{prefix}Sink{apophis_sink_id:04d}N*.ev or combined {prefix}*.sink)."
    )


def _parse_utime_from_phantom_log(log_path: Path) -> Optional[float]:
    """Seconds per PHANTOM code time unit from ``phantom.log`` (units block)."""
    if not log_path.is_file():
        return None
    text = log_path.read_text(encoding="utf-8", errors="replace")
    m = re.search(r"Time:\s*([\d.Ee+-]+)\s*s", text)
    if not m:
        return None
    try:
        utime = float(m.group(1))
    except ValueError:
        return None
    return utime if utime > 0.0 else None


def _parse_utime_from_setup_log(log_path: Path) -> Optional[float]:
    """Fallback: code time unit (seconds) printed by phantomsetup in ``setup.log``."""
    if not log_path.is_file():
        return None
    text = log_path.read_text(encoding="utf-8", errors="replace")
    m = re.search(r"Time:\s*([\d.Ee+-]+)\s*s", text)
    if not m:
        return None
    try:
        utime = float(m.group(1))
    except ValueError:
        return None
    return utime if utime > 0.0 else None


def _parse_run_utime(run_dir: Path) -> Optional[float]:
    return _parse_utime_from_phantom_log(run_dir / "phantom.log") or _parse_utime_from_setup_log(
        run_dir / "setup.log"
    )


SinkFullRow = Tuple[float, Tuple[float, float, float], Tuple[float, float, float], float]


def parse_sink_full_rows(path: Path) -> List[SinkFullRow]:
    """Parse a sink `.ev` file into (time, position, velocity, mass) rows (PHANTOM code units)."""
    rows: List[SinkFullRow] = []
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            parts = stripped.split()
            if len(parts) < 8:
                continue
            try:
                t = float(parts[0])
                pos = (float(parts[1]), float(parts[2]), float(parts[3]))
                mass = float(parts[4])
                vel = (float(parts[5]), float(parts[6]), float(parts[7]))
            except ValueError:
                continue
            rows.append((t, pos, vel, mass))
    return rows


def _parse_sink_ev_array(path: Path) -> np.ndarray:
    """Parse a sink .ev file to a float64 array of shape (N, 8): [t, x, y, z, mass, vx, vy, vz].

    Uses numpy.loadtxt for C-level text parsing — ~10–30× faster than the pure-Python
    ``parse_sink_full_rows``. Returns an empty (0, 8) array on any error.
    """
    try:
        arr = np.loadtxt(path, comments="#", dtype=np.float64)
    except Exception:
        return np.empty((0, 8), dtype=np.float64)
    if arr.ndim == 1:
        arr = arr[np.newaxis, :]
    if arr.shape[0] == 0 or arr.shape[1] < 8:
        return np.empty((0, 8), dtype=np.float64)
    return np.ascontiguousarray(arr[:, :8])


def _sink_id_from_ev_name(name: str) -> Optional[int]:
    """Sink ID (the 4-digit field) from a `<prefix>Sink<NNNN>N<dd>.ev` filename, or None."""
    m = re.search(r"Sink(\d+)N\d+\.ev$", name, re.IGNORECASE)
    return int(m.group(1)) if m else None


def extract_breakup_metrics(
    run_dir: Path,
    prefix: str,
    apophis_sink_id: int,
    *,
    _groups: Optional[Dict[str, np.ndarray]] = None,
    _time_of_key: Optional[Dict[str, float]] = None,
) -> Tuple[float, float]:
    """Peak dispersion ratio and unbound mass fraction of the Apophis rubble pile.

    Apophis particles are the equal-mass sinks with ID >= ``apophis_sink_id``. Dispersion uses
    the global mass-weighted centre of mass (unchanged). Unbound fraction is peak mass not
    belonging to the main rubble fragment — see ``_unbound_mass_main_body``. Returns ``(nan, nan)``
    when fewer than two Apophis particles are present (non-DEM / single-sink runs).

    ``_groups`` / ``_time_of_key``: pre-built from ``_apophis_time_groups`` to avoid a redundant
    full re-read of all sink files when multiple metrics are computed for the same run.
    """
    if _groups is None or _time_of_key is None:
        _groups, _time_of_key, n_sinks = _apophis_time_groups(run_dir, prefix, apophis_sink_id)
        if n_sinks < 2:
            return float("nan"), float("nan")

    rg_initial: Optional[float] = None
    disp_peak = float("nan")
    unbound_peak = float("nan")
    for key in sorted(_groups, key=lambda k: _time_of_key[k]):
        rg_initial, disp_peak, unbound_peak, _ = _process_dem_dump_frame(
            _groups[key], rg_initial, disp_peak, unbound_peak, compute_spin=False
        )

    return disp_peak, unbound_peak


def _mass_weighted_rg(pos: np.ndarray, mass: np.ndarray) -> float:
    """Mass-weighted radius of gyration for Apophis grains at one dump."""
    m_tot = float(mass.sum())
    if m_tot <= 0.0:
        return float("nan")
    rcom = (pos * mass[:, None]).sum(0) / m_tot
    dr = pos - rcom
    r2 = (dr * dr).sum(1)
    return math.sqrt(float((mass * r2).sum()) / m_tot)


def _process_dem_dump_frame(
    arr: np.ndarray,
    rg_initial: Optional[float],
    disp_peak: float,
    unbound_peak: float,
    *,
    utime: Optional[float] = None,
    spin_axis: Optional[Tuple[float, float, float]] = None,
    compute_spin: bool = True,
) -> Tuple[Optional[float], float, float, float]:
    """Incorporate one dump into breakup peaks and optionally spin period.

    Returns ``(rg_initial, disp_peak, unbound_peak, spin_period_hr)``. Spin is ``nan`` when
    ``utime`` is ``None`` or ``compute_spin`` is False. Main-body masking shares one FoF pass
    for unbound mass and spin when spin is computed.
    """
    spin_period = float("nan")
    if len(arr) < 2:
        return rg_initial, disp_peak, unbound_peak, spin_period
    pos = arr[:, :3]
    mass = arr[:, 3]
    vel = arr[:, 4:7]
    m_tot = float(mass.sum())
    if m_tot <= 0.0:
        return rg_initial, disp_peak, unbound_peak, spin_period

    rg = _mass_weighted_rg(pos, mass)
    if rg_initial is None:
        rg_initial = rg
        rg_ratio = 1.0
    elif rg_initial > 0.0:
        rg_ratio = rg / rg_initial
    else:
        rg_ratio = 1.0

    mb_all = _main_body_mask(pos, mass, rg_ratio=rg_ratio)
    unbound_mass = _unbound_mass_main_body(pos, mass, vel, main_body_mask=mb_all)

    if rg_initial and rg_initial > 0.0:
        ratio = rg / rg_initial
        if math.isnan(disp_peak) or ratio > disp_peak:
            disp_peak = ratio
    frac = unbound_mass / m_tot
    if math.isnan(unbound_peak) or frac > unbound_peak:
        unbound_peak = frac

    if compute_spin and utime is not None and utime > 0.0:
        spin_period = _spin_period_hr_bound_rubble(
            arr, utime, spin_axis=spin_axis, main_body_mask_all=mb_all
        )

    return rg_initial, disp_peak, unbound_peak, spin_period


def _parse_spin_axis_from_setup_log(log_path: Path) -> Optional[Tuple[float, float, float]]:
    """Unit spin axis ``(nx,ny,nz)`` printed by ``apply_apophis_spin`` in ``setup.log``, if present."""
    if not log_path.is_file():
        return None
    text = log_path.read_text(encoding="utf-8", errors="replace")
    m = re.search(
        r"Apophis spin axis\s*(?:\(nx,ny,nz\)\s*=\s*|=)\s*"
        r"([\d.Ee+-]+)\s+([\d.Ee+-]+)\s+([\d.Ee+-]+)",
        text,
    )
    if not m:
        return None
    try:
        axis = (float(m.group(1)), float(m.group(2)), float(m.group(3)))
    except ValueError:
        return None
    mag = math.sqrt(axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2])
    if mag <= 0.0:
        return None
    return (axis[0] / mag, axis[1] / mag, axis[2] / mag)


def _period_hr_from_omega_code(omega_code: float, utime: float) -> float:
    """Convert ``omega`` in code units to spin period (hours); matches ``apply_apophis_spin``."""
    if omega_code <= 0.0 or utime <= 0.0:
        return float("nan")
    return (2.0 * math.pi * utime) / (omega_code * 3600.0)


def _main_body_mask_fof(b_pos: np.ndarray, b_mass: np.ndarray, link_factor: float) -> np.ndarray:
    """Largest friends-of-friends cluster at a fixed link factor."""
    n = len(b_pos)
    if n < 2:
        return np.ones(n, dtype=bool)

    tree = cKDTree(b_pos)
    dists, _ = tree.query(b_pos, k=2)
    nn_dist = dists[:, 1]
    link_length = link_factor * float(np.median(nn_dist))

    parent: List[int] = list(range(n))

    def _find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i, j in tree.query_pairs(link_length):
        ri, rj = _find(int(i)), _find(int(j))
        if ri != rj:
            parent[ri] = rj

    roots = np.array([_find(i) for i in range(n)])
    unique_roots, counts = np.unique(roots, return_counts=True)
    largest_root = unique_roots[np.argmax(counts)]
    return roots == largest_root


def _main_body_mask_fof_adaptive(b_pos: np.ndarray, b_mass: np.ndarray) -> np.ndarray:
    """FoF main body with increasing link length until the largest cluster is substantial."""
    n = len(b_pos)
    if n < 2:
        return np.ones(n, dtype=bool)
    min_count = max(2, int(math.ceil(MAIN_BODY_MIN_CLUSTER_FRAC * n)))
    link_factor = MAIN_BODY_LINK_FACTOR_INIT
    mb = _main_body_mask_fof(b_pos, b_mass, link_factor)
    while int(mb.sum()) < min_count and link_factor < MAIN_BODY_LINK_FACTOR_MAX:
        link_factor = min(link_factor * 2.0, MAIN_BODY_LINK_FACTOR_MAX)
        mb = _main_body_mask_fof(b_pos, b_mass, link_factor)
    return mb


def _main_body_mask(
    b_pos: np.ndarray,
    b_mass: np.ndarray,
    *,
    rg_ratio: Optional[float] = None,
) -> np.ndarray:
    """Boolean mask of the main Apophis rubble fragment at one dump.

    When ``rg_ratio`` (current rg / initial rg) is below ``MAIN_BODY_INTACT_DISP_RATIO``,
    all grains are treated as the main body (cohesive pile with settling gaps only).
    Otherwise adaptive friends-of-friends isolates the largest spatial cluster after breakup.
    """
    n = len(b_pos)
    if n < 2:
        return np.ones(n, dtype=bool)
    if rg_ratio is not None and rg_ratio < MAIN_BODY_INTACT_DISP_RATIO:
        return np.ones(n, dtype=bool)
    return _main_body_mask_fof_adaptive(b_pos, b_mass)


def _unbound_mass_main_body(
    pos: np.ndarray,
    mass: np.ndarray,
    vel: np.ndarray,
    *,
    main_body_mask: Optional[np.ndarray] = None,
) -> float:
    """Mass not on the main Apophis fragment (largest FoF cluster), in code units.

    After breakup into two separated chunks the global-CoM energy test mis-classifies
    secondary-fragment grains; this matches the spin path (``_main_body_mask``). Mass in
    smaller spatial clusters counts as unbound from the main body. Grains still in the main
    cluster but with positive specific energy in that fragment's CoM frame
    (``eps = 0.5 v_rel^2 - M_main/r``, G=1) also count.
    """
    m_tot = float(mass.sum())
    if m_tot <= 0.0:
        return 0.0
    mb = main_body_mask if main_body_mask is not None else _main_body_mask(pos, mass)
    unbound = float(mass[~mb].sum())
    m_main = float(mass[mb].sum())
    if m_main <= 0.0:
        return unbound
    main_pos = pos[mb]
    main_vel = vel[mb]
    main_mass = mass[mb]
    rcom_m = (main_pos * main_mass[:, None]).sum(0) / m_main
    vcom_m = (main_vel * main_mass[:, None]).sum(0) / m_main
    dr_m = main_pos - rcom_m
    dv_m = main_vel - vcom_m
    r2_m = (dr_m * dr_m).sum(1)
    v2_m = (dv_m * dv_m).sum(1)
    r_m = np.sqrt(r2_m)
    safe_r_m = np.where(r_m > 0.0, r_m, 1.0)
    eps_m = np.where(r_m > 0.0, 0.5 * v2_m - m_main / safe_r_m, -1.0)
    unbound += float(main_mass[eps_m > 0.0].sum())
    return unbound


def _spin_period_hr_bound_rubble(
    arr: np.ndarray,
    utime: float,
    spin_axis: Optional[Tuple[float, float, float]] = None,
    *,
    main_body_mask_all: Optional[np.ndarray] = None,
    main_body_mask_bound: Optional[np.ndarray] = None,
) -> float:
    """Spin period (hours) of bound Apophis rubble at one dump; ``nan`` if undefined.

    ``arr``: shape (N, 7) float64 — columns [x, y, z, mass, vx, vy, vz] from ``_apophis_time_groups``.

    Inverts ``apply_apophis_spin`` (``setup_solarsystem.f90``): ``dv_i = omega_vec × (r_i - r_com)``
    with mass-weighted CoM (matches ``get_centreofmass`` / ``reset_centreofmass``). When
    ``spin_axis`` is supplied (from ``setup.log``), uses ``omega = (L·n)/(n^T I n)`` about that
    axis; otherwise uses the least-squares ``omega_vec`` and ``|omega_vec|``.
    """
    if len(arr) < 2 or utime <= 0.0:
        return float("nan")

    pos = arr[:, :3]
    mass = arr[:, 3]
    vel = arr[:, 4:7]
    m_tot = float(mass.sum())
    if m_tot <= 0.0:
        return float("nan")
    rcom = (pos * mass[:, None]).sum(0) / m_tot
    vcom = (vel * mass[:, None]).sum(0) / m_tot
    dr = pos - rcom
    dv = vel - vcom
    r = np.sqrt((dr * dr).sum(1))
    v2 = (dv * dv).sum(1)

    # Bound gate: eps = 0.5*v_rel² - M/r ≤ 0; r=0 → treat as bound.
    safe_r = np.where(r > 0.0, r, 1.0)
    eps = np.where(r > 0.0, 0.5 * v2 - m_tot / safe_r, -1.0)
    bound_mask = eps <= 0.0
    if bound_mask.sum() < 2:
        return float("nan")

    b_pos = pos[bound_mask]
    b_mass = mass[bound_mask]
    b_vel = vel[bound_mask]

    if main_body_mask_bound is not None:
        mb = main_body_mask_bound
    elif main_body_mask_all is not None:
        mb = main_body_mask_all[bound_mask]
    else:
        mb = _main_body_mask(b_pos, b_mass)
    b_pos = b_pos[mb]
    b_mass = b_mass[mb]
    b_vel = b_vel[mb]
    if len(b_pos) < 2:
        return float("nan")

    m_bound = float(b_mass.sum())
    if m_bound <= 0.0:
        return float("nan")
    b_rcom = (b_pos * b_mass[:, None]).sum(0) / m_bound
    b_vcom = (b_vel * b_mass[:, None]).sum(0) / m_bound
    b_dr = b_pos - b_rcom
    b_dv = b_vel - b_vcom

    omega_code: Optional[float] = None
    if spin_axis is not None:
        n = np.array(spin_axis, dtype=np.float64)
        r2 = (b_dr * b_dr).sum(1)
        # Inertia tensor: I = sum_k m_k (r_k² δ_ij - dr_ki dr_kj)
        I_mat = (
            np.einsum("k,k->", b_mass, r2) * np.eye(3)
            - np.einsum("k,ki,kj->ij", b_mass, b_dr, b_dr)
        )
        # Angular momentum: L = sum_k m_k (dr_k × dv_k)
        L = (b_mass[:, None] * np.cross(b_dr, b_dv)).sum(0)
        In = I_mat @ n
        n_I_n = float(n @ In)
        L_n = float(L @ n)
        if abs(n_I_n) > 0.0 and math.isfinite(n_I_n):
            w = L_n / n_I_n
            if math.isfinite(w) and abs(w) > 0.0:
                omega_code = w

    if omega_code is None:
        # Least-squares: for each bound particle the 3 cross-product equations are
        # A_k @ omega_vec = dv_k where A_k is the skew-symmetric matrix of dr_k matching
        # the Fortran rows: [0,-dz,dy], [dz,0,-dx], [-dy,dx,0] (apply_apophis_spin order).
        # AtWA = sum_k m_k * A_k^T A_k; AtWb = sum_k m_k * A_k^T dv_k.
        # A_k^T A_k = r_k² I - dr_k dr_k^T (identical to inertia contribution per particle).
        N = len(b_dr)
        A = np.zeros((N, 3, 3), dtype=np.float64)
        A[:, 0, 1] = -b_dr[:, 2]
        A[:, 0, 2] =  b_dr[:, 1]
        A[:, 1, 0] =  b_dr[:, 2]
        A[:, 1, 2] = -b_dr[:, 0]
        A[:, 2, 0] = -b_dr[:, 1]
        A[:, 2, 1] =  b_dr[:, 0]
        # AtWA[p,q] = sum_k m_k sum_i A[k,i,p]*A[k,i,q]
        AtWA = np.einsum("k,kip,kiq->pq", b_mass, A, A)
        # AtWb[p]   = sum_k m_k sum_i A[k,i,p]*dv[k,i]
        AtWb = np.einsum("k,kip,ki->p", b_mass, A, b_dv)
        try:
            omega_vec = np.linalg.solve(AtWA, AtWb)
        except np.linalg.LinAlgError:
            return float("nan")
        omega_code = float(np.sqrt((omega_vec * omega_vec).sum()))

    return _period_hr_from_omega_code(omega_code, utime)


def _parse_combined_sink_array(path: Path, apophis_sink_id: int) -> np.ndarray:
    """Parse a combined ``prefixNN.sink`` file (``icreate_sinks>0`` DEM runs).

    Returns shape (N, 8): [t, x, y, z, mass, vx, vy, vz] for sink IDs >= apophis_sink_id.
    """
    rows: List[List[float]] = []
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            parts = stripped.split()
            if len(parts) < 19:
                continue
            try:
                vals = [float(parts[i]) for i in range(19)]
            except ValueError:
                continue
            if vals[18] < float(apophis_sink_id):
                continue
            rows.append(vals[:8])
    if not rows:
        return np.empty((0, 8), dtype=np.float64)
    return np.ascontiguousarray(np.array(rows, dtype=np.float64))


def _group_sink_rows_by_time(big: np.ndarray) -> Tuple[Dict[str, np.ndarray], Dict[str, float]]:
    """Group [t,x,y,z,mass,vx,vy,vz] rows by simulation time."""
    if big.shape[0] == 0:
        return {}, {}
    times_big = big[:, 0]
    data = big[:, 1:]
    sort_idx = np.argsort(times_big, kind="stable")
    sorted_times = times_big[sort_idx]
    sorted_data = data[sort_idx]
    boundaries = np.concatenate(
        [[0], np.where(np.diff(sorted_times) != 0)[0] + 1, [len(sorted_times)]]
    )
    groups: Dict[str, np.ndarray] = {}
    time_of_key: Dict[str, float] = {}
    for i in range(len(boundaries) - 1):
        chunk = sorted_data[boundaries[i]: boundaries[i + 1]]
        t = float(sorted_times[boundaries[i]])
        key = f"{t:.9g}"
        groups[key] = chunk
        time_of_key[key] = t
    return groups, time_of_key


def _metrics_legacy_substeps() -> bool:
    """True only when ``METRICS_LEGACY_SUBSTEPS=1`` (debug/regression; not used by sweep runs)."""
    return os.environ.get("METRICS_LEGACY_SUBSTEPS", "").strip().lower() in ("1", "true", "yes")


def _use_fast_metrics_defaults() -> None:
    """Ensure sweep metric extraction uses dump-grid rows and spin window gating.

    Clears ``METRICS_LEGACY_SUBSTEPS`` so a shell export cannot slow every worker.
    Legacy substep scanning is opt-in for ``verify_metric_extraction.py --legacy`` only.
    """
    os.environ.pop("METRICS_LEGACY_SUBSTEPS", None)


def _parse_dtmax_from_in(run_dir: Path, prefix: str) -> Optional[float]:
    """``dtmax`` in PHANTOM code time units from the run ``.in`` file."""
    in_path = run_dir / f"{prefix}.in"
    if not in_path.is_file():
        return None
    text = in_path.read_text(encoding="utf-8", errors="replace")
    m = re.search(r"^\s*dtmax\s*=\s*([\d.Ee+-]+)", text, re.MULTILINE | re.IGNORECASE)
    if not m:
        return None
    try:
        dtmax = float(m.group(1))
    except ValueError:
        return None
    return dtmax if dtmax > 0.0 else None


def _subsample_sink_ev_to_dumps(arr: np.ndarray, dtmax: float) -> np.ndarray:
    """Keep dump-grid rows from one sink ``.ev`` array; fallback to last row per dump bin."""
    if arr.shape[0] == 0:
        return arr
    times = arr[:, 0]
    atol = max(1e-6 * dtmax, 1.0)
    k_grid = np.round(times / dtmax)
    dump_mask = np.isclose(times, k_grid * dtmax, rtol=0.0, atol=atol)
    if int(dump_mask.sum()) >= 3:
        return arr[dump_mask]
    bins = np.floor(times / dtmax).astype(np.int64)
    order = np.argsort(bins, kind="stable")
    sorted_bins = bins[order]
    sorted_arr = arr[order]
    boundaries = np.concatenate(
        [[0], np.where(np.diff(sorted_bins) != 0)[0] + 1, [len(sorted_bins)]]
    )
    last_rows = [sorted_arr[boundaries[i + 1] - 1] for i in range(len(boundaries) - 1)]
    if not last_rows:
        return arr
    return np.vstack(last_rows)


def _filter_complete_sink_groups(
    groups: Dict[str, np.ndarray],
    time_of_key: Dict[str, float],
) -> Tuple[Dict[str, np.ndarray], Dict[str, float], int]:
    """Keep only full snapshots from a combined ``.sink`` history.

    Combined sink files append per-dump rows; most timestamps contain every Apophis grain,
    but occasional duplicate or partial writes exist. The modal row count per time is treated
    as the complete particle count; other snapshots are dropped for breakup metrics.
    """
    if not groups:
        return {}, {}, 0
    size_counts: Dict[int, int] = {}
    for arr in groups.values():
        n = len(arr)
        size_counts[n] = size_counts.get(n, 0) + 1
    mode_size = max(size_counts, key=lambda n: size_counts[n])
    filtered_groups: Dict[str, np.ndarray] = {}
    filtered_times: Dict[str, float] = {}
    for key, arr in groups.items():
        if len(arr) != mode_size:
            continue
        filtered_groups[key] = arr
        filtered_times[key] = time_of_key[key]
    return filtered_groups, filtered_times, mode_size


def _apophis_time_groups_from_combined_sink(
    run_dir: Path, prefix: str, apophis_sink_id: int
) -> Tuple[Dict[str, np.ndarray], Dict[str, float], int]:
    """Fallback for DEM runs that write ``prefixNN.sink`` instead of per-sink ``.ev`` files."""
    sink_files = sorted(run_dir.glob(f"{prefix}*.sink"), key=_combined_sink_sort_key)
    if not sink_files:
        return {}, {}, 0
    # Prefer the highest-numbered / largest file (latest dump cycle).
    latest = sink_files[-1]
    big = _parse_combined_sink_array(latest, apophis_sink_id)
    groups, time_of_key = _group_sink_rows_by_time(big)
    return _filter_complete_sink_groups(groups, time_of_key)


def _apophis_time_groups(
    run_dir: Path, prefix: str, apophis_sink_id: int
) -> Tuple[Dict[str, np.ndarray], Dict[str, float], int]:
    """Group latest per-sink ``.ev`` rows by dump time into numpy arrays.

    Returns ``(groups, time_of_key, n_apophis_sinks)`` where ``groups[key]`` is a float64 array
    of shape ``(N_particles, 7)`` with columns ``[x, y, z, mass, vx, vy, vz]`` (time dropped —
    it is the dict key). By default only PHANTOM **dump-grid** rows are kept (not hydro substeps);
    set ``METRICS_LEGACY_SUBSTEPS=1`` to retain every unique timestamp (regression/debug only —
    sweep runs call ``_use_fast_metrics_defaults()`` first).
    """
    # Collect the latest .ev file per Apophis sink ID; cache sort keys (avoids double regex).
    sort_cache: Dict[Path, Tuple[int, str]] = {}

    def _sk(p: Path) -> Tuple[int, str]:
        if p not in sort_cache:
            sort_cache[p] = _sink_ev_dump_sort_key(p)
        return sort_cache[p]

    latest_by_id: Dict[int, Path] = {}
    for p in run_dir.glob(f"{prefix}Sink*N*.ev"):
        sid = _sink_id_from_ev_name(p.name)
        if sid is None or sid < apophis_sink_id:
            continue
        cur = latest_by_id.get(sid)
        if cur is None or _sk(p) > _sk(cur):
            latest_by_id[sid] = p

    if not latest_by_id:
        return _apophis_time_groups_from_combined_sink(run_dir, prefix, apophis_sink_id)

    # Read files in parallel. np.loadtxt holds the GIL so threads don't help; use processes.
    # Cap workers to avoid oversubscription when run_one_case is itself inside a process pool.
    paths = list(latest_by_id.values())
    max_workers = min(4, max(1, (os.cpu_count() or 4) // 2), len(paths))
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as ex:
        arrays = list(ex.map(_parse_sink_ev_array, paths, chunksize=max(1, len(paths) // (max_workers * 4))))

    dtmax: Optional[float] = None
    if not _metrics_legacy_substeps():
        dtmax = _parse_dtmax_from_in(run_dir, prefix)
        if dtmax is None:
            print(
                "[WARN] Could not read dtmax from .in; using all .ev timestamps for metrics",
                file=sys.stderr,
                flush=True,
            )

    # Deduplicate per-sink (PHANTOM may emit the same dump time twice) then concatenate.
    deduped: List[np.ndarray] = []
    for arr in arrays:
        if arr.shape[0] == 0:
            continue
        times = arr[:, 0]
        _, first_idx = np.unique(times, return_index=True)
        row = arr[np.sort(first_idx)]
        if dtmax is not None:
            row = _subsample_sink_ev_to_dumps(row, dtmax)
        deduped.append(row)

    if not deduped:
        return _apophis_time_groups_from_combined_sink(run_dir, prefix, apophis_sink_id)

    big = np.vstack(deduped)
    groups, time_of_key = _group_sink_rows_by_time(big)
    return groups, time_of_key, len(latest_by_id)


def _hours_to_code_time(hours: float, utime: float) -> float:
    """Convert physical hours to PHANTOM code time units."""
    if utime <= 0.0:
        return float("nan")
    return hours * 3600.0 / utime


def _spin_in_time_window(
    t: float,
    time_relation: str,
    t_ca: Optional[float],
    *,
    utime: Optional[float] = None,
) -> bool:
    if time_relation == "all":
        return True
    if t_ca is None:
        return False
    if time_relation == "after_ca":
        return t >= t_ca or _sink_times_close(t, t_ca)
    if time_relation == "before_ca":
        return t < t_ca and not _sink_times_close(t, t_ca)
    if utime is None or utime <= 0.0:
        return False
    if time_relation == "intrinsic_early":
        t_intrinsic_end = min(
            SPIN_INTRINSIC_MAX_CA_FRACTION * t_ca,
            _hours_to_code_time(SPIN_INTRINSIC_MAX_HOURS, utime),
        )
        t_approach_start = t_ca - _hours_to_code_time(SPIN_APPROACH_HOURS_BEFORE_CA, utime)
        return t < t_intrinsic_end and t < t_approach_start
    if time_relation == "approach_pre_ca":
        t_approach_start = t_ca - _hours_to_code_time(SPIN_APPROACH_HOURS_BEFORE_CA, utime)
        return t_approach_start <= t < t_ca and not _sink_times_close(t, t_ca)
    raise ValueError(f"unknown time_relation: {time_relation!r}")


def _dump_in_any_spin_window(
    t: float,
    *,
    apophis_only: bool,
    t_ca: Optional[float],
    utime: Optional[float],
) -> bool:
    """True when ``t`` can contribute to at least one spin metric (union of spin windows)."""
    if apophis_only:
        return _spin_in_time_window(t, "all", None)
    if t_ca is None or utime is None or utime <= 0.0:
        return False
    return (
        _spin_in_time_window(t, "intrinsic_early", t_ca, utime=utime)
        or _spin_in_time_window(t, "approach_pre_ca", t_ca, utime=utime)
        or _spin_in_time_window(t, "after_ca", t_ca)
    )


def _mean_spin_from_timed_periods(
    timed_periods: List[Tuple[float, float]],
    *,
    skip_earliest_in_window: bool,
) -> float:
    if not timed_periods:
        return float("nan")
    if skip_earliest_in_window and len(timed_periods) >= 3:
        timed_periods = sorted(timed_periods, key=lambda tp: tp[0])[1:]
    if not timed_periods:
        return float("nan")
    if len(timed_periods) > SPIN_MEAN_MAX_DUMPS:
        sorted_tp = sorted(timed_periods, key=lambda tp: tp[0])
        idx = np.linspace(0, len(sorted_tp) - 1, SPIN_MEAN_MAX_DUMPS, dtype=int)
        timed_periods = [sorted_tp[int(i)] for i in idx]
    periods = [p for _, p in timed_periods]
    return math.fsum(periods) / len(periods)


def _extract_dem_metrics_bundle(
    run_dir: Path,
    prefix: str,
    apophis_sink_id: int,
    *,
    apophis_only: bool,
    earth_sink_id: Optional[int] = None,
    t_ca: Optional[float] = None,
    _groups: Dict[str, np.ndarray],
    _time_of_key: Dict[str, float],
    _n_sinks: int,
) -> Tuple[float, float, float, float, float]:
    """Breakup peaks + intrinsic/approach/post spin in one pass over pre-built dump groups."""
    if _n_sinks < 2:
        return float("nan"), float("nan"), float("nan"), float("nan"), float("nan")

    utime = _parse_run_utime(run_dir)
    if utime is None:
        return float("nan"), float("nan"), float("nan"), float("nan"), float("nan")

    spin_axis = _parse_spin_axis_from_setup_log(run_dir / "setup.log")

    if not apophis_only and t_ca is None:
        if earth_sink_id is None:
            return float("nan"), float("nan"), float("nan"), float("nan"), float("nan")
        try:
            t_ca = _closest_approach_time(run_dir, prefix, earth_sink_id, apophis_sink_id)
        except RuntimeError:
            return float("nan"), float("nan"), float("nan"), float("nan"), float("nan")

    rg_initial: Optional[float] = None
    disp_peak = float("nan")
    unbound_peak = float("nan")
    intrinsic_periods: List[Tuple[float, float]] = []
    approach_periods: List[Tuple[float, float]] = []
    post_periods: List[Tuple[float, float]] = []
    spin_enabled = apophis_only or t_ca is not None

    for key in sorted(_groups, key=lambda k: _time_of_key[k]):
        t = _time_of_key[key]
        arr = _groups[key]
        do_spin = spin_enabled and _dump_in_any_spin_window(
            t, apophis_only=apophis_only, t_ca=t_ca, utime=utime
        )
        rg_initial, disp_peak, unbound_peak, period = _process_dem_dump_frame(
            arr,
            rg_initial,
            disp_peak,
            unbound_peak,
            utime=utime,
            spin_axis=spin_axis,
            compute_spin=do_spin,
        )
        if len(arr) < 2 or not math.isfinite(period):
            continue
        if apophis_only:
            if _spin_in_time_window(t, "all", None):
                intrinsic_periods.append((t, period))
        else:
            if _spin_in_time_window(t, "intrinsic_early", t_ca, utime=utime):
                intrinsic_periods.append((t, period))
            if _spin_in_time_window(t, "approach_pre_ca", t_ca, utime=utime):
                approach_periods.append((t, period))
            if _spin_in_time_window(t, "after_ca", t_ca):
                post_periods.append((t, period))

    intrinsic_spin = _mean_spin_from_timed_periods(
        intrinsic_periods, skip_earliest_in_window=True
    )
    approach_spin = float("nan")
    post_spin = float("nan")
    if not apophis_only:
        approach_spin = _mean_spin_from_timed_periods(
            approach_periods, skip_earliest_in_window=False
        )
        post_spin = _mean_spin_from_timed_periods(
            post_periods, skip_earliest_in_window=False
        )
    return disp_peak, unbound_peak, intrinsic_spin, approach_spin, post_spin


def _extract_mean_spin_period_hr(
    run_dir: Path,
    prefix: str,
    apophis_sink_id: int,
    *,
    earth_sink_id: Optional[int] = None,
    time_relation: str,
    skip_earliest_in_window: bool = False,
    t_ca: Optional[float] = None,
    _groups: Optional[Dict[str, np.ndarray]] = None,
    _time_of_key: Optional[Dict[str, float]] = None,
    _n_sinks: Optional[int] = None,
) -> float:
    """Mean bound-rubble spin period (hours) over a dump time window.

    ``time_relation`` is one of:
    - ``intrinsic_early``: early plateau before tidal ramp-up (min of 24 h and 30% of way to CA).
    - ``approach_pre_ca``: last ``SPIN_APPROACH_HOURS_BEFORE_CA`` strictly before closest approach.
    - ``after_ca``: dumps at or after closest approach (post-flyby tidal spin).
    - ``all``: every dump (e.g. ``apophis_only`` runs with no Earth sink).

    When ``skip_earliest_in_window`` is True and at least three dumps fall in the window, the
    earliest-time dump is omitted so the mean is not dominated by immediate post-setup transients.
    The mean uses at most ``SPIN_MEAN_MAX_DUMPS`` evenly spaced dumps in the window.

    ``t_ca``: if supplied, closest-approach time is not re-read from Earth/Apophis centroid files.

    ``_groups`` / ``_time_of_key`` / ``_n_sinks``: pre-built from ``_apophis_time_groups`` to
    avoid a redundant full re-read when multiple metrics are computed for the same run.
    """
    if _groups is None or _time_of_key is None or _n_sinks is None:
        _groups, _time_of_key, _n_sinks = _apophis_time_groups(run_dir, prefix, apophis_sink_id)
    if _n_sinks < 2:
        return float("nan")

    utime = _parse_run_utime(run_dir)
    if utime is None:
        return float("nan")

    spin_axis = _parse_spin_axis_from_setup_log(run_dir / "setup.log")

    if time_relation in ("intrinsic_early", "approach_pre_ca", "after_ca"):
        if t_ca is None:
            if earth_sink_id is None:
                return float("nan")
            try:
                t_ca = _closest_approach_time(run_dir, prefix, earth_sink_id, apophis_sink_id)
            except RuntimeError:
                return float("nan")

    timed_periods: List[Tuple[float, float]] = []
    rg_initial: Optional[float] = None
    for key in sorted(_groups, key=lambda k: _time_of_key[k]):
        t = _time_of_key[key]
        arr = _groups[key]
        if len(arr) < 2:
            continue
        in_window = _spin_in_time_window(t, time_relation, t_ca, utime=utime)
        rg_initial, _, _, period = _process_dem_dump_frame(
            arr,
            rg_initial,
            float("nan"),
            float("nan"),
            utime=utime,
            spin_axis=spin_axis,
            compute_spin=in_window,
        )
        if not in_window:
            continue
        if math.isfinite(period):
            timed_periods.append((t, period))

    return _mean_spin_from_timed_periods(
        timed_periods, skip_earliest_in_window=skip_earliest_in_window
    )


def extract_intrinsic_spin_period_hr(
    run_dir: Path,
    prefix: str,
    apophis_sink_id: int,
    earth_sink_id: int,
    *,
    apophis_only: bool = False,
) -> float:
    """Mean spin period (hours) of settled rubble before tidal ramp-up.

    For full solar-system runs, averages the early plateau (min of 24 h and 30% of way to CA,
    excluding the last 24 h before CA). For ``apophis_only`` runs, averages all dumps with the
    earliest-dump skip when three or more exist. Spin inference matches
    ``extract_post_flyby_spin_period_hr``.
    """
    if apophis_only:
        return _extract_mean_spin_period_hr(
            run_dir,
            prefix,
            apophis_sink_id,
            time_relation="all",
            skip_earliest_in_window=True,
        )
    return _extract_mean_spin_period_hr(
        run_dir,
        prefix,
        apophis_sink_id,
        earth_sink_id=earth_sink_id,
        time_relation="intrinsic_early",
        skip_earliest_in_window=True,
    )


def extract_approach_spin_period_hr(
    run_dir: Path,
    prefix: str,
    apophis_sink_id: int,
    earth_sink_id: int,
) -> float:
    """Mean spin period (hours) over the last 24 h strictly before closest approach."""
    return _extract_mean_spin_period_hr(
        run_dir,
        prefix,
        apophis_sink_id,
        earth_sink_id=earth_sink_id,
        time_relation="approach_pre_ca",
        skip_earliest_in_window=False,
    )


def extract_post_flyby_spin_period_hr(
    run_dir: Path,
    prefix: str,
    apophis_sink_id: int,
    earth_sink_id: int,
) -> float:
    """Mean spin period (hours) of bound Apophis rubble over dumps at or after closest approach.

    Spin is inferred from sink velocities (DEM does not populate ``spinx`` in ``.ev``), inverting
    the same ``omega × dr`` kick as ``set_apophis_spin``. A time-averaged value can differ from
    the setup ``apophis_spin_period`` when DEM rearranges the rubble (``I`` about the spin axis
    changes while ``L`` is nearly conserved). Returns ``nan`` when fewer than two Apophis sinks
    exist, ``utime`` cannot be read, or no valid post-CA dumps have at least two bound particles.
    """
    return _extract_mean_spin_period_hr(
        run_dir,
        prefix,
        apophis_sink_id,
        earth_sink_id=earth_sink_id,
        time_relation="after_ca",
        skip_earliest_in_window=False,
    )


def preflight(args: argparse.Namespace, base_dir: Path, output_root: Path) -> Tuple[Path, Path, Path, Path]:
    base_setup = base_dir / f"{args.prefix}.setup"
    base_input = base_dir / f"{args.prefix}.in"
    if not base_setup.is_file():
        raise FileNotFoundError(f"Missing setup file: {base_setup}")
    if not base_input.is_file():
        raise FileNotFoundError(f"Missing input file: {base_input}")

    phantom_root = Path(args.phantom_dir)
    must_exist = not args.dry_run
    phantomsetup_bin = resolve_phantom_executable(phantom_root, "phantomsetup", must_exist=must_exist)
    phantom_bin = resolve_phantom_executable(phantom_root, "phantom", must_exist=must_exist)

    output_root.mkdir(parents=True, exist_ok=True)
    return base_setup, base_input, phantomsetup_bin, phantom_bin


def write_samples_csv(path: Path, samples: List[RunSample], column_order: List[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        header = ["run_id"] + column_order
        writer.writerow(header)
        for idx, s in enumerate(samples, start=1):
            row: List[Any] = [idx]
            for col in column_order:
                if col == "mass_input_kg":
                    row.append(f"{s.mass_kg:.12g}" if s.mass_kg is not None else "")
                elif col == "use_dem":
                    row.append("T" if s.use_dem is True else "F" if s.use_dem is False else "")
                elif col == "use_shape_crop":
                    row.append("T" if s.use_shape_crop is True else "F" if s.use_shape_crop is False else "")
                elif col == "apophis_only":
                    row.append("T" if s.apophis_only is True else "F" if s.apophis_only is False else "")
                else:
                    v = getattr(s, col, None)
                    row.append(f"{float(v):.12g}" if v is not None else "")
            writer.writerow(row)


def summary_secondary_columns(col_order: List[str]) -> List[str]:
    """Columns after mass_input_kg in the summary CSV (mass is its own field)."""
    return [c for c in col_order if c != "mass_input_kg"]


def write_summary_csv(path: Path, records: List[RunRecord], param_column_order: List[str]) -> None:
    """Write full sobol_mass_outputs.csv: header plus one row per record (sorted by run_id)."""
    secondary = summary_secondary_columns(param_column_order)
    base_cols = [
        "run_id",
        "mass_input_kg",
        *secondary,
        "run_dir",
        "status",
        "closest_approach_km",
        "closest_approach_au",
        "dispersion_ratio",
        "unbound_fraction",
        "intrinsic_spin_period_hr",
        "approach_spin_period_hr",
        "post_flyby_spin_period_hr",
        "error",
    ]
    ordered = sorted(records, key=lambda r: r.run_id)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(base_cols)
        for row in ordered:
            extra = [row.param_columns.get(c, "") for c in secondary]
            writer.writerow(
                [
                    row.run_id,
                    f"{row.mass_input_kg:.12g}" if math.isfinite(row.mass_input_kg) else "",
                    *extra,
                    row.run_dir,
                    row.status,
                    f"{row.closest_approach_km:.12g}" if not math.isnan(row.closest_approach_km) else "",
                    f"{row.closest_approach_au:.12g}" if not math.isnan(row.closest_approach_au) else "",
                    f"{row.dispersion_ratio:.12g}" if not math.isnan(row.dispersion_ratio) else "",
                    f"{row.unbound_fraction:.12g}" if not math.isnan(row.unbound_fraction) else "",
                    f"{row.intrinsic_spin_period_hr:.12g}"
                    if not math.isnan(row.intrinsic_spin_period_hr)
                    else "",
                    f"{row.approach_spin_period_hr:.12g}"
                    if not math.isnan(row.approach_spin_period_hr)
                    else "",
                    f"{row.post_flyby_spin_period_hr:.12g}"
                    if not math.isnan(row.post_flyby_spin_period_hr)
                    else "",
                    row.error,
                ]
            )


def read_setup_logical(setup_path: Path, key: str, *, default: bool = False) -> bool:
    """Read a Fortran logical (T/F) from a patched or template ``.setup`` file."""
    if not setup_path.is_file():
        return default
    match = re.search(
        rf"^\s*{re.escape(key)}\s*=\s*(\S+)",
        setup_path.read_text(encoding="utf-8"),
        re.MULTILINE | re.IGNORECASE,
    )
    if not match:
        return default
    return match.group(1).strip().upper().startswith("T")


def resolve_run_apophis_only(sample: RunSample, setup_path: Path) -> bool:
    """True when this run has no Earth sink (explicit sample flag or template/setup file)."""
    if sample.apophis_only is not None:
        return sample.apophis_only
    return read_setup_logical(setup_path, "apophis_only")


def resolve_run_use_dem(sample: RunSample, setup_path: Path) -> bool:
    if sample.use_dem is not None:
        return sample.use_dem
    return read_setup_logical(setup_path, "use_dem")


def effective_apophis_sink_id(
    sample: RunSample, setup_path: Path, configured_id: int
) -> int:
    """Apophis DEM sinks start at 1 when ``apophis_only=T`` (not 11 in the full solar system)."""
    if resolve_run_apophis_only(sample, setup_path):
        if configured_id == APOPHIS_SINK_ID_DEFAULT:
            return 1
    return configured_id


def skip_closest_approach(sample: RunSample, setup_path: Optional[Path] = None) -> bool:
    """Earth is absent when apophis_only is True; default sink IDs do not apply."""
    if sample.apophis_only is True:
        return True
    if setup_path is not None:
        return resolve_run_apophis_only(sample, setup_path)
    return False


def run_one_case(
    run_id: int,
    sample: RunSample,
    base_setup: Path,
    base_input: Path,
    output_root: Path,
    prefix: str,
    phantomsetup_bin: Path,
    phantom_bin: Path,
    ref_mass_kg: Optional[float],
    dry_run: bool,
    earth_sink_id: int,
    apophis_sink_id: int,
    ephemeris_cache_dir: Optional[Path] = None,
    shape_file: Optional[str] = None,
    literature_scale_r_allowed: bool = True,
) -> RunRecord:
    run_dir = output_root / f"run_{run_id:04d}"
    run_dir.mkdir(parents=True, exist_ok=True)
    run_setup = run_dir / f"{prefix}.setup"
    run_input = run_dir / f"{prefix}.in"
    shutil.copy2(base_setup, run_setup)
    shutil.copy2(base_input, run_input)

    shape_setup_value: Optional[str] = None
    if sample.use_shape_crop and shape_file:
        shape_setup_value = stage_shape_assets(Path(shape_file), run_dir)
        _maybe_apply_literature_scale_r(
            sample, literature_scale_r_allowed=literature_scale_r_allowed
        )

    param_columns = apply_run_sample_to_setup(
        run_setup, sample, ref_mass_kg, shape_setup_value
    )

    if shape_setup_value is None:
        stage_shape_from_setup_template(run_setup, base_setup.parent, run_dir)

    if ephemeris_cache_dir is not None:
        copy_ephemeris_txt_cache(ephemeris_cache_dir, run_dir)

    apophis_only_run = resolve_run_apophis_only(sample, run_setup)
    apophis_sink_id_eff = effective_apophis_sink_id(sample, run_setup, apophis_sink_id)
    use_dem_run = resolve_run_use_dem(sample, run_setup)

    mass_for_record = float(sample.mass_kg) if sample.mass_kg is not None else float("nan")

    if dry_run:
        return RunRecord(
            run_id=run_id,
            mass_input_kg=mass_for_record,
            run_dir=str(run_dir),
            status="prepared_only",
            closest_approach_km=float("nan"),
            closest_approach_au=float("nan"),
            error="",
            param_columns=param_columns,
        )

    # Warn before starting if np_apophis is high enough that lattice overshoot will likely push
    # the actual sink count above the compiled MAXPTMASS limit (default 1000).
    maxptmass_warn = _np_apophis_maxptmass_warning(sample)
    if maxptmass_warn:
        print(f"[WARN] Run {run_id}: {maxptmass_warn}", file=sys.stderr, flush=True)

    try:
        # phantomsetup is caught separately so we can inspect setup.log for the specific
        # MAXPTMASS overflow message before falling back to a generic "command failed" error.
        # For DEM runs, pass --maxp to phantomsetup and phantom so they allocate O(np)
        # instead of the 6.3 GB default (maxp_alloc=5200000).  The temporary SPH lattice
        # during setup scales with np_apophis (~1.5×); after DEM conversion npart=0 so
        # phantom itself only needs a small ceiling.
        if use_dem_run:
            np_eff = sample.np_apophis
            if np_eff is None:
                np_eff = _read_np_apophis_from_setup(run_setup)
            use_shape = _resolve_use_shape_crop(sample, run_setup)
            if np_eff is not None and np_eff > 1:
                setup_maxp = _dem_maxp_for_np(np_eff, use_shape)
            else:
                setup_maxp = 4000 if use_shape else 2000
            maxp_flag = [f"--maxp={setup_maxp}"]
            print(
                f"[INFO] Run {run_id}: DEM --maxp={setup_maxp} "
                f"(np_apophis={np_eff}, shape_crop={use_shape})",
                flush=True,
            )
        else:
            maxp_flag = []
        setup_log = run_dir / "setup.log"
        try:
            run_command([str(phantomsetup_bin), prefix] + maxp_flag, cwd=run_dir, log_path=setup_log)
            _verify_phantomsetup_outputs(run_dir, prefix, setup_log)
            param_columns.update(_parse_setup_roche_metrics(setup_log))
        except RuntimeError:
            diag = _diagnose_setup_log(setup_log)
            raise RuntimeError(
                diag or f"phantomsetup failed (returncode != 0); see {setup_log}"
            ) from None
        if use_dem_run:
            in_text = run_input.read_text()
            in_text = replace_setup_assignment(in_text, "nfulldump", f"{1:>10}")
            run_input.write_text(in_text)
            print(f"[INFO] Run {run_id}: DEM enabled — set nfulldump=1 in {run_input.name}", flush=True)
        in_param_cols = apply_run_sample_to_in(run_input, sample)
        param_columns.update(in_param_cols)
        run_command([str(phantom_bin), f"{prefix}.in"] + maxp_flag, cwd=run_dir, log_path=run_dir / "phantom.log")

        # Dump-grid metrics only (fast path); ignore METRICS_LEGACY_SUBSTEPS from the shell.
        _use_fast_metrics_defaults()
        _groups, _time_of_key, _n_sinks = _apophis_time_groups(run_dir, prefix, apophis_sink_id_eff)

        dispersion_ratio = float("nan")
        unbound_fraction = float("nan")
        intrinsic_spin_period_hr = float("nan")
        approach_spin_period_hr = float("nan")
        post_flyby_spin_period_hr = float("nan")
        closest_km = float("nan")
        closest_au = float("nan")
        t_ca: Optional[float] = None

        if use_dem_run and _n_sinks >= 2:
            if apophis_only_run:
                (
                    dispersion_ratio,
                    unbound_fraction,
                    intrinsic_spin_period_hr,
                    _,
                    _,
                ) = _extract_dem_metrics_bundle(
                    run_dir,
                    prefix,
                    apophis_sink_id_eff,
                    apophis_only=True,
                    _groups=_groups,
                    _time_of_key=_time_of_key,
                    _n_sinks=_n_sinks,
                )
            else:
                try:
                    t_ca = _closest_approach_time(
                        run_dir, prefix, earth_sink_id, apophis_sink_id_eff
                    )
                except RuntimeError:
                    t_ca = None
                closest_km, closest_au = _resolve_closest_approach(
                    run_dir, prefix, earth_sink_id, apophis_sink_id_eff, param_columns
                )
                (
                    dispersion_ratio,
                    unbound_fraction,
                    intrinsic_spin_period_hr,
                    approach_spin_period_hr,
                    post_flyby_spin_period_hr,
                ) = _extract_dem_metrics_bundle(
                    run_dir,
                    prefix,
                    apophis_sink_id_eff,
                    apophis_only=False,
                    earth_sink_id=earth_sink_id,
                    t_ca=t_ca,
                    _groups=_groups,
                    _time_of_key=_time_of_key,
                    _n_sinks=_n_sinks,
                )
        elif apophis_only_run:
            return RunRecord(
                run_id=run_id,
                mass_input_kg=mass_for_record,
                run_dir=str(run_dir),
                status="ok",
                closest_approach_km=float("nan"),
                closest_approach_au=float("nan"),
                error="",
                dispersion_ratio=dispersion_ratio,
                unbound_fraction=unbound_fraction,
                intrinsic_spin_period_hr=intrinsic_spin_period_hr,
                approach_spin_period_hr=approach_spin_period_hr,
                post_flyby_spin_period_hr=post_flyby_spin_period_hr,
                param_columns=param_columns,
            )
        else:
            closest_km, closest_au = _resolve_closest_approach(
                run_dir, prefix, earth_sink_id, apophis_sink_id_eff, param_columns
            )
        return RunRecord(
            run_id=run_id,
            mass_input_kg=mass_for_record,
            run_dir=str(run_dir),
            status="ok",
            closest_approach_km=closest_km,
            closest_approach_au=closest_au,
            error="",
            dispersion_ratio=dispersion_ratio,
            unbound_fraction=unbound_fraction,
            intrinsic_spin_period_hr=intrinsic_spin_period_hr,
            approach_spin_period_hr=approach_spin_period_hr,
            post_flyby_spin_period_hr=post_flyby_spin_period_hr,
            param_columns=param_columns,
        )
    except Exception as exc:  # pragma: no cover - runtime path
        return RunRecord(
            run_id=run_id,
            mass_input_kg=mass_for_record,
            run_dir=str(run_dir),
            status="failed",
            closest_approach_km=float("nan"),
            closest_approach_au=float("nan"),
            error=str(exc),
            param_columns=param_columns,
        )


def _execute_run_worker(payload: RunWorkerPayload) -> RunRecord:
    return run_one_case(
        payload.run_id,
        payload.sample,
        Path(payload.base_setup),
        Path(payload.base_input),
        Path(payload.output_root),
        payload.prefix,
        Path(payload.phantomsetup_bin),
        Path(payload.phantom_bin),
        payload.ref_mass_kg,
        payload.dry_run,
        payload.earth_sink_id,
        payload.apophis_sink_id,
        Path(payload.ephemeris_cache_dir) if payload.ephemeris_cache_dir else None,
        payload.shape_file,
        payload.literature_scale_r_allowed,
    )


def _cleanup_run_dir(run_dir: Path, prefix: str) -> None:
    """Delete heavy output files after metrics are safely extracted.

    Removes binary dump files (sobol_NNNNN), per-sink .ev files, and phantom.log.
    Keeps .in, .setup, setup.log, and ephemeris .txt files for audit purposes.
    """
    removed = 0
    for pattern in (f"{prefix}_[0-9]*", f"{prefix}Sink*N*.ev", f"{prefix}*.ev", "phantom.log"):
        for f in run_dir.glob(pattern):
            if f.is_file():
                f.unlink()
                removed += 1
    print(f"[INFO] Cleaned {run_dir.name}: removed {removed} heavy files ({prefix}_* dumps, *.ev, phantom.log)", flush=True)


def _print_run_progress(result: RunRecord, samples: List[RunSample], total: int) -> None:
    idx = result.run_id
    if not (1 <= idx <= len(samples)):
        raise IndexError(f"run_id {idx} out of range for samples list of length {len(samples)}")
    sample = samples[idx - 1]
    mass_str = f"{sample.mass_kg:.6e} kg" if sample.mass_kg is not None else "(template mass)"
    np_str = f" np_apophis={sample.np_apophis}" if sample.np_apophis is not None else ""
    print(f"[INFO] Run {idx}/{total} mass={mass_str}{np_str}", flush=True)
    print(f"[INFO]   status={result.status} run_dir={result.run_dir}", flush=True)


def main() -> int:
    argv_cli, interactive = _strip_interactive_flags(sys.argv[1:])
    if interactive:
        from tui_run import SobolTUIApp

        initial_args = parse_args(argv_cli)
        wizard_argv = SobolTUIApp(initial_args).run()
        if wizard_argv is None:
            return 0
        args = parse_args(wizard_argv)
    else:
        args = parse_args(argv_cli if argv_cli else None)
    base_dir = Path(args.base_dir).resolve()
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    problem: Optional[Dict[str, object]] = None
    saltelli_meta: Optional[Dict[str, object]] = None
    X_saltelli: Optional[Any] = None
    ephemeris_cache: Optional[Path] = None

    try:
        validate_args(args)
        ephemeris_cache = resolve_ephemeris_cache_dir(base_dir, args.ephemeris_cache_dir)
        if args.ephemeris_cache_dir:
            if ephemeris_cache is None or not ephemeris_cache.is_dir():
                raise FileNotFoundError(
                    f"Ephemeris cache path is not a directory: {args.ephemeris_cache_dir}"
                )
        if ephemeris_cache is not None:
            txt_n = sum(1 for p in ephemeris_cache.glob("*.txt") if p.is_file())
            if txt_n == 0:
                print(
                    f"[WARN] No *.txt files in ephemeris cache {ephemeris_cache}; "
                    "PHANTOM will attempt live Horizons downloads.",
                    file=sys.stderr,
                )
            else:
                src = "base_dir" if not args.ephemeris_cache_dir else "--ephemeris-cache-dir"
                print(
                    f"[INFO] Ephemeris cache ({src}): {txt_n} *.txt file(s) from {ephemeris_cache}",
                    flush=True,
                )
        elif not args.dry_run:
            print(
                "[WARN] No ephemeris *.txt cache configured (set --ephemeris-cache-dir or "
                "place Horizons .txt files in --base-dir). phantomsetup may fail in isolated "
                "run folders.",
                file=sys.stderr,
            )
        batch_basename = build_batch_directory_basename(args, timestamp)
        output_root = Path(args.output_root).resolve() / batch_basename

        if getattr(args, "np_apophis_list", None) is not None:
            # One run per particle count; Sobol sampling is not used in this mode.
            samples = build_np_list_samples(args)
        elif getattr(args, "scale_earth_sep_list", None) is not None:
            samples = build_scale_earth_sep_list_samples(args)
        elif args.saltelli_n is not None:
            problem = build_salib_problem(args)
            try:
                from SALib.sample import sobol as sobol_sample
            except ImportError as exc:
                raise RuntimeError("Saltelli mode requires SALib (pip install -r requirements.txt)") from exc
            X_saltelli = sobol_sample.sample(
                problem,
                args.saltelli_n,
                calc_second_order=args.saltelli_calc_second_order,
                seed=args.seed,
            )
            samples = [run_sample_from_salib_row(X_saltelli[i], args) for i in range(len(X_saltelli))]
            ne = expected_saltelli_num_evals(
                int(problem["num_vars"]), args.saltelli_n, args.saltelli_calc_second_order
            )
            if len(samples) != ne:
                raise RuntimeError(f"internal error: expected {ne} Saltelli rows, got {len(samples)}")
            saltelli_meta = {
                "base_n": args.saltelli_n,
                "calc_second_order": args.saltelli_calc_second_order,
                "num_evals": ne,
                "num_vars": problem["num_vars"],
                "seed": args.seed,
            }
        else:
            samples = build_run_samples(args.num_samples, args)

        base_setup, base_input, phantomsetup_bin, phantom_bin = preflight(args, base_dir, output_root)

        if args.saltelli_n is not None and problem is not None and saltelli_meta is not None:
            (output_root / "saltelli_problem.json").write_text(
                json.dumps(problem, indent=2), encoding="utf-8"
            )
            (output_root / "saltelli_meta.json").write_text(
                json.dumps(saltelli_meta, indent=2), encoding="utf-8"
            )
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    col_order = sample_column_order(args)
    samples_csv = output_root / "sobol_mass_samples.csv"
    summary_csv = output_root / "sobol_mass_outputs.csv"
    write_samples_csv(samples_csv, samples, col_order)
    print(f"[INFO] Wrote Sobol sample table: {samples_csv}")
    dim_msg = f"Varying dimensions ({count_dimensions(args)}): {', '.join(col_order)}"
    if args.saltelli_n is not None and saltelli_meta is not None:
        print(
            f"[INFO] Saltelli mode: base_n={args.saltelli_n} evals={saltelli_meta['num_evals']} "
            f"second_order={args.saltelli_calc_second_order}; {dim_msg}"
        )
    else:
        print(f"[INFO] {dim_msg}")
    if args.jobs > 1:
        print(
            f"[INFO] Running {len(samples)} cases with up to {args.jobs} worker process(es) "
            "(set OMP_NUM_THREADS=1 if PHANTOM is OpenMP to limit threads per process)."
        )

    literature_scale_r_allowed = "scale_r_apophis" not in {
        p for p, _, _, _ in _active_scale_variations(args)
    }
    payloads = [
        RunWorkerPayload(
            run_id=idx,
            sample=sample,
            base_setup=str(base_setup),
            base_input=str(base_input),
            output_root=str(output_root),
            prefix=args.prefix,
            phantomsetup_bin=str(phantomsetup_bin),
            phantom_bin=str(phantom_bin),
            ref_mass_kg=args.apophis_ref_mass_kg if _mass_bounds_active(args) else None,
            dry_run=args.dry_run,
            earth_sink_id=args.sink_earth_id,
            apophis_sink_id=args.sink_apophis_id,
            ephemeris_cache_dir=str(ephemeris_cache) if ephemeris_cache is not None else None,
            shape_file=args.shape_file if args.shape_file else None,
            literature_scale_r_allowed=literature_scale_r_allowed,
        )
        for idx, sample in enumerate(samples, start=1)
    ]
    total_runs = len(samples)
    results: List[RunRecord] = []
    if args.jobs == 1:
        for p in payloads:
            result = _execute_run_worker(p)
            results.append(result)
            _print_run_progress(result, samples, total_runs)
            write_summary_csv(summary_csv, results, col_order)
            print(f"[INFO] Partial results saved: {len(results)}/{total_runs} runs → {summary_csv}", flush=True)
            if not args.no_cleanup and result.status == "ok" and result.run_dir:
                _cleanup_run_dir(Path(result.run_dir), p.prefix)
    else:
        payload_by_id = {p.run_id: p for p in payloads}
        with concurrent.futures.ProcessPoolExecutor(max_workers=args.jobs) as executor:
            futures = [executor.submit(_execute_run_worker, p) for p in payloads]
            for fut in as_completed(futures):
                result = fut.result()
                results.append(result)
                _print_run_progress(result, samples, total_runs)
                write_summary_csv(summary_csv, results, col_order)
                print(f"[INFO] Partial results saved: {len(results)}/{total_runs} runs → {summary_csv}", flush=True)
                p = payload_by_id[result.run_id]
                if not args.no_cleanup and result.status == "ok" and result.run_dir:
                    _cleanup_run_dir(Path(result.run_dir), p.prefix)

    write_summary_csv(summary_csv, results, col_order)

    sorted_results = sorted(results, key=lambda r: r.run_id)

    saltelli_y_rows: List[Dict[str, object]] = []
    if args.saltelli_n is not None:
        for result in sorted_results:
            saltelli_y_rows.append(
                {
                    "eval_index": result.run_id - 1,
                    "run_id": result.run_id,
                    "closest_approach_km": result.closest_approach_km,
                    "closest_approach_au": result.closest_approach_au,
                    "dispersion_ratio": result.dispersion_ratio,
                    "unbound_fraction": result.unbound_fraction,
                    "intrinsic_spin_period_hr": result.intrinsic_spin_period_hr,
                    "approach_spin_period_hr": result.approach_spin_period_hr,
                    "post_flyby_spin_period_hr": result.post_flyby_spin_period_hr,
                    "status": result.status,
                }
            )

    print(f"[INFO] Summary written to: {summary_csv}")

    if args.saltelli_n is not None and X_saltelli is not None and problem is not None:
        num_vars = int(problem["num_vars"])
        if X_saltelli.shape[1] != num_vars:
            raise ValueError(
                f"Saltelli design has {X_saltelli.shape[1]} columns but problem specifies "
                f"num_vars={num_vars}; design matrix and problem dict are inconsistent."
            )
        manifest_path = output_root / "saltelli_eval_manifest.csv"
        names = list(problem["names"])
        with manifest_path.open("w", newline="", encoding="utf-8") as mf:
            mw = csv.writer(mf)
            mw.writerow(["eval_index", "run_id", *names])
            for j in range(len(samples)):
                row_vals = [float(X_saltelli[j, k]) for k in range(num_vars)]
                mw.writerow([j, j + 1, *row_vals])
        print(f"[INFO] Wrote Saltelli eval manifest: {manifest_path}")

        y_path = output_root / "saltelli_Y.csv"
        y_fields = [
            "eval_index",
            "run_id",
            "closest_approach_km",
            "closest_approach_au",
            "dispersion_ratio",
            "unbound_fraction",
            "intrinsic_spin_period_hr",
            "approach_spin_period_hr",
            "post_flyby_spin_period_hr",
            "status",
        ]
        with y_path.open("w", newline="", encoding="utf-8") as yf:
            yw = csv.DictWriter(yf, fieldnames=y_fields)
            yw.writeheader()
            for row in saltelli_y_rows:
                yw.writerow(row)
        print(f"[INFO] Wrote Saltelli model outputs (evaluation order): {y_path}")
        prob_p = (output_root / "saltelli_problem.json").resolve()
        meta_p = (output_root / "saltelli_meta.json").resolve()
        print(
            "[INFO] Sobol indices: python3 Analysis/Analysis.py --method saltelli \\\n"
            f"       --sobol-problem-json {prob_p} \\\n"
            f"       --saltelli-meta-json {meta_p} \\\n"
            f"       --saltelli-y-csv {y_path.resolve()} \\\n"
            "       --saltelli-y-column closest_approach_au  "
            "# also: dispersion_ratio, unbound_fraction, intrinsic_spin_period_hr, "
            "approach_spin_period_hr, post_flyby_spin_period_hr (DEM sweeps) "
            + (
                " \\\n       --saltelli-calc-second-order"
                if args.saltelli_calc_second_order
                else ""
            )
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
