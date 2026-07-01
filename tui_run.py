#!/usr/bin/env python3
"""Textual TUI for configuring and launching the PHANTOM Sobol sweep.

Launch:
    python3 tui_run.py [any flags from run_mass_sobol_phantom.py]
    python3 run_mass_sobol_phantom.py --interactive   # same TUI, then runs in-process

CLI flags pre-fill form defaults.  Navigate with Tab / Shift-Tab, Space to
toggle checkboxes, then press [Run Sweep], [Dry Run], or [Quit] (also bound
to R, D, Q).
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from typing import List, Optional

sys.path.insert(0, str(Path(__file__).parent))

try:
    from textual import on
    from textual.app import App, ComposeResult
    from textual.binding import Binding
    from textual.containers import Horizontal, ScrollableContainer
    from textual.widgets import Button, Checkbox, Footer, Header, Input, Select, Static
except ImportError:
    sys.exit(
        "textual is not installed for this Python interpreter.\n"
        f"  Python: {sys.executable}\n"
        "Install with:\n"
        "  python3 -m pip install textual\n"
        "or:\n"
        "  pip3 install textual\n"
        "(On macOS, `pip` is often missing; use one of the above.)"
    )

import run_mass_sobol_phantom as _runner

# ── constants ────────────────────────────────────────────────────────────────

# (CSS-id-stem, attr-name, checkbox label, is_integer)
# attr-name is the argparse dest STEM: min/max are read as attr+"_min"/attr+"_max".
# Spin entries use "spin_period" etc. (argparse dest stem), not the full setup key.
_SCALE_DIMS = [
    ("scale-vel",       "scale_vel",        "scale_vel  — Apophis velocity scale",             False),
    ("scale-pos",       "scale_pos",        "scale_pos  — heliocentric position (not Earth distance)", False),
    ("scale-earth-sep", "scale_earth_sep",  "scale_earth_sep  — Earth–Apophis distance (1=ephemeris)", False),
    ("scale-r-apophis", "scale_r_apophis",  "scale_r_apophis  — Apophis radius scale",         False),
    ("scale-rho",       "scale_rho",        "scale_rho  — bulk density scale",                 False),
]

# Spin keys match setup_solarsystem.f90 (period in seconds; axis as x/y/z components).
# argparse dest stems: spin_period → --spin-period-min; spin_axis_x → --spin-axis-x-min.
_SPIN_DIMS = [
    ("spin-period",  "spin_period",    "apophis_spin_period (s)  — spin period (0 = no spin)", False),
    ("spin-axis-x",  "spin_axis_x",    "apophis_spin_axis_x  — spin axis x component",         False),
    ("spin-axis-y",  "spin_axis_y",    "apophis_spin_axis_y  — spin axis y component",         False),
    ("spin-axis-z",  "spin_axis_z",    "apophis_spin_axis_z  — spin axis z component",         False),
]

_MASS_GATE_IDS = ("mass-min-kg", "mass-max-kg", "apophis-ref-mass-kg")

# DEM contact parameters: patched into .in after phantomsetup.
# (css-id-stem, argparse-dest-stem, label, is_integer)
# dest stem follows argparse convention: --kc-min → kc_min, getattr(d, "kc_min", None)
_IN_DIMS = [
    ("kc",    "kc",    "kt_cgs  — DEM tensile spring (g/s^2 per cm gap; 0=off; DEM only)", False),
    ("ct",    "ct",    "ct_dem  — tangential damping coefficient (DEM only)",                   False),
    ("eps-n", "eps_n", "epsilon_n_dem  — normal restitution coeff [0,1] (DEM only)",            False),
    ("kn",    "kn",    "kn_cgs  — normal spring constant (dyne/cm; DEM only)",                  False),
]


# ── helper: one labelled row ─────────────────────────────────────────────────

class _Row(Horizontal):
    """Label + interactive widget on a single row, with optional right-side hint."""

    DEFAULT_CSS = """
    _Row {
        height: 3;
    }
    _Row .rl {
        width: 26;
        content-align: right middle;
        padding-right: 1;
        color: $text-muted;
    }
    _Row .rw {
        width: 1fr;
    }
    _Row .rh {
        width: 15;
        content-align: left middle;
        padding-left: 1;
        color: $text-muted;
        text-style: italic;
    }
    """

    def __init__(self, label: str, widget, hint: str = "", **kw) -> None:
        super().__init__(**kw)
        self._label = label
        self._widget = widget
        self._hint = hint

    def compose(self) -> ComposeResult:
        yield Static(self._label, classes="rl")
        self._widget.add_class("rw")
        yield self._widget
        if self._hint:
            yield Static(self._hint, classes="rh")


# ── main app ─────────────────────────────────────────────────────────────────

class SobolTUIApp(App[Optional[List[str]]]):
    """Full-screen parameter form for the PHANTOM Sobol sweep."""

    TITLE = "PHANTOM · Sobol Sweep"
    SUB_TITLE = "edit parameters → launch"

    BINDINGS = [
        Binding("r", "run_sweep", "Run"),
        Binding("d", "do_dry_run", "Dry-run"),
        Binding("q", "app_quit", "Quit"),
    ]

    DEFAULT_CSS = """
    Screen { background: $surface; }

    #scroll {
        height: 1fr;
        border: solid $primary;
        margin: 0 1;
        padding: 0 2;
    }

    .sec {
        color: $accent;
        text-style: bold;
        padding: 1 0 0 0;
    }

    #bar {
        height: 5;
        padding: 1 2;
        background: $surface-darken-1;
        border-top: solid $primary;
        align: center middle;
    }

    #status {
        width: 1fr;
        content-align: left middle;
        padding-left: 2;
    }

    .err { color: $error; }
    .ok  { color: $success; }
    .note { color: $warning; text-style: italic; padding: 0 0 0 27; }

    Button { margin: 0 1; }
    """

    def __init__(self, defaults) -> None:
        super().__init__()
        self._d = defaults

    # ── compose ──────────────────────────────────────────────────────────────

    def compose(self) -> ComposeResult:  # noqa: C901
        d = self._d

        yield Header()

        with ScrollableContainer(id="scroll"):

            # ── Paths ──────────────────────────────────────────────────────
            yield Static("Paths", classes="sec")
            yield _Row("prefix",
                       Input(str(d.prefix or _runner._DEFAULT_PREFIX), id="prefix"))
            yield _Row("base_dir",
                       Input(str(d.base_dir or _runner._DEFAULT_BASE_DIR), id="base-dir"), "dir")
            yield _Row("ephemeris_cache_dir",
                       Input(str(d.ephemeris_cache_dir or ""), id="eph-cache",
                             placeholder="(blank → copy *.txt from base_dir)"),
                       "optional dir")
            yield _Row("phantom_dir",
                       Input(str(d.phantom_dir or ""), id="phantom-dir"), "dir")
            yield _Row("output_root",
                       Input(str(d.output_root or "sobol_mass_runs"), id="output-root"),
                       "dir")

            # ── Sampling ───────────────────────────────────────────────────
            yield Static("Sampling", classes="sec")
            yield _Row("num_samples",
                       Input(str(d.num_samples if d.num_samples is not None else 8),
                             id="num-samples"),
                       "int ≥ 1")
            yield _Row("seed",
                       Input(str(d.seed if d.seed is not None else 42), id="seed"),
                       "int")
            yield _Row("saltelli_n",
                       Input(str(d.saltelli_n or ""), id="saltelli-n",
                             placeholder="(blank → use num_samples)"),
                       "int ≥ 2 or blank")
            yield _Row("saltelli_2nd_order",
                       Checkbox("also estimate 2nd-order Sobol indices",
                                value=bool(getattr(d, "saltelli_calc_second_order", False)),
                                id="saltelli-2nd"))

            # ── Mass ───────────────────────────────────────────────────────
            yield Static("Mass", classes="sec")
            mass_on = d.mass_min_kg is not None and d.mass_max_kg is not None
            yield _Row("vary mass",
                       Checkbox("vary Apophis mass as a Sobol dimension",
                                value=mass_on, id="vary-mass"))
            yield _Row("mass_min_kg",
                       Input(str(d.mass_min_kg or ""), id="mass-min-kg",
                             disabled=not mass_on, placeholder="e.g. 1e10"),
                       "kg")
            yield _Row("mass_max_kg",
                       Input(str(d.mass_max_kg or ""), id="mass-max-kg",
                             disabled=not mass_on, placeholder="e.g. 1e11"),
                       "kg")
            yield _Row("apophis_ref_mass_kg",
                       Input(str(getattr(d, "apophis_ref_mass_kg", "") or ""), id="apophis-ref-mass-kg",
                             disabled=not mass_on,
                             placeholder="e.g. 2.7e10  (mass at scale_rho=1)"),
                       "kg")

            # ── Scale bounds ───────────────────────────────────────────────
            yield Static("Scale bounds  (optional Sobol dimensions)", classes="sec")
            for css, attr, label, is_int in _SCALE_DIMS:
                lo = getattr(d, attr + "_min", None)
                hi = getattr(d, attr + "_max", None)
                active = lo is not None and hi is not None
                hint = "int" if is_int else "float"
                yield _Row(f"vary {attr}",
                           Checkbox(label, value=active, id=f"vary-{css}"))
                yield _Row(f"  {attr}_min",
                           Input(str(lo or ""), id=f"{css}-min",
                                 disabled=not active, placeholder="min"),
                           hint)
                yield _Row(f"  {attr}_max",
                           Input(str(hi or ""), id=f"{css}-max",
                                 disabled=not active, placeholder="max"),
                           hint)

            # ── Spin (setup_solarsystem.f90) ───────────────────────────────
            yield Static(
                "Spin  (optional Sobol dimensions — seconds + axis x/y/z; DEM runs only)",
                classes="sec",
            )
            yield Static(
                "Fixed spin comes from your template .setup; check a box below to sweep that parameter.",
                classes="note",
            )
            for css, attr, label, is_int in _SPIN_DIMS:
                lo = getattr(d, attr + "_min", None)
                hi = getattr(d, attr + "_max", None)
                active = lo is not None and hi is not None
                hint = "int" if is_int else "float"
                yield _Row(f"vary {attr}",
                           Checkbox(label, value=active, id=f"vary-{css}"))
                yield _Row(f"  {attr}_min",
                           Input(str(lo or ""), id=f"{css}-min",
                                 disabled=not active, placeholder="min"),
                           hint)
                yield _Row(f"  {attr}_max",
                           Input(str(hi or ""), id=f"{css}-max",
                                 disabled=not active, placeholder="max"),
                           hint)

            # ── DEM Contact parameters ─────────────────────────────────────
            yield Static("DEM Contact  (optional Sobol dimensions — .in file, DEM mode only)", classes="sec")
            for css, attr, label, is_int in _IN_DIMS:
                lo = getattr(d, attr + "_min", None)
                hi = getattr(d, attr + "_max", None)
                active = lo is not None and hi is not None
                hint = "int" if is_int else "float"
                yield _Row(f"vary {attr}",
                           Checkbox(label, value=active, id=f"vary-{css}"))
                yield _Row(f"  {attr}_min",
                           Input(str(lo or ""), id=f"{css}-min",
                                 disabled=not active, placeholder="min"),
                           hint)
                yield _Row(f"  {attr}_max",
                           Input(str(hi or ""), id=f"{css}-max",
                                 disabled=not active, placeholder="max"),
                           hint)

            # ── Timeframe ──────────────────────────────────────────────────
            yield Static("Timeframe  (fix for all runs, or sweep as Sobol dimension)", classes="sec")
            for t_css, t_attr, t_label in (
                ("tmax",  "tmax_hours",  "tmax_in — simulation end time"),
                ("dtmax", "dtmax_hours", "dtmax_in — dump interval"),
            ):
                t_lo  = getattr(d, t_attr + "_min", None)
                t_hi  = getattr(d, t_attr + "_max", None)
                t_fix = getattr(d, t_attr, None)
                t_vary = t_lo is not None and t_hi is not None
                yield _Row(f"vary {t_attr}",
                           Checkbox(f"{t_label}  (sweep as Sobol dim)",
                                    value=t_vary, id=f"vary-{t_css}"))
                yield _Row(f"  {t_attr} fixed",
                           Input(str(t_fix or ""), id=f"{t_css}-fixed",
                                 disabled=t_vary,
                                 placeholder="hours  e.g. 108 = 4.5 days  (blank = keep template)"),
                           "hr  fixed")
                yield _Row(f"  {t_attr}_min",
                           Input(str(t_lo or ""), id=f"{t_css}-min",
                                 disabled=not t_vary, placeholder="min hr"),
                           "hr")
                yield _Row(f"  {t_attr}_max",
                           Input(str(t_hi or ""), id=f"{t_css}-max",
                                 disabled=not t_vary, placeholder="max hr"),
                           "hr")

            # ── Roche / Earth distance ─────────────────────────────────────
            yield Static("Roche limit  (Earth–Apophis separation sweep)", classes="sec")
            yield Static(
                "Use apophis_only=F (Earth present). Vary scale_earth_sep, not scale_pos. "
                "Analysis plots use setup pericentre r_p/r_tidal (close approach), not initial "
                "separation at epoch. Set spin=0 and kt_cgs=0 in template for baseline. "
                "Rebuild PHANTOM after pulling setup changes (make SETUP=solarsystem setup && make).",
                classes="note",
            )
            sep_list_val = " ".join(
                str(x) for x in (getattr(d, "scale_earth_sep_list", None) or [])
            )
            yield _Row("  scale_earth_sep_list",
                       Input(sep_list_val, id="scale-earth-sep-list",
                             placeholder="e.g. 0.4 0.6 0.8 1.0 1.2 1.5  — one run per value"),
                       "Roche sweep")
            yield _Row("  spin_period_fixed",
                       Input(str(getattr(d, "spin_period_fixed", "") or ""),
                             id="spin-period-fixed",
                             placeholder="0 = no spin (recommended for Roche baseline)"),
                       "seconds")
            yield _Row("  apophis_only_fixed",
                       Select([("(from template)", ""),
                               ("False — Earth present (Roche)", "false"),
                               ("True — no Earth", "true")],
                              value=getattr(d, "apophis_only_fixed", None) or "",
                              id="apophis-only-fixed"),
                       "required for Roche")

            # ── Setup toggles ──────────────────────────────────────────────
            yield Static("Setup toggles", classes="sec")
            vary_dem = bool(getattr(d, "vary_use_dem", False))
            dem_fixed_raw = getattr(d, "use_dem_fixed", None)
            dem_fixed_val = dem_fixed_raw if dem_fixed_raw in ("true", "false") else ""
            yield _Row("vary_use_dem",
                       Checkbox("vary use_dem as a Sobol dimension  (u≥0.5 → True)",
                                value=vary_dem, id="vary-use-dem"))
            yield _Row("  use_dem_fixed",
                       Select([("(from template)", ""),
                               ("True", "true"),
                               ("False", "false")],
                              value=dem_fixed_val,
                              id="use-dem-fixed", disabled=vary_dem),
                       "if not varying")
            np_ap_val = str(getattr(d, "np_apophis", "") or "")
            np_list_raw = getattr(d, "np_apophis_list", None)
            np_list_val = " ".join(str(n) for n in np_list_raw) if np_list_raw else ""
            # np_apophis and np_apophis_list are mutually exclusive; disable whichever is empty.
            np_list_active = bool(np_list_val)
            yield _Row("  np_apophis",
                       Input(np_ap_val, id="np-apophis",
                             disabled=np_list_active,
                             placeholder="(from template)  0=none  1=sink  N=gas/DEM"),
                       "int ≥ 0")
            yield _Row("  np_apophis_list",
                       Input(np_list_val, id="np-apophis-list",
                             disabled=not np_list_active,
                             placeholder="e.g. 250 500 1000  — one run per value in this batch"),
                       "space-sep ints")
            yield Static(
                f"⚠ np_apophis ≥ {_runner._MAXPTMASS_WARN_THRESHOLD}: lattice overshoot (~2-3%) "
                f"will exceed default MAXPTMASS={_runner._MAXPTMASS_DEFAULT}. "
                "Rebuild with MAXPTMASS=2000 in the Makefile ('make setup && make').",
                classes="note",
            )
            vary_shape = bool(getattr(d, "vary_use_shape_crop", False))
            shape_fixed_raw = getattr(d, "use_shape_crop_fixed", None)
            shape_fixed_val = shape_fixed_raw if shape_fixed_raw in ("true", "false") else ""
            shape_path_active = vary_shape or shape_fixed_val == "true"
            yield _Row("vary_use_shape_crop",
                       Checkbox("vary shape cropping as a Sobol dimension  (u≥0.5 → crop)",
                                value=vary_shape, id="vary-use-shape-crop"))
            yield _Row("  use_shape_crop_fixed",
                       Select([("(from template)", ""),
                               ("True", "true"),
                               ("False", "false")],
                              value=shape_fixed_val,
                              id="use-shape-crop-fixed", disabled=vary_shape),
                       "if not varying")
            yield _Row("  shape_file",
                       Input(str(getattr(d, "shape_file", "") or ""), id="shape-file",
                             disabled=not shape_path_active,
                             placeholder="path/to/apophis.shape or apophis.obj  (required when cropping on)"),
                       "shape/OBJ path")
            yield _Row("vary_apophis_only",
                       Checkbox("vary apophis_only (Earth absent when True; CA → NaN)",
                                value=bool(getattr(d, "vary_apophis_only", False)),
                                id="vary-apophis-only"))

            # ── Sinks / post-processing ────────────────────────────────────
            yield Static("Sinks / post-processing", classes="sec")
            yield _Row("sink_earth_id",
                       Input(str(getattr(d, "sink_earth_id", _runner.EARTH_SINK_ID_DEFAULT)),
                             id="sink-earth"),
                       "int")
            yield _Row("sink_apophis_id",
                       Input(str(getattr(d, "sink_apophis_id", _runner.APOPHIS_SINK_ID_DEFAULT)),
                             id="sink-apophis"),
                       "int  (use 1 when apophis_only=T)")

            # ── Batch naming ───────────────────────────────────────────────
            yield Static("Batch naming", classes="sec")
            yield _Row("batch_label",
                       Input(str(getattr(d, "batch_label", "") or ""), id="batch-label",
                             placeholder="(auto slug from varied params)"),
                       "optional")
            yield _Row("batch_slug_max_len",
                       Input(str(getattr(d, "batch_slug_max_len",
                                         _runner._BATCH_SLUG_DEFAULT_MAX_LEN)),
                             id="slug-max-len"),
                       "int ≥ 17")

            # ── Execution ──────────────────────────────────────────────────
            yield Static("Execution", classes="sec")
            yield _Row("dry_run",
                       Checkbox("dry run  (prepare dirs only, skip PHANTOM binary)",
                                value=bool(getattr(d, "dry_run", False)),
                                id="dry-run"))
            yield _Row("jobs",
                       Input(str(getattr(d, "jobs", 1)), id="jobs"),
                       "parallel workers")

        # ── button bar ─────────────────────────────────────────────────────
        with Horizontal(id="bar"):
            yield Button("Run Sweep", id="btn-run",    variant="success")
            yield Button("Dry Run",   id="btn-dryrun", variant="warning")
            yield Button("Quit",      id="btn-quit",   variant="error")
            yield Static("Ready — edit fields above, then press Run Sweep.", id="status")

        yield Footer()

    # ── gate handlers ─────────────────────────────────────────────────────────

    @on(Checkbox.Changed)
    def _checkbox_gate(self, event: Checkbox.Changed) -> None:
        cb_id = event.checkbox.id
        active = event.value

        if cb_id == "vary-mass":
            for fid in _MASS_GATE_IDS:
                self.query_one(f"#{fid}").disabled = not active

        elif cb_id == "vary-use-dem":
            self.query_one("#use-dem-fixed").disabled = active


        elif cb_id == "vary-use-shape-crop":
            self.query_one("#use-shape-crop-fixed").disabled = active
            if active:
                self.query_one("#shape-file").disabled = False
            else:
                fixed = self._sel("use-shape-crop-fixed")
                self.query_one("#shape-file").disabled = (fixed != "true")

        elif cb_id in ("vary-tmax", "vary-dtmax"):
            stem = cb_id[len("vary-"):]   # "tmax" or "dtmax"
            self.query_one(f"#{stem}-fixed").disabled = active
            self.query_one(f"#{stem}-min").disabled   = not active
            self.query_one(f"#{stem}-max").disabled   = not active

        else:
            for css, _, _, _ in (*_SCALE_DIMS, *_SPIN_DIMS, *_IN_DIMS):
                if cb_id == f"vary-{css}":
                    self.query_one(f"#{css}-min").disabled = not active
                    self.query_one(f"#{css}-max").disabled = not active
                    break

    @on(Select.Changed)
    def _select_gate(self, event: Select.Changed) -> None:
        if event.select.id == "use-shape-crop-fixed":
            val = "" if event.value is Select.BLANK else str(event.value)
            self.query_one("#shape-file").disabled = (val != "true")

    # np_apophis and np_apophis_list are mutually exclusive: filling one disables the other.
    @on(Input.Changed, "#np-apophis")
    def _np_apophis_changed(self, event: Input.Changed) -> None:
        self.query_one("#np-apophis-list").disabled = bool(event.value.strip())

    @on(Input.Changed, "#np-apophis-list")
    def _np_apophis_list_changed(self, event: Input.Changed) -> None:
        self.query_one("#np-apophis").disabled = bool(event.value.strip())

    # ── button / action handlers ───────────────────────────────────────────────

    @on(Button.Pressed)
    def _on_button(self, event: Button.Pressed) -> None:
        bid = event.button.id
        if bid == "btn-run":
            self._launch(dry=False)
        elif bid == "btn-dryrun":
            self._launch(dry=True)
        elif bid == "btn-quit":
            self.exit(None)

    def action_run_sweep(self) -> None:
        self._launch(dry=False)

    def action_do_dry_run(self) -> None:
        self._launch(dry=True)

    def action_app_quit(self) -> None:
        self.exit(None)

    # ── internals ─────────────────────────────────────────────────────────────

    def _set_status(self, msg: str, *, error: bool = False) -> None:
        s = self.query_one("#status", Static)
        s.update(msg)
        s.remove_class("err", "ok")
        s.add_class("err" if error else "ok")

    def _launch(self, dry: bool) -> None:
        try:
            argv = self._build_argv(force_dry=dry)
        except ValueError as exc:
            self._set_status(f"Error: {exc}", error=True)
            return
        self.exit(argv)

    # ── query helpers ─────────────────────────────────────────────────────────

    def _iv(self, wid: str) -> str:
        return self.query_one(f"#{wid}", Input).value.strip()

    def _cb(self, wid: str) -> bool:
        return self.query_one(f"#{wid}", Checkbox).value

    def _sel(self, wid: str) -> str:
        val = self.query_one(f"#{wid}", Select).value
        return "" if val is Select.BLANK else str(val)

    # ── argv builder ──────────────────────────────────────────────────────────

    def _build_argv(self, force_dry: bool = False) -> List[str]:  # noqa: C901
        argv: List[str] = []

        def si(wid: str, flag: str) -> None:
            """Emit string flag if non-empty."""
            v = self._iv(wid)
            if v:
                argv.extend([flag, v])

        def ii(wid: str, flag: str) -> None:
            """Emit integer flag if non-empty; raise on bad input."""
            v = self._iv(wid)
            if not v:
                return
            try:
                int(v)
            except ValueError:
                raise ValueError(f"{flag} must be an integer (got '{v}')")
            argv.extend([flag, v])

        def fi(wid: str, flag: str) -> None:
            """Emit float flag if non-empty; raise on bad input."""
            v = self._iv(wid)
            if not v:
                return
            try:
                float(v)
            except ValueError:
                raise ValueError(f"{flag} must be a number (got '{v}')")
            argv.extend([flag, v])

        # Paths
        si("prefix",      "--prefix")
        si("base-dir",    "--base-dir")
        si("eph-cache",   "--ephemeris-cache-dir")
        si("phantom-dir", "--phantom-dir")
        si("output-root", "--output-root")

        # Sampling
        ii("num-samples", "--num-samples")
        ii("seed",        "--seed")
        sn_raw = self._iv("saltelli-n")
        if sn_raw:
            try:
                sn = int(sn_raw)
            except ValueError:
                raise ValueError(f"--saltelli-n must be an integer (got '{sn_raw}')")
            if sn < 2:
                raise ValueError("--saltelli-n must be ≥ 2")
            argv.extend(["--saltelli-n", str(sn)])
        if self._cb("saltelli-2nd"):
            argv.append("--saltelli-calc-second-order")

        # Mass
        if self._cb("vary-mass"):
            fi("mass-min-kg", "--mass-min-kg")
            fi("mass-max-kg", "--mass-max-kg")
            fi("apophis-ref-mass-kg", "--apophis-ref-mass-kg")

        # Scale bounds
        for css, attr, _, is_int in _SCALE_DIMS:
            if self._cb(f"vary-{css}"):
                flag_base = "--" + attr.replace("_", "-")
                if is_int:
                    ii(f"{css}-min", f"{flag_base}-min")
                    ii(f"{css}-max", f"{flag_base}-max")
                else:
                    fi(f"{css}-min", f"{flag_base}-min")
                    fi(f"{css}-max", f"{flag_base}-max")

        # Spin
        for css, attr, _, is_int in _SPIN_DIMS:
            if self._cb(f"vary-{css}"):
                flag_base = "--" + attr.replace("_", "-")
                if is_int:
                    ii(f"{css}-min", f"{flag_base}-min")
                    ii(f"{css}-max", f"{flag_base}-max")
                else:
                    fi(f"{css}-min", f"{flag_base}-min")
                    fi(f"{css}-max", f"{flag_base}-max")

        # DEM contact parameters
        for css, attr, _, is_int in _IN_DIMS:
            if self._cb(f"vary-{css}"):
                flag_base = "--" + attr.replace("_", "-")
                if is_int:
                    ii(f"{css}-min", f"{flag_base}-min")
                    ii(f"{css}-max", f"{flag_base}-max")
                else:
                    fi(f"{css}-min", f"{flag_base}-min")
                    fi(f"{css}-max", f"{flag_base}-max")

        # Timeframe
        for t_css, t_flag in (("tmax", "--tmax-hours"), ("dtmax", "--dtmax-hours")):
            if self._cb(f"vary-{t_css}"):
                fi(f"{t_css}-min", f"{t_flag}-min")
                fi(f"{t_css}-max", f"{t_flag}-max")
            else:
                fi(f"{t_css}-fixed", t_flag)

        # Setup toggles
        if self._cb("vary-use-dem"):
            argv.append("--vary-use-dem")
        else:
            df = self._sel("use-dem-fixed")
            if df:
                argv.extend(["--use-dem-fixed", df])
        # np_apophis and np_apophis_list are mutually exclusive.
        np_list_raw = self._iv("np-apophis-list")
        if np_list_raw:
            try:
                np_list_vals = [int(x) for x in np_list_raw.split()]
            except ValueError:
                raise ValueError(f"np_apophis_list must be space-separated integers (got '{np_list_raw}')")
            argv.extend(["--np-apophis-list"] + [str(v) for v in np_list_vals])
        else:
            ii("np-apophis", "--np-apophis")
        if self._cb("vary-use-shape-crop"):
            argv.append("--vary-use-shape-crop")
            si("shape-file", "--shape-file")
        else:
            ocf = self._sel("use-shape-crop-fixed")
            if ocf:
                argv.extend(["--use-shape-crop-fixed", ocf])
                if ocf == "true":
                    si("shape-file", "--shape-file")
        if self._cb("vary-apophis-only"):
            argv.append("--vary-apophis-only")
        else:
            aof = self._sel("apophis-only-fixed")
            if aof:
                argv.extend(["--apophis-only-fixed", aof])

        spf = self._iv("spin-period-fixed")
        if spf:
            argv.extend(["--spin-period-fixed", spf])

        sep_list_raw = self._iv("scale-earth-sep-list")
        if sep_list_raw:
            try:
                sep_vals = [float(x) for x in sep_list_raw.split()]
            except ValueError:
                raise ValueError(
                    f"scale_earth_sep_list must be space-separated floats (got '{sep_list_raw}')"
                )
            argv.extend(["--scale-earth-sep-list"] + [str(v) for v in sep_vals])

        # Sinks
        ii("sink-earth",   "--sink-earth-id")
        ii("sink-apophis", "--sink-apophis-id")

        # Batch naming
        bl = self._iv("batch-label")
        if bl:
            argv.extend(["--batch-label", bl])
        ii("slug-max-len", "--batch-slug-max-len")

        # Execution
        if force_dry or self._cb("dry-run"):
            argv.append("--dry-run")
        ii("jobs", "--jobs")

        return argv


# ── entry point ───────────────────────────────────────────────────────────────

def main() -> None:
    cli_argv = [a for a in sys.argv[1:] if a not in ("--interactive", "-i", "--tui")]
    initial_args = _runner.build_parser().parse_args(cli_argv)

    app = SobolTUIApp(initial_args)
    result = app.run()

    if result is None:
        return  # user quit without running

    runner = Path(__file__).parent / "run_mass_sobol_phantom.py"
    cmd = [sys.executable, str(runner), *result]
    print("Running:", " ".join(cmd), flush=True)
    sys.exit(subprocess.call(cmd))


if __name__ == "__main__":
    main()
