#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
One-button runner with calibration checking and optional default calibration.

Usage examples (from project root):
  python tools\press_start.py --assume-default-cal
  python tools\press_start.py --strict-cal         (requires real cal CSVs to pass)
  python tools\press_start.py --force-sim          (ignore any real data CSVs)
"""

from __future__ import annotations
import argparse, csv, datetime as dt, hashlib, json, math, os, shutil, subprocess, sys
from pathlib import Path
from statistics import mean
from typing import Dict, List, Tuple, Optional

ROOT = Path(__file__).resolve().parents[1]
PY = sys.executable

# ----------------- util -----------------
def ensure_dirs(paths: List[Path]) -> None:
    for p in paths:
        p.mkdir(parents=True, exist_ok=True)

def sha256sum(p: Path) -> str:
    h = hashlib.sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda: f.read(65536), b''):
            h.update(chunk)
    return h.hexdigest()

def copy_if_exists(src: Path, dst: Path) -> Optional[Path]:
    if src.exists():
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        return dst
    return None

def read_csv_rows(p: Path) -> List[Dict[str, str]]:
    if not p.exists():
        return []
    with p.open('r', encoding='utf-8', errors='ignore') as f:
        sample = f.read(2048)
        f.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=",;")
        except csv.Error:
            dialect = csv.excel
        reader = csv.DictReader(f, dialect=dialect)
        return [row for row in reader if row and any((v or "").strip() for v in row.values())]

def write_csv(p: Path, header: List[str], rows: List[List[object]]) -> None:
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open('w', encoding='utf-8', newline='') as f:
        w = csv.writer(f)
        w.writerow(header)
        for r in rows:
            w.writerow(r)

def log(msg: str) -> None:
    print(msg, flush=True)

# ----------------- bootstrap -----------------
def bootstrap_folders() -> None:
    ensure_dirs([
        ROOT / "docs",
        ROOT / "cal",
        ROOT / "runs" / "cable",
        ROOT / "runs" / "optical",
        ROOT / "data" / "cable",
        ROOT / "data" / "optical",
        ROOT / "paper_assets",
        ROOT / "tools",
    ])
    prereg = ROOT / "docs" / "prereg_v1.0.md"
    if not prereg.exists():
        prereg.write_text(
            "# Preregistration (v1.0)\n\n"
            "- Decision rule: 95% HDI on |gamma|; ROPE = +/- 2*sigma_hat.\n"
            "- Smoking-gun: length scaling, null dependence, cross-platform.\n"
            "- Normalization: cable vs 1.0 Vpp; optical vs integrated mainlobe.\n",
            encoding="utf-8",
        )

    # Create empty cal CSVs with headers if missing (non-blocking)
    cal_dir = ROOT / "cal"
    templates = {
        "cable_tau.csv":              ["tau_ns"],
        "optical_tau.csv":            ["tau_ns"],
        "cable_null.csv":             ["state","depth_dB"],
        "optical_null.csv":           ["state","depth_dB"],
        "cable_sigma_vs_N.csv":       ["N","sigma"],
        "optical_sigma_vs_N.csv":     ["N","sigma"],
        "optical_dark.csv":           ["rate_Hz"],
    }
    for name, header in templates.items():
        p = cal_dir / name
        if not p.exists():
            write_csv(p, header, [])

# ----------------- freeze analysis -----------------
def freeze_analysis() -> Path:
    stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    out = ROOT / "docs" / f"analysis_freeze_{stamp}"
    ensure_dirs([out])
    files = [
        ROOT / "final_analysis.py",
        ROOT / "final_analysis_optical.py",
        ROOT / "tools" / "press_start.py",
    ]
    manifest = []
    for f in files:
        if f.exists():
            shutil.copy2(f, out / f.name)
            manifest.append({"name": f.name, "sha256": sha256sum(f)})
    # snapshot current cal CSVs
    cal_src = ROOT / "cal"
    cal_dst = out / "cal"
    if cal_src.exists():
        ensure_dirs([cal_dst])
        for csvf in cal_src.glob("*.csv"):
            shutil.copy2(csvf, cal_dst / csvf.name)
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    with (out / "checksums.txt").open("w", encoding="utf-8") as f:
        for m in manifest:
            f.write(f"{m['sha256']}  {m['name']}\n")
    return out

# ----------------- default calibration -----------------
def ensure_default_calibration(cal_dir: Path,
                               tau_nom_cable_ns: float,
                               tau_nom_optical_ns: float,
                               sigma0_cable_uV: float = 500.0,
                               sigma0_opt_counts: float = 1.0) -> List[Path]:
    """
    If any required cal CSV is empty/missing, fill with nominal defaults:
      - tau: nominal ns
      - null depths: A=-40 dB, B=-20 dB
      - sigma vs N: N in [1,4,16,64,256,1024], sigma0/sqrt(N)
      - dark: flat 100 Hz x3
    Returns list of files created/filled.
    """
    created: List[Path] = []
    # tau
    for name, tau in [("cable_tau.csv", tau_nom_cable_ns),
                      ("optical_tau.csv", tau_nom_optical_ns)]:
        p = cal_dir / name
        rows = read_csv_rows(p)
        if not rows:
            write_csv(p, ["tau_ns"], [[f"{tau:.3f}"]])
            created.append(p)

    # null depths
    for name in ("cable_null.csv", "optical_null.csv"):
        p = cal_dir / name
        rows = read_csv_rows(p)
        if not rows:
            write_csv(p, ["state","depth_dB"], [["A","-40.0"],["B","-20.0"]])
            created.append(p)

    # sigma vs N
    Ns = [1,4,16,64,256,1024]
    for name, sigma0 in [("cable_sigma_vs_N.csv", sigma0_cable_uV),
                         ("optical_sigma_vs_N.csv", sigma0_opt_counts)]:
        p = cal_dir / name
        rows = read_csv_rows(p)
        if not rows:
            table = [[n, f"{sigma0/(n**0.5):.6g}"] for n in Ns]
            write_csv(p, ["N","sigma"], table)
            created.append(p)

    # dark
    p_dark = cal_dir / "optical_dark.csv"
    rows = read_csv_rows(p_dark)
    if not rows:
        write_csv(p_dark, ["rate_Hz"], [[100],[100],[100]])
        created.append(p_dark)

    if created:
        log("[cal] default calibration written for: " + ", ".join(x.name for x in created))
    else:
        log("[cal] default calibration not needed (all present).")
    return created

# ----------------- cal checks -----------------
def regress_loglog(xs: List[float], ys: List[float]) -> Tuple[float, float, float]:
    n = len(xs)
    if n < 2:
        return float('nan'), float('nan'), float('nan')
    xbar, ybar = mean(xs), mean(ys)
    sxx = sum((x - xbar) ** 2 for x in xs)
    sxy = sum((x - xbar) * (y - ybar) for x, y in zip(xs, ys))
    if sxx == 0:
        return float('nan'), float('nan'), float('nan')
    b = sxy / sxx
    a = ybar - b * xbar
    sst = sum((y - ybar) ** 2 for y in ys)
    ssr = sum((a + b * x - y) ** 2 for x, y in zip(xs, ys))
    r2 = 1.0 - (ssr / sst if sst > 0 else float('nan'))
    return b, a, r2

def parse_float_list(rows: List[Dict[str, str]], key: str) -> List[float]:
    out = []
    for r in rows:
        v = (r.get(key) or "").strip()
        if v:
            try:
                out.append(float(v))
            except ValueError:
                pass
    return out

def check_calibration(cal_dir: Path,
                      tau_nom_cable: float,
                      tau_nom_optical: float,
                      strict: bool = False) -> Dict[str, Dict[str, str]]:
    report: Dict[str, Dict[str, str]] = {}

    def add(name: str, status: str, message: str):
        report[name] = {"status": status, "message": message}

    # tau
    for label, fname, tau_nom in [
        ("cable_tau", "cable_tau.csv", tau_nom_cable),
        ("optical_tau", "optical_tau.csv", tau_nom_optical),
    ]:
        rows = read_csv_rows(cal_dir / fname)
        taus = parse_float_list(rows, "tau_ns")
        if not taus:
            add(label, "MISSING", f"{fname} not found or empty")
        else:
            t = mean(taus)
            tol = 0.02 * tau_nom
            status = "PASS" if abs(t - tau_nom) <= tol else "WARN"
            add(label, status, f"mean tau={t:.3f} ns vs nominal {tau_nom:.3f} ns (±2%)")

    # null depths
    for label, fname in [("cable_null","cable_null.csv"),
                         ("optical_null","optical_null.csv")]:
        rows = read_csv_rows(cal_dir / fname)
        if not rows:
            add(label, "MISSING", f"{fname} not found or empty")
        else:
            dep = {}
            for r in rows:
                s = (r.get("state") or "").strip().upper()
                try:
                    d = float((r.get("depth_dB") or "").strip())
                except ValueError:
                    continue
                if s:
                    dep[s] = d
            if "A" in dep and "B" in dep:
                a_ok = dep["A"] <= -40.0
                b_ok = abs(dep["B"] - (-20.0)) <= 3.0
                status = "PASS" if (a_ok and b_ok) else "WARN"
                add(label, status, f"A={dep['A']:.1f} dB (<=-40), B={dep['B']:.1f} dB (-20±3)")
            else:
                add(label, "WARN", f"need rows for states A and B in {fname}")

    # sigma vs N
    for label, fname in [("cable_sigma_vs_N","cable_sigma_vs_N.csv"),
                         ("optical_sigma_vs_N","optical_sigma_vs_N.csv")]:
        rows = read_csv_rows(cal_dir / fname)
        if not rows or sum(1 for r in rows if (r.get("N") and r.get("sigma"))) < 3:
            add(label, "MISSING", f"{fname} not found or too few points")
        else:
            try:
                Ns = [float(r["N"]) for r in rows]
                sigs= [float(r["sigma"]) for r in rows]
                xs = [math.log10(n) for n in Ns]
                ys = [math.log10(s) for s in sigs]
                slope, intercept, r2 = regress_loglog(xs, ys)
                ok = (abs(slope + 0.5) <= 0.1) and (r2 >= 0.95)
                add(label, "PASS" if ok else "WARN",
                    f"log10(sigma) vs log10(N): slope={slope:.3f} (target -0.5±0.1), R2={r2:.3f} (>=0.95)")
            except Exception:
                add(label, "WARN", f"could not parse {fname}")

    # optical dark stability
    rows = read_csv_rows(cal_dir / "optical_dark.csv")
    if not rows:
        add("optical_dark", "MISSING", "optical_dark.csv not found or empty")
    else:
        rates = parse_float_list(rows, "rate_Hz")
        if len(rates) >= 2:
            mu = mean(rates)
            span = (max(rates) - min(rates)) / (mu if mu else 1.0)
            status = "PASS" if span <= 0.20 else "WARN"
            add("optical_dark", status, f"span={span*100:.1f}% (<=20%)")
        else:
            add("optical_dark", "WARN", "need >=2 measurements")

    # write report
    pa = ROOT / "paper_assets"
    ensure_dirs([pa, ROOT / "runs"])
    (ROOT / "runs" / "cal_summary.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    lines = ["# Calibration report\n"]
    for k, v in report.items():
        lines.append(f"- {k}: {v['status']} — {v['message']}")
    (pa / "cal_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    if strict:
        bad = [k for k, v in report.items() if v["status"] in {"MISSING","WARN"}]
        if bad:
            raise SystemExit(f"[cal] STRICT mode: failing due to {', '.join(bad)}")

    return report

# ----------------- run analyses -----------------
def has_real_data(path: Path) -> Optional[Path]:
    for p in path.glob("*.csv"):
        return p
    return None

def run_cable(out_dir: Path, vpp: float, seed: int, force_sim: bool = False) -> None:
    script = ROOT / "final_analysis.py"
    cmd = [PY, str(script), "--out", str(out_dir)]
    data_csv = None if force_sim else has_real_data(ROOT / "data" / "cable")
    if data_csv:
        cmd += ["--data", str(data_csv), "--Vpp", f"{vpp:.6g}"]
        log(f"[cable] analyzing real: {data_csv.name}")
    else:
        cmd += ["--simulate", "--Vpp", f"{vpp:.6g}", "--seed", str(seed)]
        log("[cable] simulating data…")
    subprocess.check_call(cmd, cwd=ROOT)

def run_optical(out_dir: Path, seed: int, force_sim: bool = False) -> None:
    script = ROOT / "final_analysis_optical.py"
    cmd = [PY, str(script), "--out", str(out_dir)]
    data_csv = None if force_sim else has_real_data(ROOT / "data" / "optical")
    if data_csv:
        cmd += ["--data", str(data_csv), "--Vref", "mainlobe"]
        log(f"[optical] analyzing real: {data_csv.name}")
    else:
        cmd += ["--simulate", "--seed", str(seed), "--Vref", "mainlobe"]
        log("[optical] simulating data…")
    subprocess.check_call(cmd, cwd=ROOT)

# ----------------- merge + assets -----------------
def get_first(d: Dict[str,str], keys: List[str], default: str="") -> str:
    for k in keys:
        if k in d and d[k] != "":
            return d[k]
    return default

def read_csv_dict(p: Path) -> Dict[str, str]:
    if not p.exists():
        return {}
    rows = list(csv.DictReader(p.open("r", encoding="utf-8", errors="ignore")))
    return rows[0] if rows else {}

def merge_and_copy_assets() -> None:
    runs = ROOT / "runs"
    pa = ROOT / "paper_assets"
    ensure_dirs([pa])

    d1 = read_csv_dict(runs / "cable" / "results_table.csv")
    d2 = read_csv_dict(runs / "optical" / "results_table_optical.csv")

    merged = runs / "results_table_merged.csv"
    headers = ["platform","bound_abs","bound_rel","bound_dBc","prewindow_rms","rope_pm","HDI_low","HDI_high"]
    with merged.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(headers)
        if d1:
            w.writerow([
                "cable",
                get_first(d1, ["bound_abs_uV"]), d1.get("bound_rel",""), d1.get("bound_dBc",""),
                d1.get("prewindow_rms_uV",""), d1.get("ROPE_pm_uV",""),
                d1.get("HDI_low_uV",""), d1.get("HDI_high_uV",""),
            ])
        if d2:
            # be tolerant to alternate column names
            w.writerow([
                "optical",
                get_first(d2, ["bound_abs_counts","bound_abs_counts_per_bin"]),
                d2.get("bound_rel",""), d2.get("bound_dBc",""),
                get_first(d2, ["prewindow_rms_counts","prewindow_rms_counts_per_bin"]),
                get_first(d2, ["ROPE_pm_counts","ROPE_pm_counts_per_bin"]),
                get_first(d2, ["HDI_low_counts","HDI_low_counts_per_bin"]),
                get_first(d2, ["HDI_high_counts","HDI_high_counts_per_bin"]),
            ])

    # ASCII-only summary (avoid greek/µ for Windows encodings)
    lines = ["# Results summary\n"]
    if d1:
        try:
            b_abs = float(get_first(d1, ["bound_abs_uV"], "0"))
            b_dbc = float(d1.get("bound_dBc","0"))
            lines.append(f"- Cable: |gamma_adv| < {b_abs:.1f} uV ({b_dbc:.1f} dBc).")
        except Exception:
            lines.append("- Cable: see runs/cable/results_table.csv")
    if d2:
        try:
            b_abs = float(get_first(d2, ["bound_abs_counts","bound_abs_counts_per_bin"], "0"))
            b_dbc = float(d2.get("bound_dBc","0"))
            lines.append(f"- Optical: |gamma_adv| < {b_abs:.2f} counts/bin ({b_dbc:.1f} dB).")
        except Exception:
            lines.append("- Optical: see runs/optical/results_table_optical.csv")
    (pa / "results_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    # plots
    copy_if_exists(runs / "cable" / "posterior_gamma_plot.png", pa / "fig_cable_posterior.png")
    copy_if_exists(runs / "optical" / "posterior_gamma_plot_optical.png", pa / "fig_optical_posterior.png")
    # fallback name some older scripts used
    copy_if_exists(runs / "optical" / "posterior_gamma_plot.png", pa / "fig_optical_posterior.png")

# ----------------- main -----------------
def main():
    ap = argparse.ArgumentParser(description="One-button runner with calibration checks and defaults.")
    ap.add_argument("--strict-cal", action="store_true",
                    help="fail if calibration checks WARN/MISSING")
    ap.add_argument("--assume-default-cal", action="store_true",
                    help="if cal CSVs are missing/empty, write nominal defaults first")
    ap.add_argument("--force-sim", action="store_true",
                    help="ignore any real data; always simulate")
    ap.add_argument("--Vpp", type=float, default=1.0,
                    help="drive amplitude for cable sim (Vpp)")
    ap.add_argument("--seed-cable", type=int, default=1234)
    ap.add_argument("--seed-optical", type=int, default=2025)
    ap.add_argument("--tau-nom-cable", type=float, default=4.817,
                    help="nominal tau (ns) for cable platform")
    ap.add_argument("--tau-nom-optical", type=float, default=3.336,
                    help="nominal tau (ns) for optical platform")
    ap.add_argument("--sigma0-cable-uV", type=float, default=500.0,
                    help="default single-shot sigma for cable (uV) when auto-cal")
    ap.add_argument("--sigma0-optical-counts", type=float, default=1.0,
                    help="default single-shot sigma for optical (counts/bin) when auto-cal")
    args = ap.parse_args()

    log("[setup] ensuring folder structure…")
    bootstrap_folders()

    if args.assume_default_cal:
        ensure_default_calibration(
            ROOT / "cal",
            tau_nom_cable_ns=args.tau_nom_cable,
            tau_nom_optical_ns=args.tau_nom_optical,
            sigma0_cable_uV=args.sigma0_cable_uV,
            sigma0_opt_counts=args.sigma0_optical_counts,
        )

    log("[freeze] snapshotting analysis…")
    freeze_dir = freeze_analysis()
    log(f"[freeze] -> {freeze_dir.name}")

    log("[cal] running calibration checks…")
    try:
        check_calibration(
            ROOT / "cal",
            tau_nom_cable=args.tau_nom_cable,
            tau_nom_optical=args.tau_nom_optical,
            strict=args.strict_cal,
        )
    except SystemExit as e:
        log(str(e))
        sys.exit(2)

    log("[run] cable platform…")
    run_cable(ROOT / "runs" / "cable", vpp=args.Vpp, seed=args.seed_cable, force_sim=args.force_sim)
    log("[run] optical platform…")
    run_optical(ROOT / "runs" / "optical", seed=args.seed_optical, force_sim=args.force_sim)

    log("[merge] collecting outputs…")
    merge_and_copy_assets()
    log("[done] results in paper_assets/ ; merged CSV in runs/results_table_merged.csv ; cal report in paper_assets/cal_report.md")

if __name__ == "__main__":
    main()
