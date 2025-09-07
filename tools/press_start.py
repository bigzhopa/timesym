# tools/press_start.py
from __future__ import annotations
import argparse
import subprocess
import sys
from pathlib import Path
from typing import Optional, List
import csv

ROOT = Path(__file__).resolve().parents[1]

def find_script(name: str, extra_dirs: Optional[List[Path]] = None) -> Path:
    """
    Find a script by name searching common locations:
    repo root, tools/, scripts/, and any extra_dirs provided.
    """
    candidates = [
        ROOT / name,
        ROOT / "tools" / name,
        ROOT / "scripts" / name,
    ]
    if extra_dirs:
        candidates.extend([d / name for d in extra_dirs])
    for p in candidates:
        if p.exists():
            return p
    raise FileNotFoundError(
        f"Could not find {name} in " +
        ", ".join({str(p.parent) for p in candidates})
    )

def run(cmd: List[str], cwd: Optional[Path] = None) -> None:
    here = cwd if cwd is not None else ROOT
    print(f"[run] {sys.executable} {' '.join(cmd[1:])}")
    subprocess.check_call(cmd, cwd=str(here))

def find_result_csvs(search_dir: Path) -> List[Path]:
    # look for results_table*.csv anywhere under search_dir
    return sorted(search_dir.rglob("results_table*.csv"))

def merge_results(cable_csv: Optional[Path], optical_csv: Optional[Path], out_path: Path) -> None:
    """
    Merge two results_table*.csv files into one CSV with a 'platform' column.
    This does a union of headers so it's resilient if columns differ slightly.
    """
    rows = []
    header_union: List[str] = []

    def load(path: Path, platform: str):
        nonlocal header_union
        if not path or not path.exists():
            return
        with path.open(newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                row = dict(row)  # copy
                row["platform"] = platform
                rows.append(row)
                # grow header union
                for k in row.keys():
                    if k not in header_union:
                        header_union.append(k)

    load(cable_csv, "cable")
    load(optical_csv, "optical")

    if not rows:
        print("[merge] no input rows found; skipping merged CSV.")
        return

    # Ensure 'platform' is first column for readability
    if "platform" in header_union:
        header_union = ["platform"] + [h for h in header_union if h != "platform"]

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=header_union)
        writer.writeheader()
        for r in rows:
            # fill missing keys with ""
            for h in header_union:
                r.setdefault(h, "")
            writer.writerow(r)
    print(f"[merge] wrote {out_path}")

def main() -> None:
    ap = argparse.ArgumentParser(description="One-button RF+Optical run (simulate + inline-merge).")
    ap.add_argument("--out", default=str(ROOT / "runs" / "demo"),
                    help="Output directory root (default: runs/demo)")
    # knobs (defaults match your prior runs/logs)
    ap.add_argument("--strict-cal", action="store_true", help="Require calibration assets (if applicable).")
    ap.add_argument("--assume-default-cal", action="store_true", help="Proceed with default/no calibration.")
    ap.add_argument("--force-sim", action="store_true", help="Force simulation even if data present.")
    ap.add_argument("--Vpp", type=float, default=1.0, help="Reference Vpp for RF dBc normalization.")
    ap.add_argument("--seed-cable", type=int, default=1234, help="Seed for RF simulation.")
    ap.add_argument("--seed-optical", type=int, default=2025, help="Seed for optical simulation.")
    ap.add_argument("--tau-nom-cable", type=float, default=None, help="Optional RF nominal tau (ns).")
    ap.add_argument("--tau-nom-optical", type=float, default=None, help="Optional optical nominal tau (ns).")
    ap.add_argument("--sigma0-cable-uV", type=float, default=None, help="Optional RF noise override (uV).")
    ap.add_argument("--sigma0-optical-counts", type=float, default=None, help="Optional optical noise override (counts/bin).")
    ap.add_argument("--skip-merge", action="store_true", help="Skip inline merge step.")
    args = ap.parse_args()

    out_root = Path(args.out).resolve()
    out_root.mkdir(parents=True, exist_ok=True)

    # Per-platform subdirs under the chosen --out
    cable_out = out_root / "cable"
    optical_out = out_root / "optical"
    cable_out.mkdir(parents=True, exist_ok=True)
    optical_out.mkdir(parents=True, exist_ok=True)

    # Locate scripts regardless of layout
    fa_cable = find_script("final_analysis.py")
    fa_opt   = find_script("final_analysis_optical.py")

    # ---- RF cable ----
    cmd_rf = [
        sys.executable, str(fa_cable),
        "--simulate",
        "--out", str(cable_out),
        "--Vpp", str(args.Vpp),
        "--seed", str(args.seed_cable),
    ]
    if args.tau_nom_cable is not None:
        cmd_rf += ["--tau_nom_ns", str(args.tau_nom_cable)]
    if args.sigma0_cable_uV is not None:
        cmd_rf += ["--sigma0_uV", str(args.sigma0_cable_uV)]

    run(cmd_rf)

    # ---- Optical ----
    cmd_opt = [
        sys.executable, str(fa_opt),
        "--simulate",
        "--out", str(optical_out),
        "--seed", str(args.seed_optical),
    ]
    if args.tau_nom_optical is not None:
        cmd_opt += ["--tau_nom_ns", str(args.tau_nom_optical)]
    if args.sigma0_optical_counts is not None:
        cmd_opt += ["--sigma0_counts", str(args.sigma0_optical_counts)]

    run(cmd_opt)

    # ---- Inline merge (robust, no external script) ----
    if not args.skip_merge:
        cable_csvs = find_result_csvs(cable_out)
        optical_csvs = find_result_csvs(optical_out)
        cable_csv = cable_csvs[0] if cable_csvs else None
        optical_csv = optical_csvs[0] if optical_csvs else None

        if not cable_csv and not optical_csv:
            print("[merge] no results_table*.csv found under cable/ or optical/; skipping merge.")
        else:
            merged_csv = out_root / "results_table_merged.csv"
            merge_results(cable_csv, optical_csv, merged_csv)
    else:
        print("[merge] skipped by --skip-merge")

    print(f"[done] RF → {cable_out}")
    print(f"[done] Optical → {optical_out}")
    print(f"[done] Root → {out_root}")

if __name__ == "__main__":
    main()
