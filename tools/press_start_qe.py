# tools/press_start_qe.py
from __future__ import annotations
import argparse
import subprocess
import sys
from pathlib import Path
from typing import Optional, List

ROOT = Path(__file__).resolve().parents[1]

def find_script(name: str, extra_dirs: Optional[List[Path]] = None) -> Path:
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
    raise FileNotFoundError(f"Could not find {name} in expected locations (root/tools/scripts).")

def run(cmd: list[str], cwd: Optional[Path] = None) -> None:
    here = cwd if cwd is not None else ROOT
    print(f"[run] {sys.executable} {' '.join(cmd[1:])}")
    subprocess.check_call(cmd, cwd=str(here))

def main() -> None:
    ap = argparse.ArgumentParser(description="Generate + analyze a QE validation dataset.")
    ap.add_argument("--out", default=str(ROOT / "runs" / "qe_demo"),
                    help="Output directory for analysis results (default: runs/qe_demo)")
    ap.add_argument("--events", type=int, default=50000, help="Total events to simulate.")
    ap.add_argument("--gamma", type=float, default=0.03, help="Advanced-wave injection strength.")
    ap.add_argument("--seed", type=int, default=42, help="RNG seed.")
    ap.add_argument("--delay_short_ns", type=float, default=0.0, help="Short-path delay (ns).")
    ap.add_argument("--delay_long_ns", type=float, default=1.67, help="Long-path delay (ns).")
    ap.add_argument("--jitter_ps", type=float, default=50.0, help="Timing jitter (ps).")
    ap.add_argument("--coin_p", type=float, default=0.30, help="Coincidence probability.")
    ap.add_argument("--runs", type=int, default=100, help="Number of runs (blocks).")
    ap.add_argument("--idler", default="D3", help="Idler detector to analyze (default: D3).")
    ap.add_argument("--bin_ns", type=float, default=0.05, help="Histogram bin width (ns).")
    ap.add_argument("--range_ns", type=float, default=5.0, help="Histogram half-range (ns).")
    ap.add_argument("--expected_ns", type=float, default=-1.67, help="Expected advanced timing (ns).")
    args = ap.parse_args()

    out_dir = Path(args.out).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    data_dir = ROOT / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    csv_path = data_dir / "qe_demo.csv"

    gen = find_script("gen_qe_dataset.py")
    ana = find_script("analyze_qe_dataset.py")

    # ---- Generate dataset ----
    run([
        sys.executable, str(gen),
        "--out", str(csv_path),
        "--events", str(args.events),
        "--gamma", str(args.gamma),
        "--seed", str(args.seed),
        "--delay_short_ns", str(args.delay_short_ns),
        "--delay_long_ns", str(args.delay_long_ns),
        "--jitter_ps", str(args.jitter_ps),
        "--coin_p", str(args.coin_p),
        "--runs", str(args.runs),
    ])

    # ---- Analyze dataset ----
    run([
        sys.executable, str(ana),
        "--data", str(csv_path),
        "--out", str(out_dir),
        "--idler", str(args.idler),
        "--bin_ns", str(args.bin_ns),
        "--range_ns", str(args.range_ns),
        "--expected_ns", str(args.expected_ns),
    ])

    print(f"[done] QE validation assets are in {out_dir}")

if __name__ == "__main__":
    main()
