#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
One-button runner for QE validation:
  - generates synthetic data (with chosen gamma)
  - runs analysis
  - prints a compact summary
"""
import argparse, subprocess, sys
from pathlib import Path

def run(cmd):
    print("[run]", " ".join(cmd))
    subprocess.check_call(cmd)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=str, default='runs/qe_demo')
    ap.add_argument('--data', type=str, default='data/qe_demo.csv')
    ap.add_argument('--gamma', type=float, default=0.03)
    ap.add_argument('--seed', type=int, default=42)
    ap.add_argument('--delay_long_ns', type=float, default=1.67)
    ap.add_argument('--bin_ns', type=float, default=0.05)
    ap.add_argument('--range_ns', type=float, default=5.0)
    ap.add_argument('--idler', type=str, default='D3')
    args = ap.parse_args()

    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)

    # 1) generate
    run([sys.executable, 'gen_qe_dataset.py',
         '--out', args.data,
         '--events', '50000',
         '--gamma', str(args.gamma),
         '--seed', str(args.seed),
         '--delay_short_ns', '0.0',
         '--delay_long_ns', str(args.delay_long_ns),
         '--jitter_ps', '50.0',
         '--coin_p', '0.30',
         '--runs', '100'])

    # 2) analyze (expected negative lag = -delay_long)
    run([sys.executable, 'analyze_qe_dataset.py',
         '--data', args.data,
         '--out', args.out,
         '--idler', args.idler,
         '--bin_ns', str(args.bin_ns),
         '--range_ns', str(args.range_ns),
         '--expected_ns', str(-args.delay_long_ns)])

    print(f"\n[done] QE validation assets are in {args.out}")

if __name__ == '__main__':
    main()
