#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
analyze_dcqe.py — DCQE-style re-analysis using the same HDI/ROPE pipeline.

Inputs (choose one):
  A) time-tags CSV with columns:
     t_signal_ns, t_idler_ns, [eraser_state or idler_det], [run_id], [deltaL_ns]
  B) histogram CSV with columns:
     bin_center_ns, hist_A_counts, hist_B_counts
     OR bin_center_ns, hist_D1, hist_D2, hist_D3, hist_D4  (+ CLI mapping)

Outputs:
  - results_table_dcqe.csv
  - posterior_gamma_plot_dcqe.png
  - (if built from time-tags) dcqe_histograms_built.csv

Usage examples:
  # time-tags with explicit eraser_state and a single deltaL
  python tools/analyze_dcqe.py --time-tags data/dcqe_tags.csv --deltaL-ns 2.5 --out runs/dcqe

  # time-tags + mapping from detectors to states
  python tools/analyze_dcqe.py --time-tags data/dcqe_tags.csv --deltaL-ns 2.5 \
      --detmap "D1:erase,D2:erase,D3:mark,D4:mark" --out runs/dcqe

  # four-detector histograms
  python tools/analyze_dcqe.py --hists data/dcqe_hist_4dets.csv --deltaL-ns 2.5 \
      --detmap "D1:erase,D2:erase,D3:mark,D4:mark" --out runs/dcqe

  # A/B histograms already provided
  python tools/analyze_dcqe.py --hists data/dcqe_hist_AB.csv --deltaL-ns 2.5 --out runs/dcqe
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Dict, Tuple, Optional

# ---------- helpers ----------
def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)

def parse_detmap(s: Optional[str]) -> Dict[str, str]:
    # "D1:erase,D2:erase,D3:mark,D4:mark"
    m = {}
    if not s:
        return m
    for tok in s.split(","):
        k, v = tok.split(":")
        m[k.strip()] = v.strip().lower()
    return m

def gaussian_center_and_sigma(x: np.ndarray, y: np.ndarray) -> Tuple[float,float]:
    # rough moments around peak; robust enough for our purpose
    if y.sum() <= 0:
        return 0.0, (x[1]-x[0]) * 3.0
    y = y - y.min()
    y_sum = y.sum()
    mu = (x * y).sum() / y_sum
    var = ((x - mu)**2 * y).sum() / y_sum
    sigma = np.sqrt(max(var, (x[1]-x[0])**2))
    return mu, sigma

def hdi_from_samples(abs_samples: np.ndarray, mass: float=0.95) -> Tuple[float,float]:
    z = np.sort(abs_samples)
    n = len(z)
    k = int(np.floor(mass * n))
    if k < 1:
        return (0.0, float(z[-1]) if n else 0.0)
    widths = z[k:] - z[:n-k]
    j = np.argmin(widths)
    return float(z[j]), float(z[j+k])

# ---------- build A/B ----------
def build_from_time_tags(csv_path: Path,
                         detmap: Dict[str,str],
                         dt_ns: float,
                         tmin_ns: float,
                         tmax_ns: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    df = pd.read_csv(csv_path)
    if "eraser_state" not in df.columns:
        if "idler_det" not in df.columns or not detmap:
            raise ValueError("Need eraser_state or idler_det+detmap for time-tags.")
        df["eraser_state"] = df["idler_det"].map(lambda d: detmap.get(str(d), "unknown"))
    # relative delay
    df["t_rel_ns"] = df["t_signal_ns"] - df["t_idler_ns"]
    # bins
    bins = np.arange(tmin_ns, tmax_ns + dt_ns/2, dt_ns)
    centers = (bins[:-1] + bins[1:]) / 2.0
    A = np.histogram(df.loc[df["eraser_state"]=="erase", "t_rel_ns"].values, bins=bins)[0].astype(float)
    B = np.histogram(df.loc[df["eraser_state"]=="mark",  "t_rel_ns"].values, bins=bins)[0].astype(float)
    # save reconstructed histograms for audit
    out = pd.DataFrame({"bin_center_ns": centers, "hist_A_counts": A, "hist_B_counts": B})
    out.to_csv(csv_path.parent / "dcqe_histograms_built.csv", index=False)
    return centers, A, B

def build_from_histograms(csv_path: Path,
                          detmap: Dict[str,str]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    df = pd.read_csv(csv_path)
    t = df["bin_center_ns"].values
    cols = set(df.columns)
    if {"hist_A_counts","hist_B_counts"}.issubset(cols):
        A = df["hist_A_counts"].astype(float).values
        B = df["hist_B_counts"].astype(float).values
    elif {"hist_D1","hist_D2","hist_D3","hist_D4"}.issubset(cols) and detmap:
        erase_cols = [c for c,s in detmap.items() if s=="erase" and c in cols]
        mark_cols  = [c for c,s in detmap.items() if s=="mark"  and c in cols]
        if not erase_cols or not mark_cols:
            raise ValueError("detmap must map some detectors to erase and mark.")
        A = df[erase_cols].sum(axis=1).astype(float).values
        B = df[mark_cols].sum(axis=1).astype(float).values
    else:
        raise ValueError("Histogram CSV must have hist_A/B or 4 detectors with detmap.")
    return t, A, B

# ---------- core analysis ----------
def analyze_dcqe(t_ns: np.ndarray,
                 A: np.ndarray,
                 B: np.ndarray,
                 deltaL_ns: float,
                 pre_start_ns: Optional[float]=None,
                 rope_mult: float=2.0,
                 prior_mult: float=10.0,
                 sigma_floor: float=1.0) -> Dict[str, float]:
    # align main lobe near zero using A+B
    S = A + B
    t0, sig = gaussian_center_and_sigma(t_ns, S)
    t_align = t_ns - t0
    # difference and pre-window
    D = A - B
    if pre_start_ns is None:
        pre_end = -5.0 * sig
        pre_start = t_align.min()
    else:
        pre_end = min(-5.0 * sig, -abs(deltaL_ns)/2.0)  # keep conservative if provided
        pre_start = pre_start_ns
    pre_mask = (t_align >= pre_start) & (t_align < pre_end)
    sigma_hat = float(np.std(D[pre_mask], ddof=1)) if np.any(pre_mask) else 0.0
    if sigma_hat == 0.0:
        sigma_hat = float(sigma_floor)

    # test point
    t_star = -abs(deltaL_ns)  # sign convention: precursor at -Δt
    # snap to nearest bin
    j = int(np.argmin(np.abs(t_align - t_star)))
    y = float(D[j])

    # Normal–Normal posterior for gamma
    tau0 = prior_mult * sigma_hat
    post_var = 1.0 / (1.0/tau0**2 + 1.0/sigma_hat**2)
    post_mu  = post_var * (y / sigma_hat**2)
    # sample posterior and compute |gamma|
    rng = np.random.default_rng(12345)
    samples = rng.normal(post_mu, np.sqrt(post_var), size=200000)
    abs_samp = np.abs(samples)
    hdi_lo, hdi_hi = hdi_from_samples(abs_samp, 0.95)
    rope_pm = rope_mult * sigma_hat
    p_in_rope = float(np.mean(abs_samp <= rope_pm))

    return {
        "t0_ns": float(t0),
        "sigma_main_ns": float(sig),
        "t_star_ns": float(t_star),
        "y_at_tstar": y,
        "sigma_hat": sigma_hat,
        "ROPE_pm": rope_pm,
        "P_in_ROPE": p_in_rope,
        "HDI_low": hdi_lo,
        "HDI_high": hdi_hi,
    }, (t_align, D, pre_mask, samples)

def plot_posterior(samples: np.ndarray, rope_pm: float, hdi: Tuple[float,float], out_png: Path) -> None:
    abs_samp = np.abs(samples)
    plt.figure(figsize=(7.0,4.0))
    plt.hist(abs_samp, bins=200, density=True, alpha=0.8)
    plt.axvline( rope_pm, ls="--")
    plt.axvline(-rope_pm, ls="--")
    plt.axvspan(hdi[0], hdi[1], alpha=0.2)
    plt.xlabel("|gamma_adv| (counts/bin)")
    plt.ylabel("Posterior density")
    plt.title("Posterior for |gamma_adv| (DCQE re-analysis)")
    plt.tight_layout()
    plt.savefig(out_png, dpi=180)
    plt.close()

# ---------- CLI ----------
def main():
    ap = argparse.ArgumentParser(description="DCQE re-analysis (A/B → D, HDI/ROPE bound).")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--time-tags", type=str, help="CSV with time-tag coincidences")
    g.add_argument("--hists", type=str, help="CSV with histograms (A/B or 4-detector)")

    ap.add_argument("--out", type=str, default="runs/dcqe", help="output folder")
    ap.add_argument("--detmap", type=str, default="", help='e.g. "D1:erase,D2:erase,D3:mark,D4:mark"')
    ap.add_argument("--deltaL-ns", type=float, required=True, help="idler path excess as time (ns)")

    # time-tags binning window (if using --time-tags)
    ap.add_argument("--dt-ns", type=float, default=0.05, help="bin width for time-tags (ns)")
    ap.add_argument("--tmin-ns", type=float, default=-50.0)
    ap.add_argument("--tmax-ns", type=float, default= +50.0)

    # analysis knobs (match your pipeline)
    ap.add_argument("--rope-mult", type=float, default=2.0)
    ap.add_argument("--prior-mult", type=float, default=10.0)
    ap.add_argument("--sigma-floor", type=float, default=1.0)
    args = ap.parse_args()

    out_dir = Path(args.out)
    ensure_dir(out_dir)

    detmap = parse_detmap(args.detmap)

    if args.time_tags:
        t, A, B = build_from_time_tags(Path(args.time_tags), detmap,
                                       dt_ns=args.dt_ns,
                                       tmin_ns=args.tmin_ns,
                                       tmax_ns=args.tmax_ns)
    else:
        t, A, B = build_from_histograms(Path(args.hists), detmap)

    res, aux = analyze_dcqe(t, A, B, deltaL_ns=args.deltaL_ns,
                            rope_mult=args.rope_mult,
                            prior_mult=args.prior_mult,
                            sigma_floor=args.sigma_floor)
    t_align, D, pre_mask, samples = aux

    # save results
    pd.DataFrame([res]).to_csv(out_dir / "results_table_dcqe.csv", index=False)
    (out_dir / "results_table_dcqe.json").write_text(json.dumps(res, indent=2), encoding="utf-8")

    # quick posterior plot
    plot_posterior(samples, res["ROPE_pm"], (res["HDI_low"], res["HDI_high"]),
                   out_png=out_dir / "posterior_gamma_plot_dcqe.png")

    # optional: save differenced histogram for sanity
    dfD = pd.DataFrame({"t_ns": t_align, "D_counts": D, "in_pre": pre_mask.astype(int)})
    dfD.to_csv(out_dir / "dcqe_difference_trace.csv", index=False)

    # console summary (ASCII-safe)
    print("\n--- DCQE re-analysis ---")
    print(f"Aligned main-lobe center t0       : {res['t0_ns']:.3f} ns")
    print(f"Test time t* = -Δt                 : {res['t_star_ns']:.3f} ns")
    print(f"y = D(t*)                          : {res['y_at_tstar']:.3f} counts/bin")
    print(f"sigma_hat (pre-window RMS)         : {res['sigma_hat']:.3f} counts/bin")
    print(f"95% HDI(|gamma|)                   : [{res['HDI_low']:.3f}, {res['HDI_high']:.3f}] counts/bin")
    print(f"ROPE ±                             : {res['ROPE_pm']:.3f}  (P in ROPE = {res['P_in_ROPE']*100:.1f}%)")
    print("Exports: results_table_dcqe.csv/json, posterior_gamma_plot_dcqe.png, dcqe_difference_trace.csv")
    print()

if __name__ == "__main__":
    main()
