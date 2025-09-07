#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Analyze synthetic QE time-tag data:
  1) form coincidences between D0 and a chosen idler detector (default D3)
  2) build Δt = t_idler - t_D0 histograms separately for 'erase' and 'mark'
  3) difference them: H_diff = H_erase - H_mark (choice-locked)
  4) Bayesian bound on the bin at expected advanced time (negative lag)

Saves:
  - coincidence_hist.png  (Erase, Mark, and their difference)
  - posterior_gamma_plot_qe.png
  - results_table_qe.csv
"""

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def create_dirs(p: Path):
    p.mkdir(parents=True, exist_ok=True)


def extract_deltas(df: pd.DataFrame,
                   idler: str,
                   choice: str,
                   window_ns: float = 5.0) -> np.ndarray:
    """
    Compute Δt = t_idler - t_D0 for coincidences within +/- window,
    run by run, using a two-pointer sweep for efficiency.
    """
    deltas = []
    for run, g in df[df['choice_type'] == choice].groupby('experiment_run'):
        t0 = g.loc[g['detector_id'] == 'D0', 'time_tag_ns'].to_numpy()
        ti = g.loc[g['detector_id'] == idler, 'time_tag_ns'].to_numpy()
        if len(t0) == 0 or len(ti) == 0:
            continue
        t0.sort(); ti.sort()

        i = j = 0
        n0, ni = len(t0), len(ti)
        while i < n0 and j < ni:
            dt = ti[j] - t0[i]
            if abs(dt) <= window_ns:
                deltas.append(dt)
                # advance the closer side to look for the next match
                if j+1 < ni and abs(ti[j+1] - t0[i]) < abs(dt):
                    j += 1
                else:
                    i += 1
            else:
                # shift pointers to reduce |dt|
                if dt > 0:
                    i += 1
                else:
                    j += 1
    return np.asarray(deltas, dtype=float)


def hdi_of_abs(samples: np.ndarray, cred=0.95):
    s = np.sort(np.abs(samples))
    n = len(s)
    k = int(np.floor(cred * n))
    widths = s[k:] - s[:n-k]
    j = np.argmin(widths)
    return s[j], s[j + k]


def normal_posterior(y, sigma, sigma0):
    """
    Conjugate Normal-Normal posterior for a single-bin amplitude γ.
    Prior: N(0, sigma0^2); Likelihood: N(y, sigma^2).
    Returns mean, std.
    """
    sigma2 = sigma**2
    sigma02 = sigma0**2
    post_var = 1.0 / (1.0/sigma02 + 1.0/sigma2)
    post_mean = post_var * (y / sigma2)
    post_std = np.sqrt(post_var)
    return post_mean, post_std


def main():
    ap = argparse.ArgumentParser(description="Analyze QE time-tags and bound advanced precursor.")
    ap.add_argument('--data', type=str, default='data/causally_disconnected_eraser_simulation.csv')
    ap.add_argument('--out', type=str, default='runs/qe')
    ap.add_argument('--idler', type=str, default='D3', choices=['D1', 'D2', 'D3', 'D4'])
    ap.add_argument('--bin_ns', type=float, default=0.05)
    ap.add_argument('--range_ns', type=float, default=5.0)
    ap.add_argument('--expected_ns', type=float, default=-1.67, help='expected advanced lag (negative)')
    ap.add_argument('--rope_mult', type=float, default=2.0)
    ap.add_argument('--prior_mult', type=float, default=10.0)
    ap.add_argument('--sigma_floor', type=float, default=0.30, help='min RMS per bin (counts)')
    args = ap.parse_args()

    out = Path(args.out)
    create_dirs(out)

    df = pd.read_csv(args.data)

    # Build Δt lists for erase/mark
    del_erase = extract_deltas(df, idler=args.idler, choice='erase', window_ns=args.range_ns)
    del_mark  = extract_deltas(df, idler=args.idler, choice='mark',  window_ns=args.range_ns)

    # Histograms
    edges = np.arange(-args.range_ns, args.range_ns + args.bin_ns, args.bin_ns)
    ctrs = 0.5 * (edges[:-1] + edges[1:])
    H_e, _ = np.histogram(del_erase, bins=edges)
    H_m, _ = np.histogram(del_mark,  bins=edges)
    H_d = H_e - H_m  # choice-locked difference

    # Noise estimate (pre-window) from negative side excluding ±2 bins around expected
    mask_pre = (ctrs < args.expected_ns - 2*args.bin_ns)
    sigma_hat = np.std(H_d[mask_pre], ddof=1)
    sigma_hat = max(sigma_hat, args.sigma_floor)

    # Identify target bin nearest expected_ns
    k = int(np.argmin(np.abs(ctrs - args.expected_ns)))
    y = H_d[k]

    # Prior and posterior
    sigma0 = args.prior_mult * sigma_hat
    mu_post, sd_post = normal_posterior(y, sigma_hat, sigma0)

    # Monte Carlo samples for HDI on |γ|
    N = 200_000
    samples = np.random.default_rng(0).normal(mu_post, sd_post, size=N)
    hdi_lo, hdi_hi = hdi_of_abs(samples, cred=0.95)
    rope = args.rope_mult * sigma_hat
    p_in_rope = np.mean(np.abs(samples) <= rope)

    # ---- Plots ----
    # 1) Histograms
    plt.figure(figsize=(9, 4.8))
    plt.step(ctrs, H_e, where='mid', label='Erase', linewidth=1.5)
    plt.step(ctrs, H_m, where='mid', label='Mark', linewidth=1.5)
    plt.step(ctrs, H_d, where='mid', label='Erase - Mark', linewidth=1.8)
    plt.axvline(args.expected_ns, linestyle='--', alpha=0.7, label='expected advanced time')
    plt.xlabel(r'$\Delta t \equiv t_{\rm idler}-t_{D0}$ (ns)')
    plt.ylabel('counts / bin')
    plt.legend()
    plt.tight_layout()
    plt.savefig(out / 'coincidence_hist.png', dpi=140)
    plt.close()

    # 2) Posterior plot (same style as your optical)
    xs = np.linspace(mu_post - 8*sd_post, mu_post + 8*sd_post, 1000)
    pdf = np.exp(-0.5*((xs-mu_post)/sd_post)**2)/(sd_post*np.sqrt(2*np.pi))

    plt.figure(figsize=(9, 4.8))
    plt.plot(xs, pdf, label=r'Posterior PDF for $\gamma_{\rm adv}$ (counts/bin)')
    plt.axvline(-rope, ls='--', alpha=0.6, label=r'$\pm$ ROPE')
    plt.axvline(+rope, ls='--', alpha=0.6)
    # Shade 95% HDI on |γ|
    # We visualize as [0, hdi_hi] since it's an interval on |γ|
    plt.fill_betweenx([0, pdf.max()*1.02], 0, hdi_hi, alpha=0.15, label='95% HDI on $|\\gamma|$')
    plt.xlabel(r'$\gamma_{\rm adv}$ (counts/bin)')
    plt.ylabel('Density (a.u.)')
    plt.title('Posterior for $\\gamma_{\\rm adv}$ with ROPE & 95% HDI (QE)')
    plt.legend()
    plt.tight_layout()
    plt.savefig(out / 'posterior_gamma_plot_qe.png', dpi=140)
    plt.close()

    # ---- Results table ----
    res = {
        'expected_ns': args.expected_ns,
        'bin_ns': args.bin_ns,
        'y_at_expected_counts_per_bin': float(y),
        'prewindow_rms_counts_per_bin': float(sigma_hat),
        'ROPE_counts_per_bin': float(rope),
        'P_in_ROPE': float(p_in_rope),
        'HDI_low_abs_counts_per_bin': float(hdi_lo),
        'HDI_high_abs_counts_per_bin': float(hdi_hi)
    }
    pd.DataFrame([res]).to_csv(out / 'results_table_qe.csv', index=False)

    print("\n--- QE Analysis Results ---")
    print(f"Expected advanced time    : {args.expected_ns:.3f} ns")
    print(f"y at expected bin         : {y:.2f} counts/bin")
    print(f"σ̂ (pre-window RMS)       : {sigma_hat:.2f} counts/bin")
    print(f"95% HDI(|γ|)              : [{hdi_lo:.2f}, {hdi_hi:.2f}] counts/bin")
    print(f"ROPE ±                    : {rope:.2f} (P in ROPE = {100*p_in_rope:.1f}%)")
    print(f"Exports → {out}/coincidence_hist.png, {out}/posterior_gamma_plot_qe.png, {out}/results_table_qe.csv")


if __name__ == '__main__':
    main()
