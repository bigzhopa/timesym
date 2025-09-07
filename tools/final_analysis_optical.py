#!/usr/bin/env python3
"""
final_analysis_optical.py
Free-space/heralded platform — harmonized Bayesian bound on a pre-window precursor.

Outputs:
  - results_table_optical.csv
  - posterior_gamma_samples_optical.csv
  - posterior_gamma_plot_optical.png
  - params_optical.yaml/json
  - optical_histograms_simulated.csv  (when --simulate)
"""

import os, sys, json
from dataclasses import dataclass, asdict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ----------------- helpers -----------------
def to_native(obj):
    """Convert numpy types -> pure Python for YAML/JSON."""
    import numpy as _np
    if isinstance(obj, dict):
        return {k: to_native(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [to_native(v) for v in obj]
    if isinstance(obj, _np.generic):
        return obj.item()
    if isinstance(obj, _np.ndarray):
        return obj.tolist()
    return obj

def hdi_from_samples(samples, cred=0.95):
    s = np.sort(samples)
    n = len(s)
    m = max(1, int(np.floor(cred * n)))
    idx = np.argmin(s[m:] - s[:n - m])
    return s[idx], s[idx + m]

def smooth(x, w=9):
    if w <= 1: 
        return x
    w = int(w) | 1
    k = np.ones(w, dtype=float) / w
    return np.convolve(x, k, mode="same")

# YAML (optional)
try:
    import yaml
    _HAVE_YAML = True
except Exception:
    _HAVE_YAML = False

# ----------------- parameters -----------------
@dataclass
class Params:
    # timing & grid
    dt: float = 5.0e-12
    Trec: float = 50.0e-9
    t0: float = 10.0e-9
    L: float = 1.0
    c: float = 299792458.0
    sigma: float = 200.0e-12

    # stats / physics
    N_trials: int = 200_000
    p_pair: float = 0.01
    eta_det: float = 0.9
    jitter_sigma: float = 80.0e-12
    dark_rate: float = 100.0  # Hz, uniform accidental rate

    # interferometer depths (explicit)
    null_depth_dB: float = 40.0       # A: deep cancel (amplitude, dB)
    null_depth_B_dB: float = 20.0     # B: spoiled cancel (amplitude, dB)
    spoil_phase_deg: float = 5.0      # kept in case you want to use it elsewhere

    # analysis window
    pre_window_ns: float = 20.0
    pre_end_sigma: float = 5.0
    rope_mult: float = 2.0
    prior_scale_mult: float = 10.0
    Nsamples: int = 200_000
    sigma_floor_counts: float = 1.0
    min_pre_bins: int = 20

    # normalization
    Vref_mode: str = "mainpeak"       # "mainpeak" | "mainlobe" | "custom"
    Vref_custom: float = 1.0

    # io & seed
    out_dir: str = "."
    seed: int = 2025

# ----------------- simulation -----------------
def simulate_histograms(par: Params):
    rng = np.random.default_rng(par.seed)

    # time axis
    t = np.arange(0.0, par.Trec, par.dt, dtype=float)

    # pulse with jitter folded in (gaussian convolution → variance add)
    sigma_eff = np.sqrt(par.sigma**2 + par.jitter_sigma**2)
    g = np.exp(-0.5 * ((t - par.t0) / sigma_eff) ** 2)
    area = np.trapezoid(g, t)
    if area <= 0:
        area = 1.0
    g_norm = g / area  # integrates to 1 over time

    # amplitudes (explicit depths; amplitude ratio in linear units)
    D_A = 10 ** (-par.null_depth_dB  / 20.0)
    D_B = 10 ** (-par.null_depth_B_dB / 20.0)

    # expected *signal* counts per bin (Poisson mean)
    K = par.N_trials * par.p_pair * par.eta_det
    lam_true_A = K * (D_A**2) * g_norm * par.dt
    lam_true_B = K * (D_B**2) * g_norm * par.dt

    # uniform darks per bin across record
    lam_dark = par.dark_rate * par.dt * par.N_trials
    lam_A = lam_true_A + lam_dark
    lam_B = lam_true_B + lam_dark

    # sample Poisson
    A = rng.poisson(lam_A)
    B = rng.poisson(lam_B)

    return t, A.astype(float), B.astype(float)

# ----------------- analysis -----------------
def analyze(par: Params, t, A, B):
    # differenced trace
    D = A - B

    # timing
    tau = par.L / par.c
    idx_tau = int(np.argmin(np.abs(t - (par.t0 - tau))))
    y_obs = float(D[idx_tau])

    # pre-window mask
    pre_start = par.t0 - par.pre_window_ns * 1e-9
    pre_end = par.t0 - par.pre_end_sigma * par.sigma
    pre_mask = (t >= pre_start) & (t < pre_end)
    if np.count_nonzero(pre_mask) < par.min_pre_bins:
        # if not enough bins, expand to earliest
        pre_mask = (t < pre_end)
    sigma_hat = float(np.std(D[pre_mask], ddof=1))
    if not np.isfinite(sigma_hat) or sigma_hat <= 0:
        sigma_hat = par.sigma_floor_counts

    # Bayesian Normal-Normal (closed form)
    tau0 = par.prior_scale_mult * sigma_hat
    sigma2 = sigma_hat**2
    tau02 = tau0**2
    post_var = 1.0 / (1.0/tau02 + 1.0/sigma2)
    post_mean = post_var * (0.0/tau02 + y_obs/sigma2)

    # draw posterior samples
    samples = np.random.default_rng(par.seed + 1).normal(post_mean, np.sqrt(post_var), size=par.Nsamples)
    abs_samples = np.abs(samples)
    lo, hi = hdi_from_samples(abs_samples, 0.95)

    # ROPE
    rope = par.rope_mult * sigma_hat
    p_in_rope = float(np.mean(abs_samples < rope))

    # normalization
    mean_trace = 0.5 * (A + B)
    if par.Vref_mode == "mainpeak":
        Vref = float(np.max(smooth(mean_trace, w=9)))
        if Vref <= 0: Vref = 1.0
    elif par.Vref_mode == "mainlobe":
        mask_lobe = (t >= par.t0 - 3*par.sigma) & (t <= par.t0 + 3*par.sigma)
        Vref = float(max(np.sum(A[mask_lobe]), np.sum(B[mask_lobe])))
        if Vref <= 0: Vref = 1.0
    else:
        Vref = max(float(par.Vref_custom), 1e-30)

    bound_abs = hi
    bound_rel = bound_abs / Vref
    bound_dBc = 20.0 * np.log10(bound_rel)  # reported as dBc relative to chosen Vref

    res = {
        "tau_ns": 1e9 * (par.L / par.c),
        "t0_ns": 1e9 * par.t0,
        "t0_minus_tau_ns": 1e9 * (par.t0 - par.L/par.c),
        "y_obs_counts": y_obs,
        "prewindow_rms_counts": sigma_hat,
        "HDI_low_counts": lo,
        "HDI_high_counts": hi,
        "P_in_ROPE": p_in_rope,
        "bound_abs_counts_per_bin": bound_abs,
        "bound_rel": bound_rel,
        "bound_dBc": bound_dBc,
        "Vref_used": Vref,
    }
    return res, samples, D, pre_mask

# ----------------- plotting -----------------
def plot_posterior(par: Params, samples, res, out_png):
    fig, ax = plt.subplots(figsize=(8, 4.5))
    # simple KDE via histogram
    xs = np.linspace(np.min(samples)-3*np.std(samples), np.max(samples)+3*np.std(samples), 1024)
    hist, edges = np.histogram(samples, bins=200, density=True)
    centers = 0.5*(edges[:-1]+edges[1:])
    ax.plot(centers*1e0, hist, label="Posterior PDF for γ_adv (counts/bin)")
    rope = par.rope_mult * res["prewindow_rms_counts"]
    ax.axvline(-rope, ls="--", color="C0", alpha=0.7, label="± ROPE")
    ax.axvline(+rope, ls="--", color="C0", alpha=0.7)
    lo, hi = res["HDI_low_counts"], res["HDI_high_counts"]
    ax.fill_betweenx([0, hist.max()*1.05], lo, hi, color="C0", alpha=0.15, label="95% HDI on |γ_adv|")
    ax.set_xlabel("γ_adv (counts/bin)")
    ax.set_ylabel("Density (a.u.)")
    ax.set_title("Posterior for γ_adv with ROPE & 95% HDI (optical)")
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(out_png, dpi=160)
    plt.close(fig)

# ----------------- io -----------------
def save_outputs(par: Params, t, A, B, res, samples, D, pre_mask):
    os.makedirs(par.out_dir, exist_ok=True)

    # params
    p_dict = to_native(asdict(par))
    if _HAVE_YAML:
        try:
            with open(os.path.join(par.out_dir, "params_optical.yaml"), "w", encoding="utf-8") as f:
                yaml.safe_dump(p_dict, f, sort_keys=False)
        except Exception:
            with open(os.path.join(par.out_dir, "params_optical.json"), "w", encoding="utf-8") as f:
                json.dump(p_dict, f, indent=2)
    else:
        with open(os.path.join(par.out_dir, "params_optical.json"), "w", encoding="utf-8") as f:
            json.dump(p_dict, f, indent=2)

    # table
    df = pd.DataFrame([res])
    df.to_csv(os.path.join(par.out_dir, "results_table_optical.csv"), index=False)

    # posterior samples
    pd.DataFrame({"gamma_adv_samples_counts": samples}).to_csv(
        os.path.join(par.out_dir, "posterior_gamma_samples_optical.csv"), index=False
    )

    # plot
    plot_posterior(par, samples, res, os.path.join(par.out_dir, "posterior_gamma_plot_optical.png"))

    # traces (if simulated)
    traces = pd.DataFrame({
        "time_ns": 1e9*t,
        "hist_A_counts": A,
        "hist_B_counts": B,
        "D_counts": A - B,
        "pre_mask": pre_mask.astype(int),
    })
    traces.to_csv(os.path.join(par.out_dir, "optical_histograms_simulated.csv"), index=False)

# ----------------- main -----------------
def main():
    import argparse
    ap = argparse.ArgumentParser(description="Free-space/heralded Bayesian bound (harmonized, robust).")
    ap.add_argument("--data", type=str, default=None, help="CSV: time_ns,hist_A_counts,hist_B_counts")
    ap.add_argument("--simulate", action="store_true", help="Force simulation even if --data provided")
    ap.add_argument("--out", type=str, default=".", help="Output directory")

    ap.add_argument("--Vref", type=str, default="mainpeak",
                    choices=["mainpeak","mainlobe","custom"],
                    help="Normalization: mainpeak (max bin), mainlobe (∑±3σ), or custom")
    ap.add_argument("--Vref_custom", type=float, default=1.0, help="If --Vref custom, numeric reference")

    ap.add_argument("--sigma_floor_counts", type=float, default=None,
                    help="Noise floor if pre-window variance is zero (default 1.0)")
    ap.add_argument("--seed", type=int, default=2025, help="RNG seed (simulation)")

    # overrides for simulation fidelity
    ap.add_argument("--N_trials", type=int, default=None)
    ap.add_argument("--p_pair", type=float, default=None)
    ap.add_argument("--null_depth_dB", type=float, default=None)
    ap.add_argument("--null_depth_B_dB", type=float, default=None)
    ap.add_argument("--spoil_phase_deg", type=float, default=None)
    ap.add_argument("--dt_override", type=float, default=None)

    args = ap.parse_args()

    par = Params(out_dir=args.out, Vref_mode=args.Vref, Vref_custom=args.Vref_custom, seed=args.seed)

    # optional analysis knobs
    if args.sigma_floor_counts is not None: par.sigma_floor_counts = float(args.sigma_floor_counts)

    # apply overrides
    if args.N_trials is not None:        par.N_trials = int(args.N_trials)
    if args.p_pair   is not None:        par.p_pair = float(args.p_pair)
    if args.null_depth_dB is not None:   par.null_depth_dB = float(args.null_depth_dB)
    if args.null_depth_B_dB is not None: par.null_depth_B_dB = float(args.null_depth_B_dB)
    if args.spoil_phase_deg is not None: par.spoil_phase_deg = float(args.spoil_phase_deg)
    if args.dt_override is not None:     par.dt = float(args.dt_override)

    # ingest or simulate
    if (args.data is not None) and (not args.simulate) and os.path.exists(args.data):
        df = pd.read_csv(args.data)
        t = df["time_ns"].values * 1e-9
        A = df["hist_A_counts"].values.astype(float)
        B = df["hist_B_counts"].values.astype(float)
    else:
        print("[+] Simulating optical histograms…")
        t, A, B = simulate_histograms(par)

    res, samples, D, pre_mask = analyze(par, t, A, B)

    # pretty print
    print("\n--- Optical Platform: Final Results ---")
    print(f"Delay τ                : {res['tau_ns']:.3f} ns")
    print(f"y = D(t0−τ)            : {res['y_obs_counts']:.2f} counts/bin")
    print(f"σ̂ (pre-window RMS)    : {res['prewindow_rms_counts']:.2f} counts/bin")
    print(f"95% HDI(|γ_adv|)       : [{res['HDI_low_counts']:.2f}, {res['HDI_high_counts']:.2f}] counts/bin")
    rope = 2.0 * res['prewindow_rms_counts']  # display only; actual multiplier stored in params
    print(f"ROPE ±                 : {rope:.2f}  (P in ROPE = {res['P_in_ROPE']*100:.1f}%)")
    print(f"Bound (absolute)       : |γ_adv| < {res['bound_abs_counts_per_bin']:.2f} counts/bin  (95% credible upper)")
    print(f"Normalization Vref     : {res['Vref_used']:.2f} (per {par.Vref_mode})")
    print(f"Bound (relative)       : |γ_adv|/Vref < {res['bound_rel']:.3e}")
    print(f"Bound (dBc)            : {res['bound_dBc']:.1f} dB re chosen Vref\n")
    print("Exports:\n  results_table_optical.csv, posterior_gamma_samples_optical.csv,\n"
          "  posterior_gamma_plot_optical.png, params_optical.yaml/json\n"
          "  optical_histograms_simulated.csv")

    save_outputs(par, t, A, B, res, samples, D, pre_mask)

if __name__ == "__main__":
    main()
