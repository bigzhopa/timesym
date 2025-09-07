#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
final_analysis.py — Publication-grade bound for advanced-wave precursor

Features
- Accepts real oscilloscope data (CSV) or simulates traces if not provided.
- Interferometer (3 dB hybrid abstraction): deep-cancel (A) vs spoiled-null (B).
- Coax model: 1 m PTFE-like line, small VSWR at source/load, frequency-dependent
  round-trip loss (skin + dielectric).
- DAQ realism: 33 GHz low-pass, 1000× averaging.
- Statistical analysis: Bayesian Normal–Normal posterior for γ_adv at t0−τ,
  95% HDI for |γ_adv|, ROPE probability; normalized and dBc bounds.
- Exports: results_table.csv, posterior samples, posterior figure, params.yaml,
  and (for sim) traces_with_loss_and_bayes.csv.

Usage
  python final_analysis.py --data real_data.csv --out ./out --Vpp 1.0
  python final_analysis.py --simulate --out ./out

CSV schema (real data)
  time_ns, trace_A_V, trace_B_V
"""

import argparse, os, sys, json
from dataclasses import dataclass, asdict
from typing import Tuple, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import json
try:
    import yaml  # for params.yaml
    _HAVE_YAML = True
except Exception:
    _HAVE_YAML = False


# ---------------- Parameters ----------------
@dataclass
class Params:
    # Line & timing
    Z0: float = 50.0
    vp: float = 2.076e8            # m/s (PTFE ~0.692c)
    L: float = 1.0                 # m
    dt: float = 5e-12              # s (5 ps)
    Trec: float = 50e-9            # s total record
    # Source pulse (Gaussian monocycle)
    t0: float = 10e-9              # s
    sigma: float = 200e-12         # s
    Vpp: float = 1.0               # drive amplitude used for normalization
    # Interferometer
    fc: float = 5e9                # Hz (for 5° mapping)
    delta_phase_deg: float = 5.0   # state-B phase error at fc
    # DAQ realism
    Navg: int = 1000
    noise_single_shot: float = 5e-4  # 500 µV RMS per shot
    scope_bw_Hz: float = 33e9
    # Mismatches (VSWR → |Γ|)
    VSWR_source: float = 1.05
    VSWR_load: float = 1.05
    # Frequency-dependent losses (≈0.1 dB/m @1 GHz skin; ≈0.05 dB/m @10 GHz diel)
    k_skin: float = 0.01151 / (1e9 ** 0.5)  # Np/m/√Hz
    k_diel: float = 0.005755 / 1e10         # Np/m/Hz
    # Analysis windows
    pre_window_ns: float = 20.0   # pre-trigger window length
    pre_end_sigma: float = 5.0    # pre-window ends at t0 - 5σ
    rope_mult: float = 2.0        # ROPE = rope_mult * sigma_hat
    # Posterior sampling
    Nsamples: int = 200_000
    prior_scale_mult: float = 10.0  # τ0 = prior_scale_mult * noise_rms
    # I/O
    out_dir: str = "."
    seed: int = 1234


# ---------------- Helpers ----------------
def to_native(obj):
    """Recursively convert NumPy scalars/arrays to plain Python types for YAML/JSON."""
    import numpy as np
    if isinstance(obj, dict):
        return {k: to_native(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [to_native(v) for v in obj]
    if isinstance(obj, np.generic):
        return obj.item()
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    return obj


def gaussian_monocycle(t: np.ndarray, t0: float, sigma: float, Vpp: float) -> np.ndarray:
    """Derivative of Gaussian, normalized to ~Vpp peak-to-peak."""
    g = np.exp(-0.5 * ((t - t0) / sigma) ** 2)
    dg = -(t - t0) / (sigma ** 2) * g
    pkpk = dg.max() - dg.min()
    return np.zeros_like(t) if pkpk == 0 else Vpp * dg / pkpk


def splitter_3dB(v_in: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    return v_in / np.sqrt(2.0), v_in / np.sqrt(2.0)


def delay_signal(x: np.ndarray, dt: float, delay_seconds: float) -> np.ndarray:
    if delay_seconds == 0:
        return x.copy()
    shift = delay_seconds / dt
    n = np.arange(len(x)) - shift
    n0 = np.floor(n).astype(int)
    frac = n - n0
    y = np.zeros_like(x)
    m = (n0 >= 0) & (n0 + 1 < len(x))
    y[m] = (1 - frac[m]) * x[n0[m]] + frac[m] * x[n0[m] + 1]
    return y


def interferometer_output(Vs: np.ndarray, dt: float, fc: float,
                          delta_phase_deg: float, state: str) -> Tuple[np.ndarray, np.ndarray]:
    arm1, arm2 = splitter_3dB(Vs)
    arm2 = -arm2  # π phase
    if state == 'B':
        dphi = np.deg2rad(delta_phase_deg)
        dt_extra = dphi / (2 * np.pi * fc)
    else:
        dt_extra = 0.0
    arm2_d = delay_signal(arm2, dt, dt_extra)
    v_line = (arm1 + arm2_d) / np.sqrt(2.0)  # to coax (dark port when balanced)
    v_dump = (arm1 - arm2_d) / np.sqrt(2.0)  # to dump
    return v_line, v_dump


def run_line(v_in_line: np.ndarray, Nx: int, Gamma_s: float, Gamma_L: float, Nt: int) -> np.ndarray:
    """Digital waveguide: right/left traveling waves, one-cell per step."""
    vplus = np.zeros(Nx)
    vminus = np.zeros(Nx)
    Vsrc = np.zeros(Nt)
    for k in range(Nt):
        b0 = vminus[0]
        a0 = v_in_line[k] + Gamma_s * b0
        vplus[0] = a0
        aN = vplus[-1]
        vminus[-1] = Gamma_L * aN
        if Nx > 1:
            vplus[1:] = vplus[:-1]
            vminus[:-1] = vminus[1:]
        Vsrc[k] = vplus[0] + vminus[0]
    return Vsrc


def apply_scope_filter(x: np.ndarray, dt: float, fcut: float) -> np.ndarray:
    X = np.fft.rfft(x)
    f = np.fft.rfftfreq(len(x), d=dt)
    H = (f <= fcut).astype(float)
    return np.fft.irfft(X * H, n=len(x))


def hdi_from_samples(a: np.ndarray, cred_mass: float = 0.95) -> Tuple[float, float]:
    a_sorted = np.sort(a)
    n = len(a_sorted)
    k = int(np.floor(cred_mass * n))
    widths = a_sorted[k:] - a_sorted[:n - k]
    i = np.argmin(widths)
    return float(a_sorted[i]), float(a_sorted[i + k])


def rms(a: np.ndarray) -> float:
    return float(np.sqrt(np.mean(a * a)))


# ---------------- Simulation (if no real data) ----------------
def simulate_traces(par: Params, rng: np.random.Generator) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    dx = par.vp * par.dt
    Nx = int(np.ceil(par.L / dx))
    Nt = int(np.ceil(par.Trec / par.dt))
    t = np.arange(Nt) * par.dt

    Vs = gaussian_monocycle(t, par.t0, par.sigma, par.Vpp)
    vline_A, _ = interferometer_output(Vs, par.dt, par.fc, par.delta_phase_deg, 'A')
    vline_B, _ = interferometer_output(Vs, par.dt, par.fc, par.delta_phase_deg, 'B')

    Gamma_s = (par.VSWR_source - 1) / (par.VSWR_source + 1)
    Gamma_L = (par.VSWR_load - 1) / (par.VSWR_load + 1)

    V_A_nom = run_line(vline_A, Nx, Gamma_s, 0.0, Nt)
    V_A_ref = run_line(vline_A, Nx, Gamma_s, Gamma_L, Nt) - V_A_nom
    V_B_nom = run_line(vline_B, Nx, Gamma_s, 0.0, Nt)
    V_B_ref = run_line(vline_B, Nx, Gamma_s, Gamma_L, Nt) - V_B_nom

    # Frequency-dependent round-trip attenuation for reflection
    freqs = np.fft.rfftfreq(Nt, d=par.dt)
    alpha_skin = par.k_skin * np.sqrt(np.maximum(freqs, 0.0))
    alpha_diel = par.k_diel * freqs
    A_refl = np.exp(-2 * par.L * (alpha_skin + alpha_diel))

    def apply_roundtrip_loss(x):
        X = np.fft.rfft(x)
        return np.fft.irfft(X * A_refl, n=len(x))

    V_A_raw = V_A_nom + apply_roundtrip_loss(V_A_ref)
    V_B_raw = V_B_nom + apply_roundtrip_loss(V_B_ref)

    noise_rms = par.noise_single_shot / np.sqrt(par.Navg)

    def measure(V_raw):
        V = V_raw + rng.normal(0.0, noise_rms, V_raw.shape)
        return apply_scope_filter(V, par.dt, par.scope_bw_Hz)

    V_A = measure(V_A_raw)
    V_B = measure(V_B_raw)
    return t, V_A, V_B


# ---------------- Analysis ----------------
def analyze(par: Params, t: np.ndarray, V_A: np.ndarray, V_B: np.ndarray) -> dict:
    tau = par.L / par.vp
    idx_adv = int(np.clip(round(((par.t0 - tau) - t[0]) / par.dt), 0, len(t) - 1))

    D = V_A - V_B

    # Pre-window: [t0 - pre_window_ns, t0 - pre_end_sigma*sigma)
    pre_start = par.t0 - par.pre_window_ns * 1e-9
    pre_end = par.t0 - par.pre_end_sigma * par.sigma
    pre_mask = (t >= pre_start) & (t < pre_end)
    D_pre = D[pre_mask]
    sigma_hat = np.std(D_pre, ddof=1)

    # Bayesian Normal–Normal: y ~ N(γ, σ^2), prior γ ~ N(0, τ0^2) with τ0 = prior_scale_mult * noise_rms
    noise_rms = par.noise_single_shot / np.sqrt(par.Navg)
    tau0 = par.prior_scale_mult * noise_rms
    y_obs = float(D[idx_adv])
    sigma2 = float(sigma_hat ** 2)
    tau02 = float(tau0 ** 2)
    post_var = 1.0 / (1.0 / tau02 + 1.0 / sigma2)
    post_mean = post_var * (y_obs / sigma2)
    post_std = float(np.sqrt(post_var))

    # ROPE and HDI
    rope = par.rope_mult * sigma_hat
    rng = np.random.default_rng(par.seed + 1)
    samples = rng.normal(post_mean, post_std, size=par.Nsamples)
    abs_samples = np.abs(samples)
    hdi_low, hdi_high = hdi_from_samples(abs_samples, 0.95)
    p_in_rope = float(np.mean(np.abs(samples) < rope))

    # Normalized bounds
    bound_abs = hdi_high
    bound_rel = bound_abs / max(par.Vpp, 1e-30)
    bound_dBc = float(20.0 * np.log10(bound_rel + 1e-30))

    return dict(
        tau=tau, idx_adv=idx_adv, y_obs=y_obs, sigma_hat=sigma_hat,
        post_mean=post_mean, post_std=post_std,
        rope=rope, p_in_rope=p_in_rope, hdi_low=hdi_low, hdi_high=hdi_high,
        bound_abs=bound_abs, bound_rel=bound_rel, bound_dBc=bound_dBc,
        D=D, D_pre=D_pre, samples=samples
    )


def save_outputs(par: Params, t: np.ndarray, V_A: np.ndarray, V_B: np.ndarray, res: dict,
                 save_traces: bool) -> None:
    os.makedirs(par.out_dir, exist_ok=True)

    # Posterior plot
    xs = np.linspace(res["post_mean"] - 6 * res["post_std"],
                     res["post_mean"] + 6 * res["post_std"], 1000)
    pdf = (1.0 / (res["post_std"] * np.sqrt(2 * np.pi))) * np.exp(
        -0.5 * ((xs - res["post_mean"]) / res["post_std"]) ** 2
    )
    plt.figure(figsize=(7, 4.2))
    plt.plot(xs * 1e6, pdf, label="Posterior PDF for γ_adv")
    plt.axvline(res["rope"] * 1e6, linestyle="--", label="± ROPE")
    plt.axvline(-res["rope"] * 1e6, linestyle="--")
    plt.axvspan(-res["hdi_high"] * 1e6, res["hdi_high"] * 1e6, alpha=0.2, label="95% HDI on |γ_adv|")
    plt.xlabel("γ_adv (µV)"); plt.ylabel("Density (a.u.)")
    plt.title("Posterior for γ_adv with ROPE & 95% HDI")
    plt.legend(); plt.tight_layout()
    fig_path = os.path.join(par.out_dir, "posterior_gamma_plot.png")
    plt.savefig(fig_path, dpi=200)
    plt.close()

    # Results table
    results = {
        "bound_abs_uV": res["bound_abs"] * 1e6,
        "bound_rel": res["bound_rel"],
        "bound_dBc": res["bound_dBc"],
        "P_in_ROPE": res["p_in_rope"],
        "prewindow_rms_uV": res["sigma_hat"] * 1e6,
        "HDI_low_uV": res["hdi_low"] * 1e6,
        "HDI_high_uV": res["hdi_high"] * 1e6,
        "y_obs_uV": res["y_obs"] * 1e6,
        "tau_ns": res["tau"] * 1e9,
        "t0_ns": par.t0 * 1e9,
        "t0_minus_tau_ns": (par.t0 - res["tau"]) * 1e9,
        "Vpp": par.Vpp,
    }
    pd.DataFrame([results]).to_csv(os.path.join(par.out_dir, "results_table.csv"), index=False)

    # Posterior samples (downsampled)
    pd.DataFrame({"gamma_adv_V": res["samples"][:20_000]}).to_csv(
        os.path.join(par.out_dir, "posterior_gamma_samples.csv"), index=False
    )

    # Params audit trail
    if _HAVE_YAML:
        with open(os.path.join(par.out_dir, "params.yaml"), "w") as f:
            yaml.safe_dump(asdict(par), f, sort_keys=False)
    else:
        with open(os.path.join(par.out_dir, "params.json"), "w") as f:
            json.dump(asdict(par), f, indent=2)

    # Traces (if we simulated)
    if save_traces:
        pd.DataFrame({
            "time_ns": t * 1e9,
            "trace_A_V": V_A,
            "trace_B_V": V_B,
            "D_AminusB_V": res["D"]
        }).to_csv(os.path.join(par.out_dir, "traces_with_loss_and_bayes.csv"), index=False)


# ---------------- CLI ----------------
def main():
    ap = argparse.ArgumentParser(description="Publication-grade bound for advanced-wave precursor.")
    ap.add_argument("--data", type=str, default=None, help="CSV with columns: time_ns,trace_A_V,trace_B_V")
    ap.add_argument("--simulate", action="store_true", help="Force simulation even if --data is provided")
    ap.add_argument("--out", type=str, default=".", help="Output directory")
    ap.add_argument("--Vpp", type=float, default=1.0, help="Drive amplitude used for normalization (Vpp)")
    ap.add_argument("--seed", type=int, default=1234, help="RNG seed")
    args = ap.parse_args()

    par = Params(out_dir=args.out, Vpp=args.Vpp, seed=args.seed)

    rng = np.random.default_rng(par.seed)
    save_traces = False

    if args.data and not args.simulate and os.path.exists(args.data):
        df = pd.read_csv(args.data)
        for col in ("time_ns", "trace_A_V", "trace_B_V"):
            if col not in df.columns:
                sys.exit(f"[!] Missing column '{col}' in {args.data}")
        t = df["time_ns"].values * 1e-9
        V_A = df["trace_A_V"].values
        V_B = df["trace_B_V"].values
        print(f"[+] Loaded real data: {args.data}  (N={len(t)})")
    else:
        print("[+] Generating simulated data…")
        t, V_A, V_B = simulate_traces(par, rng)
        save_traces = True

    # Sanity check: ensure dt matches loaded/simulated t
    par.dt = float(np.median(np.diff(t)))
    par.Trec = float(t[-1] - t[0] + par.dt)

    # Analyze and save
    res = analyze(par, t, V_A, V_B)
    save_outputs(par, t, V_A, V_B, res, save_traces)

    # Console summary (publication-ready numbers)
    print("\n--- Final Results ---")
    print(f"Delay τ                : {res['tau']*1e9:.3f} ns")
    print(f"y = D(t0−τ)            : {res['y_obs']*1e6:.3f} µV")
    print(f"σ̂ (pre-window RMS)    : {res['sigma_hat']*1e6:.3f} µV")
    print(f"95% HDI(|γ_adv|)       : [{res['hdi_low']*1e6:.3f}, {res['hdi_high']*1e6:.3f}] µV")
    print(f"ROPE ±                 : {res['rope']*1e6:.3f} µV  (P in ROPE = {res['p_in_rope']*100:.1f}%)")
    print(f"Bound (absolute)       : |γ_adv| < {res['bound_abs']*1e6:.1f} µV  (95% credible upper)")
    print(f"Bound (relative)       : |γ_adv|/Vpp < {res['bound_rel']:.3e}")
    print(f"Bound (dBc)            : {res['bound_dBc']:.1f} dB re {par.Vpp:.3f} Vpp")
    print("\nExports:")
    print(f"  results_table.csv, posterior_gamma_samples.csv, posterior_gamma_plot.png, params.yaml/json")
    if save_traces:
        print(f"  traces_with_loss_and_bayes.csv  (simulated)")

if __name__ == "__main__":
    main()
