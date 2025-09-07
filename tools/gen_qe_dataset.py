#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate synthetic time-tagged data for a causally disconnected quantum eraser (QE)
in the style of Ma et al. Outputs a single CSV with columns:
  event_id, detector_id, time_tag_ns, experiment_run, choice_type

Key idea: for idler detector D3 (the long-path arm), we optionally inject an
"advanced component" by placing some coincidences at t_idler = t_D0 - Δt_long.
"""

import argparse
import numpy as np
import pandas as pd
from pathlib import Path

def generate_qe_time_tags(num_events=50_000,
                          gamma_adv_strength=0.03,
                          seed=42,
                          delay_short_ns=0.0,
                          delay_long_ns=1.67,
                          jitter_ps=50.0,
                          coincidence_prob=0.30,
                          nruns=100):
    rng = np.random.default_rng(seed)

    # Containers
    event_ids, dets, tt, runs, choices = [], [], [], [], []
    eid = 0

    base_time_ns = 100.0  # arbitrary offset

    for run in range(1, nruns + 1):
        choice = rng.choice(['erase', 'mark'])
        n_run = num_events // nruns

        # Poissonian D0 stream via exponential interarrival (mean 0.1 ns)
        d0_times = base_time_ns + np.cumsum(rng.exponential(0.1, n_run))

        for t0 in d0_times:
            # With some probability, create a D0–idler coincidence pair
            if rng.random() < coincidence_prob:
                # Choose idler detector given the delayed-choice branch
                if choice == 'erase':
                    # roughly balanced across D1..D4
                    det = rng.choice(['D1', 'D2', 'D3', 'D4'], p=[0.25, 0.25, 0.25, 0.25])
                else:
                    # emphasize which-path tagging in the "mark" branch
                    det = rng.choice(['D1', 'D2', 'D3', 'D4'], p=[0.40, 0.40, 0.10, 0.10])

                # Idler time relative to D0
                if det in ('D1', 'D2'):
                    tid = t0 + delay_short_ns
                else:  # D3 or D4 (long arm)
                    tid = t0 + delay_long_ns

                # Inject "advanced" component for D3 only (long-path arm),
                # with probability gamma_adv_strength.
                if det == 'D3' and rng.random() < gamma_adv_strength:
                    tid = t0 - delay_long_ns  # ADVANCED PEAK at negative lag

                # Apply Gaussian timing jitter (ps → ns)
                sigma_ns = jitter_ps * 1e-3
                t0_j = t0 + rng.normal(0.0, sigma_ns)
                tid_j = tid + rng.normal(0.0, sigma_ns)

                # Record D0
                event_ids.append(eid); eid += 1
                dets.append('D0')
                tt.append(t0_j)
                runs.append(run)
                choices.append(choice)

                # Record idler
                event_ids.append(eid); eid += 1
                dets.append(det)
                tt.append(tid_j)
                runs.append(run)
                choices.append(choice)

    df = pd.DataFrame({
        'event_id': event_ids,
        'detector_id': dets,
        'time_tag_ns': tt,
        'experiment_run': runs,
        'choice_type': choices
    }).sort_values('time_tag_ns').reset_index(drop=True)

    # Reassign ordered IDs after sort
    df['event_id'] = np.arange(len(df), dtype=int)
    return df


def main():
    ap = argparse.ArgumentParser(description="Generate QE time-tagged dataset (synthetic).")
    ap.add_argument('--out', type=str, default='data/causally_disconnected_eraser_simulation.csv')
    ap.add_argument('--events', type=int, default=50_000)
    ap.add_argument('--gamma', type=float, default=0.03, help='advanced fraction for D3')
    ap.add_argument('--seed', type=int, default=42)
    ap.add_argument('--delay_short_ns', type=float, default=0.0)
    ap.add_argument('--delay_long_ns', type=float, default=1.67)
    ap.add_argument('--jitter_ps', type=float, default=50.0)
    ap.add_argument('--coin_p', type=float, default=0.30)
    ap.add_argument('--runs', type=int, default=100)
    args = ap.parse_args()

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    df = generate_qe_time_tags(num_events=args.events,
                               gamma_adv_strength=args.gamma,
                               seed=args.seed,
                               delay_short_ns=args.delay_short_ns,
                               delay_long_ns=args.delay_long_ns,
                               jitter_ps=args.jitter_ps,
                               coincidence_prob=args.coin_p,
                               nruns=args.runs)
    df.to_csv(out, index=False)
    print(f"[qe-gen] wrote {out}  shape={df.shape}")


if __name__ == '__main__':
    main()
