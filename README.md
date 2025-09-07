# timesym

Temporal interferometer simulations and analysis (RF cable, optical) plus a QE validation dataset.
This repo keeps your original CLI tools intact under `tools/` and exposes simple Python shims under
the package `timesym` so you can call them from code.

## Install
```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
source .venv/bin/activate
pip install -e .
```

## Quick start
```bash
# RF and Optical platform
python tools\press_start.py --out runs\demo

# DCQE re-analysis
python tools\press_start_qe.py --out runs\qe_demo
```

## Quick start (CLIs preserved)
```bash
# RF cable (simulate)
python tools/final_analysis.py --simulate --out runs/cable_demo --Vpp 1.0 --seed 42

# Optical (simulate)
python tools/final_analysis_optical.py --simulate --out runs/opt_demo --seed 2025

# QE validation dataset + analysis
python tools/press_start_qe.py --out runs/qe_demo --gamma 0.03
```

## Quick start (Python API)
```python
from timesym import run_cable, run_optical
# assume you're calling from the repo root (or pass repo_root="...")
run_cable(out_dir="runs/cable_demo", simulate=True, Vpp=1.0, seed=42)
run_optical(out_dir="runs/opt_demo", simulate=True, seed=2025, Vref="mainlobe")
```

## Folders
- `tools/` — your original scripts (unchanged)
- `src/timesym/` — shims to call the tools from Python
- `runs/` — outputs land here
- `data/`, `cal/` — inputs/calibration
- `paper_assets/` — optional figure export target

## Notes
- The shims shell out to the scripts in `tools/`. If you call them from a different working
  directory, pass `repo_root="path/to/repo"`.
- Dependencies are declared in `pyproject.toml`.
