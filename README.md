# timesym — registered temporal‑interferometer analysis (RF/optical) + DCQE validation

[![CI](https://github.com/bigzhopa/timesym/actions/workflows/ci.yml/badge.svg)](https://github.com/bigzhopa/timesym/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17074673.svg)](https://doi.org/10.5281/zenodo.17074673)

Temporal‑interferometer simulations and analysis (RF cable, optical) plus a delayed‑choice quantum eraser (DCQE) validation dataset.  
Your original CLIs live under `tools/`; thin Python shims in `timesym` make them callable from code.

---

## Install

```bash
# Create & activate a virtual environment
python -m venv .venv

# Windows PowerShell
#   EITHER: .\.venv\Scripts\Activate.ps1  (may require: Set-ExecutionPolicy -Scope Process Bypass)
#   OR     : .\.venv\Scripts\activate.bat
# macOS/Linux
#   source .venv/bin/activate

# Editable install
pip install -U pip
pip install -e .
```

**Optional:** include plotting/export extras
```bash
pip install -e ".[plots]"
```

---

## Quick start (CLIs)

```bash
# RF cable (simulate)
python tools/final_analysis.py --simulate --out runs/cable_demo --Vpp 1.0 --seed 42

# Optical (simulate)
python tools/final_analysis_optical.py --simulate --out runs/opt_demo --seed 2025

# DCQE validation dataset + analysis
python tools/press_start_qe.py --out runs/qe_demo --gamma 0.03
```

> Note: a legacy `tools/press_start.py` may exist in some checkouts, but it’s optional.  
> Prefer the individual scripts above.

### Outputs
- RF: `runs/cable_demo/` → results table, posterior samples, posterior plot, parameters
- Optical: `runs/opt_demo/` → same as above, in count units
- DCQE: `runs/qe_demo/` → coincidence histogram, posterior, summary CSV

---

## Quick start (Python API)

```python
from timesym import run_cable, run_optical

# assuming you call from repo root; otherwise set repo_root="..."
run_cable(out_dir="runs/cable_demo", simulate=True, Vpp=1.0, seed=42)
run_optical(out_dir="runs/opt_demo", simulate=True, seed=2025, Vref="mainlobe")
```

---

## Repository layout

```
timesym/
├─ tools/                # original CLI scripts (kept intact)
│  ├─ final_analysis.py
│  ├─ final_analysis_optical.py
│  ├─ gen_qe_dataset.py
│  ├─ analyze_qe_dataset.py
│  ├─ press_start_qe.py
│  └─ (optional) press_start.py
├─ src/timesym/          # light Python shims
│  ├─ __init__.py        # exposes __version__
│  ├─ cable.py           # run_cable()
│  └─ optical.py         # run_optical()
├─ data/                 # (created as needed)
├─ runs/                 # outputs land here
├─ cal/                  # calibration inputs (optional)
├─ paper_assets/         # optional figure export target
├─ pyproject.toml
├─ LICENSE
└─ README.md
```

---

## Reproducibility & provenance

- All analyses use a **registered ROPE/HDI** decision at a geometry‑fixed time \(t_\star\).
- Scripts emit CSV summaries and figures; seeds are exposed on the CLI for determinism.
- Windows/macOS/Linux are supported; see CI badge for smoke‑test status.

### Unicode on Windows (console)
If you see encoding errors when printing symbols (e.g. σ̂), run with UTF‑8:
```powershell
$env:PYTHONIOENCODING="utf-8"
python tools\analyze_qe_dataset.py --help
# or switch console to UTF-8 once:
chcp 65001 > $null
```

---

## Citing

If you use this software, please cite the **Zenodo concept DOI** (covers all versions):

**Concept DOI:** https://doi.org/10.5281/zenodo.17074673

For strict reproducibility, also cite a **versioned DOI** from the “Versions” panel on Zenodo (each GitHub release archived by Zenodo gets a unique DOI).

Example (APA):
> Mishchenko, I. (2025). *timesym — registered temporal‑interferometer analysis (RF/optical) + DCQE validation* (v0.2.x) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.17074673

### CITATION.cff
A `CITATION.cff` is included so GitHub shows a “Cite this repository” box automatically. Keep `version:` and `date-released:` in sync with your releases.

---

## Releasing

```bash
# bump version in pyproject.toml (e.g., 0.2.2), optionally update CITATION.cff
git add -A && git commit -m "release: v0.2.2"
git tag v0.2.2
git push origin main --tags

# create GitHub Release (auto‑notes)
gh release create v0.2.2 --generate-notes
```
Zenodo will archive the release and attach a version‑specific DOI automatically (if GitHub‑Zenodo linking is enabled).

---

## License

This project is licensed under the **MIT License**. See [LICENSE](LICENSE).