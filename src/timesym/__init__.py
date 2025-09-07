"""timesym: shims to invoke the existing CLI tools programmatically.
These helpers shell out to the scripts kept under `tools/` in the repo.
For robustness, pass `repo_root` if your working directory is not the repo root.
"""
from .cable import run_cable
from .optical import run_optical

__all__ = ["run_cable", "run_optical"]
