import subprocess, sys
from pathlib import Path
from typing import Optional

def run_cable(out_dir: str = "runs/cable",
              simulate: bool = True,
              Vpp: float = 1.0,
              seed: Optional[int] = None,
              repo_root: Optional[str] = None):
    """
    Run the RF cable analysis pipeline by invoking tools/final_analysis.py.

    Args:
        out_dir: output directory (relative to repo_root or CWD)
        simulate: if True, pass --simulate
        Vpp: peak-to-peak normalization (volts) when simulating
        seed: RNG seed (optional)
        repo_root: path to the repo root that contains `tools/`

    Returns: CompletedProcess
    """
    rr = Path(repo_root) if repo_root else Path.cwd()
    script = rr / "tools" / "final_analysis.py"
    if not script.exists():
        raise FileNotFoundError(f"Could not find {script} — set repo_root to your repo path.")
    cmd = [sys.executable, str(script)]
    if simulate:
        cmd += ["--simulate"]
        cmd += ["--Vpp", str(Vpp)]
        if seed is not None:
            cmd += ["--seed", str(seed)]
    cmd += ["--out", out_dir]
    return subprocess.run(cmd, cwd=str(rr), check=True)
