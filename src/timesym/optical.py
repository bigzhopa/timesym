import subprocess, sys
from pathlib import Path
from typing import Optional

def run_optical(out_dir: str = "runs/optical",
                simulate: bool = True,
                seed: Optional[int] = None,
                Vref: Optional[str] = None,
                repo_root: Optional[str] = None,
                **kwargs):
    """
    Run the optical analysis pipeline by invoking tools/final_analysis_optical.py.

    Args:
        out_dir: output directory (relative to repo_root or CWD)
        simulate: if True, pass --simulate
        seed: RNG seed (optional)
        Vref: normalization ('mainpeak', 'mainlobe', or 'custom'), optional
        repo_root: path to the repo root that contains `tools/`
        **kwargs: extra CLI flags to forward, e.g. N_trials=..., p_pair=..., null_depth_dB=...

    Returns: CompletedProcess
    """
    rr = Path(repo_root) if repo_root else Path.cwd()
    script = rr / "tools" / "final_analysis_optical.py"
    if not script.exists():
        raise FileNotFoundError(f"Could not find {script} — set repo_root to your repo path.")
    cmd = [sys.executable, str(script)]
    if simulate:
        cmd += ["--simulate"]
    if seed is not None:
        cmd += ["--seed", str(seed)]
    if Vref is not None:
        cmd += ["--Vref", str(Vref)]
    # forward simple kwargs as --key value
    for k, v in kwargs.items():
        flag = f"--{k}".replace("_", "-")
        cmd += [flag, str(v)]
    cmd += ["--out", out_dir]
    return subprocess.run(cmd, cwd=str(rr), check=True)
