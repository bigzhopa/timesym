"""
timesym: registered temporal-interferometer analysis (RF/optical) + DCQE validation.
"""

# Keep the package version in sync with pyproject.toml
try:
    from importlib.metadata import version as _pkg_version, PackageNotFoundError
except Exception:  # very old Python fallback (not expected on 3.10+)
    _pkg_version = None
    PackageNotFoundError = Exception  # type: ignore

try:
    __version__ = _pkg_version("timesym") if _pkg_version else "0+unknown"
except PackageNotFoundError:
    # Package not installed (e.g., running from a checkout without `pip install -e .`)
    __version__ = "0+unknown"

__all__ = ["__version__"]
