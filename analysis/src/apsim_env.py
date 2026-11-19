"""APSIM installation discovery.

Keeps the pipeline independent of which APSIM build is installed:

  1. ``APSIM_BIN`` environment variable wins if set
     (e.g. set APSIM_BIN=C:\\...\\APSIM2026.2.7980.0\\bin).
  2. otherwise every ``APSIM<version>`` folder under the usual Windows
     install roots is scanned and the highest version is used.

``apsimNGpy`` (optional) is configured from the same resolution, so switching
builds needs no code change in steps 04-07.
"""
from __future__ import annotations

import os
import re
from dataclasses import dataclass
from pathlib import Path

SEARCH_ROOTS = [
    Path(os.environ.get("LOCALAPPDATA", "")) / "Programs",
    Path(r"C:/Program Files"),
    Path(r"C:/APSIM"),
    Path(os.environ.get("LOCALAPPDATA", "")) / "APSIM",
]

_VERSION_RE = re.compile(r"APSIM(?:NG)?[ _-]?(\d{4})\.(\d+)\.(\d+)(?:\.(\d+))?", re.I)


@dataclass(frozen=True)
class ApsimInstall:
    bin_dir: Path
    version: str
    key: tuple[int, ...]

    @property
    def models_exe(self) -> Path:
        return self.bin_dir / "Models.exe"

    @property
    def examples(self) -> Path:
        return self.bin_dir.parent / "Examples"

    @property
    def has_apsim_core(self) -> bool:
        """APSIM.Core.dll arrived with the 2026.x line; apsimNGpy 1.5+ needs it."""
        return any(self.bin_dir.glob("APSIM.[Cc]ore.dll"))


def _candidates() -> list[ApsimInstall]:
    out: list[ApsimInstall] = []
    for root in SEARCH_ROOTS:
        if not root or not root.exists():
            continue
        for d in root.iterdir():
            if not d.is_dir():
                continue
            m = _VERSION_RE.match(d.name)
            if not m:
                continue
            bin_dir = d / "bin"
            if not (bin_dir / "Models.exe").exists():
                continue
            ver = ".".join(g for g in m.groups() if g)
            key = tuple(int(g) for g in m.groups() if g)
            out.append(ApsimInstall(bin_dir=bin_dir, version=ver, key=key))
    return out


def find_all() -> list[ApsimInstall]:
    return sorted(_candidates(), key=lambda i: i.key)


def resolve(prefer_latest: bool = True) -> ApsimInstall:
    env = os.environ.get("APSIM_BIN")
    if env:
        bin_dir = Path(env)
        if not (bin_dir / "Models.exe").exists():
            raise FileNotFoundError(f"APSIM_BIN={env} does not contain Models.exe")
        m = _VERSION_RE.match(bin_dir.parent.name) or _VERSION_RE.match(bin_dir.name)
        ver = ".".join(g for g in m.groups() if g) if m else "unknown"
        key = tuple(int(g) for g in m.groups() if g) if m else (0,)
        return ApsimInstall(bin_dir=bin_dir, version=ver, key=key)

    found = _candidates()
    if not found:
        raise FileNotFoundError(
            "no APSIM installation found; set APSIM_BIN to the folder containing Models.exe")
    return found[-1] if prefer_latest else found[0]


def configure_apsimngpy(install: ApsimInstall | None = None) -> str:
    """Point apsimNGpy at the resolved build (no-op if the package is absent)."""
    install = install or resolve()
    os.environ["APSIM_BIN_PATH"] = str(install.bin_dir)
    try:
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                from apsimNGpy.config import set_apsim_bin_path  # type: ignore
            except ImportError:
                from apsimNGpy.core.config import set_apsim_bin_path  # type: ignore
            set_apsim_bin_path(str(install.bin_dir))
    except Exception as exc:  # pragma: no cover - optional dependency
        return f"apsimNGpy not configured: {type(exc).__name__}: {exc}"
    return f"apsimNGpy -> {install.bin_dir}"


if __name__ == "__main__":
    for i in find_all():
        print(f"{i.version:>22}  core={i.has_apsim_core!s:<5}  {i.bin_dir}")
    cur = resolve()
    print(f"\nresolved: {cur.version}  ({cur.bin_dir})")
    print(f"apsimNGpy compatible: {cur.has_apsim_core}")
    print(configure_apsimngpy(cur))
