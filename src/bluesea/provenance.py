"""Provenance captured at three levels, so a result can be re-derived.

1. the container digest and tool versions,
2. the git commit of the workflow and the config files it ran with,
3. for every API call, the model id, the full prompt and the raw response
   (handled by `llm.CachedClient`).

The third is what makes the LLM-judged scores defensible. A reviewer who
doubts a literature score can read exactly what was asked and exactly what
came back.
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
from pathlib import Path
from typing import Any

from . import __version__


def git_commit(repo: str | Path = ".") -> str:
    try:
        result = subprocess.run(
            ["git", "-C", str(repo), "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
        )
        commit = result.stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "unknown"

    dirty = subprocess.run(
        ["git", "-C", str(repo), "status", "--porcelain"],
        capture_output=True,
        text=True,
        check=False,
    ).stdout.strip()
    # A dirty tree means the commit alone does not identify what ran, so say so
    # rather than recording a commit that is not what produced the result.
    return f"{commit}-dirty" if dirty else commit


def file_digest(path: str | Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def collect(
    config_path: str | Path,
    weights_path: str | Path,
    tool_versions: dict[str, str] | None = None,
    repo: str | Path = ".",
) -> dict[str, Any]:
    return {
        "bluesea_version": __version__,
        "git_commit": git_commit(repo),
        "config_sha256": file_digest(config_path),
        "weights_sha256": file_digest(weights_path),
        "container": os.environ.get("APPTAINER_CONTAINER", "none"),
        "container_digest": os.environ.get("BLUESEA_CONTAINER_DIGEST", "unknown"),
        "python": platform.python_version(),
        "host": platform.node(),
        "tool_versions": tool_versions or {},
    }


def write(provenance: dict[str, Any], path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
