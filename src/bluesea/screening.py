"""Bounded pilot inputs and actual host-proteome spike controls."""

from __future__ import annotations

import hashlib
import json
import random
from pathlib import Path
from typing import Any

from . import io
from .cascade import estimate_low_complexity
from .provenance import file_digest

HOST_PREFIX = "BLUESEA_HOST_SPIKE_"


def signature(config: dict[str, Any], weights: dict[str, Any]) -> str:
    keys = ("species", "host", "cascade", "orthology", "predictors", "validation")
    payload = {key: config[key] for key in keys}
    payload["weights"] = weights
    # Curation changes invalidate pilot approval even if paths stay the same.
    for key in ("positive_controls", "negative_controls", "construct_controls"):
        path = config["validation"].get(key)
        if path:
            payload[key] = Path(path).read_text()
            if key == "construct_controls":
                payload["construct_sequence_sha256"] = {
                    c["construct_id"]: (
                        file_digest(c["sequence_file"])
                        if c.get("sequence_file") and Path(c["sequence_file"]).is_file()
                        else None
                    )
                    for c in json.loads(payload[key])["constructs"]
                }
    root = Path(config.get("workdir", "."))
    reference = root / config["paths"]["reference"]
    evidence_files = [reference / config["host"]["proteome"]]
    for species in config["species"].values():
        evidence_files.extend(
            reference / species[key]
            for key in ("proteome", "gff", "expression", "rnai", "secretome", "gpi")
            if species.get(key)
        )
    for tool in config.get("predictors", {}).values():
        evidence_files.extend(Path(path) for path in tool.get("required_files", []))
    payload["code_sha256"] = {
        path.name: file_digest(path) for path in Path(__file__).parent.glob("*.py")
    }
    payload["file_sha256"] = {
        str(path): file_digest(path) if path.is_file() else None
        for path in evidence_files
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()


def require_pilot(config: dict[str, Any], weights: dict[str, Any]) -> None:
    if config.get("screening", {}).get("mode", "pilot") != "full":
        return
    path = Path(config["screening"]["pilot_report"])
    report = json.loads(path.read_text())
    if (
        not report["gate_passed"]
        or report.get("screening_mode") != "pilot"
        or report.get("screening_signature") != signature(config, weights)
    ):
        raise ValueError(
            "Full run requires a passed pilot with matching biology, tools, "
            "controls and weights"
        )


def prepare(config: dict[str, Any], species: str, weights: dict[str, Any]) -> None:
    if config.get("screening", {}).get("mode", "pilot") not in {"pilot", "full"}:
        raise ValueError("screening.mode must be pilot or full")
    require_pilot(config, weights)
    root = Path(config.get("workdir", "."))
    reference = root / config["paths"]["reference"]
    scratch = root / config["paths"]["scratch"] / species
    cfg = config.get("screening", {})
    records = list(io.read_fasta(reference / config["species"][species]["proteome"]))
    if cfg.get("mode", "pilot") == "pilot":
        wanted = set()
        for key in ("positive_controls", "negative_controls"):
            for row in io.read_tsv(config["validation"][key]):
                wanted.update(
                    [row["protein_id"], *filter(None, row.get("isoforms", "").split(";"))]
                )
        path = config["validation"].get("construct_controls")
        if path:
            for construct in json.loads(Path(path).read_text())["constructs"]:
                if construct["species"] == species:
                    for component in construct["louse_components"]:
                        wanted.update(component["mapped_accessions"])
        controls = [r for r in records if r[0] in wanted]
        background = [r for r in records if r[0] not in wanted]
        rng = random.Random(f"{cfg.get('seed', 0)}:{species}")
        records = controls + rng.sample(
            background, min(cfg.get("random_n", 200), len(background))
        )
    spikes = []
    if species == next(iter(config["species"])):
        for pid, description, sequence in io.read_fasta(
            reference / config["host"]["proteome"]
        ):
            if len(spikes) >= config["validation"].get("host_spike_n", 3):
                break
            length = config["cascade"]["length"]
            if not length["min_aa"] <= len(sequence) <= length["max_aa"]:
                continue
            if (
                estimate_low_complexity(sequence)
                > config["cascade"]["low_complexity"]["max_fraction"]
            ):
                continue
            spike_id = HOST_PREFIX + pid
            records.append(
                (spike_id, f"benchmark-only host spike: {description}", sequence)
            )
            spikes.append(
                {
                    "protein_id": spike_id,
                    "host_accession": pid,
                    "class": "host_spike",
                    "expected_drop_step": "host_similarity",
                }
            )
        if len(spikes) < config["validation"].get("host_spike_n", 3):
            raise ValueError("Insufficient host sequences to build spike controls")
    if any(
        pid.startswith(HOST_PREFIX)
        for pid, _, _ in records
        if pid not in {s["protein_id"] for s in spikes}
    ):
        raise ValueError("Input proteome uses reserved host-spike prefix")
    scratch.mkdir(parents=True, exist_ok=True)
    (scratch / "screen.faa").write_text(
        "".join(f">{pid} {desc}\n{seq}\n" for pid, desc, seq in records)
    )
    (scratch / "host_spikes.json").write_text(json.dumps(spikes, indent=2))
