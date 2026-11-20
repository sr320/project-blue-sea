"""Command-line entry points, one per workflow stage.

Stages hand off through JSON on disk rather than in memory, because the
deterministic stages run on compute nodes without network and the API stage
runs somewhere else entirely. Every stage is independently restartable.

    bluesea prefilter  --species lsal      # cheap filters -> survivor FASTA
    bluesea assemble   --species lsal      # tool output -> records.json
    bluesea cascade    --species lsal      # apply filters, write audit log
    bluesea orthology                      # RBH between the two species
    bluesea literature --species lsal      # API stage (dry-run by default)
    bluesea rank                           # score, intersect, write ranking
    bluesea dossiers                       # render the shortlist
    bluesea benchmark                      # positive/negative/robustness
    bluesea symbiont-gate                  # go/no-go on the symbiont track
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any

from . import io, orthology, provenance, scoring, validation
from .records import ProteinRecord, from_dict


def _paths(config: dict[str, Any]) -> dict[str, Path]:
    root = Path(config.get("workdir", "."))
    return {key: root / value for key, value in config["paths"].items()}


def _records_path(paths: dict[str, Path], species: str, stage: str) -> Path:
    return paths["scratch"] / species / f"{stage}.json"


def _load_records(path: Path) -> list[ProteinRecord]:
    with open(path) as handle:
        return [from_dict(item) for item in json.load(handle)]


def _save_records(records: list[ProteinRecord], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as handle:
        json.dump([record.to_dict() for record in records], handle)


# -- prefilter -------------------------------------------------------------


def cmd_prepare(args: argparse.Namespace) -> int:
    from .screening import prepare

    prepare(io.load_config(args.config), args.species, io.load_config(args.weights))
    return 0


def cmd_verify_pilot(args: argparse.Namespace) -> int:
    from .screening import require_pilot

    require_pilot(io.load_config(args.config), io.load_config(args.weights))
    return 0


def cmd_check_gate(args: argparse.Namespace) -> int:
    report = json.loads(Path(args.report).read_text())
    if not report["gate_passed"]:
        print("Gate failed; inspect the preserved benchmark report.", file=sys.stderr)
        return 1
    return 0


def cmd_prefilter(args: argparse.Namespace) -> int:
    """Run the cheap half of the cascade and emit the survivor FASTA.

    This is the step that keeps epitope prediction and the salmon homology
    search off the ~8 proteins in 9 that cannot reach host antibody anyway.
    It needs only SignalP and DeepTMHMM, which have to run proteome-wide
    because the exposure filter is built on them.

    Every protein -- dropped or not -- is written to `prefiltered.json` with
    its verdicts, so the audit log and the negative controls still account for
    the whole proteome once the later stages have run.
    """
    from .cascade import run_cascade, steps_for

    config = io.load_config(args.config)
    paths = _paths(config)
    species_cfg = config["species"][args.species]
    cascade_cfg = config["cascade"]
    reference = paths["reference"]
    scratch = paths["scratch"] / args.species

    from .predictors import require_ids

    signalp = io.read_signalp(scratch / "signalp.tsv")
    tmhmm = io.read_deeptmhmm(scratch / "deeptmhmm.gff3")
    secretome = (
        _read_id_list(reference / species_cfg["secretome"])
        if species_cfg.get("secretome")
        else set()
    )
    gpi = (
        _read_id_list(reference / species_cfg["gpi"]) if species_cfg.get("gpi") else set()
    )

    # The GFF is the only authoritative protein-to-gene link; NCBI protein
    # FASTA headers carry no gene tag, so without it every protein looks like
    # its own gene and nothing keyed on genes can join.
    gene_map: dict[str, tuple[str, str]] = {}
    gff = species_cfg.get("gff")
    if gff and (reference / gff).exists():
        gene_map = io.read_protein_gene_map(reference / gff)
        print(f"gene map: {len(gene_map)} proteins from {gff}")
    else:
        print(
            f"WARNING: no GFF for {args.species}; gene ids fall back to protein "
            f"accessions and gene-keyed evidence will not join. See "
            f"docs/reference-data.md.",
            file=sys.stderr,
        )

    from .screening import HOST_PREFIX

    fasta = list(io.read_fasta(scratch / "screen.faa"))
    expected = {pid for pid, _, _ in fasta}
    if len(expected) != len(fasta):
        raise ValueError("Duplicate proteome identifiers")
    require_ids(set(signalp), expected, "SignalP")
    require_ids(set(tmhmm), expected, "DeepTMHMM")
    records = []
    for protein_id, description, sequence in fasta:
        if tmhmm[protein_id]["sequence"].upper() != sequence.upper():
            raise ValueError(f"DeepTMHMM sequence mismatch: {protein_id}")
        gene_id, gene_symbol = gene_map.get(
            protein_id, (_gene_of(protein_id, description), "")
        )
        record = ProteinRecord(
            protein_id=protein_id,
            gene_id=gene_id,
            gene_symbol=gene_symbol,
            species=args.species,
            sequence=sequence,
            description=description,
            benchmark_only=protein_id.startswith(HOST_PREFIX),
            expected_stages=io.species_stages(config, args.species),
            early_stages=species_cfg.get("early_stages", []),
            topology=io.topology_from_tools(protein_id, signalp, tmhmm, secretome, gpi),
        )
        run_cascade(
            record,
            cascade_cfg,
            stop_on_fail=not args.profile,
            phase="prefilter",
        )
        records.append(record)

    _save_records(records, _records_path(paths, args.species, "prefiltered"))

    # With the prefilter disabled the full proteome goes to the predictors and
    # nothing is skipped. It is the expensive way to run, and the point of it
    # is to measure what the exposure filter costs in false negatives -- pair
    # it with `cascade --profile` so the later steps still return a verdict on
    # proteins exposure rejected.
    enabled = cascade_cfg.get("prefilter", {}).get("enabled", True)
    forwarded = [r for r in records if r.survived] if enabled else list(records)

    out_fasta = scratch / "candidates.faa"
    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    with open(out_fasta, "w") as handle:
        for record in forwarded:
            handle.write(f">{record.protein_id} {record.description}\n")
            for start in range(0, len(record.sequence), 60):
                handle.write(record.sequence[start : start + 60] + "\n")

    for step_name, _, _ in steps_for("prefilter"):
        dropped = sum(1 for r in records if r.dropped_at == step_name)
        print(f"{step_name:<16} dropped {dropped:>6}")
    kept = sum(1 for r in records if r.survived)
    print(
        f"{kept}/{len(records)} proteins reach the predictors "
        f"({kept / len(records):.1%}) -> {out_fasta}"
        if records
        else f"no proteins read from {species_cfg['proteome']}"
    )

    if not forwarded:
        # An empty FASTA would produce empty predictor output, an empty
        # ranking and no error anywhere. Stop here instead.
        print(
            "\nNo proteins survived the prefilter. Check the SignalP and "
            "DeepTMHMM output before rerunning.",
            file=sys.stderr,
        )
        return 1
    return 0


# -- assemble --------------------------------------------------------------


def cmd_assemble(args: argparse.Namespace) -> int:
    """Attach the candidate-stage predictor output to the prefiltered records.

    The epitope and salmon-homology tables cover only the proteins the
    prefilter forwarded, so proteins dropped before this point simply carry no
    such evidence. They stay in the record set to keep the audit complete.
    """
    config = io.load_config(args.config)
    paths = _paths(config)
    species_cfg = config["species"][args.species]
    reference = paths["reference"]
    scratch = paths["scratch"] / args.species

    records = _load_records(_records_path(paths, args.species, "prefiltered"))

    epitopes = io.read_epitopes(scratch / "epitopes.tsv")
    host_hits = io.read_blast_tabular(scratch / "vs_host.tsv")
    stages = io.species_stages(config, args.species)
    expression = (
        io.read_expression(reference / species_cfg["expression"], stages)
        if species_cfg.get("expression")
        else {}
    )
    rnai = (
        io.read_rnai(reference / species_cfg["rnai"]) if species_cfg.get("rnai") else {}
    )

    # Mirror the prefilter's own forwarding rule, so the records enriched here
    # are exactly the ones the predictors were given. A protein can legitimately
    # be absent from both tables -- no predicted epitopes, no salmon homolog --
    # and must still pick up its expression and RNAi evidence.
    enabled = config["cascade"].get("prefilter", {}).get("enabled", True)

    enriched = joined_expression = joined_rnai = 0
    key_hits: Counter[str] = Counter()
    for record in records:
        if enabled and not record.survived:
            continue
        record.epitopes = epitopes.get(record.protein_id, [])
        hits = host_hits.get(record.protein_id, [])
        record.host_hits = io.host_alignments(hits, record.length)
        record.host_hit = io.best_host_hit(
            hits, record.length, config["cascade"]["host_similarity"]
        )
        record.assayed_stages = (
            io.expression_stages(reference / species_cfg["expression"], stages)
            if species_cfg.get("expression")
            else []
        )

        value, key = _lookup_by_gene(expression, record)
        record.expression = value or {}
        if value is not None:
            joined_expression += 1
            key_hits[f"expression:{key}"] += 1

        value, key = _lookup_by_gene(rnai, record)
        if value is not None:
            record.rnai.phenotype = value["phenotype"]
            record.rnai.source = value["source"]
            record.rnai.note = value["note"]
            joined_rnai += 1
            key_hits[f"rnai:{key}"] += 1
        enriched += 1

    out = _records_path(paths, args.species, "assembled")
    _save_records(records, out)
    print(f"assembled {enriched} candidates of {len(records)} proteins -> {out}")
    _report_joins(
        enriched,
        [("expression", expression, joined_expression), ("rnai", rnai, joined_rnai)],
    )
    if any(
        table and joined == 0
        for table, joined in [(expression, joined_expression), (rnai, joined_rnai)]
    ):
        return 1
    for label, count in sorted(key_hits.items()):
        # Which key form carried the join: a table joining only on protein_id
        # is keyed per isoform, not per gene, and will miss siblings.
        print(f"  joined via {label}: {count}")
    return 0


def _gene_keys(record: ProteinRecord) -> list[tuple[str, str]]:
    """(label, key) pairs to try when joining gene-keyed evidence.

    Which of these the expression and RNAi tables use is not knowable from
    here -- they come from LiceBase and the literature, not from the assembly
    -- so rather than guess one, try the few forms that exist and report which
    one matched. Most specific first.
    """
    keys: list[tuple[str, str]] = []
    if record.gene_id:
        keys.append(("gene_id", record.gene_id))
        if record.gene_id.startswith("GeneID:"):
            # Tables that carry the bare NCBI number rather than the prefixed
            # form are common enough to be worth one extra lookup.
            keys.append(("geneid_number", record.gene_id.split(":", 1)[1]))
    if record.gene_symbol:
        keys.append(("gene_symbol", record.gene_symbol))
    keys.append(("protein_id", record.protein_id))

    seen, unique = set(), []
    for label, key in keys:
        if key not in seen:
            seen.add(key)
            unique.append((label, key))
    return unique


def _lookup_by_gene(table: dict[str, Any], record: ProteinRecord) -> tuple[Any, str]:
    for label, key in _gene_keys(record):
        if key in table:
            return table[key], label
    return None, ""


def _report_joins(enriched: int, tables: list[tuple[str, dict[str, Any], int]]) -> None:
    """Say how much of each gene-keyed table actually matched.

    Expression and RNAi join on `gene_id`, which comes from the `[gene=...]`
    tag in the FASTA header -- a tag NCBI `protein.faa` files do not carry, so
    `_gene_of` falls back to the protein accession. If the tables are keyed on
    anything else, nothing joins: every protein scores `untested` on
    essentiality, the heaviest criterion in the rubric, and the run completes
    looking healthy. A zero join rate has to be loud.
    """
    for name, table, joined in tables:
        if not table:
            continue
        rate = joined / enriched if enriched else 0.0
        print(f"  {name}: {joined}/{enriched} joined ({rate:.1%})")
        if joined == 0:
            print(
                f"\nWARNING: no record joined to the {name} table. Its keys do "
                f"not match the gene ids taken from the proteome headers, so "
                f"this evidence is being silently dropped. See "
                f"docs/reference-data.md.",
                file=sys.stderr,
            )


def _read_id_list(path: Path) -> set[str]:
    return {line.strip() for line in path.read_text().splitlines() if line.strip()}


def _gene_of(protein_id: str, description: str) -> str:
    """NCBI-style headers carry the gene in a [gene=...] tag; fall back to the
    protein id so a missing tag never silently collapses genes together."""
    marker = "[gene="
    if marker in description:
        start = description.index(marker) + len(marker)
        return description[start : description.index("]", start)]
    return protein_id


# -- cascade ---------------------------------------------------------------


def cmd_cascade(args: argparse.Namespace) -> int:
    from .cascade import cascade_summary, run_cascade

    config = io.load_config(args.config)
    paths = _paths(config)
    records = _load_records(_records_path(paths, args.species, "assembled"))

    # The prefilter already recorded verdicts for length, low_complexity and
    # exposure; this appends the rest in cascade order, so the verdict list
    # reads exactly as it did when one process ran all seven steps.
    for record in records:
        run_cascade(
            record, config["cascade"], stop_on_fail=not args.profile, phase="post"
        )

    survivors = [record for record in records if record.survived]
    _save_records(records, _records_path(paths, args.species, "cascaded"))
    _save_records(survivors, _records_path(paths, args.species, "survivors"))

    summary = cascade_summary(records)
    io.write_tsv(
        paths["results"] / args.species / "cascade_summary.tsv",
        summary,
        ["step", "entered", "dropped", "survived", "drop_rate"],
    )
    # Per-protein audit: every verdict, so any drop traces to one rule.
    audit = paths["logs"] / args.species / "cascade_audit.jsonl"
    audit.parent.mkdir(parents=True, exist_ok=True)
    with open(audit, "w") as handle:
        for record in records:
            handle.write(
                json.dumps(
                    {
                        "protein_id": record.protein_id,
                        "dropped_at": record.dropped_at,
                        "verdicts": record.verdicts,
                    }
                )
                + "\n"
            )

    for row in summary:
        print(
            f"{row['step']:<16} entered {row['entered']:>6} "
            f"dropped {row['dropped']:>6} ({row['drop_rate']:.1%})"
        )
    print(f"{len(survivors)} survivors -> {audit}")
    return 0


# -- orthology -------------------------------------------------------------


def cmd_orthology(args: argparse.Namespace) -> int:
    config = io.load_config(args.config)
    paths = _paths(config)
    species = list(config["species"])
    if len(species) == 1:
        name = species[0]
        records = _load_records(_records_path(paths, name, "cascaded"))
        _save_records(records, _records_path(paths, name, "orthologized"))
        return 0
    if len(species) != 2:
        raise ValueError("orthology expects one or two species")
    first, second = species
    forward = io.read_blast_tabular(paths["scratch"] / f"{first}_vs_{second}.tsv")
    reverse = io.read_blast_tabular(paths["scratch"] / f"{second}_vs_{first}.tsv")
    for source, partner, fwd, rev in (
        (first, second, forward, reverse),
        (second, first, reverse, forward),
    ):
        pairs = orthology.reciprocal_best_hits(fwd, rev, config["orthology"])
        records = {
            r.protein_id: r
            for r in _load_records(_records_path(paths, source, "cascaded"))
        }
        sequences = {
            pid: seq
            for pid, _, seq in io.read_fasta(
                paths["reference"] / config["species"][partner]["proteome"]
            )
        }
        for record in records.values():
            record.orthologs = []
        orthology.attach_orthologs(
            records, pairs, partner, config["orthology"], sequences
        )
        _save_records(
            list(records.values()), _records_path(paths, source, "orthologized")
        )
        print(f"{source}: {len(pairs)} cross-species matches")
    return 0


# -- literature (API stage) ------------------------------------------------


def cmd_literature(args: argparse.Namespace) -> int:
    from .literature import assess_all, load_sources
    from .llm import client_from_config

    config = io.load_config(args.config)
    paths = _paths(config)
    client = client_from_config(
        config,
        {"cache": paths["cache"]},
        dry_run=not args.live,
    )
    records = [
        r
        for r in _load_records(_records_path(paths, args.species, "orthologized"))
        if r.survived and not r.benchmark_only
    ]
    assess_all(
        records,
        client,
        blinded=args.blinded,
        blinded_terms=config["llm"]["blinded_terms"],
        sources_by_protein=load_sources(config.get("literature_sources", [])),
    )
    _save_records(records, _records_path(paths, args.species, "assessed"))

    report = client.report()
    (paths["logs"] / f"llm_usage_{args.species}.json").parent.mkdir(
        parents=True, exist_ok=True
    )
    (paths["logs"] / f"llm_usage_{args.species}.json").write_text(
        json.dumps(report, indent=2)
    )
    print(json.dumps(report, indent=2))
    if report["dry_run"]:
        print(
            "\nDry run: nothing was called. Re-run with --live to spend.",
            file=sys.stderr,
        )
    return 0


# -- rank ------------------------------------------------------------------


def cmd_rank(args: argparse.Namespace) -> int:
    config = io.load_config(args.config)
    weights_cfg = io.load_config(args.weights)
    paths = _paths(config)

    per_species: dict[str, list[ProteinRecord]] = {}
    for species in config["species"]:
        records = _load_records(_records_path(paths, species, "assessed"))
        per_species[species] = scoring.rank(records, weights_cfg)

    # Flag the intersection before writing: a cross-species candidate is the
    # set worth reading first, and the flag has to be on the record for the
    # ranking table and the dossier to show it.
    names = list(per_species)
    if len(names) == 2:
        shared = orthology.intersect_shortlists(
            per_species[names[0]], per_species[names[1]], config["orthology"]
        )
        print(f"{len(shared)} cross-species candidates")

    all_records = [r for rows in per_species.values() for r in rows]
    representatives = scoring.gene_representatives(all_records)
    sensitivity = validation.sensitivity_report(
        representatives, weights_cfg, config["validation"]
    )
    (paths["results"]).mkdir(parents=True, exist_ok=True)
    (paths["results"] / "ranking_sensitivity.json").write_text(
        json.dumps(sensitivity, indent=2)
    )
    for record in all_records:
        if any(
            r.species == record.species
            and r.gene_id == record.gene_id
            and "single_source_dependent" in r.flags
            for r in representatives
        ):
            record.flag("single_source_dependent")
    everything: list[ProteinRecord] = []
    for species, ranked in per_species.items():
        _save_records(ranked, _records_path(paths, species, "ranked"))
        io.write_tsv(
            paths["results"] / species / "ranked.tsv",
            scoring.ranking_rows(scoring.gene_representatives(ranked)),
            scoring.RANKING_COLUMNS,
        )
        everything.extend(scoring.gene_representatives(ranked))

    combined = sorted(everything, key=lambda r: (-r.total_score, r.species, r.protein_id))
    io.write_tsv(
        paths["results"] / "ranked_all.tsv",
        scoring.ranking_rows(combined),
        scoring.RANKING_COLUMNS,
    )
    provenance.write(
        provenance.collect(args.config, args.weights),
        paths["results"] / "provenance.json",
    )
    print(f"ranked {len(combined)} candidates -> {paths['results']}/ranked_all.tsv")
    return 0


# -- dossiers --------------------------------------------------------------


def cmd_dossiers(args: argparse.Namespace) -> int:
    from .dossier import write_dossiers

    config = io.load_config(args.config)
    weights_cfg = io.load_config(args.weights)
    paths = _paths(config)

    records: list[ProteinRecord] = []
    for species in config["species"]:
        path = _records_path(paths, species, "ranked")
        records.extend(scoring.gene_representatives(_load_records(path)))
    records.sort(key=lambda r: r.total_score, reverse=True)

    trail = provenance.collect(args.config, args.weights)
    written = write_dossiers(
        records[: args.top], paths["results"] / "dossiers", weights_cfg["weights"], trail
    )
    print(f"wrote {len(written)} dossiers -> {paths['results']}/dossiers")
    return 0


# -- benchmark -------------------------------------------------------------


def cmd_benchmark(args: argparse.Namespace) -> int:
    config = io.load_config(args.config)
    weights_cfg = io.load_config(args.weights)
    paths = _paths(config)
    validation_cfg = config["validation"]

    records: list[ProteinRecord] = []
    for species in config["species"]:
        path = _records_path(paths, species, "orthologized")
        records.extend(_load_records(path))

    positives = {
        row["protein_id"]: _control_accessions(row)
        for row in io.read_tsv(validation_cfg["positive_controls"])
    }
    negatives = {
        row["protein_id"]: (row["expected_drop_step"], _control_accessions(row))
        for row in io.read_tsv(validation_cfg["negative_controls"])
    }
    for species in config["species"]:
        spike_path = paths["scratch"] / species / "host_spikes.json"
        if spike_path.exists():
            for row in json.loads(spike_path.read_text()):
                negatives[row["protein_id"]] = (
                    row["expected_drop_step"],
                    {row["protein_id"]},
                )
    random_ids = _sample_random_controls(records, validation_cfg, positives, negatives)

    report = validation.run_benchmark(
        records, weights_cfg, validation_cfg, positives, negatives, random_ids
    )
    constructs_path = validation_cfg.get("construct_controls")
    if constructs_path:
        with open(constructs_path) as handle:
            constructs = json.load(handle)
        report["experimental_constructs"] = validation.construct_report(
            records, weights_cfg, validation_cfg, constructs
        )
        if not positives:
            report["gate_passed"] = (
                report["negative_controls"]["passed"]
                and report["sensitivity_validation"]["passed"]
                and report["experimental_constructs"]["recovery_passed"]
            )
        report["gate_passed"] &= report["experimental_constructs"]["ready"]
    if validation_cfg.get("require_host_spikes", False):
        generated = [pid for pid in negatives if pid.startswith("BLUESEA_HOST_SPIKE_")]
        report["host_spikes_present"] = len(generated) >= validation_cfg.get(
            "host_spike_n", 3
        ) and all(any(r.protein_id == pid for r in records) for pid in generated)
        report["gate_passed"] &= report["host_spikes_present"]
    from .screening import signature

    report["screening_signature"] = signature(config, weights_cfg)
    report["screening_mode"] = config.get("screening", {}).get("mode", "pilot")
    report["input_stage"] = "orthologized"
    out = paths["results"] / "benchmark.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2))

    print(json.dumps({k: v for k, v in report.items() if k != "criteria"}, indent=2))
    if not report["gate_passed"]:
        # Non-zero exit stops the workflow here on purpose. A failed gate means
        # rebuild the rubric, not press on to the full run.
        print(
            "\nGATE FAILED: control completeness or antigen recovery did not pass "
            "on deterministic evidence. Rebuild the rubric before the full run.",
            file=sys.stderr,
        )
        return 0 if args.report_only else 1
    print("\nGate passed.")
    return 0


def _control_accessions(row: dict[str, str]) -> set[str]:
    """Every accession representing one control.

    A control is a gene, but the proteome carries one entry per isoform. The
    optional `isoforms` column lists the siblings, semicolon separated, so a
    control is not scored as missing when a sibling is what ranked.
    """
    accessions = {row["protein_id"]}
    extra = row.get("isoforms") or ""
    return accessions | {part.strip() for part in extra.split(";") if part.strip()}


def _sample_random_controls(
    records: list[ProteinRecord],
    validation_cfg: dict[str, Any],
    positives: dict[str, set[str]],
    negatives: dict[str, tuple[str, set[str]]],
) -> set[str]:
    import random

    # Every accession of every control, not just the canonical one: drawing a
    # control's own isoform as a "random" protein would corrupt the random set.
    reserved: set[str] = set()
    for accessions in positives.values():
        reserved |= accessions
    for _, accessions in negatives.values():
        reserved |= accessions

    rng = random.Random(validation_cfg["random_seed"])
    pool = [
        record.protein_id
        for record in records
        if record.survived and record.protein_id not in reserved
    ]
    n = min(validation_cfg["random_control_n"], len(pool))
    return set(rng.sample(pool, n)) if n else set()


# -- symbiont --------------------------------------------------------------


def cmd_symbiont_gate(args: argparse.Namespace) -> int:
    from .symbiont import MAG, gate

    config = io.load_config(args.config)
    paths = _paths(config)
    with open(args.mags) as handle:
        mags = [MAG(**item) for item in json.load(handle)]

    report = gate(mags, config["symbiont"])
    out = paths["results"] / "symbiont_gate.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2))
    print(report["conclusion"])
    # A no-go is the expected result and not an error, so this exits 0 either
    # way; the workflow reads `proceed` from the JSON.
    return 0


# -- argument parsing ------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="bluesea", description=__doc__)
    parser.add_argument("--config", default="config/config.yaml")
    parser.add_argument("--weights", default="config/weights.yaml")
    subparsers = parser.add_subparsers(dest="command", required=True)

    def add(name: str, func, species: bool = False):
        sub = subparsers.add_parser(name)
        sub.set_defaults(func=func)
        if species:
            sub.add_argument("--species", required=True)
        return sub

    add("prepare", cmd_prepare, species=True)
    add("verify-pilot", cmd_verify_pilot)
    check_gate = add("check-gate", cmd_check_gate)
    check_gate.add_argument("--report", required=True)
    prefilter = add("prefilter", cmd_prefilter, species=True)
    prefilter.add_argument(
        "--profile",
        action="store_true",
        help="run every prefilter step on every protein to see independent drop rates",
    )
    add("assemble", cmd_assemble, species=True)
    cascade = add("cascade", cmd_cascade, species=True)
    cascade.add_argument(
        "--profile",
        action="store_true",
        help=(
            "run every remaining step on every protein the predictors saw, to see "
            "independent drop rates; pair with cascade.prefilter.enabled=false to "
            "measure what the exposure filter costs"
        ),
    )
    add("orthology", cmd_orthology)
    literature = add("literature", cmd_literature, species=True)
    literature.add_argument(
        "--live", action="store_true", help="actually call the API (default: dry run)"
    )
    literature.add_argument(
        "--blinded",
        action="store_true",
        help="exclude vaccine terms; required for positive-control benchmarking",
    )
    add("rank", cmd_rank)
    dossiers = add("dossiers", cmd_dossiers)
    dossiers.add_argument("--top", type=int, default=45)
    benchmark = add("benchmark", cmd_benchmark)
    benchmark.add_argument("--report-only", action="store_true")
    symbiont = add("symbiont-gate", cmd_symbiont_gate)
    symbiont.add_argument("--mags", required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
