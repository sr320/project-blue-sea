# Reference data

What the run needs before `snakemake all` is worth starting, and which of it
the workflow can get by itself.

## Fetched automatically

```bash
snakemake --profile workflow/profiles/slurm --use-apptainer provision
```

| File | Source |
| --- | --- |
| `resources/reference/<genome>/proteins.faa` | NCBI assembly FTP |
| `resources/reference/<genome>/<genome>.dmnd` | `diamond makedb` over the above |

`<genome>` is `lsal`, `crog` and `ssal`. The FTP path is built from the
`accession` and `assembly` pair in `config.yaml`, so the two must match one
release — a mismatched pair gives a 404 rather than a silent miss.

### Why L. salmonis is a GCF accession

No `GCA_016086655.x` version carries an annotation; they are unannotated WGS
assemblies with genomic FASTA only, so none of them has a proteome to screen.
`GCF_016086655.4` is the NCBI RefSeq annotation of the same UVic assembly
line (UVic_Lsal_1.4, 19,574 protein-coding genes).

The alternatives, if the RefSeq gene models turn out to be the wrong call:

| Accession | Assembly | Genes | Annotated by |
| --- | --- | --- | --- |
| `GCF_016086655.4` | UVic_Lsal_1.4 | 19,574 | NCBI RefSeq |
| `GCA_902860245.1` | LSAL_IOA00-2.0 | 15,081 | Stirling |
| `GCA_905330665.1` | LSAL_IOA00_R2 | 15,203 | Stirling |

Switching is a two-line config change plus a re-provision.

## Built by hand

**The container.** Every rule runs inside `containers/bluesea.sif`:

```bash
bash scripts/build_container.sh
```

It is a script, not a rule, because a rule would have to run inside the image
it was building. SignalP 6 and DeepTMHMM are licence-restricted and are not
installable from conda — stage their tarballs beside `workflow/envs/bluesea.def`
before building, or accept an image in which every rule except `signalp` and
`deeptmhmm` will run.

## Supplied by hand

None of these are on NCBI, and no rule fetches them.

| File | What it is | Where it comes from |
| --- | --- | --- |
| `<genome>/stage_tpm.tsv` | gene × life-stage TPM matrix | LiceBase, or the published stage-series RNA-seq |
| `lsal/rnai_phenotypes.tsv` | knockdown phenotypes | LiceBase |
| `<genome>/secretome_confirmed.txt` | accessions with experimental secretome support | published proteomics; optional, scores as a bonus |
| `<genome>/gpi_anchored.txt` | predicted GPI anchors | optional |

Columns are in `src/bluesea/io.py` (`read_expression`, `read_rnai`). Absent
files are handled: expression and RNAi evidence that is missing scores as
unknown and raises a flag, never as zero. `crog` has `rnai: null` because no
RNAi resource has been published for it.

### The gene-ID join is broken as it stands

`cmd_assemble` joins RNAi and expression on `gene_id`, which `_gene_of` takes
from a `[gene=...]` tag in the FASTA header. **NCBI `protein.faa` files do not
carry that tag**, in either proteome:

```
>XP_040563227.1 uncharacterized protein [Lepeophtheirus salmonis]
>QQP31343.1 Hypothetical protein FKW44_024911 [Caligus rogercresseyi]
```

So `gene_id` falls back to the protein accession for every protein. If the
expression matrix and the RNAi table are keyed on anything else — LiceBase
gene ids, locus tags, `LOC` numbers — nothing joins. Every protein then scores
`untested` on essentiality, which carries the heaviest weight in the rubric,
and the run completes looking healthy.

`assemble` now prints a join rate per table and warns loudly at zero, so this
cannot pass unnoticed. Resolving it is a choice about where gene identity
comes from:

| Option | Gets you | Costs |
| --- | --- | --- |
| Key the expression and RNAi tables on protein accessions | no code change | remapping whatever LiceBase gives you |
| Parse the locus tag out of the description (`FKW44_024911`) | works for `crog` | RefSeq `lsal` descriptions have no locus tag |
| Fetch `translated_cds.faa` instead | real `[gene=...]` and `[db_xref=GeneID:...]` tags | ids are `lcl\|NC_052131.2_prot_XP_040583074.1_1`, needing their own parser |
| Build a protein→gene map from the GFF | authoritative, handles isoforms | another download and a parser |

The GFF route also settles a second problem these files raise: RefSeq
`protein.faa` holds one entry per **isoform**, so several near-identical
proteins of one gene are scored independently and can crowd the top tier.
`kuz` above appears as both `XP_040583096.1` and `XP_040583086.1`.

## Not implemented

`scratch/symbiont/mags.json`, the input to `symbiont_gate`, is the output of a
MAG assembly track that does not exist yet. The tools are pinned in
`workflow/envs/environment.yml`; the rules are not written. `snakemake all`
does not depend on it.

## The control sets

Mapped to `GCF_016086655.4` accessions by alignment, not by annotation text:

```bash
pip install -e '.[curation]'
python scripts/map_controls.py resources/reference/lsal/proteins.faa
```

Annotation text does not work here. RefSeq names P0 "large ribosomal subunit
protein uL10" and L7 "uL30", so neither is findable under the name the
literature uses, while two unrelated proteins are both named "peroxiredoxin 2"
and only alignment says which is the published antigen. The script fetches a
UniProt reference per control, aligns it, and prints identity and coverage so
every assignment can be re-checked.

A control is a **gene**, but the ranking holds proteins and RefSeq carries one
entry per isoform. The `isoforms` column lists the siblings that represent the
same control, semicolon separated, and the benchmark counts a control
recovered when any of them reaches the top tier — reporting which one did, in
`matched_accession`. Without that, a correct result is recorded as a miss
whenever a sibling is what ranked, and with only five active positives one
false miss moves the gate. Negative controls are held to the opposite
standard: every isoform present must drop at the same expected step, because
paralogs rejected by different rules is itself worth seeing.

Five of eight positives and all four louse negatives resolved at ≥89% identity
over full-length coverage. Three positives did not and are commented out
rather than guessed, so the gate measures only what actually resolves:

- **bm86-like** — best alignment to the canonical Bm86 (`P20736`) is 28.8%
  identity over 79% coverage, to "protein crumbs". That is twilight-zone
  similarity between two EGF-repeat proteins, not orthology. Whether
  *L. salmonis* has a Bm86 ortholog at all needs a phylogenetic analysis, not
  a best-hit search.
- **hpx1** — no "heme peroxidase 1" exists in this proteome or in UniProt for
  this species, only hpx2.
- **kcc-transporter** — matches 84 transporters by name, no UniProt reference.

### The host-spike negatives cannot run

The four salmon proteins were never scorable. `cmd_benchmark` builds its
record set from `config["species"]` — lsal and crog — and salmon is
`config["host"]`, with no `cascaded.json`. A salmon accession matches no
record, contributes no row, and `negative_control_report` returns
`passed: False` on an empty list.

Fixing it needs a spike-in — the named salmon proteins appended to the louse
proteome for a benchmark run, which `tests/test_cli_end_to_end.py` does
synthetically — so it is a workflow change, not an id change. Until then the
cross-reactivity filter is exercised by the cascade but unverified by the
negative controls.
