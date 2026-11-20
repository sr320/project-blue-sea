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

## The control sets are still placeholders

`resources/controls/*.tsv` use symbolic ids — `LSAL_P0`, `SSAL_ACTB` — not
accessions from any proteome. Nothing matches them, so
`positive_control_report` recovers nothing and `negative_control_report`
returns `passed: False` on an empty list: **the gate fails by construction
until these are mapped.** Map them once `provision` has put the proteomes on
disk, by searching each control's sequence against the relevant proteome.
