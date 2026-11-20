# Reference data

`snakemake -s workflow/Snakefile provision` downloads each configured NCBI
proteome and the louse GFF annotations, then builds DIAMOND databases. Assembly
accession/name pairs in `config.yaml` determine the exact FTP URLs. Treat the
reference directory as read-only after provisioning.

The GFF links each protein accession to its gene. RefSeq `GeneID:` identifiers,
GenBank locus tags, gene symbols and protein accessions are tried in that order
when joining curated expression/RNAi tables; join counts are reported. A
missing GFF causes a warning and protein-accession fallback. Snakemake requires
the configured GFF before the normal workflow runs.

Expression and RNAi sources are optional and default to `null`. If a path is
configured, that file is required and declared as a workflow dependency. A
missing configured file is an error, rather than silently becoming unknown.

| Evidence | Schema | Source |
| --- | --- | --- |
| Stage expression | `gene_id` plus species-valid stage columns; nonnegative TPM | Curated stage-series RNA-seq / LiceBase |
| RNAi | `gene_id`, `phenotype`, `source`, `note` | Curated primary phenotypes / LiceBase |
| Secretome support | One experimentally supported accession per line | Primary proteomics |
| GPI anchors | One accession per line | Recorded prediction/experimental source |

Use the stage lists in each species's config block. Unknown stage columns
fail validation. Missing columns/`NA` are unassayed rather than zero TPM.
Insufficient assayed stages cannot support an expression exclusion. Partial
stage coverage remains visible in scores and flags.

Optional accession lists are configured through each species's `secretome` and
`gpi` paths relative to the reference directory. Configured files are required,
declared workflow inputs and included in the pilot signature; retain the primary
source and prediction method with the reference curation.

The [tool installation guide](tools.md) covers licensed local predictors,
model weights, container staging and native-output adapters. Merely building
the baseline container does not install those predictors.

The [validation guide](validation.md) covers experimental constructs, inactive
legacy controls and host spikes. Control mapping is annotation-release-specific;
a new annotation requires rechecking fragments and isoforms.

`scratch/.../symbiont/mags.json` is user-supplied output from a MAG assembly
track that is not implemented. Its quality fields use percentages. The MAG
signature is exploratory and does not establish obligate symbiosis.
