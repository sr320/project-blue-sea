# Open decisions

Two of these change the build; the rest change the deliverable or the budget.

- [ ] **Standalone project, or preliminary-data arm of a specific proposal?**
      Changes the deliverable and the timeline more than anything else here.
- [x] **Both species.** *L. salmonis* and *C. rogercresseyi* run through the
      same cascade and their shortlists are intersected. Implemented:
      `config.yaml:species`, `orthology.intersect_shortlists`.
- [x] **RefSeq gene models for *L. salmonis*.** The config previously named
      `GCA_016086655.3`, which carries no annotation — no `GCA_016086655.x`
      version does, so none of them has a proteome to screen. Now
      `GCF_016086655.4`, the RefSeq annotation of the same UVic assembly line,
      19,574 genes. The Stirling annotations are the alternative and are
      closer to the 14,014 figure this README used to quote; switching is a
      two-line config change. See `docs/reference-data.md`.
- [ ] **Who runs wet-lab validation, and is there a collaborator for challenge
      trials?** Steps 1–3 of the bench sequence are plausibly in-house; step 4
      is not.
- [x] **Control sets mapped by alignment.** Five of eight positives and all
      four louse negatives now carry `GCF_016086655.4` accessions, established
      with `scripts/map_controls.py` and recorded with identity and coverage.
      Review the assignments before the first gated run.
- [ ] **Three positive controls remain unmapped**, and each needs a judgment
      rather than a search: whether *L. salmonis* has a Bm86 ortholog at all
      (the best hit is 28.8% to an unrelated EGF-repeat protein), what the
      source paper's `hpx1` corresponds to when only `hpx2` exists, and which
      of 84 transporters the `kcc-transporter` antigen is. They are commented
      out meanwhile, so the gate is scored over five.
- [ ] **Host-spike negative controls need a spike-in mechanism.** Salmon
      accessions can never match a record, because the benchmark reads only
      `config["species"]`. Until the named salmon proteins are appended to the
      louse proteome for benchmark runs, the cross-reactivity filter is
      unverified by the negative controls.
- [ ] **Total API budget ceiling for the first full run.** Currently
      `config.yaml:llm.budget_usd` sums to $400 with the symbiont arm capped
      at $60 (15%). These are planning placeholders, not a quote.
- [ ] **Public GitHub repo from the first commit, or after the first validated
      run?** Affects nothing in the code; affects when the negative results
      become citable.

## Build-affecting notes

**Verify before building further**: Hyak compute nodes do not reliably have
outbound internet access. The architecture already assumes they do not — the
API stage is separate and reads tool output from disk — so if that assumption
turns out to be wrong, nothing breaks. If it is right, nothing needs
rearranging.

**Secretome data is thin.** There is no comprehensive, stage-resolved
*L. salmonis* secretome. The pipeline predicts it computationally and treats
experimental confirmation as a scoring bonus (`exposure.class_scores.
secretome_confirmed`), never a requirement. Drop confirmed accessions into
`scratch/<species>/secretome_confirmed.txt` as they become available.

**RNAi coverage is partial.** `weights.yaml:essentiality.phenotype_scores`
maps `untested` to 0.4, above `no_phenotype` at 0.1, and flags the record.
Do not "fix" this by setting it to zero.

## Bench sequence, cheapest discriminator first

1. **RNAi knockdown in copepodids or pre-adults.** Directly tests the
   essentiality axis, which carries the most weight and the most uncertainty.
   A candidate whose knockdown produces no phenotype should probably not
   advance regardless of score.
2. **Recombinant expression and antibody generation.** Immunohistochemistry
   here is the check on the exposure axis.
3. **In vitro or ex vivo antibody challenge**, if a louse culture system
   allows it, before committing to a fish trial.
4. **Immunization and challenge trial.** The only step that tests the actual
   claim. Design to twelve months, not three — that is where every reported
   effect has decayed.
