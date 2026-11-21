# Open scientific and deployment decisions

- [ ] Provision licensed local SignalP/DeepTMHMM and BepiPred/ESM weights,
      document versions, configure required files, and run the real control pilot.
- [ ] Assemble the administered TT-P0 sequence. Leal et al. 2019
      (doi 10.1016/j.fsi.2019.06.018) has been read and supplied both T-cell
      epitopes verbatim, the fusion orientation, the vector and the strain,
      and it confirms the P0 peptide is aa 267-301 of BT077866 for TT-P0 too.
      What remains is the linker and His-tag residues, which neither paper
      states, and the P0 35-mer itself. Only the submitted manuscript was
      readable; its restriction-site assignment is internally inconsistent,
      so the published version is still worth obtaining.
- [x] Map the fragments to the configured annotation. The P0 peptide derived
      from BT077866 matches XP_040578982.1 exactly at 267-301, with the single
      clone-versus-RefSeq difference falling outside it, and both cloning
      primers confirm the boundaries independently. Both P0 components are
      now `verified`.
- [ ] Settle two things my32 still lacks: the P0-to-my32 junction, which no
      paper states, and a single residue difference at position 132 between
      the administered ADD38399 sequence and all four mapped RefSeq entries
      (99.4% identity). Its mapping stays pending until both are explained.
- [ ] Resolve two timing conflicts in the 2025 study, both recorded unresolved:
      the "2, 5, 11 months" labels match no challenge date under either the
      prime or the boost anchor, and Figure 5's 32 dpc contradicts the week
      44-to-51 interval. Arithmetic hypotheses are recorded but not applied.
- [ ] Add independent primary benchmark studies with traceable constructs,
      comparators and endpoints. Multiple constructs in one study are not
      independent replications. Two unassessed leads: Swain et al. 2020
      (doi 10.1371/journal.pone.0239827), and Carpio et al., "Novel chimeric
      antigen reduced adult sea louse Lepeophtheirus salmonis infestations
      under different prime-boost vaccination strategies", cited as submitted
      to Scientific Reports in 2018 and not yet located in published form.
      Leal et al. 2019 does *not* qualify: it measures IgM titers only, with
      no challenge or infestation endpoint, and the registry records that.
- [ ] Decide how to handle experimental components rejected by exposure using
      localization evidence. Do not change the policy solely to recover controls.
- [ ] Supply curated stage TPM/RNAi tables keyed to the configured annotation.
      Current defaults are explicit missing evidence.
- [ ] Establish wet-lab collaborators for localization, antibody function,
      recombinant expression and longitudinal challenge experiments.
- [ ] Set and document live API spending policy before enabling sourced calls.
      Current budget estimates are planning values and enforced per client,
      not an account-wide billing limit.
- [ ] Decide whether the exploratory symbiont arm merits MAG assembly and direct
      transmission/dependence experiments. Genomic signatures alone cannot
      answer those questions.

Implemented decisions include species-specific life stages, actual host
spikes, immutable stage records, native predictor adapters, explicit missing
control failures, gene-representative ranking, and sensitivity flags in dossiers.
See [validation.md](validation.md) for the current readiness policy.
