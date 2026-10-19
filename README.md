# Pan-cancer C1–C4 lactate-associated state atlas

Code supporting the C1–C4 atlas manuscript: state definitions, clinical and molecular comparisons, drug-response associations, single-cell/ICB analyses and CD8 perturbation results. The package starts from matching prepared inputs; data access follows the manuscript. Study data, participant-level tables and complete result packages are not bundled here.

## Quick start

```sh
Rscript run.R --list
Rscript run.R --module clinical --plan
Rscript run.R --task ag_survival --plan
```

Copy the appropriate `config/*.dcf` outside the checkout, complete the selected inputs and new output paths, then use `--check-config` before running. A plan does not verify data contents. See [reproduction.md](docs/reproduction.md) for input acquisition boundaries, prepared handoffs, result-to-task mapping and detailed usage.

## Structure

- `run.R`: one entry point; stages execute in independent R processes.
- `R/analysis/`: 18 analysis modules containing named stage functions.
- `R/common/`: five files for configuration, I/O, runtime, statistics and plotting.
- `config/`: eight path templates and one authoritative stage registry.
- `docs/`: reproduction guide, stage dependencies and manuscript/supplement content maps.
- `checks/run_checks.R`: syntax, interface and synthetic dispatcher checks.

| Analysis module | Configuration | Task aliases |
|---|---|---|
| `signature_acid` | `states.dcf` | `acid_signature_selection`, `signature_correlations`, `acid_ag_groups`, `acid_ph_correlations`, `acid_distributions` |
| `immune_states` | `states.dcf` | `tcga_immune_discovery`, `lm12_definition`, `c1c4_state_assignment`, `external_bulk_enrichment` |
| `clinical_survival` | `clinical.dcf` | `ag_survival`, `c1c4_clinical` |
| `immune_ecology` | `clinical.dcf` | `immune_ecology` |
| `mrna_programmes` | `molecular.dcf` | `mrna_programmes` |
| `lncrna_rppa` | `molecular.dcf` | `lncrna_rppa` |
| `boundary_concordance` | `molecular.dcf` | `boundary_gsea`, `mrna_rppa_concordance` |
| `subtype_context` | `molecular.dcf` | `subtype_context` |
| `genomics` | `genomics.dcf` | `somatic_mutations`, `tmb_msi`, `cnv_burden`, `cnv_regions` |
| `drug_response` | `drug.dcf` | `predicted_drug_response`, `drug_crossplatform`, `drug_shared_axis` |
| `programme_associations` | `drug.dcf` | `programme_feature_associations` |
| `single_cell_states` | `single_cell.dcf` | `single_cell_states` |
| `pan_pseudobulk` | `single_cell.dcf` | `single_cell_pseudobulk` |
| `single_cell_sensitivity` | `single_cell.dcf` | `bassez_sensitivity` |
| `melanoma_pseudobulk` | `single_cell.dcf` | `melanoma_pseudobulk` |
| `icb` | `icb.dcf` | `icb_three_cohort`, `liu_icb`, `gide_c1` |
| `cd8_response` | `cd8.dcf` | `cd8_qc`, `cd8_response` |
| `cd8_figures` | `cd8.dcf` | `cd8_figures` |

The 35 task aliases and 87 execution records are retained within these modules. E1 and E3 share an implementation but remain distinct executions. Original output filenames and handoff directories are preserved. The consolidation removes repeated script/CLI/configuration files while retaining the sample-selection, covariate and testing rules needed to reproduce the reported comparisons.

Package dependencies are listed in [software_dependencies.tsv](docs/software_dependencies.tsv). Coverage and retention reasons are in [manuscript_scope.tsv](docs/manuscript_scope.tsv), and the seven submitted supplement mappings are in [supplementary_inputs.tsv](docs/supplementary_inputs.tsv). Final manuscript panel assembly and original sequencing preprocessing are outside this prepared-input code package.

## Validation

```sh
Rscript checks/run_checks.R /path/outside/checkout/new_checks.tsv
```

Validation exercises code structure, filesystem metadata guards and synthetic process behavior. It does not run study analyses or claim end-to-end numerical reproduction. Pre-refactor implementations remain available in Git history at `41670141faa1902e623bebcb725ec79e33b0b563`; `config/stages.tsv` provides the source mapping.

Licensed under the included [LICENSE](LICENSE).
