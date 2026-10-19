# C1-C4 atlas analysis code

Explore lactate-associated transcription and inflammatory remodelling through the paper's analysis routes.

One dispatcher exposes 36 task configurations: 23 core analyses, 3 definition-construction tasks, 9 supporting analyses and 1 CD8 figure/table task. Ninety retained implementations are internal steps; their scientific code has not been refactored.

This code-only release provides the task interface below. Interface checks are recorded in [checks](checks/README.md); study-data numerical reproduction has not been performed.

| Paper result | Tasks |
|---|---|
| Fig. 1 | acid_ag_groups, acid_ph_correlations, acid_distributions, ag_survival |
| Fig. 2 | c1c4_state_assignment, external_bulk_enrichment, c1c4_clinical, immune_ecology |
| Fig. 3 | mrna_programmes, lncrna_rppa, boundary_gsea |
| Fig. 4 | somatic_mutations, tmb_msi, cnv_burden, cnv_regions |
| Fig. 5 | predicted_drug_response, programme_feature_associations |
| Fig. 6 | single_cell_states, single_cell_pseudobulk, melanoma_pseudobulk, icb_three_cohort |
| Fig. 7 | cd8_qc, cd8_response, cd8_figures |

Supporting tasks reproduce reported supplementary results and remain part of full-paper reproduction. Definition construction can be skipped when starting downstream from adopted frozen definitions.

```sh
Rscript run.R --list
Rscript run.R --task cd8_response --plan
Rscript run.R --task cd8_response --config /path/to/cd8_response.dcf --steps s02 --check-config
Rscript run.R --task cd8_response --config /path/to/cd8_response.dcf
```

Obtain matching inputs and edit a copied task configuration first. Example placeholders are not delivered data. Use a new output root outside this package. There is no repository-wide run-all command.

See [usage](docs/task_usage.md), [input delivery](docs/input_delivery.md), [task index](docs/task_index.tsv) and [step mapping](docs/step_index.tsv).
The existing [result mapping](https://github.com/liedog123/LM_pan_cancer_C1-C4_atlas/tree/7fc95d5e205c7dfced0ffdfa4b17945bf62d13e5/docs/result_output_mapping.md) and [environment evidence](https://github.com/liedog123/LM_pan_cancer_C1-C4_atlas/tree/7fc95d5e205c7dfced0ffdfa4b17945bf62d13e5/environment/README.md) remain the provenance sources.

MIT applies to the code. Supplementary tables provide study definitions, memberships and reported results; see the [supplementary content map](docs/supplementary_inputs.tsv). Obtain study data through the sources specified in the manuscript and set the required paths in a copied task configuration.

Last updated: 2026-10-19
