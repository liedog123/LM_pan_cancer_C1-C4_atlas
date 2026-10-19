# C1-C4 atlas analysis code

Code for the lactate-associated transcriptional atlas described in the manuscript. Start with the figure or supplementary table below; only run the corresponding tasks and their required producers.

| Paper results | Analysis route |
|---|---|
| Fig. 1 and Supplementary Table 1 | Signature selection, Acid_Score, AG grouping and survival |
| Fig. 2 and Supplementary Table 2 | LM12 construction, C1-C4 labels, clinical and immune comparisons |
| Fig. 3 and Supplementary Table 3 | mRNA, lncRNA, RPPA, adjacent states and subtype adjustment |
| Fig. 4 and Supplementary Table 4 | Somatic mutation, TMB/MSI and copy number |
| Fig. 5 and Supplementary Table 5 | Drug comparisons, PRISM concordance and programme-feature associations |
| Fig. 6, ICB and Supplementary Table 6 | Single-cell states, pseudobulk, reported sensitivities and ICB |
| Fig. 7 and Supplementary Table 7 | CD8 treatment-time models and figure source tables |

See [reproduction routes](docs/reproduction_routes.md) for task IDs, execution boundaries and output descriptions. The [scope register](docs/manuscript_scope.tsv) records why each retained task is needed. The package has 35 task configurations and 86 implementation scripts (87 step records because E1/E3 share an implementation); script count reflects separate analysis steps, not 86 independent workflows.

```sh
Rscript run.R --list
Rscript run.R --task cd8_response --plan
# Copy a task DCF outside this checkout and replace SET_ME paths before execution.
Rscript run.R --task cd8_response --config /path/to/my_task.dcf --steps s02 --check-config
```

Use [task usage](docs/task_usage.md), [task index](docs/task_index.tsv), [step index](docs/step_index.tsv) and [software requirements](docs/software_requirements.md). Every implementation begins with its purpose, task/step, input fields, output boundary and usage. Existing methodological comments are retained.

The repository starts from the analysis inputs specified in each DCF. It is not a raw-data-to-figure pipeline: upstream raw-data/QC preparation, original oncoPredict training/prediction and GLDS correction are not included. The [input guide](docs/input_delivery.md) distinguishes study matrices from definitions and reported results in the seven supplementary workbooks. CD8 deposition is pending ethics approval; use the accession in the manuscript when available.

[Checks](checks/README.md) cover source preservation, syntax, task wiring and metadata-only safety; no study-data numerical reproduction was performed during this update. Run `Rscript checks/run_checks.R /path/to/new/check_report.tsv` to repeat the interface checks without loading study data.

MIT applies to the code. Third-party datasets and gene-set resources retain their own access and licensing terms. Historical commits remain available; the current branch contains the analysis scope mapped above.
