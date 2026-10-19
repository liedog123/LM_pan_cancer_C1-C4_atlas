# Task usage

`Rscript run.R --list` lists the 35 tasks. `--task ID --plan` validates a task configuration and prints selected scripts, output paths and unresolved placeholders without reading study contents or creating outputs.

Copy the task DCF outside this checkout, replace `SET_ME` input and output paths, then use `--config /path/to/my_task.dcf`. `--steps sNN,sMM` selects steps in configuration order. Matching prepared intermediates may replace skipped producers; other tasks and absent prepared inputs are not auto-generated.

`--check-config` invokes each selected implementation's metadata-only checks. Producer outputs must already exist for downstream checks. Select the first step or provide matching existing intermediates when checking a continuation.

Ordinary execution checks and runs each selected step in a separate Rscript process, stopping at the first failure. Each step has an independent output directory. Logs record STARTED and COMPLETED; STARTED alone does not establish completion. Existing nonempty outputs are protected.

The first DCF record contains Task, Layer and Output-Root plus Input-* shared values. Later records contain Step, Script and the exact fields in [step_index.tsv](step_index.tsv). DCF comments begin with #. `@input/name` uses Input-name; `@root/task/step` resolves under Output-Root; `@sNN/file` resolves under that producer's output_dir. Relative paths resolve from the task config directory. Child configs use absolute paths. `edge_id` is a scalar; other retained fields are paths.

E1 and E3 use separate records and output directories through the same implementation. Definition-construction tasks may be skipped when using the adopted definitions. Independent Acid distributions, CNV modes, Bassez sensitivities and signature correlations are not automatically chained.

The PRISM task now contains only s03 and requires three platform inputs; old task DCFs containing internal-support fields must be replaced with the current template. Shared-axis adjustment now takes only the complete drug results, shared-patient manifest, prediction audit and prediction directory. Programme-feature preparation s05 no longer requires the pilot or definition fields of the old MOFA input route. Survival s00 prepares the original group-selected object without intermediate preview plots; the formal survival steps remain.

The CD8 export task also generates reported quadrant and trajectory source tables. Use its matching prepared metadata and expression inputs; it is not only a cosmetic plotting step. Final manuscript composite assembly is maintained with the manuscript.

The [reproduction routes](reproduction_routes.md) map tasks to the paper and state their starting boundaries. Structural/source checks do not establish numerical equality to frozen study results. Earlier interfaces remain in Git history; public history is not rewritten.

The source_config column in the step index records the original template name for provenance; current executable task configurations are the paths in task_index.tsv.
