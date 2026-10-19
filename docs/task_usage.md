# Task usage

`--list` lists the tasks. `--task ID --plan` validates structure and prints selected scripts, output paths and placeholders. It creates no outputs and does not read study contents.

`--steps sNN,sMM` selects steps in configuration order. Matching prepared intermediates may replace skipped producers. Other tasks and absent prepared inputs are not auto-generated.

`--check-config` invokes each selected implementation's metadata-only checks. Producer outputs must already exist. A downstream metadata failure before its producer runs is expected; select the first step or supply matching existing intermediates.

Ordinary execution checks and runs each selected step in a separate Rscript process, stopping at the first failure. Original guards and analysis bodies remain in force. Each step has a distinct output directory. Logs record STARTED and COMPLETED; STARTED does not mean completion.

The first DCF record contains Task, Layer and Output-Root plus Input-* shared values. Subsequent records contain Step, Script and exact retained fields. Task configs support comments beginning with #; generated child DCFs contain only valid fields.

`@input/name` uses Input-name. `@root/task/step` resolves under Output-Root. `@sNN/file` resolves under that step's output_dir. Replace a producer reference with a matching existing path when continuing from a prepared intermediate. Relative paths resolve from the task config directory; child configs use absolute paths.

`edge_id` is a scalar; other retained fields are paths. E1 and E3 use separate records/output roots through the same retained implementation. E2 remains independent. Statistical families and supplied P/FDR are preserved.

Independent Acid distributions, CNV burden modes, Bassez sensitivities and signature correlations are not chained. The lncRNA/RPPA task retains separate branches after shared preparation. The retained CellChat branch keeps rep1 and its distinct rep2-10 seeds; it is outside the current manuscript and submitted supplements.

CD8 figures use the same Output-Root convention as cd8_qc/cd8_response, or explicit prepared paths. Hypoxia trajectory and audited sample metadata remain separately prepared inputs. Descriptive table generation stays with figures and is required for complete paper reproduction.

Structural/source inspection does not establish real-data compatibility or equality to frozen results.

## Existing users

The previous numbered workflows and their 90 DCF examples remain available at commit 7fc95d5. This release uses one dispatcher and 36 task DCFs; legacy DCFs are not task configs. Copy the matching task template and transfer the exact input paths. Retained implementations are mapped in step_index.tsv. Use independent new output directories; matching frozen intermediates can replace earlier producers.

## Manuscript scope

Internal drug-model performance is not a reported result in the current manuscript. The drug_crossplatform task retains s02 because s03 reads its preparation outputs. To run s03 from existing preparations, supply matching paths for internal_support_summary and drug_identity in a copied task configuration and select --steps s03. Keep this dependency when running the original route; the submitted PRISM comparisons are mapped to Supplementary Table 5.
