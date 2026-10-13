# Shared definitions: assignment, axis validation and canonical-format export

LM12 reference inputs can be obtained by exporting the complete `Pathway_definitions` and `Gene_membership` worksheets of Supplementary Table 2 to separate TSVs. See the [supplementary-input mapping](../../docs/reference_file_candidates.md#supplementary-material-first) for fields and DCF keys. No duplicate repository dataset is required; the exported files and downstream numerical agreement have not been tested in this documentation task.

This workflow continues state assignment from the formally processed pan-cancer object and the frozen LM12 definition. It consumes existing Acid_Score values and High/Medium/Low acid groups; it does not discover the nine-gene signature, rebuild LM12, recompute acid scores or reconstruct AG thresholds.

It is not a raw-TCGA or FASTQ pipeline, complete atlas reproduction, automatic promotion to the study's canonical freeze or fixed external classifier. Clustering and IG orientation are defined within each selected cancer/cohort. Existing canonical labels remain the formal downstream reference; outputs from these entries have not been numerically validated and do not replace them.

## Workflow order

| Entry | Inputs | Main outputs | Source script in the original module |
|---|---|---|---|
| `01_assign_ig_and_c1c4_from_processed_input.R` | Processed object + frozen LM12 manifest/membership | Assignment labels, LM12 scores and feature metadata | `01_build_ag3ig2_c1c4_subtypes.R` |
| `02_validate_state_axes_from_labels.R` | 01 assignment labels + feature metadata | Axis sample scores, summaries and statistical tables | `02_validate_c1c4_subtype_axes.R` |
| `03_export_canonical_labels.R` | 01 labels/metadata + 02 axis scores | Independent full/minimal labels in the canonical format, dictionary/counts/manifest/checksums | `04_freeze_c1c4_labels.R` |

Source scripts are under `6.0pancancer_subtype_ag3ig2_c1c4/scripts/`. Run 01, then 02, then 03; source 03 plotting is outside this workflow and is not a dependency of export.

## 01: Assignment interface

Copy [the example DCF](../../config/shared_states.paths.example.dcf) to a private configuration and replace every placeholder. Paths resolve relative to the DCF directory, with existing ancestors and symlinks resolved. The 01 DCF must contain exactly one record and four fields:

| Field | Required content |
|---|---|
| `processed_data` | Readable nonempty regular Rdata file containing `pancancer_data_list`. |
| `lm12_manifest` | Readable nonempty regular TSV with the frozen LM12 feature metadata. |
| `lm12_genes` | Readable nonempty regular TSV with the frozen LM12 gene membership. |
| `output_dir` | Independent absent or empty directory, with a writable existing parent. |

```sh
Rscript workflows/shared_definitions/01_assign_ig_and_c1c4_from_processed_input.R --help
Rscript workflows/shared_definitions/01_assign_ig_and_c1c4_from_processed_input.R --config config/shared_states.local.dcf --check-config
Rscript workflows/shared_definitions/01_assign_ig_and_c1c4_from_processed_input.R --config config/shared_states.local.dcf
```

Only `--check-config` stops before input deserialization, analysis dependencies and output creation. It verifies path metadata, not schema or values. Output/input/configuration overlap, root/home output locations, inaccessible inputs and nonempty outputs are rejected; symlinks are checked through resolved paths.

`--test` actually analyzes the first three cancers in input-object order unless `--cancers=BRCA,COAD` is supplied. `--cancers` actually analyzes the requested names intersected with input names; an empty intersection fails. The original filtering logic is preserved, including use of the first `--cancers=` argument if repeated. These flags are not non-analysis checks.

## 01: Input object contract

`pancancer_data_list` is a named list keyed by cancer. Each element supplies:

- `expression_data`: genes in rows (gene identifiers in row names matching the gene table) and samples in columns, on the unlogged TPM scale. The body computes `log2(TPM + 1)`.
- `acid_scores`: a matrix with sample column names; its first row supplies the existing Acid_Score values.
- `clinical_data`: a table with `sample` and `Group`, where `Group` uses exactly `High`, `Medium` and `Low` to map to AG1, AG2 and AG3.

Samples are intersected across these three inputs. At least ten common samples per cancer are required. Duplicate clinical sample rows retain the first occurrence; clinical-table order determines input sample order. Acid_Score NA values fail. These are existing body rules, not newly introduced preprocessing.

The manifest must contain `display_order`, `feature_id`, `short_label`, `pathway`, `functional_group`, `direction_label` and `stable_gene_count`; it is ordered by `display_order`. The long gene table must contain `feature_id` and `gene`. Matching uses gene identifiers exactly as stored; no identifier conversion is added. Gene sets are unique/sorted within feature, then intersected with expression row names; zero coverage fails.

The preserved body scores all LM12 sets with `ssgseaParam(normalize = TRUE)` and `gsva`, then clusters sample score vectors within cancer using Euclidean distance, `hclust(method = "ward.D")` and `cutree(k = 2)`. It assigns IG1 to the cluster with the higher mean of the 11 non-`LM12_Antimicrobial_Peptides` features and IG2 to the other cluster. IG is an inflammatory-remodeling program state, not a direct immune-cell abundance measurement.

The fixed recoding is C1 = AG1_IG1; C2 = AG1_IG2 + AG2_IG2; C3 = AG2_IG1 + AG3_IG1; C4 = AG3_IG2. C1–C4 are joint states, not a single ordered gradient.

## 01: Outputs

The body writes directly under `output_dir`: `ag3ig2_c1c4_sample_labels.csv` (AG, IG, AGIG, C_subtype, scores and interpretation labels), `pancancer_lm12_scores_by_sample.tsv.gz`, `lm12_feature_metadata.tsv`, `lm12_gene_coverage.tsv`, `ig2_cluster_manifest_by_cancer.tsv`, `ig2_cluster_feature_medians_by_cancer.tsv`, overall/long/wide C1–C4 counts, atomic-group counts, sample counts, `ag3ig2_c1c4_subtype_build_objects.Rdata` and `run_summary.txt`. These are assignment outputs, not the later canonical frozen export.

## 02: State-axis validation from assignment labels

Use [shared_axes.paths.example.dcf](../../config/shared_axes.paths.example.dcf), with exactly `assignment_labels`, `feature_metadata` and `output_dir`. Point the first two fields to 01's `ag3ig2_c1c4_sample_labels.csv` and `lm12_feature_metadata.tsv`; use a separate new output directory. Required label columns are `Sample`, `Cancer`, `Acid_Score`, `AG`, `IG`, `AGIG`, `C_subtype`, `C_label` and all manifest feature IDs.

An ordinary run actually computes within-cancer z-scores for Acid_Score, the mean of the 11 non-antimicrobial LM12 features, the all-feature mean and each LM12 feature. A zero/undefined standard deviation produces zero z-scores under the preserved `safe_z` rule. It writes `c1c4_axis_sample_scores.tsv.gz`, counts, within-cancer/equal-weight state centroids, axis/feature summaries and statistical tables.

The fixed six contrasts are C1/C2, C3/C4, C2/C3, C1/C4, C1/C3 and C2/C4. Axis statistics use within-cancer Kruskal–Wallis tests (with epsilon-squared), pairwise Wilcoxon tests with `exact = FALSE`, and signed-rank tests of cancer-median differences. BH correction is applied separately across cancers for each axis's omnibus tests, within each cancer across its six contrasts, and across six contrasts for each axis's cross-cancer summaries. Per-feature pairwise tests are corrected within feature/cancer across six contrasts; cross-cancer feature summaries are corrected within feature across six contrasts. Complete four-state cancer medians are required for the separately written pan-cancer paired table. These are source-code rules; no statistics were run during preparation of this entry.

## 03: Export in the canonical label format

Use [canonical_labels.paths.example.dcf](../../config/canonical_labels.paths.example.dcf), with exactly `assignment_labels`, `axis_scores`, `feature_metadata` and `output_dir`. Assignment labels and feature metadata come from 01; axis scores come from 02. The original module's 03 plotting step is not required.

An ordinary run reads these tables, checks required columns, duplicate `Sample`–`Cancer` keys, complete key correspondence, agreement of AG/IG/AGIG/C_subtype, and the fixed AGIG-to-C_subtype mapping. It joins assignment metadata onto axis scores and writes full/minimal CSVs, Rdata objects, a subtype dictionary, long/wide/overall counts, an input manifest, output checksums and `run_summary.txt`, retaining the source column order.

The body retains the source constant `canonical_version = "AG3IG2_C1C4_v1_20260615"`. That string records the source export convention; it does not establish equivalence of a new export to the study's frozen canonical labels. Runtime `created_at`/freeze_date, recorded paths and resulting checksums can change. A run produces an independent export copy and does not promote it or replace the project's existing canonical labels. Ordinary execution calculates the original input/output MD5 records; `--check-config` and the checks performed for this code preparation do not hash research files.

## Running the three entries

Each configuration resolves paths from its own directory. The examples show 02 consuming 01 outputs and 03 consuming 01/02 outputs, with three independent output directories. Replace every `SET_ME` placeholder consistently. For 02 and 03, the only optional flag is `--check-config`; neither entry accepts `--test` or cancer filtering.

```sh
Rscript workflows/shared_definitions/02_validate_state_axes_from_labels.R --config config/shared_axes.local.dcf --check-config
Rscript workflows/shared_definitions/03_export_canonical_labels.R --config config/canonical_labels.local.dcf --check-config
# Ordinary execution (calculates statistics or exports files):
Rscript workflows/shared_definitions/02_validate_state_axes_from_labels.R --config config/shared_axes.local.dcf
Rscript workflows/shared_definitions/03_export_canonical_labels.R --config config/canonical_labels.local.dcf
```

## Provenance, access and verification

Source relative to the originating project:
`6.0pancancer_subtype_ag3ig2_c1c4/scripts/01_build_ag3ig2_c1c4_subtypes.R`.

Current source-script SHA-256 (observed during this code preparation on 2026-09-30):
`0ea2c683ed0f3587847a6cdcec10fc0f0b1a803c4e78f4c9b7ac92fbf5cf0a7c`.
This is a fingerprint of the currently inspected source text, not a recovered hash from the historical analysis run.

Public continuation-script SHA-256:
`8bec8fb23766440cdff6c1ba0c0ba3a53f754a27645357e1cb49743f06e51f04`.
Only script text was hashed; no research input or frozen result was hashed.

The source module records a completed formal build and canonical freeze on 2026-06-15; its canonical freeze tag is `AG3IG2_C1C4_v1_20260615`. The source 01 script itself has no explicit version tag. This continuation copies text from `message("=== 6.0 AG3 x IG2 -> C1-C4 subtype build ===")` through end of file byte for byte. Only the execution header is replaced with the public interface and safety checks.

Source: `6.0pancancer_subtype_ag3ig2_c1c4/scripts/02_validate_c1c4_subtype_axes.R`. Current source SHA-256: `6d6b7f50a58e360224127b3702838330e7304fbe64172e7d01ca27879aa29f74`; public `02_validate_state_axes_from_labels.R` SHA-256: `e561b397913a1fa6cfa85a2cf715552e4e35853d1f1c8f17b6b822b682d5cee9`. Text from `c_levels <- c("C1", "C2", "C3", "C4")` through EOF is byte-identical (16180 bytes).

Source: `6.0pancancer_subtype_ag3ig2_c1c4/scripts/04_freeze_c1c4_labels.R`. Current source SHA-256: `2fbfc6410767bab9f54bee918475ce1ae872d5bd5415d07a455adeb6a46a4929`; public `03_export_canonical_labels.R` SHA-256: `8a055aa0c8f36d349d88ee61ca05432158eefd8914c25f999d938da5a90f7af6`. Text from `canonical_version <- "AG3IG2_C1C4_v1_20260615"` through EOF is byte-identical (12975 bytes).

These two source fingerprints identify the currently inspected source text, not historical-run hashes. Only their headers were adapted. The fixed export version and runtime metadata are retained from the source.

No research object, patient labels, gene-set files, third-party source, frozen assets or internal meeting material are distributed. Input access, data release and licensing remain subject to author decisions; no license is granted here. Obtain the processed object and frozen definitions through an author-approved arrangement. Generic TCGA downloads alone do not supply these intermediates.

Checks cover R syntax, help/error behavior and metadata-only interfaces using synthetic nonempty files, plus body-byte preservation and public-text scans. No research inputs were loaded, GSVA/clustering or axis statistics were not run, canonical-format files were not exported, and schemas, values, numerical equivalence and end-to-end reproduction remain unverified. Current R 4.5.1 is the selected test environment; original historical dependency versions are unknown. See [dependencies](../../environment/README.md) and [input/access notes](../../docs/inputs_and_provenance.md). Downstream drug-response statistics continue from their separately frozen canonical labels, not automatically from this entry's output.
