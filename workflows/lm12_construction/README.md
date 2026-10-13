# LM12 construction from formal GSEA results

Downstream users can start from the adopted LM12 definitions using the [complete Table 2 worksheet exports](../../docs/reference_file_candidates.md#supplementary-material-first). Reconstructing selection with 02 still requires the complete curated candidate table and family dictionary; the final 12 pathways cannot replace the preselection universe. These acquisition decisions remain open, and no worksheet export or analysis was run in this documentation task.

This workflow continues from the study's formal per-cancer GSEA CSV results. It computes pathway recurrence votes and stable leading-edge members, then finalizes the study-specific LM12 definitions using separately supplied curated candidate/family reference tables. It does not generate the upstream differential-expression rankings or GSEA, discover the acid signature, recompute AG, or reproduce the full atlas.

The second entry retains an explicit assertion against 12 named pathways encoded in the source. It is a finalization of the study-specific LM12 result, not a general engine that discovers 12 gene sets from arbitrary GSEA inputs. New outputs are independent copies and do not replace the project's frozen LM12 definitions or canonical labels.

| Order | Entry | Inputs | Main outputs |
|---|---|---|---|
| 01 | `01_build_votes_and_stable_genes.R` | Formal per-cancer GSEA CSV directory | All/selected pathway votes, stable-member summary and long table, per-set gene files |
| 02 | `02_finalize_lm12_gene_sets.R` | Three 01 tables + curated candidate terms + family dictionary | Candidate audit tables, LM12 manifest/long gene table, functional annotations and per-feature gene files |

## Configuration and execution

[01 example](../../config/lm12_votes.paths.example.dcf) contains exactly `gsea_results_dir` and `output_dir`. [02 example](../../config/lm12_finalize.paths.example.dcf) contains exactly `all_votes`, `gene_summary`, `gene_long`, `candidate_terms`, `family_dictionary` and `output_dir`. Both DCFs must have one record; replace all `SET_ME` placeholders. Relative paths resolve from the configuration directory, including existing ancestors and symlinks.

Each output must be absent or empty and independent of inputs/configuration. Root/home output locations, overlap, inaccessible files and nonempty outputs are rejected. For 01, the original `_GSEA_full_results\\.csv$` filename rule is preserved: at least one matching readable nonempty regular file is required, including resolved symlink targets. No fixed 22-cancer count or cancer-name filter is added.

```sh
Rscript workflows/lm12_construction/01_build_votes_and_stable_genes.R --help
Rscript workflows/lm12_construction/01_build_votes_and_stable_genes.R --config config/lm12_votes.local.dcf --check-config
Rscript workflows/lm12_construction/02_finalize_lm12_gene_sets.R --config config/lm12_finalize.local.dcf --check-config
# Ordinary execution calculates votes/member sets and writes outputs:
Rscript workflows/lm12_construction/01_build_votes_and_stable_genes.R --config config/lm12_votes.local.dcf
Rscript workflows/lm12_construction/02_finalize_lm12_gene_sets.R --config config/lm12_finalize.local.dcf
```

Only `--check-config` returns before package loading, input-content reads and output creation. It checks path/file metadata, not input schemas or scientific content. Neither entry accepts `--test` or cancer-selection flags.

## 01 input and preserved rules

The directory contains CSVs named `<CANCER>_GSEA_full_results.csv`. Cancer labels are derived by removing the filename suffix. Each CSV is read by base `read.csv(check.names = FALSE)` and requires `pathway`, `padj`, `NES` and `leadingEdge`. The latter stores comma-separated gene identifiers; the source parser splits on commas, trims whitespace and drops empty strings. Supply the formal results on the original pathway names and gene identifiers; no identifier conversion is added.

The preserved body uses `padj < 0.05` and separates positive NES (`High_Acid_Positive`) from negative NES (`Low_Acid_Negative`). Votes count distinct cancers per direction/pathway. The denominator is the number of distinct filename-derived cancers actually present in the loaded GSEA rows. Original-like classes are Core at vote frequency >= 0.60, Minor at >= 0.40 and Sporadic otherwise; selected pathways have frequency >= 0.40. Within each selected direction/pathway, the first significant row per cancer is retained. Stable genes use the original leading-edge occurrence counting and a threshold of `ceiling(vote_count * 0.70)`; the parsing/counting behavior is unchanged.

Outputs retain the source names: `all_pathway_votes_original_like.tsv`, `selected_core_minor_pathways_original_like.tsv`, `original_like_gene_set_summary.tsv`, `original_like_gene_sets_long.tsv`, per-set `gene_sets/*.txt` and `run_summary.txt`.

## 02 input and preserved rules

The five file fields point to readable nonempty regular TSVs, read with `readr::read_tsv`:

| DCF field | Source-required columns |
|---|---|
| `all_votes` | `direction`, `pathway`, `vote_count`, `vote_freq`, `mean_NES`, `median_NES` |
| `gene_summary` | `direction`, `pathway`, `stable_gene_count`, `leading_edge_threshold_count` |
| `gene_long` | `pathway`, `gene`, `leading_edge_count`, `vote_count`, `leading_edge_threshold_count` |
| `candidate_terms` | `analysis_tier`, `candidate_set`, `scope_layer`, `family_id`, `local_term`, `match_type`, `note` |
| `family_dictionary` | `scope_layer`, `family_id`, `fig1_or_project_basis`, `note` |

The first three tables come from 01. The last two are separately curated study references; this workflow does not create or distribute them. It selects `analysis_tier == "primary"`, matches pathway names to candidate `local_term`, joins gene summaries by direction/pathway and uses the original family annotation keys. Original-like thresholds 0.60/0.40 and strict thresholds 0.90 (Core)/0.80 (recurrent selection) are unchanged. Strict selected pathways must match the source's 12 expected pathway names, or execution stops. Fixed feature IDs, display order and functional labels are retained. The original pathway-only member joins are also unchanged.

Outputs retain `predefined_primary_candidate_all_votes.tsv`, `selected_primary_candidates_original_like_60_40.tsv`, `selected_primary_candidates_strict_80.tsv`, `lm_imm_12_manifest.tsv`, `lm_imm_12_pathway_vote_summary.tsv`, `lm_imm_12_gene_sets_long.tsv`, `lm_imm_12_functional_annotation.tsv`, `gene_sets/*.txt` and `run_summary.txt`.

The manifest and long gene table use the formats consumed by [shared state assignment](../shared_definitions/README.md). Format compatibility has been traced from source; no new LM12 output has been generated or validated against the frozen definitions.

## Access, provenance and verification

No formal GSEA results, patient data, curated candidate/family tables, gene-set materials, third-party source or frozen assets are included. Obtain inputs through an author-approved access/release arrangement. Public intermediate access, curated-reference distribution and licensing remain undecided; no license is granted here. Generic TCGA or pathway downloads alone do not supply these exact intermediates.

Source: `5.0LM_Imm_12geneset/scripts/01_build_original_like_votes_and_stable_genes.R`. Current source SHA-256: `2a184d7958d064bb00fd3f43028972a1b4f6d00613b4d47cd9c3c972777ee545`. Public `01_build_votes_and_stable_genes.R` SHA-256: `044082924b3b03bf32e0003ad1a9d443e3b5c59b278f8fadb4fdd89a5744eba5`. Text from `padj_cutoff <- 0.05` through EOF is byte-identical (5155 bytes).

Source: `5.0LM_Imm_12geneset/scripts/02_finalize_lm_imm_12_gene_sets.R`. Current source SHA-256: `45018d48a8575e2ba67f962f4f922ef52e2444d0b87729bd17115158cca7078c`. Public `02_finalize_lm12_gene_sets.R` SHA-256: `996d14a1b957efe28ef56f2384394e16f6499cf609e79e216e46075f3be83bf2`. Text from `required_files <- c(` through EOF is byte-identical (8067 bytes).

Only execution headers are adapted; original selection, counting, joins, thresholds, expected-pathway assertion and output names are retained. These fingerprints identify the currently inspected source text, not historical-run hashes. Only script text was hashed.

Checks cover syntax, help/error behavior, synthetic metadata-only configurations, directory-file/symlink output protection, body-byte preservation, example tracking and public-text scans. No real GSEA/reference contents were read and neither analysis body was run. Real-input schema, numerical equivalence, expected output contents, dependency completeness and end-to-end reproduction remain unverified. R 4.5.1 is the selected static/interface test environment; initial historical R/package versions remain unknown. See [dependencies](../../environment/README.md) and [inputs/provenance](../../docs/inputs_and_provenance.md).
