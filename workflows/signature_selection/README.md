# Acid-signature selection and final-nine correlations

The [reference-input mapping](../../docs/reference_file_candidates.md#supplementary-material-first) lists the nine genes already in Results/public code and the minimal `Gene` TSV accepted by 02/03. No duplicate membership upload is required. This does not supply 01's full 142-candidate input or the study expression/score objects; their acquisition routes remain unresolved. The minimal-format recipe has not been executed.

This module contains three study-specific calculation entrypoints: 01 selects the signature from 142 candidates; 02 and 03 independently calculate S1A/B correlations for the formal nine genes. All retain the full within-cancer TPM/score sample intersection. These workflows do not regenerate DESeq2 tumour-normal statistics or tumour-program GSVA.

## Candidate selection (entrypoint 01)

The sections below through Source adaptation and verification describe 01 only. Its all_sample route includes tumour and normal samples in the exact TPM/score intersection. It does not assert that nine genes will be selected; 02/03 require the formal nine-gene TSV. The historical tumour-score loop is not copied as a functioning upstream workflow.

The separately obtained 142-member candidate table, saved tumour-normal DEG and tumour-score objects, and per-project TPM files are required. Existing 142 candidates are reused; no reference gene materials, research objects or patient labels are distributed. Access, release and licensing remain undecided.

## Interface and input schema

[signature_selection.paths.example.dcf](../../config/signature_selection.paths.example.dcf) has exactly candidate_genes, deg_results, tumor_scores, tpm_root and output_dir. Replace placeholders; paths resolve relative to its directory with existing ancestors/symlinks resolved. Output must be absent or empty and independent of inputs/configuration.

- candidate_genes: TSV with gene column; the source asserts 142 rows, 142 unique genes and no missing gene values, preserving candidate order.
- deg_results: Rdata containing DEGs_list, named by TCGA project. Per-project tables require row, baseMean, log2FoldChange and padj; duplicate gene rows are rejected.
- tumor_scores: Rdata containing tumor_gsva_list, with matching project keys; each score item is a one-row matrix with full sample-ID column names and no duplicate names.
- tpm_root: PROJECT/PROJECT_mrna_expr_tpm.rdata files containing mrna_expr_tpm. Matrices need unique gene/sample dimnames. Raw stored TPM scale is retained, with no added log transformation.
- output_dir: separate output location.

The ordinary run preserves the source's 21 matched-project assertion and equal DEG/score project sets. --check-config only checks file/path metadata, not these keys/counts, actual eligible-project TPM coverage or schemas. Recursive TPM metadata discovery protects resolved file targets; ordinary loading still uses the original project-specific filenames.

Use Rscript workflows/signature_selection/01_select_signature_all_sample.R --help or --config config/my_signature.dcf [--check-config]. Metadata checking exits before object/table contents, any calculation or output creation. An ordinary invocation actually recalculates correlations and selection; no --test or sample-scope option is added.

## Preserved selection and outputs

Stage 1: saved DEG statistics, finite/evaluable status, baseMean > 100, log2FoldChange >= 0.5 and padj < 0.05. Stage 2: per-cancer Pearson on the all-sample TPM/score intersection, source finite/variance checks and BH across evaluable tested candidates; r > 0.25 and FDR < 0.05. Stage 3: within-cancer Spearman pairs, complete evaluable non-self pair coverage, positive ratio >= 1 and positive-correlation geometric mean > 0.25. Stage 4: count cancers passing Stage 3, retain Count >= 6. No threshold changes, pooled-cancer correlations, pairwise FDR or assertion of nine selected genes is introduced.

Outputs retain stage0_input_inventory.tsv, stage0_sample_and_scale_inventory.tsv, Stage 1–3 status/cancer-summary tables, stage3_pairwise_correlations_long.tsv, stage4_cross_cancer_counts.tsv, gene_by_cancer_retention_long.tsv and reconstructed_flow_counts.tsv. Only all_sample rows are calculated; no dual-scope/frozen comparison is generated. Applicable source QA checks retain candidate/project assertions, complete Stage 1 grid, normal samples in each all-sample cohort, nonempty intersections, retention-derived all-sample counts and nonempty outputs. The qa subdirectory contains machine checks, session information and an English run report; reconstruction_manifest.tsv retains ordinary-run file sizes/MD5. Input/output MD5 operations occur only in ordinary execution, not during this packaging check.

An explicit format-only bridge builds project-ordered Stage 3 passing-gene lists, then reuses the historical table(unlist(...)) -> Gene/Count -> Count >= 6 sequence. It preserves table-generated row order rather than the candidate-order Stage 4 display. A consistency check matches selected genes/counts to Stage 4; it does not force a nine-gene result. Outputs final_selected_genes.rdata (named final_selected_genes data.frame with Gene/Count) and final_selected_genes.tsv are compatible in schema with [acid scoring](../acid_ag_preprocessing/README.md). Point that workflow's final_genes field to the new Rdata only as an independent continuation copy; equality to the frozen nine-gene object has not been numerically verified here and frozen inputs are not replaced.

## Source adaptation and verification

Primary source: 1.0gene_select/revisions/yanmin_review_r01/v0.1_20260831/scripts/04_reconstruct_142_gene_selection_dual_scope_v0.2_20260901.R.
Current code SHA-256: f46ee88860f5755f8632c70115ee0055ea4ad355f59c9b35f863b3d175cb78c8 (not a historical execution hash).
The export-format source is the final counting/export block of 1.0gene_select/LM_gene_select.R.

Preserved helpers include source safe_pearson/safe_spearman, isolation loader, tabular helpers and logging. Retained Stage 0–4/summary/output blocks are adapted only by private-path/DCF mapping, single-mode iteration replacing the dual-mode lists, removal of frozen/formal input/comparison and tumour-only QA parts, the explicit compatible export bridge/manifest extension, English public reporting and trailing-space normalization. Statistical function bodies, candidate/project assertions, input/sample ordering and selection calculations are unchanged. The whole script is not byte-identical to either original. Script 05 plot-ready formatting, its fixed nine-gene/flow checks and display ordering are not included.

The R01 controlled reconstruction contract records R4.5.1 for that reconstruction. Earlier initial DESeq2/tumour-GSVA generation versions remain unrecovered. Current public syntax/interface checks also use author-selected R4.5.1, separately from that existing record. The entry uses base/recommended R namespaces (stats, utils, tools); it needs no new DESeq2 or GSVA execution. Checks are limited to static retained-block comparison, parse, synthetic metadata interfaces and a tiny format-bridge example. No research contents were loaded/hashed or selection statistics run. Real input schemas, coverage, numerical equivalence and end-to-end reproduction remain unverified. This is filtering transparency, not pH/lactate measurement, mechanistic evidence or independent validation.

## Final-nine supplementary correlations (S1A/B)

[02_correlate_final_nine_tumor_program.R](02_correlate_final_nine_tumor_program.R) and [03_correlate_final_nine_gene_pairs.R](03_correlate_final_nine_gene_pairs.R) are parallel, independent calculation entrypoints from existing study intermediates. They do not select the 142 candidates or generate plots. Both require the formal Fig1A_signature_genes.tsv, stored tumour-program scores and per-project TPM. Expression/score input acquisition and permission to distribute study intermediates remain unresolved; the example configurations are path templates, not delivered inputs.

The configurations [final_nine_tumor_program.example.dcf](../../config/final_nine_tumor_program.example.dcf) and [final_nine_gene_pairs.example.dcf](../../config/final_nine_gene_pairs.example.dcf) each have exactly four fields:

- final_genes: the formal TSV, using Gene when present and otherwise the first column; unique genes retain source order and must number nine.
- tumor_scores: Rdata containing tumor_gsva_list, named by 21 projects; projects use sort(names(scores)). Each project needs full, unique sample column names. S1A uses the first score row; S1B uses the score object only to define the sample intersection.
- tpm_root: PROJECT/PROJECT_mrna_expr_tpm.rdata containing mrna_expr_tpm, with unique gene and full sample identifiers. The stored TPM scale is retained.
- output_dir: an absent or empty, independent writable directory for that entrypoint.

Both use the entire exact TPM/score sample intersection within each cancer, without an added log transform or tumour-only restriction. S1A retains the source finite-pair and variance helper, Pearson cor.test, and BH across the nine genes within each cancer. S1B retains all combn(genes, 2) pairs, the finite-pair and variance helper, Spearman cor.test(exact = FALSE), and BH across the 36 pairs within each cancer. These are separate final-nine families; the candidate-selection BH family does not substitute for either. Stars retain source FDR cutoffs of <0.001, <0.01 and <0.05.

S1A writes panel_a_final_nine_complete_pearson.tsv (189 rows; project, gene, n_complete, pearson_r, pvalue, fdr, significance, gene_order, cancer_order). S1B writes panel_b_final_nine_pairwise_spearman_complete.tsv (756 rows; project, gene1, gene2, gene_pair, n_complete, spearman_rho, pvalue, fdr, significance, pair_order, cancer_order). Each also writes machine_checks.tsv and sessionInfo.txt. Applicable source QA retains nine genes, 21 cancers, the complete grid and nonmissing correlation coefficients (plus 36 pairs for B); an unevaluable correlation remains NA in the table and fails this QA rather than being silently dropped.

Run from the package root:

```sh
Rscript workflows/signature_selection/02_correlate_final_nine_tumor_program.R --help
Rscript workflows/signature_selection/02_correlate_final_nine_tumor_program.R --config config/my_s1a.dcf --check-config
Rscript workflows/signature_selection/03_correlate_final_nine_gene_pairs.R --config config/my_s1b.dcf --check-config
# Ordinary execution after inputs are available:
Rscript workflows/signature_selection/02_correlate_final_nine_tumor_program.R --config config/my_s1a.dcf
Rscript workflows/signature_selection/03_correlate_final_nine_gene_pairs.R --config config/my_s1b.dcf
```

Paths resolve relative to the DCF directory, including existing ancestors and symlinks; unresolved symlinks are rejected. Output paths must be independent of input/configuration locations and resolved TPM targets. --help and --check-config do not load research objects or tables, calculate statistics, hash inputs or create outputs. The latter checks only readable nonempty file/path metadata and does not validate the 21 project keys, nine genes, object schemas or project-specific TPM coverage. Ordinary runs perform calculations and source QA.

The source of A is final_figures/c1c4_atlas/fig1_acid_axis/panel_sources/r01_f1s1_01_prepare_and_plot_final_nine_panels.R; the source of B is r01_f1s1_02_prepare_and_plot_pairwise_synergy_and_counts.R in that directory. Correlation helper bodies, statistical blocks and output ordering are retained; private input/output paths are mapped to DCF fields, plotting dependencies and unrelated synergy/count blocks are removed, and identifier/path guards are added. Identifier rejection is an interface adaptation and does not filter eligible samples. B's source-tail 21-cancer QA is also asserted before calculation. Source A's geometric-mean synergy and source B's panel C counts are outside these entrypoints. Whole scripts are not byte-identical to their sources. No historical assets are overwritten.

Dependencies are base/recommended R (stats and utils). Original execution R/package versions for these two panel-source scripts have not been recovered. Public checks use the author-selected R4.5.1; this does not establish the historical generation environment. Verification is limited to parse/help, synthetic metadata boundary cases, source-block comparisons and small synthetic helper/statistical checks. No research input contents or hashes were inspected, and no real-input calculation was run. Input schema compatibility, numerical equivalence to the formal outputs and end-to-end reproduction remain unverified.
