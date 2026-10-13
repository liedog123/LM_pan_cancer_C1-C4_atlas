# External bulk Acid-linked LM12 direction validation

Use the [supplementary-input mapping](../../docs/reference_file_candidates.md#supplementary-material-first) for complete Table 2 LM12 worksheet-to-TSV exports and the ordered nine-gene `final_selected_genes` RData interface. The workbook is not read directly by these entries. Expression/metadata inputs remain separate; the conversion recipes have not been executed or numerically validated.

Two entries continue from separately obtained GEO expression files, prepared metadata and fixed signature definitions. They support cohort-specific Acid-linked LM12 enrichment directions, not external AG/IG or C1-C4 subtype assignment, a fixed TCGA classifier, or independent causal/mechanistic validation. They do not generate metadata, download raw data or reproduce Figure plotting. Formal Fig2C uses source 06 gene-level GSEA, not the earlier source 05 score-direction analysis.

| Entry | Continuation | Outputs | Original project source |
|---|---|---|---|
| 01_build_external_scores.R | Expression + prepared annotation + final acid genes + LM12 | Three cohort score matrices; combined external_score_matrix.tsv; expression-preparation and gene/signature-coverage QC | 5.0LM_Imm_12geneset/scripts/04_build_external_geo_score_matrix.R |
| 02_run_lm12_gene_level_gsea.R | Same expression + 01 cohort scores + LM12 | Primary/sensitivity groups, coverage, GSEA tables, gzipped full moderated-t rankings, RDS objects, combined results/QA/session information | 5.0LM_Imm_12geneset/scripts/06_external_lm12_gene_level_gsea.R |

The cohort processing order remains GSE40419/LUAD, GSE127165/LSCC, GSE96058/BRCA. Formal Fig2C displays LUAD; BRCA/LSCC are corresponding supplementary cohorts. The formal manifest records within-cohort moderated-t GSEA and BH interpretation across current 12 LM12 pathways, with no cross-cohort second BH. This packaging retains fgsea-returned padj; it does not add a new explicit n=12 adjustment.

## Inputs and exact filenames

Copy [external_scores.paths.example.dcf](../../config/external_scores.paths.example.dcf) or [external_gsea.paths.example.dcf](../../config/external_gsea.paths.example.dcf), replace placeholders and choose independent absent/empty output directories. Relative paths resolve from the DCF directory, including existing ancestors and symlinks.

01 requires exactly expression_dir, metadata_dir, final_genes, lm12_manifest, lm12_genes, output_dir. 02 requires exactly expression_dir, score_dir, lm12_manifest, lm12_genes, output_dir. expression_dir contains:

- GSE40419_LC-87_RPKM_expression.txt.gz: tab-delimited; gene column and six leading annotation columns, with sample expression columns after them.
- GSE127165_114_samples_RNA_seq_FPKM.txt.gz: tab-delimited; first column Ensembl IDs, remaining columns sample expression.
- GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz: comma-delimited; first column gene symbols, remaining sample columns. This is the author's transformed matrix, not assumed TPM.

metadata_dir contains external_cohort_manifest.tsv, GSE40419_sample_annotation.tsv, GSE127165_sample_annotation.tsv and GSE96058_expression_metadata_manifest.tsv. These are prepared source-03 intermediates, not raw GEO series matrices. Annotation needs matching Sample_ID, Cohort, Cancer and Sample_Type; Sample_Type uses Tumor/Normal for source 06 filtering. GSE96058 requires Is_Replicate interpretable by as.logical; Primary_Sample_ID is retained where available. Pair_ID/Has_Clinical_Metadata and other source-selected fields remain when present. The cohort manifest is read by source 04 but has no downstream column use in its score calculation. Metadata creation and clinical-source acquisition are outside these entries.

final_genes is Rdata containing named final_selected_genes with Gene column. lm12_manifest/lm12_genes are frozen definition TSVs: 01 uses feature_id/pathway/short_label and long feature_id/pathway/gene; 02 additionally requires display_order, direction_label and functional_group. 02 retains the source assertion of 12 unique features ordered 1..12. score_dir contains GSE40419_score_matrix.tsv, GSE127165_score_matrix.tsv and GSE96058_score_matrix.tsv from 01, with annotation columns plus Acid_Score and LM12 score columns. 02 uses Acid_Score/annotation to define groups and recomputes expression rankings; its GSEA is not performed on LM12 score values.

No expression data, patient annotation, gene sets or third-party resources are included. Matching processed expression and metadata must be obtained separately through author-approved access/release arrangements. Exact resource releases and redistribution permissions remain unresolved; generic accession downloads are not guaranteed to match these intermediates.

## Running and preserved calculations

Use Rscript workflows/external_bulk/01_build_external_scores.R --config config/my_external_scores.dcf [--check-config], followed by Rscript workflows/external_bulk/02_run_lm12_gene_level_gsea.R --config config/my_external_gsea.dcf [--check-config]. Both support --help. --check-config checks fixed-file/path metadata only: no gzip, input-content reads, analysis dependencies or output directories. Ordinary runs actually calculate ssGSEA or gene-level GSEA. No --test/cohort-filter option is added.

01 preserves trimming, PKM2-to-PKM aliasing, empty-symbol removal and duplicate-symbol means. LSCC strips Ensembl versions and maps via org.Hs.eg.db/AnnotationDbi mapIds with multiVals="first". All cohort expression uses log2(pmax(expression,0)+1); RPKM/FPKM/author-transformed sources are not relabeled TPM. Acid plus all 12 LM12 signatures are scored jointly by ssgseaParam(normalize=TRUE), preserving that normalization background and all supplied expression samples before annotation joins.

02 preserves its cohort-specific expression loaders, alias/duplicate processing and the same clipping/log transformation. Within available tumour samples with finite Acid_Score, it retains first distinct Sample_ID, orders descending Acid_Score then Sample_ID, assigns the first ceiling(n/2) High_Acid and the remainder Low_Acid. BRCA primary removes flagged technical replicates before grouping; the formal input had 136 repeats and 3273 primary samples, but the implementation filters flags rather than hardcoding a count. Its all_samples_sensitivity run is preserved.

The original model is model.matrix(~ group_factor), limma lmFit/eBayes and topTable coef=2, not voom. Finite full-gene moderated t values are ordered descending, with Gene tie ordering and duplicate-gene rejection. fgseaMultilevel uses minSize=10, maxSize=500 and eps=0. Original seed 20260804, cohort/sensitivity order, 12 finite-NES checks and combined 36-row/unique cohort-feature QA remain. No model, rank filtering, direction mapping or correction family is changed.

01 outputs *_score_matrix.tsv, external_score_matrix.tsv, external_expression_preparation_qc.tsv, acid_gene_coverage_detail.tsv and signature_coverage_qc.tsv. 02 outputs per-cohort *_sample_groups.tsv, *_coverage.tsv, *_results.tsv, *_ranked_t.tsv.gz and *_objects.rds, plus external_lm_imm_12_gene_level_gsea_results.tsv, *_qa.tsv, *_sensitivity_results.tsv and *_sessionInfo.txt. New outputs are independent copies; they do not replace formal scores, frozen GSEA or Figure source data.

## Source and verification boundary

Retained blocks begin at gene_alias_map and continue to EOF. 01 only replaces load(acid_gene_file) with isolated retrieval of final_selected_genes; its original exists check remains. 02's retained body is unchanged apart from trailing spaces/tabs. Headers replace private paths with DCF fields, preserve fixed filenames/seed and add early metadata/safe-output checks and delayed dependencies. Whole scripts are not byte-identical to their original headers.

Current small-script SHA-256 fingerprints, not historical execution hashes:

- Source 04: b8df522ef197597a2d29b28d95246f94bbf5b07c41d5f59db5c9805b73960846
- Source 06: 5e9dfd7d816e3308aa10a2d3461a7a4aa3af4ee457f16ee1d9e62b0f33fb16d8

The existing Methods provenance/session record identifies R4.5.1 for source 06 gene-level GSEA; this does not establish earlier source 04 score-generation GSVA/annotation versions, which remain unrecovered. Current public parsing/synthetic-interface checks separately use author-selected R4.5.1. Source defaults and org.Hs.eg.db mapping can depend on version; no portable historical dependency lockfile is supplied.

Only retained-block comparison, parse and synthetic metadata interfaces are checked here. No research file is read, decompressed or hashed, and no scores, limma, fgsea, sensitivity analysis or plot is run. Real schemas, sample correspondence, numerical equivalence, analysis dependency completeness and end-to-end reproduction remain unverified. Licensing/release remain pending. The formal evidence boundary is also recorded in final_figures/c1c4_atlas/fig2_c1c4_state_immune_ecology/source_data/Fig2C_external_bulk_lm12_gsea_manifest.tsv in the originating project.
