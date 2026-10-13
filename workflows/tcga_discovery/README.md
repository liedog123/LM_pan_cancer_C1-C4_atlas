# TCGA discovery: High/Low rankings and ImmPort GSEA

These two entries continue from the study's processed per-cancer object, not raw TCGA acquisition or acid-signature discovery. Independent outputs feed the existing [LM12 construction](../lm12_construction/README.md) workflow; they do not replace frozen results, LM12 definitions or canonical labels.

| Entry | Inputs | Outputs | Original project source |
|---|---|---|---|
| 01_build_high_low_deg.R | processed_data | full_degs_for_gsea/*_DEG_full_results.csv; significant_degs/*_DEG_significant.csv when nonempty | 4.group_deg/1.0pancancer_mrna_deg.R |
| 02_run_immport_gsea.R | deg_results_dir; immport_gmt | individual_full_results/*_GSEA_full_results.csv; individual_sig_results/*_GSEA_sig_results.csv when nonempty; all_gsea_results.rds | 4.group_deg/2.0pancancer_mrna_gsea.R |

## Inputs and configuration

Copy [discovery_de.paths.example.dcf](../../config/discovery_de.paths.example.dcf) or [discovery_gsea.paths.example.dcf](../../config/discovery_gsea.paths.example.dcf), replace placeholders and use a separate absent or empty output directory. DCF-relative paths resolve from the configuration file directory, including symlinks/existing ancestors. 01 requires exactly processed_data/output_dir; 02 requires exactly deg_results_dir/immport_gmt/output_dir.

The processed Rdata must contain named pancancer_data_list: each cancer supplies expression_data (gene-row, sample-column TPM matrix with gene/sample names) and clinical_data with sample and Group. The source selects High/Low, indexes expression columns in clinical sample order and skips incomplete/insufficient items; no additional sample cleaning is added. AcidScore/AG are upstream and are not recomputed here.

02 reads every file selected by the original regular expression "_DEG_full_results.csv$" (the period remains unescaped), deriving cancer names by the same substitution. CSV inputs need gene and t columns; ranks are moderated t from 01, not logFC. Only NA statistics, NA names and empty names are removed, followed by descending sort; no Inf/duplicate-gene filter was introduced. GMT must be readable by clusterProfiler::read.gmt: tab-separated pathway name, description and gene members yielding term/gene columns. The exact study release, access and redistribution permission have not been established here. Obtain matching materials separately: no GMT, patient data or third-party materials are included. Generic downloads are not guaranteed to reproduce study inputs.

Use Rscript workflows/tcga_discovery/01_build_high_low_deg.R --config config/my_de.dcf [--check-config], then Rscript workflows/tcga_discovery/02_run_immport_gsea.R --config config/my_gsea.dcf [--check-config]. Both support --help. --check-config validates metadata only, before dependencies, content reads or directory creation. It does not validate schemas, sample correspondence or numbers. Ordinary runs actually calculate DEG or GSEA; no --test/cancer-filter option is added.

## Preserved calculations and continuation

01 preserves historical TPM+1 input to voom (voom(expr_matrix + 1, design, plot = FALSE)), limma model/High-minus-Low contrast, eBayes and topTable. Packaging does not substitute counts or logTPM. Significant DEG uses adj.P.Val < 0.05 and abs(logFC) > 1.0. Full output retains topTable columns plus gene, including t.

02 prefilters pathway sizes to [5,500]. For Jaccard > 0.5 it removes the longer pathway name (ties remove the second), exactly as the source. fgsea uses minSize=10, maxSize=500, nPermSimple=100000. Significant output uses padj < 0.05 and descending NES; CSV leadingEdge is comma-separated, whereas the aggregate RDS retains per-cancer fgsea tables with list leading edges. The source has no set.seed; this remains unchanged. Randomness and unknown historical versions preclude exact numeric reproduction guarantees.

Point LM12-construction 01 gsea_results_dir to 02's individual_full_results. Historical GSEA voting/condensation/plotting are excluded; current LM12 construction remains separate and requires curated references.

## Extraction and verification boundary

01 retains from "cancer_types <- names(pancancer_data_list)" through EOF. The header replaces private paths, extracts only pancancer_data_list into an isolated environment and creates independent output directories. 02 is a cross-section extraction: original GSEA_PADJ_THRESHOLD assignment, section 3 preprocessing, section 4 through saveRDS, and section 5's empty-result stop check. Unused historical voting thresholds/dependencies and subsequent voting/condensation/plotting are not copied. Only trailing spaces/tabs are normalized in retained blocks; entire source files are not byte-identical to these entries.

Current source SHA-256 fingerprints (small code files, not historical run fingerprints):

- 4.group_deg/1.0pancancer_mrna_deg.R: b7cc567bdfbd83290e003ae59462e493896e10208541aaea2ac8add2e9f8f81b
- 4.group_deg/2.0pancancer_mrna_gsea.R: b5e4fec9402ca4b3f6ec43337c12e10423f1628acbacb875106aa7f46770db9d

The second source calls itself version 3.2; this does not recover its execution environment. The first has no explicit version tag. R 4.5.1 is the author-selected current static/interface-test environment, not a historical version. Checks are limited to retained-block comparison, parsing and synthetic metadata interfaces. No study object, CSV/GMT content was loaded; no DEG/fgsea/plot/actual export ran. Real schemas, dependency completeness, numerical/end-to-end reproduction remain unverified. Input access, licensing and release remain undecided.
