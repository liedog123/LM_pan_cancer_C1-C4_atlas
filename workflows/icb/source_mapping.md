# ICB source mapping

All paths below are relative to the originating project; SHA-256 identifies source code only. No study input was hashed.

## 01_prepare_riaz_gide.R

Source: `13.0c1c4_biological_context_and_icb/scripts/04_prepare_B.R` (111 lines including shebang); SHA-256 `6f06c811ced1f7e0f9da00a9e4efdcfcfbf08b7c13333590892c89d9c55807c9`.

The complete 110-line body after the shebang is retained below the header boundary. Reversing only this exact whitelist restores it byte-for-byte:

- `original private project-root assignment` → `# Project paths are supplied by the checked DCF.`
- `out <- file.path(root,"13.0c1c4_biological_context_and_icb/results")` → `out <- unname(cfg_paths["output_dir"])`
- `04_prepare_B.R` → `01_prepare_riaz_gide.R`
- `original private ICB input-directory assignment` → `# Explicit cohort files replace the private input directory.`
- `file.path(root,"5.0LM_Imm_12geneset/results/02_final_lm_imm_12_gene_sets/lm_imm_12_gene_sets_long.tsv")` → `unname(cfg_paths["lm12_members"])`
- `file.path(root,"7.3scRNA_c1c4_state_closure/02_program_catalog/program_nondefinition_gene_membership.tsv.gz")` → `unname(cfg_paths["program_membership"])`
- `file.path(root,"10.0paired_ici_dynamic_remodeling/00_readiness/gene_symbol_alias_manifest.tsv")` → `unname(cfg_paths["gene_alias"])`
- `file.path(root,"10.0paired_ici_dynamic_remodeling/00_readiness/GSE91061_prior_ipilimumab_mapping.tsv")` → `unname(cfg_paths["prior_ipilimumab"])`
- `file.path(dat,paste0(prefix,".Response.tsv"))` → `unname(cfg_paths[paste0(tolower(cohort),"_metadata")])`
- `file.path(dat,paste0(prefix,".Response.Rds"))` → `unname(cfg_paths[paste0(tolower(cohort),"_expression")])`
- `file.path(dat,"RCC-Braun_2020.Response.tsv")` → `unname(cfg_paths["braun_hold_metadata"])`
- `Riaz/Gide baseline states and all three-cohort source inclusion lists written` → `Riaz/Gide baseline states and historical Braun HOLD inclusion list written`

The new base-R header supplies strict DCF/CLI checks, nonempty regular-file metadata checks, symlink-aware resolution, input/config/package-root overlap protection, writable empty independent output validation, help/check-config early exit and ordinary-only output creation. It introduces no scoring/model/statistical logic.

## 02_analyze_riaz_gide.R

Source: `13.0c1c4_biological_context_and_icb/scripts/05_analyze_B.R` (145 lines including shebang); SHA-256 `147d4baffe2fc049cfaf5576638f0ecae81fa369f0cd758da1aa6f4ed8415e36`.

The complete 144-line body after the shebang is retained below the header boundary. Reversing only this exact whitelist restores it byte-for-byte:

- `original private project-root assignment` → `# Project paths are supplied by the checked DCF.`
- `out <- file.path(root,"13.0c1c4_biological_context_and_icb/results")` → `out <- unname(cfg_paths["output_dir"])`
- `05_analyze_B.R` → `02_analyze_riaz_gide.R`
- `file.path(out,"B_baseline_patient_inclusion_and_scores.tsv")` → `unname(cfg_paths["patient_scores"])`
- `three-cohort fixed response slots; exact tests, Firth sensitivity and programme exploration` → `Riaz/Gide analyses with historical unresolved Braun slots; exact tests, Firth sensitivity and programme exploration`

The new base-R header supplies strict DCF/CLI checks, nonempty regular-file metadata checks, symlink-aware resolution, input/config/package-root overlap protection, writable empty independent output validation, help/check-config early exit and ordinary-only output creation. It introduces no scoring/model/statistical logic.

The historical Braun unresolved source stage is preserved, not promoted to the resolved final integration. See README for the source 01→02 handoff and 03–05 integration route.


ICB 01-02 verification (2026-10-03): two R4.5.1 (2025-06-13) parses, two help checks and 78 synthetic invalid-schema CLI/metadata/path cases passed (01: 43; 02: 35; ten success and 68 expected rejection). Metadata successes created no output directories; nonempty-output sentinels remained unchanged. Reversing the documented whitelist restores all 111/145 source lines including shebang byte-for-byte. Source code alone was fingerprinted. All four new R/DCF files are trackable; new-file private-path scans and 66 entries/14 modules were checked. Actual stdout/exit records and synthetic fixtures are maintainer-held in interface_final_checks.log, interface_final_summary.json, static_checks.json and package_checks.log. No study object/table/metadata/schema/manifest/cache/content or research checksum, analysis dependency, scoring/clustering/Firth self-check/model/bootstrap/statistical body or figure was loaded/executed; existing 64 entries were not retested. Input acquisition/release, real schemas/dependency compatibility and numerical/end-to-end reproduction remain unverified.

## 03_prepare_braun_official.R

Source: `13.0c1c4_biological_context_and_icb/braun_source_resolution/scripts/01_prepare_braun_official.R`; 205 complete lines including shebang. Source-code SHA-256: `bd72fe7ec20240999bbc08a05f5e007fd7fe4452e64b5b3523ef3a37fe845acd`. The full source body is retained below the header boundary; inverse adaptation restores it byte-for-byte. Whitelist:

- `original private absolute-path assignment` → `# Private project root replaced by explicit DCF inputs.`
- `module <- file.path(root,"13.0c1c4_biological_context_and_icb")` → `# Module paths are explicit DCF fields.`
- `resolution <- file.path(module,"braun_source_resolution")` → `# Resolution output is independently configured.`
- `input_file <- file.path(module,"inputs","Braun_2020_official_S1_S4_MOESM2.xlsx")` → `input_file <- unname(cfg_paths["official_workbook"])`
- `out <- file.path(resolution,"results")` → `out <- unname(cfg_paths["output_dir"])`
- `01_prepare_braun_official.R` → `03_prepare_braun_official.R`
- `file.path(root,"5.0LM_Imm_12geneset/results/02_final_lm_imm_12_gene_sets/lm_imm_12_gene_sets_long.tsv")` → `unname(cfg_paths["lm12_members"])`
- `file.path(root,"7.3scRNA_c1c4_state_closure/02_program_catalog/program_nondefinition_gene_membership.tsv.gz")` → `unname(cfg_paths["program_membership"])`
- `file.path(root,"10.0paired_ici_dynamic_remodeling/00_readiness/gene_symbol_alias_manifest.tsv")` → `unname(cfg_paths["gene_alias"])`
- `original private absolute-path assignment` → `unname(cfg_paths["local_clinical_crosswalk"])`
- `system2("sha256sum",input_file,stdout=TRUE,stderr=TRUE)` → `system2("sha256sum",shQuote(input_file),stdout=TRUE,stderr=TRUE)`

Metadata checks use the inherited strict DCF/CLI, regular-file and symlink/overlap/output guards. No scientific logic, defaults, eligibility or status is changed.

## 04_analyze_braun_official.R

Source: `13.0c1c4_biological_context_and_icb/braun_source_resolution/scripts/02_analyze_braun_official.R`; 160 complete lines including shebang. Source-code SHA-256: `1c2dc3f67e2cd442beb7819c5647cd7da07627e0e2e14897d4f0b26efbef3c52`. The full source body is retained below the header boundary; inverse adaptation restores it byte-for-byte. Whitelist:

- `original private absolute-path assignment` → `# Private project root replaced by explicit DCF inputs.`
- `module <- file.path(root,"13.0c1c4_biological_context_and_icb")` → `# Module paths are explicit DCF fields.`
- `out <- file.path(module,"braun_source_resolution","results")` → `out <- unname(cfg_paths["output_dir"])`
- `02_analyze_braun_official.R` → `04_analyze_braun_official.R`
- `file.path(out,"B_Braun_patient_inclusion_and_scores.tsv")` → `unname(cfg_paths["patient_scores"])`
- `file.path(out,"B_Braun_member_coverage.tsv")` → `unname(cfg_paths["member_coverage"])`

Metadata checks use the inherited strict DCF/CLI, regular-file and symlink/overlap/output guards. No scientific logic, defaults, eligibility or status is changed.

## 05_integrate_three_cohort.R

Source: `13.0c1c4_biological_context_and_icb/scripts/09_integrate_Braun_three_cohort_v2.R`; 126 complete lines including shebang. Source-code SHA-256: `f084836e6a35eec5569ce8c2dfedf4602f316ee1291a980f3eab6b0d3afded77`. The full source body is retained below the header boundary; inverse adaptation restores it byte-for-byte. Whitelist:

- `original private absolute-path assignment` → `# Private project root replaced by explicit DCF inputs.`
- `base_res <- file.path(root,"13.0c1c4_biological_context_and_icb/results")` → `# Riaz/Gide inputs route to their separate preparation/analysis directories.`
- `braun_res <- file.path(root,"13.0c1c4_biological_context_and_icb/braun_source_resolution/results")` → `# Braun inputs route to their separate preparation/analysis directories.`
- `out_dir <- file.path(base_res,"B_three_cohort_v2_20260924")` → `out_dir <- unname(cfg_paths["output_dir"])`
- `file.path(braun_res, required)` → `icb_input("braun", required)`
- `file.path(base_res, old_name)` → `icb_input("riaz_gide", old_name)`
- `file.path(braun_res, braun_name)` → `icb_input("braun", braun_name)`
- `file.path(base_res,"B_baseline_patient_inclusion_and_scores.tsv")` → `icb_input("riaz_gide","B_baseline_patient_inclusion_and_scores.tsv")`
- `file.path(braun_res,"B_Braun_patient_inclusion_and_scores.tsv")` → `icb_input("braun","B_Braun_patient_inclusion_and_scores.tsv")`
- `file.path(base_res,"B_member_coverage.tsv")` → `icb_input("riaz_gide","B_member_coverage.tsv")`
- `file.path(braun_res,"B_Braun_member_coverage.tsv")` → `icb_input("braun","B_Braun_member_coverage.tsv")`
- `file.path(braun_res,nm)` → `icb_input("braun",nm,optional=TRUE)`

Metadata checks use the inherited strict DCF/CLI, regular-file and symlink/overlap/output guards. No scientific logic, defaults, eligibility or status is changed.

03 additionally quotes the configured workbook argument to sha256sum, retaining the original expected hash and ordinary-only check. 04 separates both patient and coverage reads from its output. 05's new fixed-basename icb_input helper routes every source read to one of the four configured stage directories; required16 file checks precede content/dependency loading, and the original required8 Braun/source body checks remain. Three optional provenance files route from Braun preparation with source missing/empty skip semantics. 05 preserves source absent-output/nonrecursive-creation behavior; help/config does not create its output. Exact private source text is confined to maintainer task artifacts, not this public document.


ICB 03–05 verification (2026-10-03): three R4.5.1 (2025-06-13) parses, three help checks, 55 synthetic invalid-schema CLI/metadata/path cases (03: 18; 04: 18; 05: 19; twelve successes and 43 expected rejection) and eleven isolated base-R fixed-basename router cases passed. Metadata successes created no outputs; nonempty-output sentinels were preserved. All 205/160/126 complete source lines compare identically after reversing only the documented path/identifier/shQuote whitelist. The workbook argument quoting was statically checked; no workbook hash command was run. Six new R/DCF files are trackable, private-path checks and 69 entries/14 modules passed. Actual stdout/exit records and fresh mktemp fixtures are maintainer-held in interface_checks.log, interface_summary.json, helper_checks.log, static_checks.json and package_checks.log. No study workbook/object/table/metadata/schema/manifest/cache/content or research checksum, analysis dependency, scoring/clustering/Firth/bootstrap/BH/statistical body or figure was loaded/executed; existing 66 entries were not retested. Input acquisition/release, real schemas/dependency/tool compatibility and numerical/end-to-end reproduction remain unverified.

## 06_prepare_liu_exploratory.R

Source: `13.0c1c4_biological_context_and_icb/scripts/07_prepare_Liu_exploratory.R`; 313 full lines including shebang. Source-code SHA-256: `4dfc3aac6cf26c910f962cfcbb024013972d16371e421ec2fd3b01b0930a01ae`. Complete-body inverse comparison is byte-for-byte; only the following path/identifier replacements are allowed:

- `original private project-root assignment` → `# Private project root replaced by explicit DCF inputs.`
- `module <- file.path(root, "13.0c1c4_biological_context_and_icb")` → `# Module paths are supplied by the checked DCF.`
- `source_dir <- file.path(module, "exploratory_Liu2019_source")` → `# Prepared source files are independently configured.`
- `out_dir <- file.path(module, "results/exploratory_Liu2019_tigeR_103_119")` → `out_dir <- unname(cfg_paths["output_dir"])`
- `source_rda <- file.path(source_dir, "MEL_phs000452.rda")` → `source_rda <- unname(cfg_paths["source_rda"])`
- `source_tpm <- file.path(source_dir, "41591_2019_654_MOESM3_ESM.txt")` → `source_tpm <- unname(cfg_paths["official_s2_tpm"])`
- `source_clin <- file.path(source_dir, "41591_2019_654_MOESM4_ESM.xlsx")` → `source_clin <- unname(cfg_paths["clinical_workbook"])`
- `file.path(root, "10.0paired_ici_dynamic_remodeling/00_readiness/gene_symbol_alias_manifest.tsv")` → `unname(cfg_paths["gene_alias"])`
- `file.path(root, "5.0LM_Imm_12geneset/results/02_final_lm_imm_12_gene_sets/lm_imm_12_gene_sets_long.tsv")` → `unname(cfg_paths["lm12_members"])`
- `file.path(root, "7.3scRNA_c1c4_state_closure/02_program_catalog/program_nondefinition_gene_membership.tsv.gz")` → `unname(cfg_paths["program_membership"])`
- `07_prepare_Liu_tigeR_FPKM_outcome_blind.R` → `06_prepare_liu_exploratory.R`

The inherited strict DCF/CLI/regular-file/symlink/overlap/empty-output guards exit before dependencies/study reads. The full source body, functions, models, seeds, source assertions, non-estimable statuses and output checks remain; no source() helper is added or required.

## 07_analyze_liu_exploratory.R

Source: `13.0c1c4_biological_context_and_icb/scripts/08_analyze_Liu_exploratory.R`; 186 full lines including shebang. Source-code SHA-256: `a7da4dd14c9e73a5bf72549e77de12041696dc3991b4d6fbc0a7d431fc281c67`. Complete-body inverse comparison is byte-for-byte; only the following path/identifier replacements are allowed:

- `original private project-root assignment` → `# Private project root replaced by explicit DCF inputs.`
- `module <- file.path(root,"13.0c1c4_biological_context_and_icb")` → `# Module paths are supplied by the checked DCF.`
- `source_clin <- file.path(module,"exploratory_Liu2019_source/41591_2019_654_MOESM4_ESM.xlsx")` → `source_clin <- unname(cfg_paths["clinical_workbook"])`
- `res_dir <- file.path(module,"results/exploratory_Liu2019_tigeR_103_119")` → `res_dir <- unname(cfg_paths["output_dir"])`
- `file.path(res_dir,"Liu_103_119_outcome_blind_state_assignments.tsv")` → `unname(cfg_paths["outcome_blind_assignments"])`
- `08_analyze_Liu_exploratory.R` → `07_analyze_liu_exploratory.R`

The inherited strict DCF/CLI/regular-file/symlink/overlap/empty-output guards exit before dependencies/study reads. The full source body, functions, models, seeds, source assertions, non-estimable statuses and output checks remain; no source() helper is added or required.

## 08_run_gide_c1_sensitivity.R

Source: `13.0c1c4_biological_context_and_icb/exploratory_c1_regimen_axis/01_explore_gide_c1.R`; 323 full lines including shebang. Source-code SHA-256: `4b393b5340e7a0e37799308a4022e9670565e3c0ff0a287fa01a622599e97f9e`. Complete-body inverse comparison is byte-for-byte; only the following path/identifier replacements are allowed:

- `original private project-root assignment` → `# Private project root replaced by explicit DCF inputs.`
- `mod <- file.path(root, "13.0c1c4_biological_context_and_icb")` → `# Module paths are supplied by the checked DCF.`
- `src <- file.path(mod, "results/B_baseline_patient_inclusion_and_scores.tsv")` → `src <- unname(cfg_paths["patient_scores"])`
- `out <- file.path(mod, "exploratory_c1_regimen_axis")` → `out <- unname(cfg_paths["output_dir"])`
- `file.path(mod, "results/A_frozen_program_members.tsv")` → `unname(cfg_paths["frozen_program_members"])`
- `file.path(root,
  "5.0LM_Imm_12geneset/results/02_final_lm_imm_12_gene_sets/lm_imm_12_gene_sets_long.tsv")` → `unname(cfg_paths["lm12_members"])`
- `01_explore_gide_c1.R` → `08_run_gide_c1_sensitivity.R`

The inherited strict DCF/CLI/regular-file/symlink/overlap/empty-output guards exit before dependencies/study reads. The full source body, functions, models, seeds, source assertions, non-estimable statuses and output checks remain; no source() helper is added or required.

06's six prepared sources retain source RDA isolation/full-gene TPM conversion, exact Patient### controls and independent 103/119 assignment. 07's sole assignment read is explicit and separate from its output. 08 explicitly reads patient scores plus transcriptome_proteome 14's A_frozen_program_members.tsv and LM12; results/figures subdirectories and eight-table/two-diagnostic exports remain. Post-hoc Liu/Gide do not enter 05's fixed three-cohort BH; necessary final exports and input/release/availability remain for unified delivery, with CD8 counts next.


ICB 06–08 verification (2026-10-03): three R4.5.1 (2025-06-13) parses, three help checks and 54 synthetic invalid-schema CLI/metadata/path cases passed (18 each; twelve successes, 42 expected rejection). Metadata successes created no outputs; nonempty-output sentinels were preserved. All 313/186/323 full source lines including shebang compare byte-for-byte after reversing only the documented path/identifier whitelist. Source code alone was fingerprinted. Six new R/DCF files are trackable; narrow new-file private-path checks and 72 entries/14 modules passed. No external source() helper is present. Actual stdout/exit records and fresh mktemp fixtures are maintainer-held in interface_checks.log, interface_summary.json, static_checks.json and package_checks.log. No study RDA/workbook/object/table/metadata/schema/manifest/cache/content or research checksum, analysis dependency, source-concordance/scoring/clustering/Firth/bootstrap/BH/statistical body or diagnostic figure was loaded/executed; existing 69 entries were not retested. Input acquisition/release, real schemas/dependency/device compatibility and numerical/end-to-end reproduction remain unverified.
