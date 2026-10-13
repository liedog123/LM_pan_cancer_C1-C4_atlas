# Patient-level C1-C4 clinical outcomes

This entry continues from an upstream representative-patient manifest, separately frozen canonical labels and the formal clinical object. It is the complete Stage 1 clinical-model continuation, not AG survival, a risk-prediction model, readiness/sample selection, or complete downstream atlas reproduction. New outputs are independent copies and do not replace frozen results.

## Invocation and inputs

Run from the package root:

    Rscript workflows/clinical_outcomes/01_run_stage1_clinical_models.R --config config/clinical_outcomes.paths.example.dcf --check-config
    Rscript workflows/clinical_outcomes/01_run_stage1_clinical_models.R --config /path/to/config.dcf

Replace example paths before use. The strict DCF record has exactly four fields. Relative paths resolve from the DCF directory through existing ancestors and symlinks. Inputs must be readable nonempty regular files; output must be absent or empty, writable and independent of input parents and configuration. Root/home output targets are rejected. --help and --check-config load no analysis dependencies; the latter inspects path metadata only, reads no research contents and creates no output directories. An ordinary invocation actually performs all statistical models and the retained plots.

| Field | Required content |
|---|---|
| patient_manifest | TSV with Patient_ID, Sample, Cancer, C_subtype, AGIG, Acid_Score, LM12_inflammatory_axis, discordant_multi_sample_state, age, gender_clean, race_clean, stage4, TumorPurity, OS, OS_time, PFI, PFI_time, DSS, DSS_time, DFI, DFI_time. |
| canonical_labels | CSV with unique Sample and Cancer, C_subtype, AGIG; source checks manifest coverage and label agreement. |
| clinical_data | Rdata containing pancancer_data_list; cancer-named items expose clinical_data with sample, _PATIENT, ajcc_pathologic_tumor_stage and clinical_stage for source stage reconstruction. |
| output_dir | Independent absent/empty directory; results and figures are created below it only in an ordinary run. |

Patient_ID must already be unique. Representative-sample selection belongs to upstream 00_readiness and is not repeated here. Purity is already present in the manifest; no additional purity file is read. The manifest's stage4 is retained as stage4_stage0_interface, then source logic rebuilds stage4 from the clinical object. The source already loads clinical Rdata into an isolated environment; this is unchanged. Metadata checks do not validate any of these schemas or sample-selection properties.

## Preserved analysis and outputs

OS is primary; PFI, DSS and DFI are secondary. All seven source designs are retained: C1-C4 omnibus, AG1_IG1 versus AG2_IG1, AG2_IG1 versus AG2_IG2, AG2_IG2 versus AG3_IG2, C1 versus C4, continuous Acid and continuous inflammatory axis. These AGIG within-axis contrasts are not a separate AG-only survival workflow. Continuous axes are standardized within cancer.

The complete representative-patient analysis and exclusion of discordant multi-sample patients retain source event/sample gates, unadjusted and age/sex/stage/purity/race model branches, Efron Cox fits, evaluability/failure reporting and proportional-hazards audit. Source BH families remain intact: cancer effects by endpoint/design/model/analysis/effect; omnibus tests by endpoint/design/model/analysis; meta-results across designs within endpoint/model/analysis; clinical associations by feature/design/method/analysis. No OS-only or displayed-contrast subset is used to redefine multiplicity.

Custom REML random-effects pooling, heterogeneity and prediction intervals, leave-one-cancer-out refits, and age/sex/stage/race/purity associations are retained without algorithm changes. Clinical tests include Kruskal-Wallis, Wilcoxon, Spearman, chi-square or simulated Fisher tests under source rules. This entry does not add final evidence classification or readiness/finalization stages.

results/ contains the version-suffixed patient manifest, cancer Cox effects/omnibus tests/model statuses, random-effects meta-analysis, leave-one-cancer-out results, clinical associations, endpoint/state distributions, two OS plot-source tables and figure manifest. figures/ contains the original OS meta-forest and cancer-effect heatmap PDF/TIFF exports. The source suffix v0.1_20260722 is retained; it is not a new freeze designation. Retaining plots does not establish reproduction of frozen publication figures.

## Provenance and scope of checks

Source relative to the originating project:
8.0pancancer_c1c4_bulk_downstream_rebuild/01_clinical_outcomes/scripts/01_run_stage1_clinical_models_v0.1_20260722.R

Current source SHA-256: 0f9c6f6edc14a6d1848260284b151a0c1df43ff9ea810827ce88e69ed5992ff7

Text from version_tag <- "v0.1_20260722" through EOF is byte-identical. Only the execution header is adapted: private paths become DCF fields, CLI/metadata/output protection is added, result/figure paths move below independent output, and script_name uses the public filename for logs. Original seed 20260722 is retained for an ordinary run. This fingerprint describes the current source text, not a recovered historical run hash.

Current R4.5.1 checks cover parsing, eight synthetic metadata CLI cases and static source comparison only. No research object/table was loaded or hashed, and no clinical test, Cox/REML/LOO calculation or plot was run. Dependency completeness, real schemas, numerical reproduction and end-to-end equivalence remain unverified. Existing Stage 1 documentation records R4.5.1; historical package versions and complete execution provenance have not been recovered here.

Inputs, patient labels and frozen assets are not distributed. Access, release and licensing require an author decision; this workspace grants no license.
