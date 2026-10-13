# Reference inputs and supplementary-table mapping

Public release: code/configuration only. References described as locally packaged or supplied below belong to internal preparation and are not distributed here. Software license and destination are recorded in [release scope](../docs/release_scope.md).

Updated: 2026-10-03

This inventory includes nineteen author-approved unchanged local copies in [references](../references/README.md), connected across twenty-three distinct example DCFs. The eight added bulk/immune/drug tables retain all fields, NA values, original TCPA tokens and statistical metadata. Other candidates remain source locators; packaging does not grant a license or external redistribution permission.

“Definition” describes a file's analytical role, not a guarantee about all its contents or rights. Study-curated mappings and third-party gene-set memberships are distinguished below. Source identities and consumer keys are recorded; current download availability, licenses and numerical agreement were not tested.

## Supplementary material first

Use the paper's supplementary materials when they already supply the adopted definition. Do not duplicate a table in the code repository solely because a script expects a separate file. The code package must still identify the exact worksheet, required format and configuration field.

### Supplementary Table 2 to LM12 inputs

The current author-confirmed workbook is `Supplementary_Table_2_fig2_author_confirmed_v1.0_20260920.xlsx`, supplied through the manuscript's supplementary-material channel when released. No public download URL/accession is established in this workspace yet.

| Worksheet | Suggested local TSV name | Public consumer and DCF field |
|---|---|---|
| `Pathway_definitions` | `lm12_pathway_definitions.tsv` | shared_definitions 01 and external_bulk 01/02: `lm12_manifest`; single_cell 04: `lm12_manifest` |
| `Gene_membership` | `lm12_gene_membership.tsv` | shared_definitions 01 and external_bulk 01/02: `lm12_genes`; single_cell 04: `lm12_gene_long`; transcriptome_proteome 10/11 and cd8_perturbation 03/04: `lm12_membership`; immune_ecology 01: `lm12_genes`; icb 01/03/06/08: `lm12_members` |

Export each entire worksheet as a separate UTF-8, tab-delimited text file with its original header, no row-number column and no added title rows. Keep all columns, rows, identifiers, row order, missing values and numeric values unchanged. Do not deduplicate memberships across pathways, filter for significance, rename feature IDs or rerun pathway selection. Use a new local input directory; leave the workbook and study originals unchanged. Set the relevant DCF field to the resulting TSV path, not to the XLSX file. Paths resolve relative to the configuration file as documented by each entry.

The most demanding checked manifest consumer requires `display_order, feature_id, short_label, pathway, direction_label, feature_class, functional_group, stable_gene_count`. Membership consumers collectively require `feature_id, gene, pathway, functional_group`. The inspected worksheet headers contain all of these fields; keeping the complete worksheets also preserves their other selection metadata. Manifest and membership worksheets are not interchangeable.

The recorded workbook source for `Pathway_definitions` is `lm_imm_12_pathway_vote_summary.tsv`, not a file inferred from the DCF field name. Source finalization writes the same data frame to that file and `lm_imm_12_manifest.tsv`. The `Gene_membership` worksheet comes from `lm_imm_12_gene_sets_long.tsv`. Existing export provenance is recorded at `final_figures/c1c4_atlas/manuscript/revisions/yanmin_review_r01/fig2_section_review_20260918/fig2_supplementary_statistics_source_manifest_review_v0.1_20260918.tsv`.

This is a source/field-level connection, not a newly tested export or a numerical reproduction claim. TSVs do not replace downstream scoring bundles, CD8 bridge manifests, the prepared expression inputs or the complete upstream discovery/candidate inputs.

### Nine-gene signature: published definition, two file formats

The Results text explicitly lists `ATP6AP1, GAPDH, HIF1A, HK2, PFKP, PGK1, PKM, SLC2A1, SLC9A7`. The same ordered list is already embedded in public icb 01/03/06. The membership information therefore does not require a separate duplicate publication. Supplementary Table 1, which contains tumour–normal and AG survival statistics, is not the source of a separate nine-gene worksheet.

- For signature_selection 02/03, prepare a one-column TSV headed `Gene`, one symbol per row in the order above, and use its path for `final_genes`. This is a minimal membership interface, not a regenerated copy of the four-column formal figure-source table. An approved existing `Fig1A_signature_genes.tsv` remains another source.
- For acid_ag_preprocessing 01/02/05, external_bulk 01 and immune_ecology 01, the current consumer expects an RData file containing a data frame named `final_selected_genes` with a character `Gene` column. In R, construct that data frame from the same ordered nine symbols and use `save(final_selected_genes, file = <new RData path>)`. Configure `final_genes`, or immune_ecology's `acid_genes`, with that RData path. A TSV or RDS path is not interchangeable with this interface.
- These consumers use the `Gene` column, not a historical selection `Count`. A minimal membership-only object does not reproduce the original candidate-selection statistics or historical RData contents. Do not invent counts or overwrite the original object. This recipe has not been executed in the documentation task.

### Remaining definition inputs not supplied by those worksheets

For reconstruction of selection, D1 (142 candidates), D3 (complete curated pathway candidates) and D4 (family dictionary) are locally packaged in references/ under their original basenames. Original source paths are below; external release/access is still undecided. D2 is an optional nine-gene interface, not a mandatory separate upload.

The final 12 pathways cannot replace the complete candidate/family inputs: lm12_construction 02 first selects `analysis_tier == "primary"` from the full candidate table, joins all pathway votes to `local_term`, and only then applies the recurrence rules. Reconstructing candidates from the final 12 would omit unselected candidates and change the audit/selection scope.

D1 needs `gene` with the adopted 142 entries and order. D3 needs `analysis_tier, candidate_set, scope_layer, family_id, local_term, match_type, note`; D4 needs `scope_layer, family_id, fig1_or_project_basis, note`. The current supplementary worksheets do not provide these complete preselection interfaces. The three unchanged local copies are now supplied; external release remains a separate decision.

| ID | Exact existing source file | Consumer / role | Basis and remaining check |
|---|---|---|---|
| D1 | `1.0gene_select/revisions/yanmin_review_r01/v0.1_20260831/inputs/candidate_acid_lactate_genes_142_v0.1_20260831.tsv` | signature_selection 01, `candidate_genes` | 142 rows and three columns: source_row, gene, function_annotation. The candidate-set source is author-confirmed Wang et al. (2024), doi:10.1016/j.celrep.2024.113796. The freeze imports Function from the source workbook; original attribution of each function annotation remains unverified. Preserve the full table and adopted order; locally packaged unchanged; external release/access remains undecided. |
| D2 | `final_figures/c1c4_atlas/fig1_acid_axis/source_data/Fig1A_signature_genes.tsv` | signature_selection 02/03, `final_genes` | The formal nine-gene panel source reads this TSV and asserts nine genes. Columns include `Gene_Order, Gene, Score_Method, Status`. Keep separate from the RData interface below. |
| D3 | `5.0LM_Imm_12geneset/references/curated_candidate_pathway_terms.tsv` | lm12_construction 02, `candidate_terms` | Locally packaged unchanged: 69 rows and seven columns. Study-curated candidate/tier mappings over local ImmPort terms; retain the complete universe, IDs, classifications and notes. note participates in distinct() and is carried into candidate outputs, so do not silently remove or rewrite it. |
| D4 | `5.0LM_Imm_12geneset/references/candidate_pathway_family_dictionary.tsv` | lm12_construction 02, `family_dictionary` | 17 rows and seven columns: twelve paper_matched families and five acid_tme_extension families. Preserve the full mapping and rationale; fig1_or_project_basis/note become result family_label/family_description. Locally packaged unchanged; external release/access remains undecided. |

All four source paths exist. D1/D3/D4 local copies were compared byte-for-byte with their sources using cmp. D2 can use the published membership recipe. No source was modified and nothing was externally uploaded.

The local ImmPort term/GMT provenance and study-assigned family/tier classifications are distinct. A subsequent methodological review points to Xia et al., *A method for predicting drugs that can boost the efficacy of immune checkpoint blockade*, Nature Immunology 25, 659–670 (2024), doi:10.1038/s41590-024-01789-x, as the reference behind the Fig1/paper terminology; this is not the original table-creation log or row-by-row attribution. Later ImmPort/Reactome/ImmuCellAI and C7 catalogues are separate resources. Read-only inspection of D1/D3/D4 found no patient identifiers or private paths. Remaining decisions concern exact ImmPort release/acquisition, original curation records where available, and authorized access/publication routes—not a separate paper for every pathway.

## Existing shared-definition sources and alternatives

| ID | Exact existing source | Consumer / issue |
|---|---|---|
| D5 | `1.0gene_select/final_selected_genes.rdata` | acid/AG and external bulk use `final_genes`; other audits use `acid_genes`. Source code expects `final_selected_genes` with a `Gene` column. The historical file exists but was not opened. The membership-only RData recipe above offers a format-specific alternative for the listed consumers without distributing this original object. |
| D6 | `5.0LM_Imm_12geneset/results/02_final_lm_imm_12_gene_sets/lm_imm_12_manifest.tsv` | Shared/external consumers use `lm12_manifest`. The module README designates this final result directory. The file combines feature/pathway annotations with selection statistics. Shared assignment explicitly consumes `display_order, feature_id, short_label, pathway, functional_group, direction_label, stable_gene_count`; use the complete Table 2 Pathway_definitions export for the mapped consumers, retaining these columns. |
| D7 | `5.0LM_Imm_12geneset/results/02_final_lm_imm_12_gene_sets/lm_imm_12_gene_sets_long.tsv` | Shared/external consumers use `lm12_genes`; RNA workflows also use `lm12_membership`. The file combines `feature_id/pathway/gene` memberships with leading-edge/vote fields. The Table 2 Gene_membership export supplies the mapped membership interface. Do not recompute or alter the frozen membership. |
| D8 | `5.0LM_Imm_12geneset/exploratory/02_acid_state_specific_immune_ecology/v0.1_20260725/references/candidate_immune_gene_set_registry.tsv` and `candidate_immune_gene_sets_long.tsv` in the same directory | RNA local-GSEA consumers use `immune_registry` and `immune_membership`. Both are locally packaged unchanged. This branch's README designates the registry as its Stage A authority and records 106 ImmPort, 12 Reactome and one ImmuCellAI set. Study-selected metadata and external-derived memberships require exact source/permission treatment; they are not LM12 or a generic latest GMT. |

D6/D7 combine definitions with statistical metadata. A complete Table 2 worksheet export retains the existing published fields and avoids a second uploaded dataset. No export was performed here.

## Definitions already encoded in public scripts

The formal 12-member pH panel is already explicit in `workflows/acid_ag_preprocessing/04_correlate_formal_ph_panel.R`; no external gene-list configuration is consumed and no duplicate list file is required. It is distinct from the nine-gene scoring signature. Source provenance and analytical role remain documented in that workflow.

CD8 03 also encodes its fixed nine MSigDB selections; CD8 04 encodes its 18-marker/four-module definition. Do not duplicate these as mandatory input files. Their generated manifests can contain resource digests, detection summaries or other derived fields, and are not automatically pure reference tables.

## Other candidate sources: check supplementary coverage before copying

For compactness, the following prefixes are exact project-relative directories:

- `T = 8.0pancancer_c1c4_bulk_downstream_rebuild/03_transcriptome_proteome/`
- `S = 7.3scRNA_c1c4_state_closure/`
- `G = 7.1scRNA_lactate_ecosystem_rebuild/01_input_contract/`

The files below exist according to targeted metadata checks. Analytical roles and field examples are inferred from consumers/builders unless explicitly noted; the three adopted programme tables, including the compressed FULL/NONDEFINITION memberships, were subsequently read only for their bounded input-interface/content check; other compressed candidate memberships were not opened. These are source locators, not a requirement to bundle every file. Apply the same supplementary-first rule before considering a separate release; detailed interface checks cover the LM12/nine-gene routes and the locally packaged Fig. 6 selections and programme references, plus the packaged exclusion/contract/alias tables; other rows remain source locators.

| Group | Exact path (expand the prefix above) | Consumer and preparation boundary |
|---|---|---|
| Cross-layer definition registry | `T/readiness/results/stage3_definition_gene_registry_v0.1_20260723.tsv` | RNA `definition_registry`; drug 05 `definition`. Locally packaged unchanged (497 rows / four columns). Study-combined definition genes, consumer `gene_symbol`; source retains the 497-gene assertion. Attribution of component memberships remains to be recorded. |
| lncRNA annotation | `T/readiness/results/lncrna_grch38_113_annotation_registry_v0.2_20260723.tsv` | RNA 04/06 `lncrna_registry`. Locally packaged unchanged (16882 rows / eighteen columns). Prepared Ensembl/GRCh38.113 mapping and eligibility, not an expression matrix. External annotation provenance/acquisition needs documentation. |
| RPPA dictionary | `T/readiness/results/rppa_antibody_dictionary_frozen_v0.2_20260723.tsv` | RNA `antibody_dictionary`; drug 05 `rppa_dict`. Locally packaged unchanged (258 rows / twenty-four columns). Platform/antibody/gene/site mappings; consumer schemas differ. Review all required fields and their external sources. |
| TCPA pathway annotation | `T/resources/tcpa_official_pathway_annotation/akbani2014_supp_table13_official_pathway_members.tsv` | RNA `tcpa_table13`; drug 05 `tcpa`. Locally packaged unchanged (75 rows / three columns). Structured derivative of a third-party supplementary table, with source-specific correction logic. Acquisition and redistribution must be addressed separately from study-owned definitions. |
| Drug-association programmes | `T/results/09_mrna_c1c4_extreme_display/v0.2_20260730/mrna_c1c4_extreme_display_genes_v0.2_20260730.tsv` and `mrna_c1c4_extreme_cluster_registry_v0.2_20260730.tsv` in that directory | drug 05 `program_genes/program_clusters`; 06 `program_clusters`. The gene table is a result-derived membership/display interface; the registry provides programme labels/boundaries. RNA 16 can generate the gene-table interface but does not produce the cluster registry. Both exact adopted tables are locally packaged unchanged; Table 3 Programme_members is a different RNA 14 input/output interface. |
| Programme catalogue | `S/config/program_manifest_v0.1_20260722.tsv`; `S/02_program_catalog/program_catalog_frozen_v0.1.tsv` | The first is the builder's study selection framework; the second is now locally packaged unchanged for single_cell 12/14/15. Its 23 rows and 38 columns include marker/GMT source metadata, gene-set sizes and overlap summaries; underlying resource terms remain separate. |
| Programme memberships | `S/02_program_catalog/program_gene_membership.tsv.gz`; `S/02_program_catalog/program_nondefinition_gene_membership.tsv.gz` | Both are now locally packaged unchanged (2814/2533 rows, 16 columns each, 23 programmes). single_cell 14 uses FULL/NONDEFINITION; 15, RNA14 and ICB01/03/06 use NONDEFINITION. Gene-level outputs combine study markers, external sets and study-universe/definition flags; keep their distinct semantics. |
| Definition-exclusion references | `G/definition_exclusion_genes.tsv`; `G/definition_exclusion_membership.tsv.gz` | The gene table is locally supplied unchanged (497 rows, seven columns) for 06/08 exclusion_genes and 12 definition_genes; gene and exclusion flags annotate Acid/original/refined LM12 overlap, not count removal. The compressed membership is unbundled builder provenance and is not consumed by these public entries. |
| Melanoma fixed contract | `S/config/stage4_melanoma_program_validation_contract_v0.1_20260722.tsv` | single_cell 14 `contract`. Locally supplied unchanged: 21 rows, three columns, unique contract_item. SC14 checks contract_item/frozen_value only for the adopted contract_version; rationale and relative source_object are historical provenance, not patient observations. Its old C3 comparisons do not control current SC15 GSEA contrasts. |
| ICB symbol alias and fixed programme subset | `10.0paired_ici_dynamic_remodeling/00_readiness/gene_symbol_alias_manifest.tsv`; `13.0c1c4_biological_context_and_icb/results/A_frozen_program_members.tsv` | ICB preparation `gene_alias` and ICB 08 `frozen_program_members`. The alias is locally supplied unchanged (one row, four columns): alias_symbol/canonical_symbol map PKM2 to PKM; rule_scope/source are provenance, not a clinical crosswalk; the programme/gene table is produced by RNA 14 from NONDEFINITION memberships. Historical directory naming does not make the stopped dynamic-ICI analysis an active workflow. |
| Optional display references | `final_figures/shared_assets/frozen_pan_cancer_order_by_integrated_ssgsea.csv`; `final_figures/c1c4_atlas/shared_assets/frozen_c1c4_state_palette.tsv` | genomics 04 `cancer_order/state_palette`. Existing plot entry references, not new mandatory plotting work. The order is study-derived. |

Authority is the source path recorded by the existing [RNA/RPPA](../workflows/transcriptome_proteome/README.md), [drug](../workflows/drug_response/README.md), [single-cell](../workflows/single_cell/README.md) and [ICB](../workflows/icb/README.md) documentation and their source mappings, rather than a filename-date sort. In particular, drug programme references are explicit in `11.0c1c4_multiomics_immune_drug_bridge/scripts/03_build_stage1_shared_numeric_inputs_v0.1_20260810.R`; programme construction is explicit in `S/scripts/01_build_program_catalog.R`.

## Fig. 6 selection references locally packaged

The historical C/D references were recovered from `final_figures/c1c4_atlas/fig6_pancancer_scrna/archive/pre_r01_freeze_20260917/source_data/` and copied unchanged into [references](../references/README.md), with cmp confirmation. The C selection is also byte-identical to the adopted R01 `fig6_c_compact_functional_gene_selection_candidate.tsv`; candidate naming does not invalidate its frozen authority.

| Packaged file | Preserved interface | Example configuration |
|---|---|---|
| `fig6_c_compact_functional_gene_selection.tsv` | 36 unique genes, eight original columns; marker/compartment/programme/order identifiers | single_cell 13 `marker_panel` |
| `fig6_d_pathway_selection.tsv` | 15 unique IDs, native `comparison, pathway, pathway_display, pathway_order` columns | single_cell 12 `fig6_pathway_source` and 15 `fig6_pathways` |

The historical D selection and 60-row statistical source have identical ID sets, first-occurrence ID order and per-ID display labels. Therefore the original four-column selection serves both 12/15 without reordering, projection or R-code changes. The mixed GSEA statistical table was not copied. Display-name repeats (Hypoxia/EMT across different IDs) must not be collapsed. Original source mapping follows the archived freeze writer and manifest; The 23 September R01 D update changed display labels, not the statistical inputs. GMT acquisition, matched patient/count/model inputs and external distribution permissions remain separate decisions.

## External references and mixed results kept separate

Exact documented local sources include `../LM_pan_cancer_v1/9.0GDSC/drug_annotation_master.tsv` (legacy annotation), `../data/GDSC_data/compounds_annotation.csv` (official compound annotation), `../data/TCGA_Data/Pan_Cancer/xena_tcga/probeMap_gencode.v23.annotation.gene.probemap` and `../LM_pan_cancer_v1/data/cytoBand.txt.gz`. They exist, but are not automatically study-owned materials. Parent-relative paths are provenance locators only, not portable acquisition instructions.

MSigDB GMTs, local ImmPort GMT and `CellChatDB.human` remain resource-acquisition/version questions; no raw GMT or database payload is bundled here. The three locally packaged programme tables contain study-prepared marker/MSigDB-derived references and do not grant external redistribution rights. PRISM coverage, programme-state centres, programme scores, patient/cell annotations, clinical crosswalks and expression/count matrices are outside this focused definition-interface task.

## Scoped completion and later decisions

Use Supplementary Table 2 for mapped LM12 interfaces and the published nine-gene list for minimal TSV/RData preparation. D1/D3/D4, Fig. 6 C/D selections and the three adopted programme references are now locally packaged unchanged; do not bulk-copy the other candidates. External publication/licensing and other prepared-input channels remain author decisions in the [release register](input_access_and_permissions.md#release-decision-register).

No additional data upload or reference/prepared-input packaging is required to close the current code-organization stage. Missing inputs remain documented dependencies for execution; access, provider terms, publication and CD8 data submission are separate decisions. See the [input-access inventory](input_access_and_permissions.md) and [coverage map](figure_coverage.md). Nineteen unchanged reference copies and their paths across twenty-three distinct example configurations are now packaged; no workbook, analysis entrypoint or manuscript was changed. The package remains 90 entries across 15 modules; documented format routes are not complete input delivery or numerical reproduction.

## Programme references locally packaged

The unchanged program_catalog_frozen_v0.1.tsv, program_gene_membership.tsv.gz and program_nondefinition_gene_membership.tsv.gz copies in references/ originate from S/02_program_catalog/. All 38 catalogue columns and both 16-column memberships, order, missing values, Chinese interpretation/claim text and compressed bytes are retained. The bounded read-only check found 23 programmes in each table, no duplicate membership program-gene keys or missing required IDs, exact NONDEFINITION correspondence to FULL where gene_in_pseudobulk & !overlaps_any_definition, and matching catalogue member counts/overlap fractions. Field/value review found no patient-level identifiers, private absolute paths or credentials. Seven DCFs now connect ten keys; other prepared prerequisites remain independent.

The unbundled builder resolves marker members after gene_present filtering and canonical members from exact-version GMTs, then adds definition/universe flags to FULL. Only NONDEFINITION applies the indicated subset. This packaging reuses adopted intermediates rather than rerunning that builder. Table 6 Melanoma_gene_membership is SC15's canonical-plus-project fixed GSEA library, not a substitute for these three input interfaces. Local inclusion is authorized; third-party provenance/terms and external distribution remain separate decisions.

## Bulk, immune-library and drug-programme interfaces locally packaged

Eight exact adopted tables are now copied unchanged and connected through twenty-three keys in twelve DCFs. The gene/cluster drug tables retain mixed statistical and definition fields. Detailed actual schemas, sources and preservation boundaries are in [references](../references/README.md). The 119-set long table retains one NA membership: ImmuCellAI Cytotoxic has fourteen valid symbols plus one NA, while source_gene_count=15 records vector length. A read-only check on 2026-10-06 reconstructed rank names from the formal model tables for E1/E2/E3, covering FULL_RANK and NON_DEFINITION_RANK (120 nonempty ranks across 22 cancer inputs). No rank names were missing, empty or the literal string NA, and the sole missing member was excluded before pathway construction. This closes the specific NA-input concern for those historical runs, not numerical reproduction. Source code still lacks an explicit valid-name guard for arbitrary future inputs; the reference file is unchanged.

TCPA mapping remains an unresolved historical branch difference: RNA 03 maps the retained source token HER3PY1298 to HER3PY1289, whereas drug 05 excludes HER3PY1298 and uses six of the seven original RTK members. The historical "printing error" judgment has not been independently confirmed externally; NOT_MEASURED_IN_CURRENT_RPPA_PANEL is a legacy script label, not verified platform absence. RNA 03 applies the alias after single-feature modelling and panel-wide BH, as a pathway annotation for Fig. 3D; it does not compute a seven-member RTK score. Drug 05 computes the six-member RTK score. Fig. 5C has no displayed RTK row, but its first five RPPA relationships share the complete within-cancer RPPA_PATHWAY BH family with RTK (up to nine programmes x ten pathways, restricted to estimated tests). Changing RTK alone would not directly change those five relationships' rho or raw P values, but could change their FDR, significance markers and significant-cancer counts. The CNV and six immune rows use separate BH families. A 2026-10-06 R 4.5.1 local sensitivity comparison added HER3_pY1289 with positive sign to the historical six-member RTK score. The historical baseline passed numerical checks. Seven RTK tests and one non-RTK test crossed FDR < 0.05; RTK significant tests changed from 13 to 10 and non-RTK from 136 to 135. Fig. 5C's first five rows retained all 100 evaluable rho values, significance markers and significant-cancer counts (11/9/8/5/5), although FDR values changed. RPPA orthogonal-candidate membership was unchanged. This is a conditional mapping scenario, not an adopted correction or full package reproduction; TrackB/MOFA effects were not evaluated and formal assets remain unchanged.
