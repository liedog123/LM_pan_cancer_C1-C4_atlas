# C1–C4 pan-cancer atlas

## Purpose and scope

This public repository contains **90 R entry points in 15 workflow modules**, example DCF configurations and source/environment documentation. It supports scoped continuations from matching study-prepared inputs. Complete statistical families and historical implementations are retained. See [release scope](docs/release_scope.md) for the distribution boundary.

Verification is limited to the recorded static/syntax/interface checks. Real-data end-to-end reproduction and complete input availability have not been established. See the [result/output map](docs/result_output_mapping.md) for Fig. 1–7 / Supplementary Tables 1–7, complete-result companions and prepared-export boundaries, and the [coverage map](docs/figure_coverage.md) for included computation.

## Choose a module

Start with the relevant workflow README rather than running all entries in numerical order. Entry numbers are local to each module.

| Module | Entries | Starting point / role |
|---|---:|---|
| [Signature selection](workflows/signature_selection/README.md) | 3 | Prepared candidates/DEG/scores/TPM; selection and separate final-nine correlation families |
| [Acid/AG preprocessing](workflows/acid_ag_preprocessing/README.md) | 6 | Prepared expression/clinical inputs; scoring, AG groups, exports and separate integrated/validation branches |
| [TCGA discovery](workflows/tcga_discovery/README.md) | 2 | Processed High/Low inputs; full DEG rankings and per-cancer ImmPort GSEA |
| [LM12 construction](workflows/lm12_construction/README.md) | 2 | Formal discovery GSEA plus curated candidates/families; recurrence and adopted definitions |
| [Shared definitions](workflows/shared_definitions/README.md) | 3 | Processed inputs and frozen LM12; IG/C1–C4 labels, validation and independent canonical-format export |
| [External bulk](workflows/external_bulk/README.md) | 2 | Prepared GEO expression/metadata; joint scores and gene-level GSEA, not external C1–C4 assignment |
| [AG survival](workflows/ag_survival/README.md) | 4 | Prepared High/Low clinical inputs; KM, Cox, complete-family BH and formal table export |
| [Clinical outcomes](workflows/clinical_outcomes/README.md) | 1 | Prepared patient inputs; Stage 1 clinical associations, Cox, pooling and diagnostics |
| [Immune ecology](workflows/immune_ecology/README.md) | 2 | Prepared patient/deconvolution inputs; Stage 2 comparisons, sensitivities and evidence finalization |
| [Transcriptome/proteome](workflows/transcriptome_proteome/README.md) | 16 | Prepared mRNA/lncRNA/RPPA/readiness inputs; models, recurrence, ORA/GSEA and cross-layer/subtype continuations |
| [Genomics](workflows/genomics/README.md) | 5 | Accepted callable/readiness and event/CNV matrices; mutation, TMB/MSI, burden and gene-to-region branches |
| [Drug response](workflows/drug_response/README.md) | 6 | Corrected predictions or matching prepared support/bundles; C1/C4, OOF, PRISM, shared-PC1 and programme associations |
| [Single cell](workflows/single_cell/README.md) | 18 | Prepared Seurat/audit/G1 inputs; scoring, pseudobulk, Bassez sensitivity, melanoma and CellChat branches |
| [ICB](workflows/icb/README.md) | 8 | Harmonized cohort inputs; fixed Riaz/Gide/Braun analyses and separate post-hoc Liu/Gide branches |
| [CD8 perturbation](workflows/cd8_perturbation/README.md) | 12 | Matched counts/design plus required references/bridge inputs; QC, DE, GSEA and source exports |

## Configuration and execution

Read the chosen workflow's README and source mapping, select its matching example DCF in `config/`, and obtain the exact required inputs. Relative input paths are resolved from the DCF directory, not the shell's working directory. The `../references/` paths identify expected filenames; obtain matching reference inputs separately and place them at those locations or configure another matching location. Placeholder paths are not delivered data.

For example, from this workspace root, drug-response 05 supports the following actual interface:

```sh
cp config/drug_shared_numeric.example.dcf config/drug_shared_numeric.local.dcf
# Edit the local copy: supply every required input and choose a new output_dir.
Rscript workflows/drug_response/05_build_shared_numeric_inputs.R --help
Rscript workflows/drug_response/05_build_shared_numeric_inputs.R --config config/drug_shared_numeric.local.dcf --check-config
```

For this entry, `--check-config` checks configuration/file metadata and output safety, not research contents or numerical correctness. Removing `--check-config` runs the analysis and requires the documented inputs/dependencies; that run has not been performed as part of packaging. Consult each entry's help and workflow README: this example is not a universal CLI for all 90 scripts.

Use a new, absent or empty **independent output directory for each stage**, outside inputs, the config directory and this code workspace. Do not continue in a nonempty historical result directory. Matching adopted intermediates may allow downstream continuation without rerunning upstream stages; they do not automatically replace frozen definitions or canonical labels.

### Important branch distinctions

Detailed input/output names and conditional order are in the [result map](docs/result_output_mapping.md) and module READMEs.

- Acid/AG 05 uses integrated expression and existing AG groups; 06 independently consumes historical `full_validation_data`, not 05's output. Its historical joint scoring background is not current LM12.
- Transcriptome/proteome 04 → 03 supplies RPPA models to recurrence/TCPA; lncRNA continues 04 → 05 → 06 → 07 → 08. mRNA 01 → 02 → 16 supplies 13's representative input, but 13's other layers remain separately required. For boundary routes, 09 supplies model inputs to 10 (E2) and two separate runs of 11 (E1/E3), then 12 combines those outputs; subtype context is separate 14 → 15. C1/C4 is formal and OVR is historical support.
- Drug 01 starts after prediction/GLDS correction; 02 supplies OOF support and 03 can instead consume matching existing support tables. Shared-PC1 04 is separate from 05 → 06 programme associations; matching numerical bundles can feed 06 directly. Predictions/associations are not clinical efficacy.
- Genomics 03, 04 and 05 are independent branches; 05 consumes prepared CNV inputs, not their statistical outputs.
- Single-cell scoring follows 04 → 05 → (01/02) → 03; prepared G1 contracts bridge into 06–11. Bassez 12/13, melanoma 14 → 15 and CellChat 16 → 17 → 18 are distinct branches. 18 imports 17's rep1 objects and computes rep2–10, retaining complete forty-network summaries. Melanoma outputs also support CD8 03; pancancer labels are not melanoma substitutes. Current manuscript analyses do not use pySCENIC.
- ICB follows (01 → 02) + (03 → 04) → 05 for the fixed three-cohort family. Liu 06 → 07 and Gide 08 are separate post-hoc families, never pooled into the fixed BH family.
- CD8 01 QC and 02 DE independently start from the same counts/design; review 01 before 02, then continue 02 → 03 → 04 with separate outputs. Optional 10/11 need independently prepared trajectory/audited-metadata inputs. Their absence does not invalidate the counts-only interface.

## Inputs: included, external and not supplied

The [reference inventory](references/README.md) documents nineteen internally prepared reference files and their consuming configuration fields. Their payloads are not included in this public release; provenance and external acquisition/distribution decisions remain explicit. Example DCFs retain the expected filenames.

The [input-access inventory](docs/input_access_and_permissions.md) distinguishes documented public providers/accessions, conditional packaged generation routes and study-prepared materials without an established delivery route. Generic TCGA/GEO/GDSC/PRISM downloads do not replace matching study derivatives. Use [Supplementary Table 2's mapped worksheets and the published nine-gene definitions](docs/reference_file_candidates.md#supplementary-material-first) where applicable; the actual supplementary release channel remains to be established.

Patient labels, expression/count matrices, processed patient/cell objects, raw GMT payloads, trained models, third-party source code and frozen figure assets are **not supplied**. Upstream raw-data/QC/metadata preparation, initial tumour–normal DEG/tumour-program generation, original oncoPredict prediction/training, GLDS correction and MOFA training are outside these scoped entries. Unsupplied prepared inputs remain explicit dependencies, not a requirement to add data before this code-organization stage can close. CD8 data submission and other access/release decisions are separate workstreams.

### Known reference boundaries

The [TCPA and 119-set NA notes](references/README.md) are required reading before using those references. RNA 03 applies HER3PY1298 → HER3PY1289 while drug 05 excludes the retained source token; historical "printing error" judgment is not independently externally confirmed, and the legacy NOT_MEASURED label does not establish platform absence. RNA uses the alias for post-statistical Fig. 3D annotation, not a seven-member pathway score; drug RTK scoring uses six of seven original members. Fig. 5C does not display RTK, but its first five RPPA relationships share the complete within-cancer RPPA_PATHWAY BH family. Their rho/raw P values would not change directly after an RTK-only mapping change, whereas FDR, significance markers and significant-cancer counts could change. The CNV and six immune rows use independent BH families. A 2026-10-06 R 4.5.1 local sensitivity comparison added HER3_pY1289 with positive sign to the historical six-member RTK score. The historical baseline passed numerical checks. Seven RTK tests and one non-RTK test crossed FDR < 0.05; RTK significant tests changed from 13 to 10 and non-RTK from 136 to 135. Fig. 5C's first five rows retained all 100 evaluable rho values, significance markers and significant-cancer counts (11/9/8/5/5), although FDR values changed. RPPA orthogonal-candidate membership was unchanged. This is a conditional mapping scenario, not an adopted correction or full package reproduction; TrackB/MOFA effects were not evaluated and formal assets remain unchanged.

The 119-set library retains one NA: ImmuCellAI Cytotoxic has 14 valid symbols plus NA, and source_gene_count=15 is vector length. A read-only check on 2026-10-06 reconstructed rank names from the formal model tables for E1/E2/E3, covering FULL_RANK and NON_DEFINITION_RANK (120 nonempty ranks across 22 cancer inputs). No rank names were missing, empty or the literal string NA, and the sole missing member was excluded before pathway construction. This closes the specific NA-input concern for those historical runs, not numerical reproduction. Source code still lacks an explicit valid-name guard for arbitrary future inputs; the reference file is unchanged. RNA 16 can generate the drug gene-table interface but not the cluster registry; Table 3 Programme_members is a different RNA 14 interface, not a substitute.

## Environment and verification level

See the [environment/dependency guide](environment/README.md) for module requirements and evidence levels. R 4.5.1 is the current packaging/check environment, not a retrospective claim about every historical run. Dependencies and external runtime resources remain module-specific; complete installation and ordinary-run compatibility are unverified.

Existing evidence consists of syntax/help checks, synthetic metadata/CLI safety cases and static retained-source comparisons, with detailed adaptations and check records in each workflow's `source_mapping.md` and the environment guide. Stage 5 additionally passed formal-path file-metadata checks without loading research rows or predictions. Reference-copy comparisons and configured-path checks passed. Those records are reused, not rerun during README closeout.

The closeout check is limited to current reader-document relative links, configured `../references/` targets, obvious private-path/credential patterns in package text and file-type inventory. It is not a real-data schema check, dependency execution, model/statistical run, figure redraw or numerical reproduction.

## Provenance and publication status

[Inputs and provenance](docs/inputs_and_provenance.md), module source mappings and the [result map](docs/result_output_mapping.md) trace source logic and adopted outputs. Source-project paths in those records are provenance locators, not resources promised inside a standalone release; interface adaptation is not a claim of byte-identical whole original scripts.

This code repository is public at https://github.com/liedog123/LM_pan_cancer_C1-C4_atlas. Software and accompanying documentation are licensed under [MIT](LICENSE). Third-party resources and study data are governed separately; this repository supplies no reference-table payloads or study-data accession. See [release scope](docs/release_scope.md) and the [input inventory](docs/input_access_and_permissions.md).
