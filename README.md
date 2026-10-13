# C1–C4 pan-cancer atlas

Core analysis code for *Joint states of lactate-associated transcription and inflammatory remodelling across human cancers*.

The repository contains 90 R entries in 15 modules, with one example DCF per entry. Start with the research question below, then follow the corresponding route in the [workflow guide](docs/workflow_guide.md). Entry numbers are local to each module.

## Start here

| Research question | Workflow route |
|---|---|
| Define the acid signature, LM12 and joint states | [Signature and states](docs/workflow_guide.md#signature-and-states) |
| Examine external bulk, survival and immune ecology | [Bulk validation and clinical context](docs/workflow_guide.md#bulk-validation-and-clinical-context) |
| Compare RNA, protein and genomic layers | [Multiomics](docs/workflow_guide.md#multiomics) |
| Examine predicted drug vulnerability and programme associations | [Predicted drug response](docs/workflow_guide.md#predicted-drug-response) |
| Resolve cell compartments and independent melanoma evidence | [Cellular context](docs/workflow_guide.md#cellular-context) |
| Examine immunotherapy associations | [Immunotherapy context](docs/workflow_guide.md#immunotherapy-context) |
| Examine experimental CD8 time responses | [CD8 time course](docs/workflow_guide.md#cd8-time-course) |

## Use a workflow

1. Read its module README and obtain the exact matching inputs.
2. Copy its example DCF to a `*.local.dcf` file and replace the placeholders. Relative paths resolve from the DCF directory.
3. Review the entry's `--help` and documented configuration check. Use a separate absent or empty output directory for each stage.

For example:

```sh
cp config/drug_shared_numeric.example.dcf config/drug_shared_numeric.local.dcf
# Edit the local DCF: provide every required input and an independent output_dir.
Rscript workflows/drug_response/05_build_shared_numeric_inputs.R --help
Rscript workflows/drug_response/05_build_shared_numeric_inputs.R --config config/drug_shared_numeric.local.dcf --check-config
```

For this entry, `--check-config` checks file/path metadata and creates no analysis outputs. Ordinary execution requires the documented inputs and dependencies. Consult the chosen entry; this example does not define a universal CLI.

## Inputs and verification

This is an input-dependent code release. Study-prepared expression/counts, patient/cell objects and the nineteen reference-table payloads are not supplied. Expected DCF filenames are not download channels. Use the [input inventory](docs/input_access_and_permissions.md), [reference requirements](references/README.md) and [supplementary-input mapping](docs/reference_file_candidates.md#supplementary-material-first).

Existing checks cover limited syntax, static and synthetic interfaces; end-to-end real-data reproduction and complete runtime compatibility are unverified. The [environment guide](environment/README.md) records dependencies and evidence. R 4.5.1 is the packaging/check environment, not every historical generating environment.

## Further information

- [Workflow guide and all 90 entry roles](docs/workflow_guide.md)
- [Result/output map for Fig. 1–7 and Supplementary Tables 1–7](docs/result_output_mapping.md)
- [Implementation notes and reference boundaries](references/README.md)
- [Release scope](docs/release_scope.md)

Software and accompanying documentation are licensed under [MIT](LICENSE). Third-party resources and study data retain their separate terms.
