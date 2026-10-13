# Interface verification

Environment: author-selected remote R4.5.1. Date: 2026-10-13.

- 90 retained implementations are byte-identical to the public 7fc95d5 code snapshot; see source_copy_check.tsv.
- 36 task configurations cover 91 execution records and 90 distinct implementations. E1 and E3 require two separate records for the same implementation.
- 61 interface checks passed: dispatcher parsing/help/list, all 36 task plans, invalid CLI/configuration rejection, genuine CD8 metadata-only checks with a fake non-RData file, and synthetic process-order/link/failure/output-protection checks.
- Spaces, apostrophes and DCF-folded paths are supported. Child configs use absolute paths without line wrapping.
- Synthetic failures stop later steps and never receive a COMPLETED event; nonempty outputs remain unchanged.

See interface_checks.tsv for results. Local test paths are replaced with descriptive directory tokens. Fixtures and development scripts remain in the project temporary workspace, outside this candidate.

No study object, expression/count matrix, FASTQ, statistical model, GSEA body or frozen result was loaded or recomputed. Structural checks do not establish real schemas, output-value equivalence, all analysis dependencies or end-to-end reproduction.

The nineteen reference tables are not bundled; only their expected supplementary-material filenames are listed.
