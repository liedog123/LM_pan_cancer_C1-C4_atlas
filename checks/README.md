# Code verification

Environment: R 4.5.1. Date: 2026-10-19.

- 35 task configurations, 87 execution records and 86 distinct implementation scripts; E1 and E3 share one implementation.
- 315 checks passed in [interface_checks.tsv](interface_checks.tsv): all implementation/entry syntax, all implementation headers and help modes, all 35 task plans, task/index coverage, obsolete-task rejection, the three extracted implementations metadata-only checks, invalid fields, input/output overlap and nonempty-output protection.
- [source_scope_check.tsv](source_scope_check.tsv) records 83 unchanged scientific bodies after stripping only their opening comment blocks, plus three scoped extractions.
- Ten [retained-logic checks](retained_logic_check.tsv) compare the extracted source blocks: the original PRISM five-drug comparison and BH family, shared-axis 18-cancer/three-drug computations, and survival group selection and serialized object.

Repeat interface checks with `Rscript checks/run_checks.R /path/to/new/report.tsv`. The report path must be new and outside this checkout. Synthetic metadata fixtures are created in a temporary directory and removed. No analysis packages or study objects are loaded; the checks do not execute fitting, enrichment, scientific statistical bodies or plotting.

The source comparison baseline is public-code commit c0790c516fd18a0f859fdad59ae961a7d51f0a75. Historical check records remain in Git history and the author internal archive. This update does not establish real study-input schema compatibility, numerical equality to frozen results or end-to-end reproduction. The starting-input and upstream-processing boundaries are stated in the [reproduction guide](../docs/reproduction_routes.md).
