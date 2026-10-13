# Release scope

Repository: https://github.com/liedog123/LM_pan_cancer_C1-C4_atlas
Prepared: 13 October 2026.

The public release contains 90 unchanged R entry points in 15 modules, example DCF configurations, workflow/source-mapping documentation and environment guidance. Software and accompanying documentation are provided under the MIT license in the root LICENSE file.

Input availability is separate from code availability. The nineteen reference filenames catalogued in [the reference inventory](../references/README.md) describe the exact inputs used during internal preparation; those table payloads are not distributed in this release. DCF paths identify expected locations and do not establish that a file is provided. Obtain matching inputs through documented providers, published supplementary materials when available, or a future explicit study-data release. No unspecified "upon request" route is promised.

Raw FASTQ, counts, sample annotations, patient/cell objects, manuscript files, frozen figures and statistical results are outside this code release. GEO preparation remains a separate workstream with no accession at this release date. Third-party datasets, gene sets and runtime resources remain subject to their providers' terms; this software license does not license those resources.

The verification level remains the prior limited static/syntax/interface evidence. Publication does not establish end-to-end real-data reproduction or supply missing study-prepared inputs. Historical implementation boundaries, including the six-member RTK calculation and its documented sensitivity assessment, are retained.

SOURCE_MANIFEST.tsv identifies each unchanged copied file and each release-adjusted documentation file by source-relative path and SHA-256. No original analysis script, research data or frozen result was modified.
