# Software requirements

The interface and syntax checks use R 4.5.1. Required analysis packages are indexed per implementation in [software_dependencies.tsv](software_dependencies.tsv), extracted from library/require and namespace references in source. Install the packages needed for the selected route in a separate R environment; the repository does not vendor standard method implementations.

Key analysis dependencies include survival, limma, edgeR, fgsea, AUCell, Seurat-related interfaces, data.table and the plotting/file-reader packages used by the corresponding scripts. R package versions should follow the paper's recorded Methods and original source environment where available. This repository does not provide a recovered historical lockfile or claim that all listed packages were installed or executed in this update. Dynamic or external-system dependencies require inspection of the selected script.

PRISM concordance calls the external `sha256sum` utility to verify its fixed provider snapshots. Use a system providing that executable. Historical random seeds, statistical families and source thresholds are retained in the scientific code. The standard-method packages themselves need not be copied into this repository; paper-specific selection and analysis code remains included.
