# Purpose: Shared statistics helpers; exact matching implementations only.
# Helpers are rebound to the calling stage scope to retain original lexical lookup.
# 仅合并实现完全相同的函数；作用域保留在各独立阶段。

shared_seed_from_keys_001 <- function(...) {
    key <- paste(..., collapse = "|")
    as.integer((20260806 + sum(utf8ToInt(key)))%%.Machine$integer.max)
}

shared_safe_ci_030 <- function(r, n) if (n > 0L) binom.test(r, n)$conf.int else c(NA_real_, NA_real_)

shared_extract_scoring_matrix_035 <- function(seurat_obj, assay_name) {
    assay_obj <- seurat_obj@assays[[assay_name]]
    if (is.null(assay_obj)) {
        stop("Assay not found: ", assay_name)
    }
    counts_mat <- assay_obj@counts
    if (!is.null(counts_mat) && nrow(counts_mat) > 0L && ncol(counts_mat) > 0L) {
        return(list(matrix = counts_mat, slot_used = "counts"))
    }
    data_mat <- assay_obj@data
    if (!is.null(data_mat) && nrow(data_mat) > 0L && ncol(data_mat) > 0L) {
        return(list(matrix = data_mat, slot_used = "data"))
    }
    stop("No usable expression matrix found in assay: ", assay_name)
}

shared_norm_sample_037 <- function(x) substr(gsub("\\.", "-", as.character(x)), 1L, 15L)

shared_norm_symbol_039 <- function(x) toupper(trimws(as.character(x)))

shared_wilcox_summary_056 <- function(x, y) {
    test <- suppressWarnings(wilcox.test(x, y, alternative = "two.sided", exact = FALSE, conf.int = TRUE, conf.level = 0.95))
    combined <- c(x, y)
    ranks <- rank(combined, ties.method = "average")
    n_x <- length(x)
    n_y <- length(y)
    u_x <- sum(ranks[seq_len(n_x)]) - n_x * (n_x + 1)/2
    cliff_delta <- 2 * u_x/(n_x * n_y) - 1
    c(p_value = unname(test$p.value), hl_shift_c1_minus_c4 = if (length(test$estimate)) unname(test$estimate) else NA_real_, hl_ci_low = if (!is.null(test$conf.int)) unname(test$conf.int[1L]) else NA_real_, hl_ci_high = if (!is.null(test$conf.int)) unname(test$conf.int[2L]) else NA_real_, cliff_delta_c1_minus_c4 = cliff_delta)
}
