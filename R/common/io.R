# Purpose: Shared io helpers; exact matching implementations only.
# Helpers are rebound to the calling stage scope to retain original lexical lookup.
# 仅合并实现完全相同的函数；作用域保留在各独立阶段。

shared_firth_002 <- function(X, y, index) {
    X <- as.matrix(X)
    storage.mode(X) <- "double"
    y <- as.numeric(y)
    if (qr(X)$rank < ncol(X))
        stop("Firth design is rank deficient")
    ll <- function(b) {
        eta <- drop(X %*% b)
        pr <- plogis(eta)
        w <- pmax(pr * (1 - pr), 1e-15)
        inf <- crossprod(X * sqrt(w))
        ld <- determinant(inf, logarithm = TRUE)
        if (ld$sign <= 0 || !is.finite(ld$modulus))
            return(-1e+50)
        sum(y * eta - (pmax(eta, 0) + log1p(exp(-abs(eta))))) + 0.5 * as.numeric(ld$modulus)
    }
    full <- optim(rep(0, ncol(X)), function(b) -ll(b), method = "BFGS", control = list(maxit = 1000, reltol = 1e-10))
    if (full$convergence != 0)
        stop("Firth full optimization failed")
    bh <- full$par
    top <- ll(bh)
    profile_ll <- function(v) {
        rem <- setdiff(seq_along(bh), index)
        fun <- function(z) {
            b <- bh
            b[index] <- v
            b[rem] <- z
            -ll(b)
        }
        op <- optim(bh[rem], fun, method = "BFGS", control = list(maxit = 1000, reltol = 1e-09))
        if (op$convergence != 0)
            stop("Firth profile optimization failed")
        -op$value
    }
    target <- top - qchisq(0.95, 1)/2
    bound <- function(side) {
        step <- 0.5
        b <- bh[index] + side * step
        while (profile_ll(b) > target && step < 32) {
            step <- step * 2
            b <- bh[index] + side * step
        }
        if (profile_ll(b) > target)
            return(if (side < 0) -Inf else Inf)
        if (side < 0)
            uniroot(function(z) profile_ll(z) - target, c(b, bh[index]), tol = 1e-05)$root
        else uniroot(function(z) profile_ll(z) - target, c(bh[index], b), tol = 1e-05)$root
    }
    p <- pchisq(2 * (top - profile_ll(0)), df = 1, lower.tail = FALSE)
    c(OR = exp(bh[index]), lo = exp(bound(-1)), hi = exp(bound(1)), p = p)
}

shared_stop_if_003 <- function(cond, msg) if (isTRUE(cond)) stop(msg, call. = FALSE)

shared_stop_if_004 <- function(condition, message) {
    if (isTRUE(condition))
        stop(message, call. = FALSE)
}

shared_select_one_per_patient_005 <- function(dat) {
    dd <- copy(dat)
    dd[, `:=`(pre_priority, fifelse(pre_post == "Pre", 0L, 1L))]
    setorder(dd, patient_key, pre_priority, -n_cells, -library_size, sample_unit_key)
    dd[, .SD[1L], by = patient_key]
}

shared_write_tsv_006 <- function(df, path, compress = FALSE) {
    con <- if (compress)
        gzfile(path, open = "wt")
    else file(path, open = "wt")
    on.exit(close(con), add = TRUE)
    write.table(df, file = con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

shared_qa_add_007 <- function(id, passed, note) {
    qa_rows[[length(qa_rows) + 1L]] <<- data.table(check_id = id, passed = isTRUE(passed), note = note)
}

shared_align_matrix_to_patient_017 <- function(mat, patient_cancer) {
    sample_id <- norm_sample(colnames(mat))
    keep_first <- !duplicated(sample_id)
    mat <- mat[, keep_first, drop = FALSE]
    colnames(mat) <- sample_id[keep_first]
    common <- intersect(patient_cancer$Sample, colnames(mat))
    mat <- mat[, common, drop = FALSE]
    meta <- patient_cancer[match(common, Sample)]
    if (!identical(meta$Sample, colnames(mat)))
        stop("\346\240\267\346\234\254\345\257\271\351\275\220\345\244\261\350\264\245")
    list(mat = mat, meta = meta)
}

shared_md5_or_na_020 <- function(path) {
    if (!file.exists(path) || dir.exists(path))
        return(NA_character_)
    unname(tools::md5sum(path))
}

shared_assert_nonempty_021 <- function(path) {
    if (!file.exists(path)) {
        stop("Missing input: ", path)
    }
    size <- file.info(path)$size
    if (is.na(size) || size <= 0) {
        stop("Empty input: ", path)
    }
    invisible(size)
}

shared_sha256_022 <- function(path) {
    out <- system2("sha256sum", shQuote(path), stdout = TRUE)
    strsplit(out[[1]], "[[:space:]]+")[[1]][1]
}

shared_read_gmt_023 <- function(path, collection, evidence_role) {
    lines <- readLines(path, warn = FALSE)
    parts <- strsplit(lines[nzchar(trimws(lines))], "\t", fixed = TRUE)
    if (any(lengths(parts) < 3L))
        fail(paste(collection, "GMT\345\255\230\345\234\250\345\260\221\344\272\2163\345\210\227\347\232\204\350\256\260\345\275\225"))
    registry <- rbindlist(lapply(seq_along(parts), function(i) data.table(Set_ID = paste0(collection, "::", trimws(parts[[i]][1])), Display_Name = trimws(parts[[i]][1]), Collection = collection, Functional_Domain = "GENERAL_CONTEXT", Evidence_Role = evidence_role, Source_Gene_Count = length(unique(norm_symbol(parts[[i]][-c(1, 2)]))))))
    membership <- rbindlist(lapply(seq_along(parts), function(i) {
        genes <- unique(norm_symbol(parts[[i]][-c(1, 2)]))
        genes <- genes[nzchar(genes)]
        data.table(Set_ID = paste0(collection, "::", trimws(parts[[i]][1])), Gene = genes)
    }))
    list(registry = registry, membership = membership)
}

shared_load_object_024 <- function(path, expected) {
    e <- new.env(parent = emptyenv())
    loaded <- load(path, envir = e)
    if (!expected %in% loaded)
        stop("Missing object ", expected, " in ", path)
    e[[expected]]
}

shared_load_named_025 <- function(path, name) {
    env <- new.env(parent = emptyenv())
    load(path, envir = env)
    if (!exists(name, envir = env, inherits = FALSE))
        stop("Missing named object: ", name)
    get(name, envir = env, inherits = FALSE)
}

shared_load_named_object_026 <- function(path, object_name) {
    input_env <- new.env(parent = emptyenv())
    load(path, envir = input_env)
    if (!exists(object_name, envir = input_env, inherits = FALSE)) {
        stop("Missing input object: ", object_name, call. = FALSE)
    }
    get(object_name, envir = input_env, inherits = FALSE)
}

shared_load_named_object_027 <- function(path, object_name) {
    input_env <- new.env(parent = emptyenv())
    load(path, envir = input_env)
    if (!exists(object_name, envir = input_env, inherits = FALSE)) {
        stop("Required object not found: ", object_name, " in ", path, call. = FALSE)
    }
    get(object_name, envir = input_env, inherits = FALSE)
}

shared_get_tcga_data_tpm_029 <- function(project) {
    project_dir <- file.path(tpm_root, project)
    if (!dir.exists(project_dir)) {
        stop("\351\241\271\347\233\256\346\226\207\344\273\266\345\244\271\344\270\215\345\255\230\345\234\250\357\274\232", project_dir)
    }
    message("\345\274\200\345\247\213\345\244\204\347\220\206\351\241\271\347\233\256 ", project)
    tpm_file <- file.path(project_dir, paste0(project, "_mrna_expr_tpm.rdata"))
    if (file.exists(tpm_file)) {
        mrna_expr_tpm <- load_named_object(tpm_file, "mrna_expr_tpm")
    }
    else {
        stop("TPM \350\241\250\350\276\276\346\225\260\346\215\256\346\226\207\344\273\266\344\270\215\345\255\230\345\234\250\357\274\232", tpm_file)
    }
    return(mrna_expr_tpm)
}

shared_relocate_input_031 <- function(recorded, input_dir) {
    if (length(recorded) != 1L || is.na(recorded) || !nzchar(recorded))
        stop("Invalid recorded input path.")
    leaf <- basename(recorded)
    if (leaf %in% c(".", "..") || !nzchar(leaf))
        stop("Invalid input basename.")
    candidate <- resolve_path(file.path(input_dir, leaf), dirname(config_path))
    if (!is_within(candidate, input_dir) || identical(candidate, input_dir))
        stop("Input escapes configured directory.")
    if (!file_test("-f", candidate) || file.access(candidate, 4L) != 0L || is.na(file.info(candidate)$size) || file.info(candidate)$size <= 0)
        stop("Selected input must be a readable nonempty regular file.")
    candidate
}

shared_relocate_input_032 <- function(recorded, input_dir, allow_missing = FALSE) {
    if (length(recorded) != 1L)
        stop("Invalid recorded input path.")
    if (is.na(recorded)) {
        if (isTRUE(allow_missing))
            return(NA_character_)
        stop("Invalid recorded input path.")
    }
    if (!nzchar(recorded))
        stop("Invalid recorded input path.")
    leaf <- basename(recorded)
    if (leaf %in% c(".", "..") || !nzchar(leaf))
        stop("Invalid input basename.")
    candidate <- resolve_path(file.path(input_dir, leaf), dirname(config_path))
    if (!is_within(candidate, input_dir) || identical(candidate, input_dir))
        stop("Input escapes configured directory.")
    if (isTRUE(allow_missing) && !file.exists(candidate) && !dir.exists(candidate))
        return(candidate)
    if (!file_test("-f", candidate) || file.access(candidate, 4L) != 0L || is.na(file.info(candidate)$size) || file.info(candidate)$size <= 0)
        stop("Selected input must be a readable nonempty regular file.")
    candidate
}

shared_summarize_cross_cancer_033 <- function(results, feature_cols, sig_col = "significant", effect_col = "effect_c1_minus_c4") {
    x <- copy(results)
    x[, `:=`(direction, fifelse(get(effect_col) > 0, "C1_higher", fifelse(get(effect_col) < 0, "C4_higher", "tie")))]
    x[, `:=`(sig_direction, fifelse(get(sig_col) & direction == "C1_higher", "C1_higher", fifelse(get(sig_col) & direction == "C4_higher", "C4_higher", "not_significant")))]
    out <- x[, .(n_testable_cancers = uniqueN(Cancer), n_sig_total = sum(get(sig_col), na.rm = TRUE), n_sig_C1_higher = sum(sig_direction == "C1_higher", na.rm = TRUE), n_sig_C4_higher = sum(sig_direction == "C4_higher", na.rm = TRUE), median_effect_c1_minus_c4 = median(get(effect_col), na.rm = TRUE), sig_cancers_C1_higher = paste(sort(Cancer[sig_direction == "C1_higher"]), collapse = ";"), sig_cancers_C4_higher = paste(sort(Cancer[sig_direction == "C4_higher"]), collapse = ";")), by = feature_cols]
    out[, `:=`(dominant_n_sig, pmax(n_sig_C1_higher, n_sig_C4_higher))]
    out[, `:=`(dominant_direction, fifelse(n_sig_C1_higher > n_sig_C4_higher, "C1_higher", fifelse(n_sig_C4_higher > n_sig_C1_higher, "C4_higher", "tie")))]
    out[, `:=`(context_dependent, n_sig_C1_higher > 0L & n_sig_C4_higher > 0L)]
    out[, `:=`(recurrent_ge3, dominant_n_sig >= 3L)]
    out
}

shared_read_col_034 <- function(rng) readxl::read_xlsx(source_clin, sheet = "Supplemental Table 1", range = rng, col_names = FALSE)

shared_safe_name_036 <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

shared_clean_numeric_038 <- function(x) suppressWarnings(as.numeric(as.character(x)))

shared_parse_filter_values_040 <- function(x) {
    vals <- safe_char(unlist(strsplit(x, ",")))
    vals[vals != "" & vals != "NA"]
}

shared_safe_filename_041 <- function(x) {
    x %>% str_replace_all("[^A-Za-z0-9_.-]+", "_") %>% str_replace_all("_+", "_") %>% str_replace_all("^_|_$", "")
}

shared_safe_char_042 <- function(x) {
    x <- as.character(x)
    x[is.na(x)] <- "NA"
    trimws(x)
}

shared_collapse_or_none_043 <- function(x) {
    x <- unique(safe_char(x))
    x <- x[x != "" & x != "NA"]
    if (length(x) == 0L) {
        return("none")
    }
    paste(sort(x), collapse = ",")
}

shared_safe_char_044 <- function(x) {
    y <- as.character(x)
    y[is.na(y)] <- "NA"
    trimws(y)
}

shared_round_or_na_045 <- function(x, digits = 6L) {
    x <- as.numeric(x)
    x[!is.finite(x)] <- NA_real_
    round(x, digits = digits)
}

shared_write_tsv_046 <- function(x, path) fwrite(as.data.table(x), path, sep = "\t", quote = TRUE, na = "NA")

shared_write_tsv_047 <- function(x, path) fwrite(x, path, sep = "\t", quote = TRUE, na = "NA")

shared_write_tsv_048 <- function(x, path) write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

shared_write_tsv_049 <- function(x, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    fwrite(as.data.table(x), path, sep = "\t", quote = FALSE, na = "NA")
}

shared_write_tsv_050 <- function(x, path) {
    fwrite(as.data.table(x), path, sep = "\t", quote = TRUE, na = "NA")
}

shared_write_new_051 <- function(x, path) {
    if (file.exists(path) && Sys.getenv("BRAUN_FORCE", "false") != "true")
        stop("Refusing to overwrite existing output: ", path)
    fwrite(x, path, sep = "\t", na = "NA", quote = FALSE)
}

shared_write_tsv_052 <- function(x, path) {
    write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

shared_write_tsv_053 <- function(x, path, row.names = FALSE) {
    utils::write.table(x, file = path, sep = "\t", quote = FALSE, row.names = row.names, na = "NA", col.names = TRUE)
}

shared_assert_cols_054 <- function(x, required, object_name) {
    missing_cols <- setdiff(required, names(x))
    if (length(missing_cols) > 0L) {
        stop(object_name, " missing columns: ", paste(missing_cols, collapse = ", "))
    }
    invisible(TRUE)
}

shared__notin__055 <- function(x, table) !(x %in% table)
