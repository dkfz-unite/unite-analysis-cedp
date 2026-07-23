library(testthat)
#source(file.path(getwd(),"helpers", "preprocessing.R")) # run test from "run" directory

# ─── min_det_prob_imputation ───────────────────────────────────────────────────

test_that("min_det_prob_imputation deterministic imputes with row quantile", {
  # Row has values 10, 20, 30 and one NA; 1st-percentile ≈ 10
  mat <- matrix(c(10, 20, 30, NA), nrow = 1)
  result <- min_det_prob_imputation(mat, pct = 0.01, method = "deterministic")
  expected <- unname(quantile(c(10, 20, 30), probs = 0.01, na.rm = TRUE))
  expect_equal(as.numeric(result[1, 4]), expected)
})

test_that("min_det_prob_imputation global_min uses same fill value for all rows", {
  mat <- matrix(c(10, 20, NA, 100, 200, NA), nrow = 2, byrow = TRUE)
  result <- min_det_prob_imputation(mat, global_min = TRUE, method = "deterministic")
  expect_equal(as.numeric(result[1, 3]), as.numeric(result[2, 3]))
})

test_that("min_det_prob_imputation row-specific min differs across rows", {
  # Row 1: low values (1–3); Row 2: high values (100–300)
  # Imputed values should differ between rows
  mat <- matrix(c(1, 2, 3, NA, 100, 200, 300, NA), nrow = 2, byrow = TRUE)
  result <- min_det_prob_imputation(mat, global_min = FALSE, method = "deterministic")
  expect_false(as.numeric(result[1, 4]) == as.numeric(result[2, 4]))
})


test_that("min_det_prob_imputation errors on non-matrix / non-data.frame input", {
  expect_error(min_det_prob_imputation(c(1, 2, 3)), "matrix or data frame")
})

test_that("min_det_prob_imputation errors on invalid method", {
  mat <- matrix(c(1, 2, 3, NA), nrow = 1)
  expect_error(min_det_prob_imputation(mat, method = "invalid"), "Invalid method")
})


test_that("min_det_prob_imputation stochastic imputed values differ across runs", {
  mat <- matrix(c(1, 2, NA, 4, 5, NA), nrow = 2)
  set.seed(1); r1 <- min_det_prob_imputation(mat, method = "stochastic")
  set.seed(2); r2 <- min_det_prob_imputation(mat, method = "stochastic")
  expect_false(identical(r1, r2))
})

# ─── normalize_data ────────────────────────────────────────────────────────────

test_that("normalize_data (median) produces zero row medians", {
  mat <- data.frame(a = c(1, 4), b = c(3, 6), c = c(5, 8))
  result <- normalize_data(mat, method = "median")
  row_medians <- apply(result, 1, median)
  expect_equal(unname(row_medians), c(0, 0))
})

test_that("normalize_data (quantile) produces expected quantile values", {
  mat <- data.frame(a = c(1, 4), b = c(3, 6), c = c(5, 8))
  result <- normalize_data(mat, method = "quantile")
  # After quantile normalization, the sorted values in each row should be the same
  sorted_rows <- t(apply(result, 1, sort))
  expect_equal(sorted_rows[1, ], sorted_rows[2, ])

})


test_that("normalize_data errors on unknown method", {
  mat <- data.frame(a = 1:3, b = 4:6)
  expect_error(normalize_data(mat, method = "unknown"), "Unknown normalization method")
})



# ─── impute_data ───────────────────────────────────────────────────────────────
# NOTE: impute_data passes `sigma` to min_det_prob_imputation but `sigma` is not
# a parameter of impute_data, causing "object 'sigma' not found" errors.
# The tests below document the *intended* behaviour; they will pass once that
# bug is fixed.

test_that("impute_data (no stratification, MinDet) removes all NAs", {
  mat <- data.frame(a = c(1, NA, 3), b = c(NA, 5, 6), c = c(7, 8, NA))
  result <- impute_data(mat, stratify_by_batch = FALSE, method = "MinDet")
  expect_false(any(is.na(result)))
})

test_that("impute_data errors when stratify=TRUE and batch_vector is NULL", {
  mat <- data.frame(a = 1:3, b = 4:6)
  expect_error(
    impute_data(mat, stratify_by_batch = TRUE, batch_vector = NULL),
    "batch_vector must be provided"
  )
})

test_that("impute_data warns and falls back when a batch has <3 samples with non-NA values for any proteins", {
  mat <- data.frame(
    a = c(1, NA, 3, 4, 5),
    b = c(NA, 5, 6, 7, 8),
    c = c(7, 8, NA, 9, 10)
  )
  batch <- c("A", "A", "B", "B", "B")   # batch A has only 2 samples with non-NA values for any proteins
  expect_warning(
    impute_data(mat, stratify_by_batch = TRUE, batch_vector = batch, method = "MinDet"),
    "fewer than 3 non-NA values"
  )
})

test_that("impute_data (no stratification) preserves dimensions", {
  mat <- data.frame(matrix(c(runif(10), NA, NA), nrow = 3))
  result <- impute_data(mat, stratify_by_batch = FALSE, method = "MinDet")
  expect_equal(dim(result), dim(mat))
})

# ─── filter_proteins ───────────────────────────────────────────────────────────
# data: samples as rows, features (proteins) as columns

test_that("filter_proteins keeps columns meeting overall non-NA threshold", {
  # col 1: 4/4 present; col 2: 2/4 (50%); col 3: 1/4 (25%)
  mat <- matrix(c(
    1,  5,  NA,
    2,  NA, NA,
    3,  7,  NA,
    4,  NA, 9
  ), nrow = 4, byrow = TRUE, dimnames = list(NULL, c("V1", "V2", "V3")))
  result <- filter_proteins(mat, c("A","A","B","B"), min_non_na_fraction = 0.5)
  expect_equal(ncol(result), 2)

  # check that the correct column was removed (the one with only 25% non-NA)
  expect_equal(colnames(result), c("V1", "V2"))

})

test_that("filter_proteins removes columns below overall non-NA threshold", {
  mat <- matrix(c(1, 2, 3, 4,   # col 1: 100%
                  NA, NA, NA, 9), # col 2: 25%
                nrow = 4)
  result <- filter_proteins(mat, c("A","A","B","B"), min_non_na_fraction = 0.5)
  expect_equal(ncol(result), 1)
})

test_that("filter_proteins overall mode: result is independent of class_labels", {
  mat <- matrix(c(1, 2, 3, 4,
                  NA, NA, NA, NA), nrow = 4)
  r1 <- filter_proteins(mat, c("A","A","B","B"), min_frac_in_one_class = FALSE, min_non_na_fraction = 0.5)
  r2 <- filter_proteins(mat, c("X","Y","X","Y"), min_frac_in_one_class = FALSE, min_non_na_fraction = 0.5)
  expect_equal(ncol(r1), ncol(r2))
})

test_that("filter_proteins with min_non_na_fraction=0 keeps all columns", {
  mat <- matrix(c(NA, NA, 1, 2), nrow = 2)  # 2 cols
  result <- filter_proteins(mat, c("A","B"), min_frac_in_one_class = FALSE, min_non_na_fraction = 0)
  expect_equal(ncol(result), 2)
})

test_that("filter_proteins with min_non_na_fraction=1 requires fully complete columns", {
  mat <- matrix(c(1, 1,    # col 1: complete
                  2, NA),  # col 2: one NA
                nrow = 2)
  result <- filter_proteins(mat, c("A","B"), min_frac_in_one_class = FALSE, min_non_na_fraction = 1)
  expect_equal(ncol(result), 1)
})

test_that("filter_proteins class mode keeps column meeting threshold in only one class", {
  # col 1: 100% in class A (rows 1-2), 0% in class B (rows 3-4)
  mat <- matrix(c(1, 2, NA, NA), nrow = 4, ncol = 1)
  result <- filter_proteins(mat, c("A","A","B","B"), min_frac_in_one_class = TRUE, min_non_na_fraction = 0.5)
  expect_equal(ncol(result), 1)

  mat <- matrix(c(NA, NA, 3, 4), nrow = 4, ncol = 1)
  result <- filter_proteins(mat, c("A","A","B","B"), min_frac_in_one_class = TRUE, min_non_na_fraction = 0.5)
  expect_equal(ncol(result), 1)
})


test_that("filter_proteins class mode removes column not meeting threshold in any class", {
  # col 1: 50% in A (1/2), 0% in B — threshold 0.6 so fails both
  # col 2: all present
  mat <- matrix(c(1, NA, NA, NA,
                  1., 2., 3., 4.), nrow = 4, dimnames = list(NULL, c("V1", "V2")))
  result <- filter_proteins(mat, c("A","A","B","B"), min_frac_in_one_class = TRUE, min_non_na_fraction = 0.6)
  expect_equal(ncol(result), 1)
  expect_equal(colnames(result), "V2")
})

test_that("filter_proteins class mode keeps column that fails overall but passes per-class", {
  # col has 2/6 non-NA overall (33%) but 2/2 in class A (100%)
  mat <- matrix(c(1, 2, NA, NA, NA, NA), nrow = 6, ncol = 1)
  class_labels <- c("A","A","B","B","B","B")
  per_class <- filter_proteins(mat, class_labels, min_frac_in_one_class = TRUE,  min_non_na_fraction = 0.5)
  expect_equal(ncol(per_class), 1)
})

test_that("filter_proteins errors when no features pass the filter", {
  mat <- matrix(c(NA, NA, NA, NA), nrow = 4, ncol = 1)
  expect_error(filter_proteins(mat, c("A","A","B","B"), min_frac_in_one_class = FALSE, min_non_na_fraction = 0.5),
               "No features passed the filtering criteria")
})

test_that("filter_proteins preserves row and column names", {
  mat <- matrix(c(1, 2, 3, 4,
                  NA, NA, NA, NA), nrow = 4,
                dimnames = list(c("s1","s2","s3","s4"), c("prot1","prot2")))
  result <- filter_proteins(mat, c("A","A","B","B"), min_non_na_fraction = 0.5)
  expect_equal(colnames(result), "prot1")
  expect_equal(rownames(result), c("s1","s2","s3","s4"))
})

test_that("filter_proteins returns all columns when all pass threshold", {
  mat <- matrix(runif(12), nrow = 4)  # 3 proteins, no NAs
  result <- filter_proteins(mat, c("A","B","A","B"), min_non_na_fraction = 0.5)
  expect_equal(ncol(result), 3)
})

# ─── batch_correct ─────────────────────────────────────────────────────────────
# data: samples as rows, features (proteins) as columns

make_batch_data <- function(n_samples_per_batch = 5, n_features = 20, batch_effect = 5, seed = 42) {
  set.seed(seed)
  n <- n_samples_per_batch * 2
  mat <- matrix(rnorm(n * n_features, mean = 10, sd = 1), nrow = n, ncol = n_features)
  # add a clear batch effect to the second batch
  batch_vector <- c(rep("A", n_samples_per_batch), rep("B", n_samples_per_batch))
  mat[batch_vector == "B", ] <- mat[batch_vector == "B", ] + batch_effect
  rownames(mat) <- paste0("s", seq_len(n))
  colnames(mat) <- paste0("prot", seq_len(n_features))
  list(data = as.data.frame(mat), batch = batch_vector)
}

test_that("batch_correct (limma) returns a data.frame", {
  d <- make_batch_data()
  result <- batch_correct(d$data, d$batch, method = "limma")
  expect_s3_class(result, "data.frame")
})

test_that("batch_correct (limma) preserves dimensions", {
  d <- make_batch_data()
  result <- batch_correct(d$data, d$batch, method = "limma")
  expect_equal(dim(result), dim(d$data))
})

test_that("batch_correct (limma) preserves row and column names", {
  d <- make_batch_data()
  result <- batch_correct(d$data, d$batch, method = "limma")
  expect_equal(rownames(result), rownames(d$data))
  expect_equal(colnames(result), colnames(d$data))
})

test_that("batch_correct (limma) reduces batch mean differences", {
  d <- make_batch_data(batch_effect = 10)
  batch_A <- d$batch == "A"
  batch_B <- d$batch == "B"

  mean_diff_before <- mean(abs(
    colMeans(d$data[batch_B, ]) - colMeans(d$data[batch_A, ])
  ))

  result <- batch_correct(d$data, d$batch, method = "limma")

  mean_diff_after <- mean(abs(
    colMeans(result[batch_B, ]) - colMeans(result[batch_A, ])
  ))

  expect_lt(mean_diff_after, mean_diff_before)
})

test_that("batch_correct (limma) accepts character batch_vector (coerces to factor)", {
  d <- make_batch_data()
  expect_no_error(batch_correct(d$data, as.character(d$batch), method = "limma"))
})

test_that("batch_correct (limma) accepts factor batch_vector", {
  d <- make_batch_data()
  expect_no_error(batch_correct(d$data, as.factor(d$batch), method = "limma"))
})

test_that("batch_correct (combat) returns a data.frame", {
  d <- make_batch_data()
  result <- batch_correct(d$data, d$batch, method = "combat")
  expect_s3_class(result, "data.frame")
})

test_that("batch_correct (combat) preserves dimensions", {
  d <- make_batch_data()
  result <- batch_correct(d$data, d$batch, method = "combat")
  expect_equal(dim(result), dim(d$data))
})

test_that("batch_correct (combat) preserves row and column names", {
  d <- make_batch_data()
  result <- batch_correct(d$data, d$batch, method = "combat")
  expect_equal(rownames(result), rownames(d$data))
  expect_equal(colnames(result), colnames(d$data))
})

test_that("batch_correct (combat) reduces batch mean differences", {
  d <- make_batch_data(batch_effect = 10)
  batch_A <- d$batch == "A"
  batch_B <- d$batch == "B"

  mean_diff_before <- mean(abs(
    colMeans(d$data[batch_B, ]) - colMeans(d$data[batch_A, ])
  ))

  result <- batch_correct(d$data, d$batch, method = "combat")

  mean_diff_after <- mean(abs(
    colMeans(result[batch_B, ]) - colMeans(result[batch_A, ])
  ))

  expect_lt(mean_diff_after, mean_diff_before)
})

test_that("batch_correct errors on unknown method", {
  d <- make_batch_data()
  expect_error(batch_correct(d$data, d$batch, method = "unknown"), "Unknown batch correction method")
})


# ─── proteomic_data_preprocessing ─────────────────────────────────────────────

test_that("proteomic_data_preprocessing preserves dimensions", {
  set.seed(42)
  mat <- data.frame(matrix(runif(20, 1, 100), nrow = 4))
  result <- proteomic_data_preprocessing(mat,
    normalization_method = "median",
    imputation_method = "MinDet",
    stratify_imputation_by_batch = FALSE,
    batch_vector = NULL,
    log_offset = 1,
    batch_correct_method = "limma",
    min_non_na_fraction = 0.5,
    min_frac_in_one_class = FALSE
  )
  expect_equal(dim(result), dim(mat))
})

test_that("proteomic_data_preprocessing removes all NAs (including zeros)", {
  set.seed(42)
  mat <- data.frame(matrix(c(runif(15, 1, 100), 0, NA, 0, NA, NA), nrow = 4))
  result <- proteomic_data_preprocessing(mat,
    normalization_method = "median",
    imputation_method = "MinDet",
    stratify_imputation_by_batch = FALSE,
    batch_vector = NULL,
    log_offset = 1,
    batch_correct_method = "limma",
    min_non_na_fraction = 0.5,
    min_frac_in_one_class = FALSE
  )
  expect_false(any(is.na(result)))
})

test_that("proteomic_data_preprocessing preserves row and column names", {
  set.seed(1)
  mat <- data.frame(matrix(runif(9, 1, 10), nrow = 3))
  rownames(mat) <- c("s1", "s2", "s3")
  result <- proteomic_data_preprocessing(mat,
    normalization_method = "median",
    imputation_method = "MinDet",
    stratify_imputation_by_batch = FALSE,
    batch_vector = NULL,
    log_offset = 1,
    batch_correct_method = "limma",
    min_non_na_fraction = 0.5,
    min_frac_in_one_class = FALSE
  )
  expect_equal(rownames(result), c("s1", "s2", "s3"))
  expect_equal(colnames(result), colnames(mat))
})

test_that("test argument case insensitivity", {
  set.seed(1)
  mat <- data.frame(matrix(runif(9, 1, 10), nrow = 3))
  ##### lowercase
  expect_no_error(proteomic_data_preprocessing(mat,
    normalization_method = "median",
    imputation_method = "mindet",
    stratify_imputation_by_batch = FALSE,
    batch_vector = NULL,
    log_offset = 1,
    batch_correct_method = "limma",
    min_non_na_fraction = 0.5,
    min_frac_in_one_class = FALSE
  ))
  expect_no_error(proteomic_data_preprocessing(mat,
    normalization_method = "quantile",
    imputation_method = "minprob",
    stratify_imputation_by_batch = FALSE,
    batch_vector = NULL,
    log_offset = 1,
    batch_correct_method = "combat",
    min_non_na_fraction = 0.5,
    min_frac_in_one_class = FALSE
  ))

  #### UPPERCASE
  expect_no_error(proteomic_data_preprocessing(mat,
    normalization_method = "MEDIAN",
    imputation_method = "MINDET",
    stratify_imputation_by_batch = FALSE,
    batch_vector = NULL,
    log_offset = 1,
    batch_correct_method = "LIMMA",
    min_non_na_fraction = 0.5,  
    min_frac_in_one_class = FALSE
  ))
  expect_no_error(proteomic_data_preprocessing(mat,
    normalization_method = "QUANTILE",
    imputation_method = "MINPROB",
    stratify_imputation_by_batch = FALSE,
    batch_vector = NULL,
    log_offset = 1,
    batch_correct_method = "COMBAT",
    min_non_na_fraction = 0.5,
    min_frac_in_one_class = FALSE
  ))

})
