library(testthat)
source(file.path(getwd(), "helpers", "linear_modelling.R"))

# ─── Shared fixtures ───────────────────────────────────────────────────────────
# 2 groups with a clear mean difference (0 vs 5) and one covariate
set.seed(42)
fixture <- data.frame(
    outcome   = c(rnorm(10, mean = 0), rnorm(10, mean = 5)),
    condition = factor(rep(c("A", "B"), each = 10)),
    age       = rnorm(20)
)
full_formula    <- outcome ~ condition + age
reduced_formula <- outcome ~ age

# ─── .build_formulas ─────────────────────────────────────────────────────────

test_that(".build_formulas returns list with correct names", {
    result <- .build_formulas(NULL)
    expect_named(result, c("cov_terms", "full_formula", "reduced_formula"))
})

test_that(".build_formulas (no covariates)", {
    result <- .build_formulas(NULL)
    expect_equal(deparse(result$full_formula), "outcome ~ condition")
    expect_equal(deparse(result$reduced_formula), "outcome ~ 1")
})


test_that(".build_formulas (one covariate)", {
    result <- .build_formulas(data.frame(age = numeric(5)))
    expect_equal(deparse(result$full_formula), "outcome ~ condition + age")
    expect_equal(deparse(result$reduced_formula), "outcome ~ age")
})


test_that(".build_formulas (multiple covariates) all terms appear in both formulas", {
    result <- .build_formulas(data.frame(age = numeric(5), sex = numeric(5)))
    full_vars    <- all.vars(result$full_formula)
    reduced_vars <- all.vars(result$reduced_formula)
    expect_true(all(c("age", "sex") %in% full_vars))
    expect_true(all(c("age", "sex") %in% reduced_vars))
})

# ─── .fit_lm ─────────────────────────────────────────────────────────────────

test_that(".fit_lm returns list with correct names", {
    fit <- .fit_lm(full_formula, reduced_formula, fixture)
    expect_named(fit, c("model", "reduced_model", "overall_test", "emm", "write_table_func"))
})

test_that(".fit_lm overall_test is an anova object", {
    fit <- .fit_lm(full_formula, reduced_formula, fixture)
    expect_s3_class(fit$overall_test, "anova")
})



# ─── .fit_rfit ───────────────────────────────────────────────────────────────

test_that(".fit_rfit returns list with correct names", {
    fit <- .fit_rfit(full_formula, reduced_formula, fixture)
    expect_named(fit, c("model", "reduced_model", "overall_test", "emm", "write_table_func"))
})



test_that(".fit_rfit emm is an emmGrid object", {
    fit <- .fit_rfit(full_formula, reduced_formula, fixture)
    expect_s4_class(fit$emm, "emmGrid")
})


# ─── .fit_model ──────────────────────────────────────────────────────────────


test_that(".fit_model dispatches to lm and returns correct structure", {
    fit <- .fit_model(full_formula, reduced_formula, fixture, method = "lm")
    expect_s3_class(fit$model, "lm")
})

test_that(".fit_model dispatches to rfit and returns correct structure", {
    fit <- .fit_model(full_formula, reduced_formula, fixture, method = "rfit")
    expect_s3_class(fit$model, "rfit")
})

test_that(".fit_model errors on invalid method", {
    expect_error(
        .fit_model(full_formula, reduced_formula, fixture, method = "invalid")
    )
})

# ─── get_covariates ───────────────────────────────────────────────────────────

test_that("get_covariates returns a data.frame with a factor batch column for valid input", {
    result <- get_covariates(c("A", "A", "B", "B"))
    expect_s3_class(result, "data.frame")
    expect_s3_class(result$batch, "factor")
})

test_that("get_covariates returns NULL when batch contains NA", {
    expect_null(get_covariates(c("A", NA, "B", "B")))
})

test_that("get_covariates returns NULL when batch contains empty string", {
    expect_null(get_covariates(c("A", "", "B", "B")))
})

test_that("get_covariates returns NULL when batch has only one unique value", {
    expect_null(get_covariates(c("A", "A", "A")))
})

test_that("get_covariates preserves batch levels in the returned factor", {
    result <- get_covariates(c("A", "B", "C", "A"))
    expect_setequal(levels(result$batch), c("A", "B", "C"))
})

# ─── write_contrasts integration ─────────────────────────────────────────────

test_that("write_contrasts creates a file", {
    fit  <- .fit_lm(full_formula, reduced_formula, fixture)
    prs  <- .pairwise_contrasts(fit$emm, adjust = "tukey")
    path <- tempfile(fileext = ".tsv")
    write_contrasts(prs, path)
    expect_true(file.exists(path))
})

test_that("write_contrasts writes # header lines from mesg attribute", {
    fit  <- .fit_lm(full_formula, reduced_formula, fixture)
    prs  <- .pairwise_contrasts(fit$emm, adjust = "tukey")
    path <- tempfile(fileext = ".tsv")
    write_contrasts(prs, path)
    header <- grep("^#", readLines(path), value = TRUE)
    expect_true(length(header) > 0)
    expect_true(any(grepl(paste(attr(prs, "mesg"), collapse = ""), header, fixed = FALSE)))
})

test_that("write_contrasts data section is readable as a table", {
    fit  <- .fit_lm(full_formula, reduced_formula, fixture)
    prs  <- .pairwise_contrasts(fit$emm, adjust = "tukey")
    path <- tempfile(fileext = ".tsv")
    write_contrasts(prs, path)
    lines <- readLines(path)
    tbl   <- read.table(text = paste(grep("^[^#]", lines, value = TRUE), collapse = "\n"),
                        sep = "\t", header = TRUE, check.names = FALSE)
    expect_s3_class(tbl, "data.frame")
    expect_gt(nrow(tbl), 0)
})

# ─── write_lm_test integration ────────────────────────────────────────────────

test_that("write_lm_test creates a file", {
    fit  <- .fit_lm(full_formula, reduced_formula, fixture)
    path <- tempfile(fileext = ".tsv")
    write_lm_test(fit$overall_test, path)
    expect_true(file.exists(path))
})

test_that("write_lm_test header lines contain both formula strings", {
    fit  <- .fit_lm(full_formula, reduced_formula, fixture)
    path <- tempfile(fileext = ".tsv")
    write_lm_test(fit$overall_test, path)
    header <- grep("^#", readLines(path), value = TRUE)
    expect_true(any(grepl(deparse1(full_formula),    header, fixed = TRUE)))
    expect_true(any(grepl(deparse1(reduced_formula), header, fixed = TRUE)))
})

test_that("write_lm_test data section is readable as a table", {
    fit  <- .fit_lm(full_formula, reduced_formula, fixture)
    path <- tempfile(fileext = ".tsv")
    write_lm_test(fit$overall_test, path)
    lines <- readLines(path)
    tbl   <- read.table(text = paste(grep("^[^#]", lines, value = TRUE), collapse = "\n"),
                        sep = "\t", header = TRUE)
    expect_s3_class(tbl, "data.frame")
    expect_gt(nrow(tbl), 0)
})

# ─── write_rfit_test integration ──────────────────────────────────────────────

test_that("write_rfit_test creates a file", {
    fit  <- .fit_rfit(full_formula, reduced_formula, fixture)
    path <- tempfile(fileext = ".tsv")
    write_rfit_test(fit$overall_test, path)
    expect_true(file.exists(path))
})

test_that("write_rfit_test header contains Rfit::drop.test and both formula strings", {
    fit  <- .fit_rfit(full_formula, reduced_formula, fixture)
    path <- tempfile(fileext = ".tsv")
    write_rfit_test(fit$overall_test, path)
    header <- grep("^#", readLines(path), value = TRUE)
    expect_true(any(grepl("Rfit::drop.test",         header, fixed = TRUE)))
    expect_true(any(grepl(deparse1(full_formula),    header, fixed = TRUE)))
    expect_true(any(grepl(deparse1(reduced_formula), header, fixed = TRUE)))
})

test_that("write_rfit_test data section has expected columns", {
    fit  <- .fit_rfit(full_formula, reduced_formula, fixture)
    path <- tempfile(fileext = ".tsv")
    write_rfit_test(fit$overall_test, path)
    lines <- readLines(path)
    tbl   <- read.table(text = paste(grep("^[^#]", lines, value = TRUE), collapse = "\n"),
                        sep = "\t", header = TRUE, check.names = FALSE)
    expect_true(all(c("Res.Df", "Df", "RD", "F", "Pr(>F)") %in% names(tbl)))
})
