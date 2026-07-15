library(testthat)
#source(file.path(getwd(), "helpers", "linear_modelling.R"))

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

# ─── fit_model: condition has a single level ─────────────────────────────────

set.seed(42)
single_condition_fixture <- data.frame(
    outcome   = rnorm(20),
    condition = factor(rep("A", 20)),
    age       = rnorm(20)
)

test_that(".fit_model_single_condition returns list with correct names", {
    fit <- .fit_model_single_condition(
        outcome    = single_condition_fixture$outcome,
        condition  = single_condition_fixture$condition,
        covariates = data.frame(age = single_condition_fixture$age),
        model_type = "lm"
    )
    expect_named(fit, c("reduced_model", "overall_test", "contrasts", "write_table_func"))
})

test_that("fit_model does not error when condition has one level (lm)", {
    expect_no_error(
        fit_model(
            outcome    = single_condition_fixture$outcome,
            condition  = single_condition_fixture$condition,
            covariates = data.frame(age = single_condition_fixture$age),
            model_type = "lm"
        )
    )
})

test_that("fit_model does not error when condition has one level (rfit)", {
    expect_no_error(
        fit_model(
            outcome    = single_condition_fixture$outcome,
            condition  = single_condition_fixture$condition,
            covariates = data.frame(age = single_condition_fixture$age),
            model_type = "rfit"
        )
    )
})

test_that("fit_model does not error when condition has one level and there are no covariates", {
    expect_no_error(
        fit_model(
            outcome    = single_condition_fixture$outcome,
            condition  = single_condition_fixture$condition,
            covariates = NULL,
            model_type = "lm"
        )
    )
})

test_that("fit_model returns valid covariate-adjusted values for every sample when condition has one level", {
    result <- fit_model(
        outcome    = single_condition_fixture$outcome,
        condition  = single_condition_fixture$condition,
        covariates = data.frame(age = single_condition_fixture$age),
        model_type = "lm",
        return_covariate_adjusted = TRUE
    )
    expect_equal(nrow(result$values), nrow(single_condition_fixture))
    expect_false(anyNA(result$values$value))
})

test_that("fit_model overall_test and contrasts carry an explanatory note when condition has one level", {
    result <- fit_model(
        outcome    = single_condition_fixture$outcome,
        condition  = single_condition_fixture$condition,
        covariates = data.frame(age = single_condition_fixture$age),
        model_type = "lm"
    )
    expect_match(attr(result$overall_test, "heading"), "condition has only 1 unique value")
    expect_match(attr(result$contrasts, "mesg"), "condition has only 1 unique value")
})

test_that("fit_model results for single-level condition can still be written to disk", {
    result <- fit_model(
        outcome    = single_condition_fixture$outcome,
        condition  = single_condition_fixture$condition,
        covariates = data.frame(age = single_condition_fixture$age),
        model_type = "lm"
    )
    dir <- tempfile()
    dir.create(dir)
    expect_no_error(write_model_results(result, dir))
    expect_true(file.exists(file.path(dir, "overall_test.tsv")))
    expect_true(file.exists(file.path(dir, "contrasts.tsv")))
    expect_true(file.exists(file.path(dir, "values.tsv")))
})

# ─── fit_model: one condition group has a single observation (>=2 groups total) ─

set.seed(42)
singleton_group_fixture <- data.frame(
    outcome   = rnorm(20),
    condition = factor(c("A", rep("B", 10), rep("C", 9))),
    age       = rnorm(20)
)

test_that("fit_model does not error when one condition group has a single observation (lm)", {
    expect_no_error(
        fit_model(
            outcome    = singleton_group_fixture$outcome,
            condition  = singleton_group_fixture$condition,
            covariates = data.frame(age = singleton_group_fixture$age),
            model_type = "lm"
        )
    )
})

test_that("fit_model does not error when one condition group has a single observation (rfit)", {
    expect_no_error(
        fit_model(
            outcome    = singleton_group_fixture$outcome,
            condition  = singleton_group_fixture$condition,
            covariates = data.frame(age = singleton_group_fixture$age),
            model_type = "rfit"
        )
    )
})

test_that("fit_model returns valid structures when one condition group has a single observation", {
    result <- fit_model(
        outcome    = singleton_group_fixture$outcome,
        condition  = singleton_group_fixture$condition,
        covariates = data.frame(age = singleton_group_fixture$age),
        model_type = "lm"
    )
    expect_gt(nrow(as.data.frame(result$contrasts)), 0)
    expect_equal(nrow(result$values), nrow(singleton_group_fixture))
    expect_false(anyNA(result$values$value))
})

test_that("fit_model results with a singleton condition group can still be written to disk", {
    result <- fit_model(
        outcome    = singleton_group_fixture$outcome,
        condition  = singleton_group_fixture$condition,
        covariates = data.frame(age = singleton_group_fixture$age),
        model_type = "lm"
    )
    dir <- tempfile()
    dir.create(dir)
    expect_no_error(write_model_results(result, dir))
    expect_true(file.exists(file.path(dir, "overall_test.tsv")))
    expect_true(file.exists(file.path(dir, "contrasts.tsv")))
    expect_true(file.exists(file.path(dir, "values.tsv")))
})

# ─── fit_model: only one condition value, and only one observation in total ────

single_observation_fixture <- data.frame(
    outcome   = 3.7,
    condition = factor("A"),
    age       = 0.5
)

test_that("fit_model does not error when there is only one observation in total (lm, with covariate)", {
    expect_no_error(
        fit_model(
            outcome    = single_observation_fixture$outcome,
            condition  = single_observation_fixture$condition,
            covariates = data.frame(age = single_observation_fixture$age),
            model_type = "lm"
        )
    )
})

test_that("fit_model does not error when there is only one observation in total (lm, no covariates)", {
    expect_no_error(
        fit_model(
            outcome    = single_observation_fixture$outcome,
            condition  = single_observation_fixture$condition,
            covariates = NULL,
            model_type = "lm"
        )
    )
})

test_that("fit_model does not error when there is only one observation in total (rfit, with covariate)", {
    expect_no_error(
        fit_model(
            outcome    = single_observation_fixture$outcome,
            condition  = single_observation_fixture$condition,
            covariates = data.frame(age = single_observation_fixture$age),
            model_type = "rfit"
        )
    )
})

test_that("fit_model does not error when there is only one observation in total (rfit, no covariates)", {
    expect_no_error(
        fit_model(
            outcome    = single_observation_fixture$outcome,
            condition  = single_observation_fixture$condition,
            covariates = NULL,
            model_type = "rfit"
        )
    )
})

test_that("fit_model returns valid placeholder structures for a single total observation (with covariate)", {
    result <- fit_model(
        outcome    = single_observation_fixture$outcome,
        condition  = single_observation_fixture$condition,
        covariates = data.frame(age = single_observation_fixture$age),
        model_type = "lm"
    )
    expect_s3_class(result$overall_test, "data.frame")
    expect_s3_class(result$contrasts, "data.frame")
    expect_match(attr(result$overall_test, "heading"), "condition has only 1 unique value")
    expect_equal(nrow(result$values), 1)
    expect_false(anyNA(result$values$value))
})

test_that("fit_model returns valid placeholder structures for a single total observation (no covariates)", {
    result <- fit_model(
        outcome    = single_observation_fixture$outcome,
        condition  = single_observation_fixture$condition,
        covariates = NULL,
        model_type = "lm"
    )
    expect_s3_class(result$overall_test, "data.frame")
    expect_s3_class(result$contrasts, "data.frame")
    expect_match(attr(result$overall_test, "heading"), "condition has only 1 unique value")
    expect_equal(nrow(result$values), 1)
    expect_false(anyNA(result$values$value))
})

test_that("fit_model results for a single total observation can still be written to disk (with covariate)", {
    result <- fit_model(
        outcome    = single_observation_fixture$outcome,
        condition  = single_observation_fixture$condition,
        covariates = data.frame(age = single_observation_fixture$age),
        model_type = "lm"
    )
    dir <- tempfile()
    dir.create(dir)
    expect_no_error(write_model_results(result, dir))
    expect_true(file.exists(file.path(dir, "overall_test.tsv")))
    expect_true(file.exists(file.path(dir, "contrasts.tsv")))
    expect_true(file.exists(file.path(dir, "values.tsv")))
})

test_that("fit_model results for a single total observation can still be written to disk (no covariates)", {
    result <- fit_model(
        outcome    = single_observation_fixture$outcome,
        condition  = single_observation_fixture$condition,
        covariates = NULL,
        model_type = "lm"
    )
    dir <- tempfile()
    dir.create(dir)
    expect_no_error(write_model_results(result, dir))
    expect_true(file.exists(file.path(dir, "overall_test.tsv")))
    expect_true(file.exists(file.path(dir, "contrasts.tsv")))
    expect_true(file.exists(file.path(dir, "values.tsv")))
})

# with 1 observation the reduced model has 0 residual df, so residuals are
# exactly 0 and the covariate-adjusted value reduces to the raw outcome
# regardless of model_type or covariates -- parameterized over both below,
# since testthat has no built-in equivalent of pytest's parametrize
for (method in c("lm", "rfit")) {
    for (has_covariate in c(TRUE, FALSE)) {
        for (adjusted in c(TRUE, FALSE)) {
            covs <- if (has_covariate) data.frame(age = single_observation_fixture$age) else NULL
            test_that(sprintf(
                "fit_model values are unchanged from outcome for a single total observation (%s, %s covariate, adjusted=%s)",
                method, if (has_covariate) "with" else "no", adjusted
            ), {
                result <- fit_model(
                    outcome    = single_observation_fixture$outcome,
                    condition  = single_observation_fixture$condition,
                    covariates = covs,
                    model_type = method,
                    return_covariate_adjusted = adjusted
                )
                expect_equal(result$values$value, single_observation_fixture$outcome)
            })
        }
    }
}

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

# ─── write_lm_test integration with writing output files ────────────────────────────────────────────────

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

# ─── write_rfit_test integration with output files ──────────────────────────────────────────────

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
    expect_true(all(c("Df", "RD", "F", "p") %in% names(tbl)))
})
