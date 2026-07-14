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
    expect_named(fit, c("model", "reduced_model", "overall_test", "emm"))
})

test_that(".fit_lm overall_test has correct columns", {
    fit <- .fit_lm(full_formula, reduced_formula, fixture)
    expect_named(fit$overall_test, c("statistic", "df", "p_value"))
})

test_that(".fit_lm overall_test$df is formatted as 'numerator, denominator'", {
    fit <- .fit_lm(full_formula, reduced_formula, fixture)
    expect_match(fit$overall_test$df, "^\\d+, \\d+$")
})

test_that(".fit_lm overall_test$p_value is in [0, 1]", {
    fit <- .fit_lm(full_formula, reduced_formula, fixture)
    expect_gte(fit$overall_test$p_value, 0)
    expect_lte(fit$overall_test$p_value, 1)
})

test_that(".fit_lm detects clear group difference", {
    fit <- .fit_lm(full_formula, reduced_formula, fixture)
    expect_lt(fit$overall_test$p_value, 0.05)
})

# ─── .fit_rfit ───────────────────────────────────────────────────────────────

test_that(".fit_rfit returns list with correct names", {
    fit <- .fit_rfit(full_formula, reduced_formula, fixture)
    expect_named(fit, c("model", "reduced_model", "overall_test", "emm"))
})

test_that(".fit_rfit overall_test has correct columns", {
    fit <- .fit_rfit(full_formula, reduced_formula, fixture)
    expect_named(fit$overall_test, c("statistic", "df", "p_value"))
})

test_that(".fit_rfit overall_test$df is formatted as 'numerator, denominator'", {
    fit <- .fit_rfit(full_formula, reduced_formula, fixture)
    expect_match(fit$overall_test$df, "^\\d+, \\d+$")
})

test_that(".fit_rfit overall_test$p_value is in [0, 1]", {
    fit <- .fit_rfit(full_formula, reduced_formula, fixture)
    expect_gte(fit$overall_test$p_value, 0)
    expect_lte(fit$overall_test$p_value, 1)
})

test_that(".fit_rfit detects clear group difference", {
    fit <- .fit_rfit(full_formula, reduced_formula, fixture)
    expect_lt(fit$overall_test$p_value, 0.05)
})

test_that(".fit_rfit emm is an emmGrid object", {
    fit <- .fit_rfit(full_formula, reduced_formula, fixture)
    expect_s4_class(fit$emm, "emmGrid")
})

# ─── .fit_gls_unequalvar ─────────────────────────────────────────────────────

test_that(".fit_gls_unequalvar returns list with correct names", {
    fit <- .fit_gls_unequalvar(full_formula, reduced_formula, fixture)
    expect_named(fit, c("model", "reduced_model", "overall_test", "emm"))
})

test_that(".fit_gls_unequalvar overall_test has correct columns", {
    fit <- .fit_gls_unequalvar(full_formula, reduced_formula, fixture)
    expect_named(fit$overall_test, c("statistic", "df", "p_value"))
})

test_that(".fit_gls_unequalvar overall_test$df is a single integer string (LRT df)", {
    fit <- .fit_gls_unequalvar(full_formula, reduced_formula, fixture)
    expect_match(fit$overall_test$df, "^\\d+$")
})

test_that(".fit_gls_unequalvar overall_test$p_value is in [0, 1]", {
    fit <- .fit_gls_unequalvar(full_formula, reduced_formula, fixture)
    expect_gte(fit$overall_test$p_value, 0)
    expect_lte(fit$overall_test$p_value, 1)
})

test_that(".fit_gls_unequalvar detects clear group difference", {
    fit <- .fit_gls_unequalvar(full_formula, reduced_formula, fixture)
    expect_lt(fit$overall_test$p_value, 0.05)
})

test_that(".fit_gls_unequalvar returned model is REML fit (not ML)", {
    fit <- .fit_gls_unequalvar(full_formula, reduced_formula, fixture)
    # call$method is NULL when REML is used (the default), "ML" when explicitly set
    expect_false(identical(fit$model$call$method, "ML"))
})

# ─── .safe_gls ───────────────────────────────────────────────────────────────

test_that(".safe_gls returns a gls object on well-behaved data", {
    fit <- .safe_gls(full_formula, data = fixture,
                     weights = nlme::varIdent(form = ~1 | condition))
    expect_s3_class(fit, "gls")
})

test_that(".safe_gls errors when a group has zero within-group variance", {
    # varIdent must estimate a scale factor of 0 for the constant group,
    # which is a boundary value — nlme emits a convergence warning
    zero_var_data <- data.frame(
        outcome   = c(rep(0, 10), rnorm(10, mean = 5)),
        condition = factor(rep(c("A", "B"), each = 10))
    )
    expect_error(
        .safe_gls(outcome ~ condition, data = zero_var_data,
                  weights = nlme::varIdent(form = ~1 | condition))
    )
})

test_that(".safe_gls errors when apVar signals variance estimation failure", {
    # mock nlme::gls to return a model whose apVar is a failure string
    # (the path that .safe_gls catches when nlme converges silently but incorrectly)
    local_mocked_bindings(
        gls = function(...) structure(
            list(apVar = "Non-positive definite approximate variance-covariance"),
            class = "gls"
        ),
        .package = "nlme"
    )
    expect_error(
        .safe_gls(outcome ~ condition, data = fixture,
                  weights = nlme::varIdent(form = ~1 | condition)),
        "variance estimation failed"
    )
})

# ─── .fit_model ──────────────────────────────────────────────────────────────


test_that(".fit_model dispatches to lm and returns correct structure", {
    fit <- .fit_model(full_formula, reduced_formula, fixture, method = "lm")
    expect_named(fit, c("model", "reduced_model", "overall_test", "emm"))
    expect_s3_class(fit$model, "lm")
})

test_that(".fit_model dispatches to rfit and returns correct structure", {
    fit <- .fit_model(full_formula, reduced_formula, fixture, method = "rfit")
    expect_named(fit, c("model", "reduced_model", "overall_test", "emm"))
    expect_s3_class(fit$model, "rfit")
})

test_that(".fit_model dispatches to gls_unequalvar and returns correct structure", {
    fit <- .fit_model(full_formula, reduced_formula, fixture, method = "gls_unequalvar")
    expect_named(fit, c("model", "reduced_model", "overall_test", "emm"))
    expect_s3_class(fit$model, "gls")
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
