# functions to return the overall models

.fit_lm <- function(full_formula, reduced_formula, model_data) {
    model         <- lm(full_formula, data = model_data)
    reduced_model <- lm(reduced_formula, data = model_data)
    av <- anova(reduced_model, model, test = "F")
    overall_test <- data.frame(
        statistic = av[2, "F"],
        df        = paste0(av[2, "Df"], ", ", av[2, "Res.Df"]),
        p_value   = av[2, "Pr(>F)"]
    )

    emm <- emmeans::emmeans(model, ~condition)
    return(list(model = model, reduced_model = reduced_model, overall_test = overall_test, emm = emm))
}

.fit_rfit <- function(full_formula, reduced_formula, model_data) {
    model         <- Rfit::rfit(full_formula, data = model_data)
    reduced_model <- Rfit::rfit(reduced_formula, data = model_data)

    dr <- Rfit::drop.test(model, reduced_model)
    overall_test <- data.frame(
        statistic = dr$F[1],
        df        = paste0(dr$df1, ", ", dr$df2),
        p_value   = dr$p.value[1]
    )

    df_resid <- nrow(model_data) - length(coef(model))
    emm <- emmeans::emmeans(
        emmeans::qdrg(full_formula, data = model_data,
                      coef = coef(model), vcov = vcov(model), df = df_resid),
        ~ condition
    )
    return(list(model = model, reduced_model = reduced_model, overall_test = overall_test, emm = emm))
}

# Wrapper around nlme::gls that converts convergence warnings to errors and
# checks the post-fit apVar for silent failures.
.safe_gls <- function(...) {
    model <- withCallingHandlers(
        nlme::gls(...),
        warning = function(w) {
            if (grepl("converge|iteration limit|singular", conditionMessage(w), ignore.case = TRUE)) {
                stop(paste("GLS convergence failure:", conditionMessage(w)), call. = FALSE)
            }
        }
    )
    # apVar is a matrix on success; nlme sets it to a character string on failure
    if (inherits(model$apVar, "character")) {
        stop(paste("GLS variance estimation failed:", model$apVar), call. = FALSE)
    }
    model
}

.fit_gls_unequalvar <- function(full_formula, reduced_formula, model_data) {
    # fit once with "maximum likelihood" method (this is needed for the likelihood ratio test)
    # comparing the full and reduced models
    # F-test is not possible with unequal variances
    model_ml   <- .safe_gls(full_formula, data = model_data,
                             weights = nlme::varIdent(form = ~1 | condition), method = "ML")

    # although 'condition' is not present in the reduced model we keep the weights specification to account for the
    # unequal variances in the conditions. This should give better estimates of the covariate effects than if
    # we did not
    reduced_ml <- .safe_gls(reduced_formula, data = model_data,
                             weights = nlme::varIdent(form = ~1 | condition), method = "ML")

    lrt    <- anova(reduced_ml, model_ml)
    chi_df <- lrt[2, "df"] - lrt[1, "df"]
    overall_test <- data.frame(
        statistic = lrt[2, "L.Ratio"],
        df        = as.character(chi_df),
        p_value   = lrt[2, "p-value"]
    )
    # REML (default "method") gives better estimates of the variance parameters and CIs
    # compared to ML
    model         <- .safe_gls(full_formula, data = model_data,
                                weights = nlme::varIdent(form = ~1 | condition))

    # the same logic applies here as above for why we keep the uneven-condition variance specification
    reduced_model <- .safe_gls(reduced_formula, data = model_data,
                                weights = nlme::varIdent(form = ~1 | condition))
    browser()
    emm <- emmeans::emmeans(model, ~ condition, data=model_data)
    return(list(model = model, reduced_model = reduced_model, overall_test = overall_test, emm = emm))
}

.build_formulas <- function(covariates) {
    cov_terms       <- if (!is.null(covariates)) paste(names(covariates), collapse = " + ") else NULL
    full_formula    <- as.formula(paste("outcome ~ condition", if (!is.null(cov_terms)) paste("+", cov_terms) else ""))
    reduced_formula <- as.formula(paste("outcome ~", if (is.null(cov_terms)) "1" else cov_terms))
    return(list(cov_terms = cov_terms, full_formula = full_formula, reduced_formula = reduced_formula))
}

.fit_model <- function(full_formula, reduced_formula, model_data, method) {
        method <- match.arg(method, c("lm", "rfit", "gls_unequalvar"))
        fit <- switch(method,
        lm             = .fit_lm(full_formula, reduced_formula, model_data),
        rfit           = .fit_rfit(full_formula, reduced_formula, model_data)
    )
    return(fit)
}

.get_model_data <- function (condition, outcome, covariates=NULL){
    condition  <- as.factor(condition)
    model_data <- data.frame(outcome = outcome, condition = condition)
    if (!is.null(covariates)) {
        covariates <- as.data.frame(covariates)
        model_data <- cbind(model_data, covariates)
    }
    return(model_data)

}

.pairwise_contrasts <- function(emm, adjust) {
    prs <- summary(pairs(emm), adjust = adjust, infer = c(TRUE, TRUE))
    contrasts_df <- data.frame(
        contrast    = as.character(prs$contrast),
        estimate    = prs$estimate,
        SE          = prs$SE,
        CI_lower = prs$lower.CL,
        CI_upper = prs$upper.CL,
        p_value = prs$p.value,
    )
    return(contrasts_df)
}

.get_plot_values <- function(model, outcome, return_covariate_adjusted) {
   values <- if (return_covariate_adjusted) {
        residuals(fit$model) + mean(outcome, na.rm = TRUE)
    } else {
        outcome
    }
    return(values)
}

get_covariates <- function(batch) {
   if (any(is.na(batch)) || any(batch == "") || length(unique(batch)) < 2) {
        batch_vector <- NULL
} else {
  batch_vector <- as.factor(batch)
}
covariates <- if (!is.null(batch_vector)) data.frame(batch = batch_vector) else NULL

}


#' Fit a linear model and return overall condition test, pairwise contrasts, and outcome values
#'
#' @param outcome Numeric vector of outcome values.
#' @param condition Factor giving the conditioning variable of interest (2 or more levels).
#' @param covariates Optional data.frame of adjustment covariates. Column names are
#'   used directly in the model formula. NULL (default) fits an unadjusted model.
#' @param method Character string selecting the modelling approach. One of:
#'   \describe{
#'     \item{"lm"}{Ordinary least squares via \code{lm}. Overall test uses a
#'       partial F-test from \code{drop1}.}
#'     \item{"rfit"}{Rank-based regression via \code{Rfit::rfit}, robust to
#'       outliers. Overall test uses \code{Rfit::drop.test}.}
#'     \item{"gls_unequalvar"}{GLS with per-condition heteroscedastic errors via
#'       \code{nlme::gls} and \code{varIdent}. Overall test uses a likelihood
#'       ratio test (models fit by ML); emmeans and residuals use the REML fit.}
#'   }
#' @param return_covariate_adjusted Logical. If \code{TRUE}, the returned \code{values}
#'   are covariate-adjusted: residuals from the covariate-only (reduced) model
#'   re-centred on the grand mean, so condition differences are preserved but
#'   covariate effects are removed. If \code{FALSE} (default), raw outcome
#'   values are returned.
#'
#' @return A named list with three elements:
#'   \describe{
#'     \item{overall_test}{data.frame with columns \code{statistic}, \code{df},
#'       and \code{p_value} for the test of the condition effect.}
#'     \item{contrasts}{data.frame with columns \code{contrast}, \code{estimate},
#'       \code{SE}, \code{CI_lower}, \code{CI_upper}, \code{p_value}, and
#'       for all pairwise condition comparisons (p valuyes are tukey adjusted)
#'     \item{values}{data.frame with columns \code{condition} and \code{value}
#'       containing the outcome values (adjusted or unadjusted) per sample.}
#'   }
fit_model <- function(outcome, condition, covariates = NULL, model_type = "lm", return_covariate_adjusted = FALSE) {

    model_data <- .get_model_data(outcome=outcome,
                                  condition=condition,
                                  covariates=covariates)


    formulas        <- .build_formulas(covariates)
    full_formula    <- formulas$full_formula
    reduced_formula <- formulas$reduced_formula
    fit <- .fit_model(full_formula=full_formula,
                      reduced_formula=reduced_formula,
                      model_type=model_type)
  
    # ---- pairwise contrasts ----
    contrasts_df <- .pairwise_contrasts(emm=fit$emm,
                                        adjust="tukey")
    
    # return either covariate adjusted or unadjusted values
    # residualising against the reduced model
    values <- .get_plot_values(model=fit$reduced_model,
                               outcome=outcome,
                               return_covariate_adjusted=return_covariate_adjusted)

    list(
        overall_test = fit$overall_test,
        contrasts    = contrasts_df,
        values       = data.frame(condition = condition, value = values)
    )
}

#' Write model results to TSV files
#'
#' @param results Named list returned by \code{fit_model}.
#' @param output_dir Directory to write files into.
#' @param prefix Optional string prepended to each filename (e.g. a feature name).
write_model_results <- function(results, output_dir, prefix = NULL) {
    fname <- function(name) {
        parts <- c(prefix, name)
        file.path(output_dir, paste0(paste(parts, collapse = "_"), ".tsv"))
    }
    write.table(results$overall_test, fname("overall_test"), sep = "\t", row.names = FALSE, quote = FALSE)
    write.table(results$contrasts,    fname("contrasts"),    sep = "\t", row.names = FALSE, quote = FALSE)
    write.table(results$values,       fname("values"),       sep = "\t", row.names = FALSE, quote = FALSE)
    invisible(NULL)
}
