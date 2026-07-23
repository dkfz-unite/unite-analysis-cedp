# functions to write the overall models

.write_table_with_header <- function(table, header, path, row.names=TRUE) {
    con <- file(path, open = "w")
    on.exit(close(con))
    writeLines(paste0("# ", strsplit(header, "\n")[[1]]), con)
    write.table(table, con, sep = "\t", quote = FALSE, row.names = row.names, na = "")
}

write_contrasts <- function(prs, path) {
    header <- paste(attr(prs, "mesg"), collapse = "\n")
    .write_table_with_header(prs, header, path, row.names = FALSE)
}

write_lm_test <- function(av, path) {
    header <- paste(attr(av, "heading"), collapse = "\n")
    .write_table_with_header(av, header, path, row.names=FALSE)
}

write_rfit_drop_test <- function(dr, path) {
    formulas <- attr(dr, "formulas")
    header <- paste0(
        "Rfit::drop.test",
        "\n  M1: ", deparse1(formulas$reduced),
        "\n  M2: ", deparse1(formulas$full)
    )
    tbl <- data.frame(
        Df       = paste(dr$df1,dr$df2,sep=","),
        RD       = dr$RD,
        F        = dr$F[1],
        p = dr$p.value[1]
    )
    .write_table_with_header(tbl, header, path, row.names=FALSE)
}

write_rfit_summary <- function(dr, path) {
    formulas <- attr(dr, "formulas")
    header <- paste0(
        "Rfit::summary",
        "\n  M1: ", deparse1(formulas$reduced),
        "\n  M2: ", deparse1(formulas$full)
    )
    tbl <- data.frame(
        Df = paste(attr(dr, "df1"), attr(dr, "df2"), sep = ","),
        F  = dr$dropstat,
        p  = dr$droppval
    )
    .write_table_with_header(tbl, header, path, row.names = FALSE)
}

write_skipped_test <- function(tbl, path) {
    header <- attr(tbl, "heading")
    .write_table_with_header(tbl, header, path, row.names = FALSE)
}

# functions to fit each model type
.fit_lm <- function(full_formula, reduced_formula, model_data) {
    model         <- lm(full_formula, data = model_data)
    reduced_model <- lm(reduced_formula, data = model_data)
    av <- anova(reduced_model, model, test = "F")
    emm <- emmeans::emmeans(model, ~condition)
    return(list(model = model, reduced_model = reduced_model, overall_test = av, emm = emm, write_table_func=write_lm_test))
}

.is_intercept_only <- function(formula) {
    length(attr(terms(formula), "term.labels")) == 0
}

.fit_rfit <- function(full_formula, reduced_formula, model_data) {
    model         <- Rfit::rfit(full_formula, data = model_data)
    if (.is_intercept_only(reduced_formula)) {
        # if the reduced formula contains only an intercept
        # we cannot use drop-test (Rfit cannot fit an intercept-only model)
        # but we can get the significance from summary of the full model
        dr <- summary(model)
        # add the degrees of freedom to the object
        pfull <- length(attr(terms(full_formula), "term.labels")) + 1
        pred  <- length(attr(terms(reduced_formula), "term.labels"))
        df1 <- pfull - pred
        df2 <- nrow(model_data) - pfull
        attr(dr, "df1") <- df1
        attr(dr, "df2") <- df2
        attr(dr, "formulas") <- list(full = full_formula, reduced = reduced_formula)
        write_output_func <- write_rfit_summary
        # we still need to return a valid reduced model object
        # for generating covariate-adjusted values (which, in the absence of covariates, will simply be equal to input values)
        # for this purpose it does not matter at all how the model is fitted, since it will have no predictive power in any case
        # so we can fall back to lm for this case
        reduced_model <- lm(reduced_formula, data=model_data)
    } else {
        # make a drop test between the full and reduced model
        reduced_model <- Rfit::rfit(reduced_formula, data = model_data)
        dr <- Rfit::drop.test(model, reduced_model)
        attr(dr, "formulas") <- list(full = full_formula, reduced = reduced_formula)
        write_output_func <- write_rfit_drop_test
    }

    df_resid <- nrow(model_data) - length(coef(model))
    emm <- emmeans::emmeans(
        emmeans::qdrg(full_formula, data = model_data,
                      coef = coef(model), vcov = vcov(model), df = df_resid),
        ~ condition
    )
    return(list(model = model, reduced_model = reduced_model, overall_test = dr, emm = emm, write_table_func=write_output_func))
}

# build the formula objects
.build_formulas <- function(covariates) {
    cov_terms       <- if (!is.null(covariates)) paste(names(covariates), collapse = " + ") else NULL
    full_formula    <- as.formula(paste("outcome ~ condition", if (!is.null(cov_terms)) paste("+", cov_terms) else ""))
    reduced_formula <- as.formula(paste("outcome ~", if (is.null(cov_terms)) "1" else cov_terms))
    return(list(cov_terms = cov_terms, full_formula = full_formula, reduced_formula = reduced_formula))
}

# dispatch correct model fitting function
.fit_model <- function(full_formula, reduced_formula, model_data, method) {
        method <- match.arg(method, c("lm", "rfit"))
        fit <- switch(method,
        lm             = .fit_lm(full_formula, reduced_formula, model_data),
        rfit           = .fit_rfit(full_formula, reduced_formula, model_data)
    )
    return(fit)
}

# condition has fewer than 2 unique values: no comparison is possible, so skip
# the overall test and pairwise contrasts, but still fit the covariate-only
# (reduced) model so covariate-adjusted values can be returned
.fit_model_single_condition <- function(outcome, condition, covariates, model_type) {
    method       <- match.arg(model_type, c("lm", "rfit"))
    model_data   <- .get_model_data(outcome = outcome, condition = condition, covariates = covariates)
    reduced_formula <- .build_formulas(covariates)$reduced_formula
    # Rfit::rfit cannot fit an intercept-only model (errors with "x cannot only
    # contain an intercept", which is what reduced_formula is when there are no
    # covariates) and errors on a single observation regardless of formula
    # ("subscript out of bounds"); fall back to lm in either case
    use_rfit <- method == "rfit" && !is.null(covariates) && nrow(model_data) >= 2
    reduced_model <- if (use_rfit) {
        Rfit::rfit(reduced_formula, data = model_data)
    } else {
        lm(reduced_formula, data = model_data)
    }

    message <- sprintf(
        "condition has only %d unique value(s); overall test and pairwise contrasts were not computed.",
        length(unique(condition))
    )
    overall_test <- data.frame(note = message)
    attr(overall_test, "heading") <- message

    contrasts <- data.frame(note = message)
    attr(contrasts, "mesg") <- message

    return(list(reduced_model = reduced_model, overall_test = overall_test,
                contrasts = contrasts, write_table_func = write_skipped_test))
}

# build data frame of data for modelling
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
    # contrasts_df <- data.frame(
    #     contrast    = as.character(prs$contrast),
    #     estimate    = prs$estimate,
    #     SE          = prs$SE,
    #     CI_lower = prs$lower.CL,
    #     CI_upper = prs$upper.CL,
    #     p_value = prs$p.value
    # )
    # # keep info
    # comment(contrasts_df) <- paste(attr(prs, "mesg"), collapse = "\n")
    return(prs)
}

.get_plot_values <- function(model, outcome, return_covariate_adjusted) {
   values <- if (return_covariate_adjusted) {
        residuals(model) + mean(outcome, na.rm = TRUE)
    } else {
        outcome
    }
    return(values)
}

get_outcome <- function(data, feature) {
    if (!feature %in% colnames(data))
        stop("Feature '", feature, "' not found in data.")
    as.numeric(data[, feature])
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
#' @param sample Optional vector of sample identifiers, one per observation.
#'   If supplied, it is included as the first column of \code{values}. NULL
#'   (default) omits this column.
#' @param method Character string selecting the modelling approach. One of:
#'   \describe{
#'     \item{"lm"}{Ordinary least squares via \code{lm}. Overall test uses a
#'       partial F-test from \code{drop1}.}
#'     \item{"rfit"}{Rank-based regression via \code{Rfit::rfit}, robust to
#'       outliers. Overall test uses \code{Rfit::drop.test}.}
#'   }
#' @param return_covariate_adjusted Logical. If \code{TRUE}, the returned \code{values}
#'   are covariate-adjusted: residuals from the covariate-only (reduced) model
#'   re-centred on the grand mean, so condition differences are preserved but
#'   covariate effects are removed. If \code{FALSE} (default), raw outcome
#'   values are returned.
#'
#' @details If \code{condition} has fewer than 2 unique values, no comparison is
#'   possible: the overall test and pairwise contrasts are skipped and replaced
#'   with a one-row placeholder explaining why, but \code{values} is still
#'   computed from the covariate-only model so covariate-adjusted values are
#'   always returned.
#'
#' @return A named list with three elements:
#'   \describe{
#'     \item{overall_test}{data.frame with columns \code{statistic}, \code{df},
#'       and \code{p_value} for the test of the condition effect.}
#'     \item{contrasts}{data.frame with columns \code{contrast}, \code{estimate},
#'       \code{SE}, \code{CI_lower}, \code{CI_upper}, \code{p_value}, and
#'       for all pairwise condition comparisons (p valuyes are tukey adjusted)
#'     \item{values}{data.frame with columns \code{condition} and \code{value}
#'       (preceded by \code{sample} if supplied) containing the outcome
#'       values (adjusted or unadjusted) per sample.}
#'   }
fit_model <- function(outcome, condition, covariates = NULL, model_type = "lm", return_covariate_adjusted = FALSE, sample = NULL) {

    if (length(unique(condition)) < 2) {
        # with only one condition we cannot test its effect, 
        # but for downstream we still need a reduced model
        # the below returns 'dummy' data structures for the outputs that cannot
        # be computed
        fit <- .fit_model_single_condition(outcome = outcome,
                                            condition = condition,
                                            covariates = covariates,
                                            model_type = model_type)
        contrasts <- fit$contrasts
    } else {
        model_data <- .get_model_data(outcome=outcome,
                                      condition=condition,
                                      covariates=covariates)

        formulas        <- .build_formulas(covariates)
        full_formula    <- formulas$full_formula
        reduced_formula <- formulas$reduced_formula
        fit <- .fit_model(full_formula=full_formula,
                          model_data=model_data,
                          reduced_formula=reduced_formula,
                          method=model_type)

        # ---- pairwise contrasts ----
        contrasts <- .pairwise_contrasts(emm=fit$emm,
                                            adjust="tukey")
    }

    # return either covariate adjusted or unadjusted values
    # residualising against the reduced model
    values <- .get_plot_values(model=fit$reduced_model,
                               outcome=outcome,
                               return_covariate_adjusted=return_covariate_adjusted)

    values_df <- data.frame(condition = condition, value = values)
    if (!is.null(sample)) {
        values_df <- cbind(sample = sample, values_df)
    }

    list(
        write_table_func = fit$write_table_func,
        overall_test = fit$overall_test,
        contrasts    = contrasts,
        values       = values_df
    )
}

#' Write model results to TSV files
#'
#' @param results Named list returned by \code{fit_model}.
#' @param output_dir Directory to write files into.
#' @param prefix Optional string prepended to each filename (e.g. a feature name).
write_model_results <- function(results, output_dir, prefix = NULL, write_test_fn = write_lm_test) {
    fname <- function(name) {
        parts <- c(prefix, name)
        file.path(output_dir, paste0(paste(parts, collapse = "_"), ".tsv"))
    }
    # execute custom writer for each model type
    results$write_table_func(results$overall_test, fname("overall_test"))
    # execute generic writers for contrasts and values
    write_contrasts(results$contrasts, fname("contrasts"))
    write.table(results$values,    fname("values"),    sep = "\t", row.names = FALSE, quote = FALSE, na="")
    invisible(NULL)
}
