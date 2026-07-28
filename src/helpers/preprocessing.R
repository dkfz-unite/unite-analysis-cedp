library(limma)
library(sva)


filter_proteins <- function(data,class_labels, min_frac_in_one_class=FALSE, min_non_na_fraction = 0.5) {
    # data: a numeric matrix with samples as rows and features (proteins) as columns
    # class_labels: a vector of class labels corresponding to each sample (row) in data
    # min_frac_in_one_class: if TRUE, keep proteins that are present in at least min_non_na_fraction of samples in at least one class;
    # if FALSE, keep proteins that are present in at least min_non_na_fraction of samples overall
    # min_non_na_fraction: minimum fraction of non-NA values required to keep a feature
    
    if (min_frac_in_one_class) {
        # Check non-NA fraction within each class
        classes <- unique(class_labels)
        keep_features <- sapply(1:ncol(data), function(i) {
            any(sapply(classes, function(cls) {
                cls_indices <- which(class_labels == cls)
                sum(!is.na(data[cls_indices, i])) >= (length(cls_indices) * min_non_na_fraction)
            }))
        })
    } else {
        # count number of features with non-NA values across all samples
        non_na_counts <- colSums(!is.na(data))
        # Check non-NA fraction across all samples
        total_samples <- nrow(data)
        keep_features <- non_na_counts >= (total_samples * min_non_na_fraction)
    }
    filtered_data <- data[, keep_features, drop = FALSE]

    # raise a sensible error if no features pass the filter
    if (ncol(filtered_data) == 0) {
        stop("No features passed the filtering criteria. Consider lowering the min_non_na_fraction threshold.")
    }

    return(filtered_data)
}

impute_data <- function(data, stratify_by_batch=TRUE, batch_vector=NULL, method = "MinDet", global_min = FALSE, pct=0.01) {
    # performs imputtaion of a dataset with missing values (NA) using either the MinDet or MinProb method, optionally stratified by batch
    method <- tolower(as.character(method))
    if (method == "mindet") {
        method <- "deterministic"
    } else if (method == "minprob") {
        method <- "stochastic"
    }


    # check if there are fewer than 3 non-NA values in any given batch, if so fall back to no stratification
    if (stratify_by_batch && !is.null(batch_vector)) {
        for (batch in unique(batch_vector)) {
            batch_indices <- which(batch_vector == batch)
            batch_data <- data[batch_indices, ]
            non_na_counts <- colSums(!is.na(batch_data))
            if (any(non_na_counts < 3)) {
                warning(paste("Batch", batch, "has fewer than 3 non-NA values for some features. Falling back to no stratification for imputation."))
                stratify_by_batch <- FALSE
                break
            }
        }
    }

    
    if (stratify_by_batch) {
        if (is.null(batch_vector)) {
            stop("batch_vector must be provided when stratify_by_batch is TRUE")
        }
        # Ensure batch_vector is a factor
        batch_vector <- as.factor(batch_vector)
        
        # Create a copy of the data to store imputed values
        imputed_data <- data
        
        # Loop through each batch and apply imputation separately
        for (batch in levels(batch_vector)) {
            batch_indices <- which(batch_vector == batch)
            # passing null here will cause the function to use the median of the per-protein standard deviations as the sigma for stochastic imputation, which is a reasonable default
            imputed_data[batch_indices, ] <- min_det_prob_imputation(data[batch_indices, ], global_min=global_min, pct=pct, method=method, sigma=NULL)
        }
    } else {
        # Apply imputation to the entire dataset without stratification
        # passing null here will cause the function to use the median of the per-protein standard deviations as the sigma for stochastic imputation, which is a reasonable default
        imputed_data <- min_det_prob_imputation(data, global_min=global_min, pct=pct, method=method, sigma=NULL) 
    }
    rownames(imputed_data) <- rownames(data)
    colnames(imputed_data) <- colnames(data)
    return(imputed_data)
}


min_det_prob_imputation <- function(data, global_min = FALSE, pct=0.01, method= "deterministic", sigma=NULL) {
    # data: a numeric matrix with missing values (NA)
    # global_min: if TRUE, use the same minimum value for all rows; if FALSE, use row-specific minimums
    # pct: the percentile to use for determining the minimum value (e.g., 0.01 for 1st percentile)
    # method: "deterministic" or "stochastic"
    # sigma: standard deviation for stochastic method, if NULL it will be set to the median of row  standard deviations
    if (!is.matrix(data) && !is.data.frame(data)) {
        stop("Input data must be a matrix or data frame")
    }
    

    data <- as.data.frame(data)
    
    # for each sample (row), find the 'minimum' value
    if (global_min) {
        # always apply the same global minimum value
        min_value <- quantile(data, probs = pct, na.rm = TRUE)
        # repeat min for each column
        mins <- rep(min_value, nrow(data))
    }
    else {
        # apply sample (row)-specific minimum values
        mins <- apply(data, 1, function(x) quantile(x, probs = pct, na.rm = TRUE))
    }

    method <- tolower(as.character(method))

    if (method == "deterministic") {
        # replace NAs for each sample (row) with the corresponding minimum value
        for (j in seq_len(nrow(data))) {
            na_indices <- which(is.na(data[j, ]))
            data[j, na_indices] <- mins[j]
        }
    } else if (method == "stochastic") {
        # replace NAs with random values drawn from N(min, sigma)
        # get the median of the per-protein (column) standard deviations

        if (is.null(sigma)) {
            col_sds <- apply(data, 2, sd, na.rm = TRUE)
            sigma <- median(col_sds, na.rm = TRUE)
        }
        print(paste("Using sigma =", sigma, "for stochastic imputation"))
        print(paste("Using mins:", paste(round(mins, 3), collapse = ", ")))
        for (j in seq_len(nrow(data))) {
            na_indices <- which(is.na(data[j, ]))
            n_na <- length(na_indices)
            if (n_na > 0) {
                imputed_values <- rnorm(n_na, mean = mins[j], sd = sigma)
                data[j, na_indices] <- imputed_values
            }
        }
    } else {
        stop(paste("Invalid method.", method, "Choose 'deterministic' or 'stochastic'."))
    }
    rownames(data) <- rownames(data)
    colnames(data) <- colnames(data)
  return(data)
}

normalize_data <- function(data, method = "median") {
    method <- tolower(as.character(method))
    if (method == "quantile") {
        # transpose rows to columns for quantile normalization, to normalize rows (samples)
        normalized_data <- as.data.frame(t(normalizeQuantiles(t(as.matrix(data)))))
    } else if (method == "median") {
        normalized_result <- apply(data, 1, function(x) x - median(x, na.rm = TRUE))
        normalized_data <- as.data.frame(t(normalized_result))
    } else {
        stop(paste("Unknown normalization method:", method))
    }
    rownames(normalized_data) <- rownames(data)
    colnames(normalized_data) <- colnames(data)
    return(normalized_data)
}

batch_correct <- function(data, batch_vector, method = "ComBat") {
    # convert batch_vector to factor if it's not already
    batch_vector <- as.factor(batch_vector)
    method <- tolower(as.character(method))
     if (method == "combat") {    

        data_matrix <- t(as.matrix(data))
        # launch figure device to capture prior plots
        corrected_matrix <- ComBat(
        dat = data_matrix, 
        batch = batch_vector,
        mod = NULL,  
        par.prior = TRUE,
        prior.plots = FALSE
        )
        # Convert back to data frame with original structure
        corrected_data <- as.data.frame(t(corrected_matrix))
  } else if (method == "limma") {
    
        # limma expects features as rows, samples as columns
        data_matrix <- t(as.matrix(data))
        
        # Apply limma removeBatchEffect
        corrected_matrix <- removeBatchEffect(
        x = data_matrix,
        batch = batch_vector
        )
        
        # Convert back to data frame
        corrected_data <- as.data.frame(t(corrected_matrix))
    
  } else {
    stop(paste("Unknown batch correction method:", method))
  }
    rownames(corrected_data) <- rownames(data)
    colnames(corrected_data) <- colnames(data)
    return(corrected_data)

}
get_required <- function(lst, key) {
  # access an element from a named list, throwing an error if the key is not present
  if (!key %in% names(lst)) stop(paste0("Required option '", key, "' not found"))
  return(lst[[key]])
}

replace_required <- function(lst, key, value) {
  # replace an element in a named list, throwing an error if the key is not present
  if (!key %in% names(lst)) stop(paste0("Required option '", key, "' not found"))
  lst[key] <- list(value)
  return(lst)           
}

preprocess_data <- function(data, batch_vector, class_labels, options) {
    print("Preprocessing data with the following options:")
    print(options)
    
    # check if batch_vector has nay NA or all empty, if so set to NULL for downstream functions
    if (any(is.na(batch_vector)) || any(batch_vector == "")) {
        batch_vector <- NULL
    } else {
        batch_vector <- as.factor(batch_vector)
    }

    data <- proteomic_data_preprocessing(data=data,
        normalization_method=get_required(options, "normalization_method"),
        imputation_method=get_required(options, "imputation_method"),
        stratify_imputation_by_batch=get_required(options, "stratify_imputation_by_batch"),
        batch_vector=batch_vector,
        class_labels=as.factor(class_labels),
        log_offset=get_required(options, "normalization_log_offset"),
        batch_correct_method=get_required(options, "batch_correction_method"),
        min_non_na_fraction=get_required(options, "min_non_missing_fraction"),
        min_frac_in_one_class=get_required(options, "require_min_fraction_one_class")
    )
    return(data)

}


proteomic_data_preprocessing <- function(data, normalization_method, imputation_method, stratify_imputation_by_batch, batch_vector, class_labels, log_offset, batch_correct_method, min_non_na_fraction, min_frac_in_one_class) {
        # convert any zeros to NA for imputation
        data[data == 0] <- NA

        data <- filter_proteins(data, class_labels = class_labels, min_frac_in_one_class = min_frac_in_one_class, min_non_na_fraction = min_non_na_fraction)

        #### log transformation
        data <- log2(data + log_offset)
        #### normalization
        data <- normalize_data(data, method = normalization_method)
        #### imputation
        data <- impute_data(data, method = imputation_method, stratify_by_batch = stratify_imputation_by_batch, batch_vector = batch_vector)
        #### batch correction
        if (!is.null(batch_correct_method)){ 
            if (!is.null(batch_vector)) {
                data <- batch_correct(data, batch_vector, method = batch_correct_method)
            }
            else {
                warning("batch_vector is NULL, skipping batch correction")
            }
        }
    return(data)
}