library(limma)
library(readr)
library(dplyr)
library(tibble)
library(jsonlite)
source(file.path(getwd(), "helpers", "preprocessing.R"))
source(file.path(getwd(), "helpers", "linear_modelling.R"))
# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
workdir <- args[1]
metadata_file <- file.path(workdir, "metadata.tsv")
results_file <- file.path(workdir, "results.tsv")
data_file <- file.path(workdir, "data.tsv")
options_file <- file.path(workdir, "options.json")

# Read data and metadata
data <- read_tsv(data_file)
metadata <- read_tsv(metadata_file)
options <- fromJSON(options_file)
print(names(options))
# Preprocess data (log, normalize, impute, batch correct)
rownames(data_matrix) <- data[[1]] # set feature names as rownames
data_matrix <- as.data.frame(data[,-1]) # remove feature column for processing



# reorder metadata rows to match data matrix columns
metadata_matrix <- metadata_matrix[match(colnames(data_matrix), rownames(metadata_matrix)), ,drop = FALSE]

# transpose data_matrix as proteomic_data_preprocessing expects samples as rows, features as columns
data_matrix <- t(data_matrix)

# process options - batch correct method must be NULL ... this will be account5ed for in the model
options <- replace_required(options, "batch_correction_method", NULL)

# preprocess data
processed_data <- preprocess_data(data=data_matrix, 
                              batch_vector=metadata_matrix$batch,
                              class_labels=metadata_matrix$condition, 
                              options = options)

# get the relevant feature values
outcome <- processed_data[get_required(options, "feature"),]
# get the condition values from secod column
condtion <- metadata_matrix$condition

# get covariates (if valid batch values it will be a data frame with a single factor column 'batch')
# votherwise will be null
covs <- get_covariates(metadata_matrix$batch)
results <- fit_model(outcome=outcome, 
            condition=condition,
            covariates = covs,
            model_type = get_required(options,"model_type"),
            return_covariate_adjusted = TRUE)