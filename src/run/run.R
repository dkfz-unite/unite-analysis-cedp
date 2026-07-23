library(limma)
library(readr)
library(dplyr)
library(tibble)
library(jsonlite)
source(file.path(getwd(), "helpers", "preprocessing.R"))
source(file.path(getwd(), "helpers", "linear_modelling.R"))
# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
inputdir <- file.path(args[1], "input")
outputdir <- file.path(args[1], "output")

metadata_file <- file.path(inputdir, "metadata.tsv")
data_file <- file.path(inputdir, "data.tsv")
options_file <- file.path(inputdir, "options.json")

# Read data and metadata
data <- read_tsv(data_file)
metadata <- read_tsv(metadata_file)
options <- fromJSON(options_file)
print(names(options))


data_matrix <- as.data.frame(data[,-1]) # remove feature column for processing
rownames(data_matrix) <- data[[1]] # set feature names as rownames

metadata_matrix <- as.data.frame(metadata[, -1]) # assuming first column is sample names
rownames(metadata_matrix) <- metadata[[1]] # set sample names as rownames

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
outcome <- get_outcome(processed_data, get_required(options, "feature"))

# get covariates (if valid batch values it will be a data frame with a single factor column 'batch')
# votherwise will be null
covs <- get_covariates(metadata_matrix$batch)
results <- fit_model(outcome=outcome,
            condition=metadata_matrix$condition,
            covariates = covs,
            model_type = get_required(options,"model_type"),
            return_covariate_adjusted = TRUE,
            sample = rownames(metadata_matrix))

write_model_results(results=results,
                    output_dir=outputdir)