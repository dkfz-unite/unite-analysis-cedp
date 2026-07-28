# Condition expression distribution for proteins

## General
R Limma based differential protein expression analysis application wrapped with web API.

The service performs the following steps:
1. Reads the input data, metadata and options files.
2. Preprocesses the data (filters features, normalizes data, and imputes missing values) see `options configuration` below.
3. For a specific protein and categorical 'condition' tests for an effect of condition on protein expression and performs pairwise contrasts among the conditions.

## Configuration
To configure the application, change environment variables as required in [commands](https://github.com/dkfz-unite/unite-commands/blob/main/README.md#configuration) web service:
- `UNITE_COMMAND` - command to run the analysis package (`Rscript`).
- `UNITE_COMMAND_ARGUMENTS` - command arguments (`run.R {data}/{proc}`).
- `UNITE_SOURCE_PATH` - location of the source code in docker container (`/src`).
- `UNITE_DATA_PATH` - location of the data in docker container (`/mnt/data`).
- `UNITE_PROCESS_LIMIT` - maximum number of concurrent jobs (`1` - process is heavy and uses a lot of CPU).

## Installation

### Docker Compose
The easiest way to install the application is to use docker-compose:
- Environment configuration and installation scripts: https://github.com/dkfz-unite/unite-environment
- Analysis service configuration and installation scripts: https://github.com/dkfz-unite/unite-environment/tree/main/applications/unite-analysis-cedp

### Docker
[Dockerfile](Dockerfile) is used to build an image of the application.
To build an image run the following command:
```
docker build -t unite.analysis.cedp:latest .
```

All application components should run in the same docker network.
To create common docker network if not yet available run the following command:
```bash
docker network create unite
```

To run application in docker run the following command:
```bash
docker run \
--name unite.analysis.cedp \
--restart unless-stopped \
--net unite \
--net-alias cedp.analysis.unite.net \
-p 127.0.0.1:5310:80 \
-e ASPNETCORE_ENVIRONMENT=Release \
-e UNITE_COMMAND=Rscript \
-v ./data:/mnt/data:rw \
-d \
unite.analysis.cedp:latest
```

## Usage
- Place the data files `data.tsv`, `metadata.tsv` and `options.json` in the `./data/{proc}/input` directory on the host machine.
- Send a POST request to the `localhost:5310/api/run?key=[key]` endpoint, where `[key]` is the process key.
- Analysis will run the command `Rscript` with the arguments `run.R {data}/{proc}` where `{proc}` is the process key.
  - All entries of `{data}` will be replaced with the path to the data location in docker container (In the example `./data` on the host machine will be mounted to `/mnt/data` in container).
  - All entries of `{proc}` will be replaced with the process key.
- Analysis will try to find the files `data.tsv`,`metadata.tsv` and `options.json` in the `{proc}/input` subfolder of the data location and use them as input.
- Analysis will save the results to the `{proc}/output` folder.

### Data format
Data file `{proc}/data.tsv` should be in the following format:
```tsv
feature sample1 sample2 sample3 sample4
protein1 10 20 30 40
protein2 15 25 35 45
protein3 20 30 40 50
protein4 25 35 45 55
```  

Where:
- `feature` - identifier of the feature (protein). Should be first column.
- `sample1`, `sample2`, `sample3`, `sample4` - names of the samples.
- Values in the table are raw protein intensity values (not normalised, not filtered).

### Metadata format
Metadata file `{proc}/metadata.tsv` should be in the following format:
```tsv
sample condition batch
sample1 A 1
sample2 A 2
sample3 B 1
sample4 B 2
sample5 C 1
sample6 C 2
```

Where:
- `sample` - name of the sample. Should be first column.
- `condition` - condition of the sample (e.g. control or treatment). Should be second column.
- `batch` - optional batch variable for batch correction. If there is no batch variable, this column should be present but empty 
- Should be at least two samples for each condition.

### Options configuration
The preprocessing of the data is configurable.

- `normalization_method`: ["median", "quantile"] data are first log2 normalized then either median centered (if "median") or quantile normalized using `preprocessCore`.
- `normalization_log_offset`:[float>0] log2 transformation applies an offset `log2(x + offset)` 
- `imputation_method`: ["mindet", "minprob"] missing values (NA or 0) are imputed using either:
  - "mindet" algorithm where each missing value is imputed with the corresponding "minimum" (actually 1st percentile) of intensity values for its sample.
  - "minprob" algorithm where imputed values are sampled from a Gaussian distribution, centered on the "minimum" (see above). The standard deviation of this distribution is set to the median of the standard deviations of the distributions of all proteins.
- `stratify_imputation_by_batch`: [true, false] if true the imputation will be done separately for each batch. In the event that any proteins have less than three non-missing samples in any batch, it will fall back to non-stratified processing.
- `batch_correction_method`: ["combat", "limma", null] if there is a batch variable will perform batch correction with `comBat` (par.prior=True) function of the `sva` package or `removeBatchEffects` from `limma`. **This option will be ignored for this analysis; batch will instead be given as a covariate to the linear model**
- `min_non_missing_fraction`:[float (0<x<=1.)] the minimum proportion of non-missing values of a protein required for it be retained for analysis
- `require_min_fraction_one_class`: [true,false] if true a protein must exceed `min_non_missing_fraction` in only one class, to retained.
- `feature`: string the name of the feature to be analysed, must correspond to a value in the `feature` column of `data.tsv`
-  `model_type`: ["lm", "rfit"] the type of linear model to fit to estimate the effects. See 'Analysis Details' below.
```json
{
  "model_type": "lm",
  "feature": "protein1",
  "normalization_method": "median",
  "normalization_log_offset": 1,
  "imputation_method": "mindet",
  "stratify_imputation_by_batch": false,
  "batch_correction_method": null,
  "min_non_missing_fraction": 0.5,
  "require_min_fraction_one_class": false
}
```
### Analysis Details
#### Model fitting
This script normalises and imputes the data, using all features then fits the full linear model:
```
y ~ condition + [batch]
```
and the reduced model without the effect of condition

```
y ~ [batch]
```
where y is the normalized and imputed values for the selected feature.

Depending on the `model_type` the full model is fit with Ordinary Least-Squares regression (`model_type=lm`) or the robust rank-based regression implemented in the R-package [`Rfit::rfit`](https://doi.org/10.32614/CRAN.package.Rfit) (`model_type=rfit`). In general the reduced model is fit the same way except when `model_type=rfit` and the reduced model is intercept-only. As `Rfit::rfit` cannot fit an intercept-only model we use `lm` in this case.

#### Effect of condition
the results of the overall test of an effect of condition are found in `{proc}/output/overall_test.tsv`
To test for the effect of condition, in general a drop-test is performed comparing the full and reduced models, generating the following output.

`model_type=lm`
```tsv
# Analysis of Variance Table
# 
# Model 1: outcome ~ batch
# Model 2: outcome ~ condition + batch
Res.Df	RSS	Df	Sum of Sq	F	Pr(>F)
15	11.0562204002151				
13	8.63490169926106	2	2.42131870095402	1.82266945291896	0.200553209750671

```
`model_type=rfit`
```tsv

# Rfit::drop.test
#   M1: outcome ~ batch
#   M2: outcome ~ condition + batch
Df	RD	F	p
2,13	1.24018844092793	1.32342051279903	0.299818505107304
```

As stated above Rfit cannot fit an intercept-only model so in the case the reduced model is intercept-only a drop-test cannot be performed, however as there is only one variable in the full model the significance can be derived from the 'summary' object of the full model resulting in the following output

```tsv
# Rfit::summary
#   M1: outcome ~ 1
#   M2: outcome ~ condition
Df	F	p
2,19	1.90593630147862	0.17750720023566
```
#### Follow-up contrasts
Follow-up contrasts are written to `{proc}/output/contrasts.tsv`. These are performed using the [`emmeans`](https://doi.org/10.32614/CRAN.package.emmeans) package and generate the following output table.

```tsv
# Results are averaged over the levels of: batch
# Confidence level used: 0.95
# Conf-level adjustment: tukey method for comparing a family of 3 estimates
# P value adjustment: tukey method for comparing a family of 3 estimates
contrast	estimate	SE	df	lower.CL	upper.CL	t.ratio	p.value
glioblastoma_IDH_wildtype_mesenchymal_type - glioblastoma_IDH_wildtype_RTK1_type	-0.630001592730843	0.515864081520734	13	-1.99210826221771	0.732105076756029	-1.22125500746949	0.462098816787775
glioblastoma_IDH_wildtype_mesenchymal_type - meningioma_benign	0.202643691415377	0.544029851236387	13	-1.23383292192924	1.63912030475999	0.372486346024653	0.926794995899045
glioblastoma_IDH_wildtype_RTK1_type - meningioma_benign	0.83264528414622	0.457131380492942	13	-0.374381381415783	2.03967194970822	1.82145728706777	0.201274523256679
```
#### Output "Values"
`values.tsv` contains the output for downstream visualisation. These are the residuals of the 'reduced' model with mean(y) added, and are in the following format.
```tsv
sample	condition	value
6	meningioma_benign	-0.877498457477708
430	meningioma_benign	-1.49666758469419
2	meningioma_benign	-1.84531037250345
345	meningioma_benign	-1.19302100201005
167	meningioma_benign	-2.90372102659362
128	meningioma_benign	-0.324048877428288
143	meningioma_benign	-1.59681559594367
```

#### If there is only one condition
If there is only one condition present in the input data, no tests will be run, but 'dummy' output files for the overall tests and contrasts will be generated. 

"values.tsv" will be output. This will simply be the normalised and imputed values of the selected feature