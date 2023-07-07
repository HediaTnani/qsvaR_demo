# qSVAR_Demo
Hedia Tnani

## qsvaR

Surrogate Variable Analysis (SVA) is a statistical method used to
identify and account for sources of noise in high-dimensional data sets,
particularly those arising from genomics and transcriptomics studies.
It’s often used in the context of gene expression data, where hidden
sources of variability can distort the relationships between genes and
the biological conditions under investigation.

SVA seeks to identify these hidden sources of variation (surrogate
variables) so that they can be included in subsequent statistical
modeling processes. These surrogate variables serve as proxies for
unknown, unmodeled, or unmeasured sources of noise, capturing the
unmodeled variability in the data. After identifying these variables,
they can be controlled for in downstream analyses, improving the ability
to detect true biological signals.

This is a basic example which shows how to obtain the quality surrogate
variables (qSVs) for `BrainSeq Phase 2`.

## BSP2 Dataset

RiboZero RNA-seq data for 900 samples across both the dorsolateral
prefrontal cortex (DLPFC) and the hippocampus (HIPPO) for 551
individuals (286 affected by schizophrenia disorder: SCZD).
Visit [here](https://eqtl.brainseq.org/phase2/) for more details.

## Load the data

For illustrative purposes, we'll use the BSP2 data at the transcript
level which reside in a `RangedSummarizedExperiment` (RSE) object
called `rse_tx`.

``` r
## Load the container package for this type of data
library("SummarizedExperiment")
library("tidySummarizedExperiment")

## Load the rse_tx object
load(url("https://s3.us-east-2.amazonaws.com/libd-brainseq2/rse_tx_unfiltered.Rdata"), verbose=T)
# Loading objects:
#   rse_tx


## General overview of the object
rse_tx
# A SummarizedExperiment-tibble abstraction:
#   178,283,700 × 86
# Features=198093 | Samples=900 | Assays=tpm
#    .feature        .sample   tpm SAMPLE_ID FQCbasicStats perBaseQual perTileQual
#    <chr>           <chr>   <dbl> <list>    <list>        <list>      <list>     
#  1 ENST0000045632… R10424   0    <chr [1]> <chr [1]>     <chr [1]>   <chr [1]>  
#  2 ENST0000045030… R10424   0    <chr [1]> <chr [1]>     <chr [1]>   <chr [1]>  
#  3 ENST0000048814… R10424   2.36 <chr [1]> <chr [1]>     <chr [1]>   <chr [1]>  
#  4 ENST0000061921… R10424   0    <chr [1]> <chr [1]>     <chr [1]>   <chr [1]>  
#  5 ENST0000047335… R10424   0    <chr [1]> <chr [1]>     <chr [1]>   <chr [1]>  
#  6 ENST0000046928… R10424   0    <chr [1]> <chr [1]>     <chr [1]>   <chr [1]>  
#  7 ENST0000060709… R10424   0    <chr [1]> <chr [1]>     <chr [1]>   <chr [1]>  
#  8 ENST0000041732… R10424   0    <chr [1]> <chr [1]>     <chr [1]>   <chr [1]>  
#  9 ENST0000046146… R10424   2.96 <chr [1]> <chr [1]>     <chr [1]>   <chr [1]>  
# 10 ENST0000060685… R10424   0    <chr [1]> <chr [1]>     <chr [1]>   <chr [1]>  
# ℹ 40 more rows
# ℹ 79 more variables: perSeqQual <list>, perBaseContent <list>,
#   GCcontent <list>, Ncontent <list>, SeqLengthDist <list>,
#   SeqDuplication <list>, OverrepSeqs <list>, AdapterContent <list>,
#   KmerContent <list>, percentGC_R1 <list>, phred100_R1 <list>,
#   phredGT30_R1 <list>, phredGT35_R1 <list>, Adapter88_R1 <list>,
#   percentGC_R2 <list>, phred100_R2 <list>, phredGT30_R2 <list>, …
# ℹ Use `print(n = ...)` to see more rows
```

This is a `SummarizedExperiment` object but it is evaluated as a
`tibble`. So it is fully compatible both with SummarizedExperiment and
`tidyverse` APIs.

## SummarizedExperiment

The `SummarizedExperiment` package contains two classes:
`SummarizedExperiment` and `RangedSummarizedExperiment`.

`RangedSummarizedExperiment` is a subclass of the SummarizedExperiment
class, meaning it inherits all the features and functionalities of its
parent class. The primary difference is that the rows of a
RangedSummarizedExperiment object represent genomic ranges of interest,
rather than a DataFrame of features. Here’s what you need to know about
RangedSummarizedExperiment:

-   `Genomic Range Representation`: RangedSummarizedExperiment objects
    use a `GRanges` or `GRangesList` object to represent genomic ranges.
    These objects store information about genomic coordinates, genomic
    features, and annotations for the ranges of interest. The
    `rowRanges()` function allows you to access and manipulate the
    genomic range information within a `RangedSummarizedExperiment`
    object.

    The following graphic displays the class geometry and highlights the
    vertical (column) and horizontal (row) relationships.

![](https://f1000research.s3.amazonaws.com/manuscripts/19158/eedc82a4-df66-4228-9d52-90b482b3df47_figure1.gif)

``` r
ncol(assays(rse_tx)$tpm) == nrow(colData(rse_tx))
# [1] TRUE
```

## Data overview

The dataset `rse_tx` contains the following assays:

``` r
assays(rse_tx)
# List of length 1
# names(1): tpm
```

`tpm`: normalized read counts of the 19,8093 genes across 900 samples .

``` r
library(tidyverse)

# Transform 'tpm' data to a long format with 'Transcript', 'Sample', and 'TPM' columns
tpm_tb <- log2(assays(rse_tx)$tpm + 0.5) %>% as.data.frame() %>%
  rownames_to_column(var = "Transcript") %>%
  pivot_longer(
    cols = -Transcript,
    names_to = "Sample",
    values_to = "TPM"
  )

# Extract 'RNum', 'Region', and 'Dx' from 'rse_tx', convert to tibble
colData_tb <- colData(rse_tx) %>% 
  as.data.frame() %>% 
  dplyr::select(RNum, Region, Dx) %>% as_tibble()

# Join 'tpm_tb' with 'colData_tb' based on 'Sample' and 'RNum'
tpm_tb <- tpm_tb %>%
  dplyr::left_join(colData_tb, by = c("Sample" = "RNum"))

# Calculate mean TPM for each Dx and Region
mean_tpm <- tpm_tb %>%
  group_by(Dx, Region) %>%
  summarise(
    mean_TPM = mean(TPM, na.rm = TRUE),
    .groups = "drop" # This drops the grouping
  )

#  Dx      Region mean_TPM
#   <chr>   <chr>     <dbl>
# 1 Control DLPFC    0.0600
# 2 Control HIPPO   -0.0385
# 3 Schizo  DLPFC    0.0850
# 4 Schizo  HIPPO   -0.0567

# Load ggplot2 for plotting
library(ggplot2)

# Create boxplot for first 100,000 rows, without outliers, split by Dx and Region
tpm_tb %>% slice(1:100000) %>%
  ggplot(aes(x = Dx, y = TPM, fill = Dx)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~ Region) +
  theme_minimal() +
  labs(title = "Boxplot of TPM for Each Dx by Region", x = "Dx", y = "TPM")

# Create a histogram of TPM for each Dx
tpm_tb %>% slice(1:100000) %>%
  ggplot(aes(x = TPM, fill = Dx)) +
  geom_histogram(alpha = 0.5, bins = 30, position = "identity") +
  facet_wrap(~ Dx) +
  theme_minimal() +
  labs(title = "Histogram of TPM for Each Dx", x = "TPM", y = "Count")

# Create a density plot of TPM for each Dx
tpm_tb %>% slice(1:100000) %>%
  ggplot(aes(x = TPM, fill = Dx)) +
  geom_density(alpha = 0.5) +
  
  facet_wrap(~ Dx) +
  theme_minimal() +
  labs(title = "Density Plot of TPM for Each Dx", x = "TPM", y = "Density")
```

## Data normalization

`Transcripts Per Million (TPM)` is a normalization method for RNA-seq
data that facilitates comparisons both within and between samples.
`It adjusts for gene length and sequencing depth, making it a helpful measure for relative gene expression`.

Here is the formula for calculating TPM:

1.  **RPK (Reads Per Kilobase):** For each gene, divide the number of
    reads by the length of the gene in kilobases. This normalizes for
    gene length.

    **`RPK = (Number of reads mapped to a gene / Length of the gene in kilobases)`**

2.  **RPKM or FPKM (Reads Per Kilobase of transcript per Million mapped
    reads):** Divide the RPK values by the per million scaling factor.
    This normalizes for sequencing depth.

    **`RPKM/FPKM = RPK / (Total number of mapped reads in the experiment / 1,000,000)`**

3.  **TPM (Transcripts Per Million):** Finally, to compute TPM, the
    RPKM/FPKM values are divided by the sum of all RPKM/FPKM values in
    the sample and then multiplied by 1,000,000.

    **`TPM = (RPKM/FPKM for a gene / Sum of RPKM/FPKM for all genes) x 1,000,000`**

Here’s how you can calculate TPM in R given a data frame **`counts`** of
raw counts and a vector **`lengths`** of gene lengths:

``` r
 # Assuming counts is a matrix with raw counts and lengths is a vector with gene lengths

# Step 1: Calculate RPK
rpk <- counts / (lengths / 1000)

# Step 2: Calculate the per-million scaling factor
scaling_factor <- colSums(rpk) / 1e6

# Step 3: Calculate TPM
tpm <- rpk / scaling_factor
```

There are many ways to access the `tpm` data. Below some of them:

``` r
# Load the magrittr library for use of the pipe operator
library(magrittr)

# Extract the first element from the list of assays in the SummarizedExperiment object
# This assumes that the first element in the list of assays is the one you are interested in
assays(rse_tx) %>% extract2(1) 

# Another way of extracting the first element from the list of assays in the SummarizedExperiment object
# Here, 'tpm' is assumed to be the name of the first assay
assays(rse_tx)$tpm 

# Load the purrr and dplyr libraries for functional programming and data manipulation, respectively

library(purrr)
library(dplyr)

# Yet another way of extracting the first element from the list of assays in the SummarizedExperiment object
# The simplify() function tries to simplify the result to a vector or matrix if possible
# The first() function extracts the first element of the simplified result
assays(rse_tx) %>%
simplify() %>%
first()
```

### Sample data

`SAMPLE_ID` : is the name of the sample.  
`ERCCsumLogErr` : a summary statistic quantifying overall difference of
expected and actual ERCC concentrations for one sample. For more
about *ERCC* check their product page
at <https://www.thermofisher.com/order/catalog/product/4456740>.
trimmed  
numReads numMapped  
numUnmapped  
`overallMapRate`: the decimal fraction of reads which successfully
mapped to the reference genome (i.e. *numMapped* / *numReads*).
concordMapRate  
`totalMapped`: the number of reads which successfully mapped to the
canonical sequences in the reference genome (excluding mitochondrial
chromosomes).  
`mitoMapped`: the number of reads which successfully mapped to the
mitochondrial chromosome. `mitoRate` : the decimal fraction of reads
which mapped to the mitochondrial chromosome, of those which map at all
(i.e. *mitoMapped* / (*totalMapped* + *mitoMapped*))  
`totalAssignedGene` : the decimal fraction of reads assigned
unambiguously to a gene (including mitochondrial genes),
with `featureCounts` (Liao et al. 2014), of those in total.  
`rRNA_rate` : the decimal fraction of reads assigned to a gene whose
type is 'rRNA', of those assigned to any gene. `RNum` :`BrNum` : stands
for Brain Number. It likely refers to a unique identifier assigned to
each brain (or brain sample) examined in the study.  
`Region` : refers to a specific region of the brain from which the
sample was taken. `RIN` : refers to RNA Integrity Number. It’s a measure
of the quality of the RNA samples used in the study. A higher RIN
indicates higher quality, and thus more reliable, RNA.  
`Age` : refers to the age of the donors from which the samples were
taken.  
`Sex` : refers to the sex of the donors from which the samples were
taken. if the sample comes from a female (F) or a male (M) donor. `Race`
: refers to the race of the donors from which the samples were taken.  
`Dx` : refers to the diagnosis of the donors from which the samples were
taken.

## Variance Partition Analysis

Evaluating the correlation among sample variables is an imperative step
in our analysis. The presence of highly correlated variables can
significantly impact the robustness of our model, creating unstable
estimates of variance fractions. This instability is detrimental as it
can obfuscate the identification of the true contributors to expression
variation.

When two or more variables are highly correlated, they essentially
contain redundant information. This redundancy can lead to what is known
as multicollinearity, a situation that can greatly destabilize our
model, potentially leading to misleading or erroneous results.
Therefore, assessing correlation enables us to prevent
**multicollinearity**, ensuring a more accurate and reliable analysis.

In essence, our objective is not just to recognize the factors
contributing to expression variation, but to accurately measure their
individual impacts. Understanding these relationships will allow us to
develop a more accurate and reliable predictive model. Let’s first
unlist `RIN`, `mitoMapped` and `totalAssignedGene` .

``` r
library(purrr)
colData(rse_tx)$RIN = as.numeric(unlist(colData(rse_tx)$RIN))
colData(rse_tx)$mitoMapped = unlist(map(colData(rse_tx)$mitoMapped,1))
colData(rse_tx)$totalAssignedGene = unlist(map(colData(rse_tx)$totalAssignedGene,1))
```

We’ll be using the **`canCorPairs`** function from the
**`variancePartition`** which performs
`canonical correlation analysis (CCA)` between pairs of variables.

`Canonical correlation` analysis is a method used to identify and
measure the associations among two sets of variables. It’s used to find
the correlation between two sets of multidimensional variables,
identifying the dimensions that provide the highest correlation between
the two sets.

In the context of the **`variancePartition`** package, the
**`canCorPairs`** function computes the canonical correlations between
each pair of variables specified in the formula argument.

``` r
library(variancePartition)
library(ComplexHeatmap)
library(colorRamp2)
formula <- ~Dx + RIN+ Age + Sex + Race + mitoMapped+ totalAssignedGene
correlation_matrix <- canCorPairs(formula, colData(rse_tx))
# Create the heatmap
# Create a color function from green to purple
color_func <- colorRamp2(c(-1, 0, 1), c("darkgreen", "white", "purple"))

# Create the heatmap
Heatmap(correlation_matrix,
        name = "Correlation",      # Title for the color key
        show_row_names = TRUE,      # Show row names
        show_column_names = TRUE,   # Show column names
        cluster_rows = FALSE,       # Do not cluster rows
        cluster_columns = FALSE,    # Do not cluster columns
        cell_fun = function(j, i, x, y, width, height, fill) {
                    grid.text(sprintf("%.2f", correlation_matrix[i, j]), x, y,
                              gp = gpar(fontsize = 10))
                    },
        col = color_func, # Use the custom color function
        heatmap_legend_param = list(title = "Correlation", # legend title
                                    color_bar = "continuous", # type of color bar
                                    legend_direction = "horizontal" # direction of legend
                                    )
)
```

The **`cell_fun`** function customizes the display of each cell in the
heatmap by adding the correlation value as text inside the cell with a
specific format and font size.

From the heatmap we can see that there’s a strong positive correlation
between `mitoMapped` and `totalAssignedGene`. One of the strategies to
address `multicollinearity` is removing some of the highly correlated
independent variables.

## Fit model and extract fraction of variance explained

`variancePartition` fits two types of models:

1.  linear mixed model where all `categorical` variables are modeled as
    `random` effects and all `continuous` variables are `fixed` effects.
    The function `lmer` from `lme4` is used to fit this model.

2.  `fixed` effected model, where all variables are modeled as `fixed`
    effects. The function `lm` is used to fit this model.

``` r
library(tidyverse)
library(matrixStats) # for rowVars function

txExpr = log2(assays(rse_tx)$tpm + 0.5)

# Identify the genes with zero variance
txExpr <- txExpr %>%
  as.matrix() %>%
  {.[rowVars(.) != 0, ]}

form <- ~ (1|Dx) + RIN + Age + (1|Sex) + (1|Race) + totalAssignedGene
varPart <- fitExtractVarPartModel(txExpr, form, colData(rse_tx))
```

## qsvaR: Quality Surrogate Variable Analysis

### Degradation transcripts

After loading the `rse_tx` object we subset for the transcripts
associated with degradation. In `qsvaR` there are three models to choose
from. Here the names `"cell_component"`, `"top1500"`, and `"standard"`
refer to models that were determined to be effective in removing
degradation effects.

### qsvaR models

-   **standard model**

The `"standard"` model involves taking the union of the top 1000
transcripts associated with degradation from the interaction model and
the main effect model.

-   **top1500 model**

The `"top1500"` model is the same as the `"standard"` model except the
union of the top 1500 genes associated with degradation is selected.

-   **cell_component model**

The `"cell_component"` model, involved deconvolution of the degradation
matrix to determine the proportion of cell types within our studied
tissue. These proportions were then added to our `model.matrix()` and
the union of the top 1000 transcripts in the interaction model, the main
effect model, and the cell proportions model were used to generate this
model of qSVs. It’s the most effective of our models.

## Run qSVA

In this example we will choose `"cell_component"` when using the
`getDegTx()` and `select_transcripts()` functions. `getDegTx()` is used
to select a specific set of transcripts from a given experiment, and
provides a warning if the mean expression level of these transcripts is
low.

``` r
library(qsvaR)
## Next we get the degraded transcripts for qSVA from the "cell_component"
## model
DegTx <- getDegTx(rse_tx, type = "cell_component")
dim(DegTx)
# [1] 2976  900
```

The `getPCs` applies PCA using the **`prcomp`** function.

``` r
## Now we can compute the Principal Components (PCs) of the degraded
## transcripts
pcTx <- getPCs(DegTx, "tpm")
```

We use a model matrix that accounts for relevant variables. We include
variables such as Age, Sex, Race, and Region in this model. Our dataset,
**`DegTx`**, a RangedSummarizedExperiment object, is used as the
**`rse_tx`** input.

The **`model.matrix()`** function takes a formula and a data frame (or
similar object) and returns the design matrix.

In this case, the formula is **`~ Dx + Age + Sex + Race + Region`**,
which specifies a model where the response variable (typically
represented on the left side of the tilde **`~`**, but it’s omitted here
as the function is used to prepare predictors for further analysis) is
being modeled as a function of the predictors **`Dx`**, **`Age`**,
**`Sex`**, **`Race`**, and **`Region`**. These predictors can be
anything - demographic factors, experimental conditions, etc.

**`model.matrix()`** will automatically create dummy variables for each
level of the factor. Dummy variables are binary (0/1) variables that
indicate the presence of a categorical level for each observation.

**Intercept or not intercept?**

For a single explanatory variable, which we simply call `var`, a design
matrix can be coded in two ways:

-   `model.matrix(~var)` : with intercept term (parametrized for
    a *mean-reference model* for factor variables)

-   `model.matrix(~0+var)` without intercept (parametrized for a *means
    model* for factor variables).

    **Note:** this also makes it easier to define *contrast groups*.

One of the most fundamental concepts in the coding of design matrices is
to understand when one should include an intercept term, when not to,
and how it affects the underlying model.  
If variable `var` is a **factor**, then the two models with and without
the intercept term are **equivalent**, but if variable is a
**covariate** the then two models are fundamentally different.

The `model.matrix` command takes the variables in the `formula` format:

`model.matrix([target] ~ [predictor / features], data = [data source])`

**`data = colData(rse_tx)`**: This argument tells R which dataset to use
for creating the design matrix. In this case, **`colData(rse_tx)`**
seems to be a function call that extracts some kind of column data from
an object **`rse_tx`**.

``` r
## Using a simple statistical model we determine the number of PCs needed (k)
mod <- model.matrix(~ Dx + Age + Sex + Race + Region,
    data = colData(rse_tx)
)

# check the names of the columns in the mod matrix
colnames(mod)
# [1] "(Intercept)" "DxSchizo"    "Age"         "SexM"        "RaceAS"     
# [6] "RaceCAUC"    "RaceHISP"    "RegionHIPPO"
```

**`~ Dx + Age + Sex + Race + Region`**: This part is called a formula.
The tilde (**`~`**) can be read as “is modeled as”. The variables on the
right side of the tilde are the predictor variables (or features). In
this case, **`Dx`**, **`Age`**, **`Sex`**, **`Race`**, and **`Region`**
are your predictors.

The **`k_qsvs()`** function is employed to determine the required number
of Principal Components (PCs) needed to capture the underlying variation
in our data.

``` r
# Get the number of qsvs
set.seed(20230703)
k <- k_qsvs(DegTx, mod, "tpm")
print(k)
# [1] 34
```

You can simplify the process by utilizing our integrated function,
**`qSVA`**. This wrapper function conveniently combines all previously
mentioned functions, streamlining the analysis into a single step.

``` r
## Example use of the wrapper function qSVA()
set.seed(20230703)
qsva_pcs_cc <- qSVA(rse_tx, type = "cell_component", mod = mod, assayname = "tpm")
dim(qsva_pcs_cc)
#> [1] 900  34
```

We can try other models such as the `standard` or `"top1500"` .

``` r
qsva_pcs_standard <- qsvaR::qSVA(rse_tx, type = "standard", mod = mod, assayname = "tpm")
dim(qsva_pcs_standard)
# [1] 900  19
```

## Differential Expression using limma

*`limma`* fits a so-called linear model; examples of linear models are
(1) linear regression, (2) multiple linear regression and (3) analysis
of variance. The limma User's Guide from the [limma
webpage](http://bioconductor.org/packages/limma).

In the context of **`limma`**, each gene’s expression level is modeled
as a linear combination of the experimental factors, with some added
noise or error term.

For example, consider a simple model for one gene:

*e**x**p**r**e**s**s**i**o**n* = *i**n**t**e**r**c**e**p**t* + *b**e**t**a* \* *t**r**e**a**t**m**e**n**t* + *e**r**r**o**r*

Here:

-   **`expression`** is the expression level of the gene.

-   **`intercept`** is the base level of gene expression (i.e., in the
    control group).

-   **`beta`** is a coefficient representing the effect size of the
    treatment.

-   **`treatment`** is a variable representing whether the sample is a
    control or treatment (often coded as 0 for control, 1 for
    treatment).

-   **`error`** is the noise or unexplained variability in gene
    expression.

``` r
library(limma)
library(qsvaR)

## Add the qSVs to our statistical model
mod_qSVA <- cbind(
    mod,
    qsva_pcs_cc
)
```

**`mod_qSVA`** includes both our original model and the `qSVs`. The
`qSVs` adjust for the degradation effect.

``` r
## Extract the transcript expression values and put them in the
## log2(TPM + 1) scale
txExprs <- log2(assays(rse_tx)$tpm + 0.5)

## Run the standard linear model for differential expression
fitTx <- lmFit(txExprs, mod_qSVA)
eBTx <- eBayes(fitTx)
```

In this part above of the code, we are running a linear model for
differential expression analysis using the **`lmFit()`** function from
the **`limma`** package in R. The function takes two main arguments: the
matrix of expression values (**`txExprs`**) and the design matrix
(**`mod_qSVA`**). The design matrix specifies the linear model that
you’re fitting to the data. In this case, **`mod_qSVA`** likely includes
the variables of interest and the qSVs.

The **`lmFit()`** function fits a separate linear model to each gene (or
transcript, in this case), and the result (**`fitTx`**) is an object
that contains the fitted model coefficients and other statistics for
each gene.

After fitting the model, **`eBayes()`** function is used to compute
empirical Bayes statistics for differential expression. The function
applies a shrinkage technique to the standard errors, improving the
stability of the estimates, particularly when the number of samples is
small. The result (**`eBTx`**) is an object containing the moderated
t-statistics, log-odds of differential expression, and other quantities
for each gene.

``` r
colnames(eBTx)
#  [1] "(Intercept)" "DxSchizo"    "Age"         "SexM"        "RaceAS"     
#  [6] "RaceCAUC"    "RaceHISP"    "RegionHIPPO" "qSV1"        "qSV2"       
# [11] "qSV3"        "qSV4"        "qSV5"        "qSV6"        "qSV7"       
# [16] "qSV8"        "qSV9"        "qSV10"       "qSV11"       "qSV12"      
# [21] "qSV13"       "qSV14"       "qSV15"       "qSV16"       "qSV17"      
# [26] "qSV18"       "qSV19"       "qSV20"       "qSV21"       "qSV22"      
# [31] "qSV23"       "qSV24"       "qSV25"       "qSV26"       "qSV27"      
# [36] "qSV28"       "qSV29"       "qSV30"       "qSV31"       "qSV32"      
# [41] "qSV33"       "qSV34"  
```

Let’s do the differential expression analysis.

``` r
## Extract the differential expression results
sigTx <- topTable(eBTx,
    coef = "DxSchizo",
    p.value = 1, number = nrow(rse_tx)
)
```

**`coef = "DxSchizo"`** specifies the model coefficient of interest. The
**`topTable`** function will return results for this specific term.

Let’s explore the top results.

``` r
## Explore the top results
head(sigTx)
```

Here’s a breakdown of what each column in this output represents:

1.  **logFC (Log Fold Change)**: This represents the log2 fold-change
    for the differential expression between two conditions. A positive
    value indicates higher expression in the second condition, while a
    negative value indicates higher expression in the first condition.
    For example, a logFC of -0.07856532 for ENST00000553142.5 means this
    transcript has lower expression in the second condition compared to
    the first.

2.  **AveExpr (Average Expression)**: This is the average log2
    expression level of the gene across all samples. It provides a
    general idea of how abundantly the gene is expressed.

3.  **t (t-statistic)**: This is the t-statistic from the differential
    expression test. It is a measure of how different the expression
    levels of the gene are between the two conditions, in units of
    standard error. A higher absolute value of the t-statistic means the
    gene’s expression is more significantly different between the
    conditions.

4.  **P.Value (p-value)**: This is the p-value from the differential
    expression test. It quantifies the statistical evidence against the
    null hypothesis that the gene’s expression levels are the same in
    the two conditions. A lower p-value suggests stronger evidence
    against the null hypothesis.

5.  **adj.P.Val (Adjusted p-value)**: This is the p-value adjusted for
    multiple testing using the method of Benjamini and Hochberg, also
    known as the false discovery rate (FDR). It helps control for false
    positives when testing multiple hypotheses (genes) at once.
    Typically, genes with an adjusted p-value (or FDR) less than 0.05 or
    0.01 are considered statistically significantly differentially
    expressed.

6.  **B (B-statistic or log-odds of differential expression)**: This is
    the B-statistic from the empirical Bayes analysis, which can be
    interpreted as the log-odds that the gene is differentially
    expressed. A higher B-statistic means the gene is more likely to be
    differentially expressed.

Each row in this table corresponds to a different transcript (denoted by
the ENST IDs). These are likely the top differentially expressed
transcripts as ranked by the t-statistic, p-value, adjusted p-value, or
B-statistic.

## DEqual plots

This R function, **`DEqual`**, is designed to compare degradation
t-statistics with the t-statistics from a differential expression (DE)
analysis. It takes as input a data frame (**`DE`**) containing
differential expression results, such as those outputted by the
**`limma`** package’s **`topTable`** function. The function identifies
the transcripts that are common between the **`DE`** dataframe and
**`degradation_tstats`** from the **`qsvaR`** package. The function
creates a new dataframe, **`common_data`**, which contains the
degradation t-statistics and DE t-statistics for the common transcripts.

``` r
## Generate a DEqual() plot using the model results with qSVs
DEqual(sigTx)
```

The color of each bin represents the number of data points (in this
case, transcripts) within that bin. This is set by the
**`scale_fill_continuous`** function, which maps the continuous number
of data points in each bin to a continuous color scale.

The **`type = "viridis"`** argument specifies the use of the viridis
color scale. The viridis color scale ranges from dark purple (for low
values) through blue, green, yellow, and up to bright yellow (for high
values).

In this plot:

-   **Dark purple to blue bins**: These colors represent bins with a low
    number of data points. If a bin is dark purple or blue, it means
    there are relatively few transcripts with those particular DE and
    degradation t-statistics.

-   **Green to yellow bins**: These colors represent bins with a
    moderate number of data points. If a bin is green or yellow, it
    means there are a moderate number of transcripts with those
    particular DE and degradation t-statistics.

**Bright yellow bins**: These colors represent bins with a high number
of data points. If a bin is bright yellow, it means there are many
transcripts with those particular DE and degradation t-statistics.

For comparison, here is the `DEqual()` plot for the model without qSVs.

``` r
## Generate a DEqual() plot using the model results without qSVs
DEqual(topTable(eBayes(lmFit(txExprs, mod)), coef = "DxSchizo", p.value = 1, number = nrow(rse_tx)))
```

The results from the differential expression analysis without qSVA
normalization are shown in the two DEqual plots. The contrast between
the two plots is evident and reveals the effectiveness of degradation
effect removal.

In the first plot, the correlation stands at -0.014, suggesting a
minimal relationship between degradation and the differential expression
results. This indicates that the degradation effects have been
effectively eliminated from the data, leading to a cleaner
interpretation of the gene expression differences.

Conversely, the second plot reveals a correlation of 0.5 with the
degradation experiment, even after adjusting for several common
variables. Such a high correlation underscores the fact that
degradation’s signal remains significantly intertwined with our data.
This can potentially introduce confounding factors in our differential
expression analysis between the schizophrenia (SCZD) and neurotypical
control groups.

In essence, the high correlation in the second plot could be masking the
true differential gene expression between the case-control groups by
mingling the degradation signal into the expression data. This
emphasizes the need for adequate normalization, such as qSVA, to reduce
such confounding factors and ensure a more accurate representation of
the biological differences in gene expression.

## References

1.  https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/designmatrices.html
2.  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7873980/
3.  https://bioconductor.org/help/course-materials/2019/BSS2019/04_Practical_CoreApproachesInBioconductor.html
4.  https://stemangiola.github.io/tidySummarizedExperiment/
