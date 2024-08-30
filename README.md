# TriCA

Trinary chart-reviewed phenotype integrated cost-effective augmented estimation

# Outline
1. Description
2. TriCA Overview
3. Example Dataset

# Description

This README is for the journal peer review of the TriCA paper, which introduces a method for cost-effective, augmented estimation in association studies. The TriCA method is particularly useful when 'undecided' cases arise during manual chart reviews. It optimally combines binary algorithm-derived phenotypes for the entire cohort with trinary chart-reviewed phenotypes from a small subset, selected through outcome-dependent sampling. This approach offers unbiased estimates with greater efficiency compared to existing methods.

To demonstrate its effectiveness (as discussed in the SIMULATION STUDY section of the paper), we use simulated data from 3,000 patients to analyze the association of two continuous covariates and one treatment indicator.

**Statement of significance**

| Summary | Description |
| ------- | ----------- |
| Problem | No methods have been proposed to include ‘undecided’ records from manual chart review phenotypes when identifying risk factors in association studies, particularly in rare disease scenarios|
| What is Already Known | Electronic health records are a valuable resource for identifying risk factors through association studies. While phenotyping algorithms are efficient for obtaining clinical outcomes, they can be error prone. Manual chart review, considered the gold standard, provides unbiased estimates but is labor-intensive and limited to a small subset of patients, potentially introducing ‘undecided’ cases. Existing methods often discard these indeterminate cases, which can reduce the efficiency of estimates, particularly in rare event conditions.|
| What this Paper Adds | We develop an augmented estimator, TriCA, that optimally combines the algorithm-derived binary phenotypes with the chart-review trinary phenotypes selected through a biased sampling strategy. By incorporating the undecided cases from manual chart review, TriCA provides unbiased estimates with higher statistical efficiency compared to existing methods. | 

# TriCA Overview

<img src="Visual Abstract.png" alt="isual Abstract for TriCA method" width="600">


# Example Dataset

file: [data_outcome_dependent_sampled.csv](https://github.com/Penncil/SSL/blob/master/data_outcome_dependent_sampled.csv) and [data_uniform_sampled.csv](https://github.com/Penncil/SSL/blob/master/data_uniform_sampled.csv)



## Full Dataset

- 3000 rows
- 6 columns: (Y, S, X)
  - Y: outcome/true phenotype, categorical data with 3 levels.
    - 1: No, 2: Yes, 3: Unknown. (1 is the reference level)
    - p(Y=2) ~ 5%
    - generate from X with parameter beta = (-3.8, 1, 1, 1, 0.5, -0.4, 0.6, -1.6)
  - S: surrogate phenotype, categorical data with 2 levels.
    - 1: No, 2: Yes. (1 is the reference level)
    - generated with p(s=1|y=1)=0.9, p(s=2|y=2)=0.9, p(s=1|y=3)=0.6
  - X = (X1, X2, X3) covariates. 
    - X1, X2 ~ N(0,1)
    - X3 ~ BER(0.5): indicate treatment/control



## Validation Dataset

- Indicator : the `y` column of full dataset is set to be `NA` if it's not included in the validation.
- the validation dataset are from the SAME full dataset



### 1 : uniform

- Idea : sampling uniformly from full dataset

- n_vu = 600 : number of observations in uniform validation set
- rho = n_v/ n = 1/5



### 2 : outcome dependent

- Idea : uniformly select n1 samples from the S-positive(S=2) patients and n0 samples the S-negative(S=1) patients to construct V
- n_v = 600
  - n0 = 300 :  number of samples from S0 (S=1)
  - n1 = 300 : number of samples from S1 (S=2)
- W: weights [calculated when estimating `Aug.Biased`]
  - for S = 1, w = n0 / ns0 (# of sample from S=1 / # of obs with S=1) where n0 = 300
  - for S = 2, w = n1 / ns1 (# of sample from S=2 / # of obs with S=2) where n1 = 300




