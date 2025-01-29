# TriCA

Trinary chart-reviewed phenotype integrated cost-effective augmented estimation

# Outline
1. Description
2. TriCA Overview
3. Example

# Description

This README is for the journal peer review of the TriCA paper, which introduces a method for cost-effective, augmented estimation in association studies. The TriCA method is particularly useful when 'undecided' cases arise during manual chart reviews. It optimally combines binary algorithm-derived phenotypes for the entire cohort with trinary chart-reviewed phenotypes from a small subset, selected through outcome-dependent sampling. This approach offers unbiased estimates with greater efficiency compared to existing methods.

**Statement of significance**

| Summary | Description |
| ------- | ----------- |
| Problem | No methods have been proposed to include ‘undecided’ records from manual chart review phenotypes when identifying risk factors in association studies, particularly in rare disease scenarios|
| What is Already Known | Electronic health records are a valuable resource for identifying risk factors through association studies. While phenotyping algorithms are efficient for obtaining clinical outcomes, they can be error prone. Manual chart review, considered the gold standard, provides unbiased estimates but is labor-intensive and limited to a small subset of patients, potentially introducing ‘undecided’ cases. Existing methods often discard these indeterminate cases, which can reduce the efficiency of estimates, particularly in rare event conditions.|
| What this Paper Adds | We develop an augmented estimator, TriCA, that optimally combines the algorithm-derived binary phenotypes with the chart-review trinary phenotypes selected through a biased sampling strategy. By incorporating the undecided cases from manual chart review, TriCA provides unbiased estimates with higher statistical efficiency compared to existing methods. | 

# TriCA Overview

<img src="Visual Abstract.png" alt="isual Abstract for TriCA method" width="1000">

# Example

file: [Example.R](https://github.com/Penncil/SSL/blob/master/Example.R) 

 

## Dataset

### Full dataset

```R
BT = c(-2.75,1.00,1.00,-0.70,1.70,-0.20) 
df = fn_dataGenF(beta = BT,
                 N = 3000, # sample size
                 pe = 0.4, # P(X2), exposure
                 pa = 0.9, # assuming sp=se
                 seed = 2025) 
```

columns: (Y, S, X)

- Y: outcome/true phenotype, categorical data with 3 levels.
  - 0: No, 1: Yes, 2: Unknown. (0 is the reference level) 
  - p(Y=1) related to the `beta`
- S: surrogate phenotype, categorical data with 2 levels.
  - 0: No, 1: Yes. (0 is the reference level)
  - generated with p(s=1|y=1)=`pa`, p(s=2|y=2)=`pa`, p(s=1|y=2)=0.4
- X = (X1, X2) covariates. 
  - X1 ~ BER(0.5)
  - X2 ~ BER(`pe`): indicate treatment/control



### Validation Dataset

- Indicator : the `y` column of full dataset is set to be `NA` if it's not included in the validation.
- the validation dataset are from the SAME full dataset



#### 1\ Uniform sampling

Idea : sampling uniformly from full dataset

```R
du = fn_dataGenU(Fd = df, # full dataset
                 n = 600, # number of samples in the validation st
                 seed = 2025,
                 threshold=5 # min num of y = 1 (yes) in validation set
                )
```



#### 2\ Biased sampling (Outcome-dependent sampling)

Idea : uniformly select n1 samples from the S-positive(S=1) patients and n0 samples the S-negative(S=0) patients to construct V

```R
dv = fn_dataGenV(Fd = df, # full dataset
                 n0 = 300, # number of samples from S0 (S=0)
                 n1 = 300, # number of samples from S1 (S=1)
                 seed = 2025,
                 threshold=5 # min num of y = 1 (yes) in validation set
                )
```



## Estimators

```R
# oracle: Y~X on F
Roral = fn_oracle(df,c('x1','x2'),'y')

# random sampling (method 1 & 3)
# method 1: Y~X on U
# method 3: augmented on U
Rrand = fn_random(du,c('x1','x2'),'y','s') 

# biased sampling (method 2 & TriCA)
# method 2: Y~X on V
# TriCA: augmented on V
Rbias = fn_bias(dv,c('x1','x2'),'y','s')   
```

