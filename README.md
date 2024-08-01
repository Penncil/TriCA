# TriCA

Trinary chart-reviewed phenotype integrated cost-effective augmented estimation

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




