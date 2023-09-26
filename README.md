# SSL
semi-supervised learning

# Example Dataset
file: `example_full_data.csv`

- 3000 rows
- 7 columns: (Y, S, X, W)
  - Y: outcome/true phenotype, categorical data with 3 levels.
    - 1: No, 2: Yes, 3: Unknown. (1 is the reference level)
    - p(Y=1) ~ 5%
    - generate from X with parameter beta = (-3.8, 1, 1, 1, 0.5, -0.4, 0.6, -1.6)
  - S: surrogate phenotype, categorical data with 2 levels.
    - 1: No, 2: Yes. (1 is the reference level)
    - generated with p(s=1|y=1)=0.9, p(s=2|y=2)=0.9, p(s=1|y=3)=0.6
  - X = (1, X1, X2, X3) covariates. 
    - X1, X2 ~ N(0,1)
    - X3 ~ BER(0.5): indicate treatment/control
  - W: weights
    - for S = 1, w = n0 / ns0 (# of sample from S=1 / # of obs with S=1) where n0 = 300
    - for S = 2, w = n1 / ns1 (# of sample from S=2 / # of obs with S=2) where n1 = 300


