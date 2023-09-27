rm(list=ls())
library(nnet)
######## 0. Data Generation  ######## 
# 1. settings
seed = 1000
set.seed(seed)

n = 3000 # number of observations s
n0 = 300 # number of samples from S0 (S=1)
n1 = 300 # number of samples from S1 (S=2)
n_v <- n0 + n1

n_vu = 600 # number of observations in validation set
rho = n_vu/n

K = 3 # number of level of Y (1-No,2-Yes,3-Not sure)
k = K-1 # 2

M = 3 # number of covariates (without 1) 
m = M+1 # number of coefficient for covariates

# 2. parameter (k*m)
beta <- c(-3.8, 1,1,1,
          0.5,-0.4,0.6,-1.6) # prevalence = 0.05
beta_matrix <- matrix(beta,nrow = k,ncol = m,byrow = TRUE)

# 3. covariates X (n*m)
x1 <- rnorm(n)
x2 <- rnorm(n)
x3 <- rbinom(n,1,0.5)
X <- cbind(x1,x2,x3)
Xwith1 <- cbind(1,X)

# 4. probability P (n*K)
P <- matrix(nrow = n,ncol = K)
for(i in 1:k){
  P[,(i+1)] <- exp(Xwith1%*%beta_matrix[i,])
}
P <- P/(rowSums(P,na.rm=TRUE) + 1)
P[,1] <- 1-rowSums(P,na.rm=TRUE)

# 5. response variable Y (n*K), y (n*1)
Y <- matrix(0,nrow = n,ncol = K)
for(i in 1:n){
  Y[i,] <- c(rmultinom(1,1,P[i,]))
}
y <- apply(Y,1,which.is.max)
# prop.table(table(y))

# 6. surrogate S (n*1)
# q (K*1) : 
# p(s=0|y=k) for each k -> p(s=1|y=k)
q_01 <- c(0.9,1-0.9,0.6) # first try
# q_01 <- c(0.75,1-0.8,0.6)
# p(s=1|y=k) for each k -> p(s=2|y=k)
q_12 <- 1 - q_01 
# p(s=1) for each i
p_s1 <- rep(NA,n)
for(i in 1:K){
  p_s1[y==i] <- q_12[i]
}
S <- rbinom(n,1,p_s1) + 1 # (0,1) -> (1,2)
table(S,y) # check

# 7. weight W (n*1)
# n_s0 <- sum(S==1)
# n_s1 <- sum(S==2)
# h0 <- n0/n_s0
# h1 <- n1/n_s1
# W <- ifelse(S==1,h0,h1)

# 5. format & combine
data <- data.frame(cbind(y,S,X,W))
data$x3 <- factor(data$x3)

data$y <- factor(data$y) 
data$y <- relevel(data$y,ref=1)

data$S <- factor(data$S)
data$S <- relevel(data$S,ref=1)

# 6. validation data
id_s0 <- as.numeric(rownames(data[data$S==1,]))
id_s1 <- as.numeric(rownames(data[data$S==2,]))
id_v <- c(sample(id_s0,n0,replace = F),
          sample(id_s1,n1,replace = F))

data_fv <- data
data_fv[-id_v,'y'] <- NA
data_fv$W <- NULL
summary(data_fv)
write.csv(data_fv,"data_outcome_dependent_sampled.csv",row.names = FALSE)


id_vu <- sample(1:n,n_v,replace = F)
data_fvu <- data
data_fvu$W <- NULL
data_fvu[-id_vu,'y'] <- NAs
write.csv(data_fvu,"data_uniform_sampled.csv",row.names = FALSE)
