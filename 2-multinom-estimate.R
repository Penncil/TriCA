# rm(list=ls())
# -------------- Library  ---------------
library(geex)
library(mvtnorm)
library(sandwich)
library(Matrix)
library(nnet)
library(dplyr)
library(matlib)
require(foreign)
require(reshape2)
library(mltools)
library(data.table)
library(stringr)

# -------------- Functions  ---------------


Aug.Unif.estimate <- function(data, 
                              formula_S = S~x1+x2+x3,
                              formula_Y = y~x1+x2+x3){
  #### 0. Parameters & Validation Data ####
  # data = data_fvu
  # formula_S = S~x1+x2+x3
  # formula_Y = y~x1+x2+x3
  
  n = nrow(data)
  data_vu <- data[!is.na(data$y),]
  n_vu <- nrow(data_vu)
  rho = n_vu / n
  
  K = nlevels(data$y) 
  k = K-1
  
  
  #### 1.1 [Aug-Unif] #### 
  # 1.1.1 (svu) S~X in V-unif -> gamma_v_u & Si-Pi
  model_svu <- glm(formula_S,data = data_vu,family="binomial",x=TRUE,y=TRUE)
  sum_svu <- summary(model_svu)
  gamma_vu <- sum_svu$coefficients[,1]
  # gamma_vu_sd <-  sum_svu$coefficients[,2]
  
  p_svu <- model_svu$fitted.values
  X_vu <- model_svu$x
  m <- ncol(X_vu) # number of coefficient
  S_vu <- model_svu$y
  R_svu <-S_vu -  p_svu
  
  
  # 1.1.2 [Surrog] (sf) S~X in F -> gamma_f & gamma_f_I -> ph2_u
  model_sf <- glm(formula_S,data = data,family="binomial",x=TRUE,y=TRUE)
  sum_sf <- summary(model_sf)
  gamma_f <- sum_sf$coefficients[,1]
  gamma_f_sd <- sum_sf$coefficients[,2]
  
  p_sf <- model_sf$fitted.values
  P_sf <- diag(p_sf * (1-p_sf), nrow = n, ncol = n)
  
  X_f <- model_sf$x
  S_f <- model_sf$y
  
  gamma_f_I <- (t(X_f) %*% P_sf %*% X_f) / n
  ph2_u <- R_svu * X_vu %*% inv(gamma_f_I)
  
  # 1.1.3 [Ori-Unif] (yvu) Y~X in V-unif -> beta_v_u -> ph1_u
  model_yvu <- multinom(formula_Y,data = data_vu,Hess = TRUE,trace=FALSE)
  sum_yvu <- summary(model_yvu)
  
  beta_vu <- sum_yvu$coefficients
  beta_vu_sd <- sum_yvu$standard.errors
  
  beta_vu_I <- sum_yvu$Hessian/n_vu
  p_yvu <- model_yvu$fitted.values
  Y_vu <- as.matrix(one_hot(as.data.table(data_vu$y)))
  colnames(Y_vu) <- str_replace(colnames(Y_vu),"^V1_","")
  R_yvu <- Y_vu[,-1] - p_yvu[,-1]
  K_yvu <- matrix(NA, nrow = n_vu, ncol = m*k)
  for(i in 1:n_vu){
    K_yvu[i,] <- as.vector(R_yvu[i,]) %x% as.vector(X_vu[i,])
  }
  ph1_u <- K_yvu %*% inv(beta_vu_I)
  
  # 1.1.4 Augment
  Sigma_u_11 <- 1/n_vu * t(ph1_u) %*% ph1_u 
  Sigma_u_12 <- (1-rho)/n_vu * t(ph1_u) %*% ph2_u
  Sigma_u_22 <- (1-rho)/n_vu * t(ph2_u) %*% ph2_u
  
  beta_au <- c(t(beta_vu)) - Sigma_u_12 %*% solve(Sigma_u_22,(gamma_vu - gamma_f))
  beta_au_var <- diag(Sigma_u_11 -Sigma_u_12 %*% inv(Sigma_u_22) %*% t(Sigma_u_12)) /n_vu
  beta_au_sd <- sqrt(beta_au_var)
  
  beta_au_matrix <- matrix(beta_au,ncol =m, nrow= k,byrow = TRUE)
  beta_au_sd_matrix <- matrix(beta_au_sd,ncol =m,nrow = k,byrow = TRUE)
  
  rownames(beta_au_matrix) = rownames(beta_au_sd_matrix) <- rownames(beta_vu)
  colnames(beta_au_matrix) = colnames(beta_au_sd_matrix) <- colnames(beta_vu)
  
  return(
    list(
      "coefficients" = beta_au_matrix,
      "standard.errors" = beta_au_sd_matrix
    )
  )
}



Aug.Bias.estimate <- function(data, 
                              n0=300,
                              n1=300,
                              formula_S = S~x1+x2+x3,
                              formula_Y = y~x1+x2+x3){
  #### 0. Parameters & Validation Data ####
  # data = data_fv
  # n0=300
  # n1=300
  # formula_S = S~x1+x2+x3
  # formula_Y = y~x1+x2+x3
  
  n_s0 <- sum(data$S==1)
  n_s1 <- sum(data$S==2)
  h0 <- n0/n_s0
  h1 <- n1/n_s1
  data$W <- ifelse(data$S==1,h0,h1)
  
  n = nrow(data)
  data_v <- data[!is.na(data$y),]
  n_v <- nrow(data_v)
  # assertthat::are_equal(n_v,n0+n1)
  
  K = nlevels(data$y) 
  k = K-1
  

  
  
  #### 1.2 [Aug-Bias | OSCA] #### 
  # 1.2.1 (sv) S~X in V -> gamma_v -> ph2
  model_sv <- glm(formula_S,data = data_v,family="binomial",x=TRUE,y=TRUE)
  sum_sv <- summary(model_sv)
  gamma_v <- sum_sv$coefficients[,1]
  # gamma_v_sd <-  sum_sv$coefficients[,2]
  
  p_sv <- model_sv$fitted.values
  P_sv <- diag(p_sv*(1-p_sv),nrow = n_v,ncol = n_v)
  X_v <- model_sv$x
  m <- ncol(X_v)
  
  S_v <- model_sv$y
  R_sv <-S_v -  p_sv
  
  gamma_v_I <- (t(X_v) %*% P_sv %*% X_v)/n_v
  ph2 <- R_sv * X_v %*% inv(gamma_v_I)
  
  
  # 1.2.2 (sfw) S~X in F-weight -> gamma_f_w  -> ph3
  W <- data$W
  W_v <- data_v$W
  
  model_sfw <- glm(formula_S,data = data,family="binomial",weights = W,x=TRUE,y=TRUE)
  sum_sfw <- summary(model_sfw)
  gamma_fw <- sum_sfw$coefficients[,1]
  gamma_fw_sd <- sum_sfw$coefficients[,2]
  
  p_sfw <- model_sfw$fitted.values
  P_sfw <- diag(W * p_sfw*(1-p_sfw),nrow = n,ncol = n)
  X_f <- model_sfw$x
  S_f <- model_sfw$y
  R_sfw <- S_f - p_sfw
  R_svw <- R_sfw[!is.na(data$y)] # id_v
  
  gamma_fw_I <- t(X_f) %*% P_sfw %*% X_f / n
  ph3 <- W_v * R_svw * X_v %*% inv(gamma_fw_I)
  
  # 1.2.3 [Ori-Bias] (yv) Y~X in V -> beta_v -> ph1
  model_yv <- multinom(formula_Y,data = data_v,Hess = TRUE,trace=FALSE)
  sum_yv <- summary(model_yv)
  
  beta_v <- sum_yv$coefficients
  beta_v_sd <- sum_yv$standard.errors
  beta_v_I <- sum_yv$Hessian/n_v
  
  p_yv <- model_yv$fitted.values
  Y_v <- as.matrix(one_hot(as.data.table(data_v$y)))
  colnames(Y_v) <- str_replace(colnames(Y_v),"^V1_","")
  R_yv <- Y_v[,-1] - p_yv[,-1]
  K_yv <- matrix(NA, nrow = n_v, ncol = m*k)
  for(i in 1:n_v){
    K_yv[i,] <- as.vector(R_yv[i,]) %x% as.vector(X_v[i,])
  }
  
  ph1 <- K_yv %*% inv(beta_v_I)
  
  # 1.2.4 Augment
  
  E11 <- t(ph1) %*% ph1 /n_v
  E22 <- t(ph2) %*% ph2 /n_v
  E33 <- t(ph3) %*% ph3 /n_v
  E12 <- t(ph1) %*% ph2 /n_v
  E13 <- t(ph1) %*% ph3 /n_v
  E23 <- t(ph2) %*% ph3 /n_v
  
  Sigma_11 <- E11/n_v
  Sigma_22 <- E22/n_v + E33/n - E23/n*2
  Sigma_12 <- E12/n_v - E13/n
  
  beta_a <-  c(t(beta_v)) - Sigma_12 %*% solve(Sigma_22,(gamma_v - gamma_fw))
  beta_a_var <- diag(Sigma_11 - Sigma_12%*% inv(Sigma_22) %*% t(Sigma_12))
  beta_a_sd <- sqrt(beta_a_var)
  
  beta_a_matrix <- matrix(beta_a,ncol =m, nrow= k,byrow = TRUE)
  beta_a_sd_matrix <- matrix(beta_a_sd,ncol =m,nrow = k,byrow = TRUE)
  
  rownames(beta_a_matrix) = rownames(beta_a_sd_matrix) <- rownames(beta_v)
  colnames(beta_a_matrix) = colnames(beta_a_sd_matrix) <- colnames(beta_v)
  
  return(
    list(
      "coefficients" = beta_a_matrix,
      "standard.errors" = beta_a_sd_matrix
    )
  )
  
}

# -------------- Example - Aug.Unif  ---------------

data_fvu <- read.csv("data_uniform_sampled.csv")
# table(data_fvu$y,data_fvu$S)

data_fvu$x3 <- factor(data_fvu$x3)

data_fvu$y <- factor(data_fvu$y) 
data_fvu$y <- relevel(data_fvu$y,ref=1)

data_fvu$S <- factor(data_fvu$S)
data_fvu$S <- relevel(data_fvu$S,ref=1)

Aug.Unif.estimate(data_fvu)

# -------------- Example - Aug.Bias  ---------------
data_fv <- read.csv("data_outcome_dependent_sampled.csv")
# table(data_fv$y,data_fv$S)

data_fv$x3 <- factor(data_fv$x3)

data_fv$y <- factor(data_fv$y) 
data_fv$y <- relevel(data_fv$y,ref=1)

data_fv$S <- factor(data_fv$S)
data_fv$S <- relevel(data_fv$S,ref=1)

Aug.Bias.estimate(data_fv)

