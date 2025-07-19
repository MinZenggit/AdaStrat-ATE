################################################################################
#                                                                              #
#                    SIMULATION STUDY FOR SAMPLING DESIGNS                       #
#                                                                              #
################################################################################

# This script conducts a simulation study to compare the performance of various
# sampling designs for estimating the Average Treatment Effect (ATE).

##################
# We consider 8 different methods: u, b, m1, m2, z, zm, f, fm, o
#####SimRan######
# u: Simple Random sampling design (SimRan)
#    - A baseline method where all units have an equal probability of selection.
#####ORACLE######
# b: the Oracle sampling design (ORACLE)
#    - An ideal, infeasible method that uses the true, unknown influence functions
#      to determine optimal sampling probabilities. Serves as a theoretical benchmark.

######AdaStrat################
# m1, m2, z, zm, f, fm are adaptive strata methods (AdaStrat)
# These methods use a two-phase design. A pilot sample is drawn to estimate
# properties of the population, which then inform the sampling strategy for the
# main sample.

# m1, m2: Stratified by error-prone influence functions.
#   - Strata are created based on influence functions estimated from incomplete data.
# z: Stratified using nnet
#   - A neural network is trained on the pilot sample to predict optimal strata for
#     the remaining units.
# zm: Stratified using nnet (error-prone influence functions also as a predictor)
#   - Same as 'z', but the error-prone influence function is included as a
#     predictor in the neural network.
# f: Stratified using random forest
#   - A random forest is trained on the pilot sample to predict optimal strata.
# fm: Stratified using random forest (error-prone influence functions also as a predictor)
#   - Same as 'f', but the error-prone influence function is included as a
#     predictor in the random forest.

# The error-prone influence functions were calculated using the
# incomplete data (x_c, treatment, y) to estimate the ATE.

######FixStrat#################
# o: Stratified by outcome (FixStrat)
#   - A fixed stratification scheme based on the observed outcome variable.

# Source helper functions for ATE estimation.
# esATE.R contains the main function to estimate ATE and influence functions.
# esATE_w2.R contains a weighted version for survey samples.
source("codes_ate_influence_functions/esATE.R")
source("codes_ate_influence_functions/esATE_w2.R")


#' Calculate Oracle Probabilities
#'
#' This function computes optimal sampling probabilities based on a vector of
#' influence function values. It implements an algorithm to cap probabilities
#' to prevent any single unit from having an excessively high chance of selection,
#' ensuring that no probability exceeds 1/sn.
#'
#' @param phi_a A numeric vector of sorted, absolute influence function values.
#' @param sn The desired sample size.
#' @param alpha A small smoothing parameter (default 0.05) to ensure all
#'   probabilities are non-zero.
#' @return A numeric vector of sampling probabilities.
getpi <- function(phi_a, sn, alpha = 0.05){
  ## phi_a is sorted abs(phi)
  pi_a = phi_a/sum(phi_a)
  if(max(pi_a)<=(1/sn)){
    # If all probabilities are below the cap, apply smoothing and return.
    pi_a = (1-alpha)*pi_a + alpha * mean(pi_a)
    return(pi_a)
  }else{
    # If some probabilities exceed the cap, find the optimal threshold H.
    N = length(phi_a)
    H <- re <- c()
    k=0
    for(g in 0:(sn-1)){
      k=g+1
      H[k] = sum(phi_a[1:(N-k)])/(sn-k)
    }
    for(k in 1:sn){
      (H[k] >= phi_a[N-k])&(H[k] <= phi_a[N-k+1]) -> re[k]
    }
    g = which(re==T)
    H[g] -> H
    # Cap the probabilities at the threshold H.
    pi = pmin(H, phi_a)
    pi = pi/sum(pi)
    # Ensure no probability exceeds the 1/sn limit and apply smoothing.
    pi <- ifelse(pi>(1/sn), 1/sn, pi)
    pi = (1-alpha)*pi + alpha * mean(pi)
    return(pi)
  }
}


#' Calculate Optimal Stratified Sampling Probabilities
#'
#' Implements Neyman allocation for stratified sampling. It calculates the
#' optimal sampling fraction for each stratum based on its estimated standard
#' deviation and size.
#'
#' @param t A numeric vector of stratum-level standard deviations (or a value
#'   proportional to it  by influence fluences).
#' @param N The total population size.
#' @param n The desired sample size.
#' @param K The number of strata.
#' @param p A numeric vector of stratum proportions in the population.
#' @param alpha A small smoothing parameter (default 0.05).
#' @return A numeric vector of stratum-specific sampling probabilities.
getpi_stra <- function(t, N, n, K, p, alpha = 0.05){
  ## t is proportional to stratum standard deviation
  pi_a = t/sum(t*p)*n/N
  if(max(pi_a)<=(1)){
    # If all sampling fractions are valid (<=1), apply smoothing.
    pi_a = (1-alpha)*pi_a + alpha * mean(pi_a)
    return(pi_a)
  }else{
    # If any sampling fraction exceeds 1, find a threshold to cap them.
    H <- re <- c()
    for(g in 1:(n*K/N-1)){
      H[g] = sum(t[1:(K-g)]*p[1:(K-g)])/(n/N-sum(p[(K-g+1): K]))
    }
    for(g in 1:(n*K/N-1)){
      (H[g] <= t[K-g+1]& H[g] >= t[K-g]) -> re[g]
    }
    g = which(re==T)[1]
    H[g] -> H
    # Cap the allocation and apply smoothing.
    pi = t/H
    pi <- ifelse(pi>1, 1, pi)
    pi = (1-alpha)*pi + alpha * mean(pi)
    return(pi)
  }
}


#' Generate Simulation Data
#'
#' Creates a dataset for the simulation with specified population size.
#' The data generating process includes a covariate (X1), an unobserved
#' confounder (U1), potential outcomes (Y0, Y1), and a treatment assignment (Z).
#'
#' @param N The population size.
#' @return A data frame containing the generated population data.
Generatedt <- function(N){
  id <- 1:N
  X1 = runif(N, 0, 2)
  U1 = 0.5 + 0.5*X1 - 2*sin(X1) + 2*sign(sin(5*X1)) + rnorm(N, sd = 2)
  Y0 = - X1 - U1 + rnorm(N, sd = 2)
  Y1 = - X1 + 4*U1 + rnorm(N, sd = 2)
  e = plogis(1 - 0.5*X1-0.5*U1)
  Z <- rbinom(N, 1, e)
  Y = Z*Y1 + (1-Z)*Y0
  # ### additive form measurement (not used in main analysis)
  Ym_a = 2*Y + X1 + rnorm(N)
  # ### multiple form measurement (not used in main analysis)
  Ym_m = Y * rnorm(N)
  ### nonlinear form measurement (not used in main analysis)
  Ym_s = plogis(Y/5 + rnorm(N, 0, 0.1))
  dt1 <- cbind(id, X1, U1, Y0, Y1, e, Z, Y, Ym_s, Ym_a, Ym_m) %>% as.data.frame()
  return(dt1)
}


#' Estimate ATE from a Single-Phase Sample
#'
#' This function takes a dataset, applies a sampling design based on provided
#' probabilities, and returns an estimate of the ATE using a weighted estimator.
#'
#' @param dt The full data frame from which to sample.
#' @param pi A numeric vector of sampling probabilities for each unit in `dt`.
#' @param method The ATE estimation method to be passed to `esATE_w`.
#' @return A numeric vector containing the ATE estimate and its variance.
get_tau_hat <- function(dt, pi, method){
  dt$pi = pi
  N = nrow(dt)
  # Simulate the sampling process.
  dt$V = rbinom(N, 1, dt$pi)
  dt_Z= filter(dt, dt$V==1)
  # Estimate ATE using inverse probability weighting.
  esATE_w(cbind(dt_Z$X1, dt_Z$U1), dt_Z$Y, dt_Z$Z, w = 1/dt_Z$pi, method, N) -> tau_hat
  return(tau_hat)
}


#' Estimate ATE from a Two-Phase Sample
#'
#' This function estimates the ATE by combining results from a pilot sample
#' and a main sample. It estimates the ATE in each sample separately and then
#' combines them using an inverse-variance weighted average.
#'
#' @param dt1_2 The data frame for the main sampling phase.
#' @param dt2 The data frame for the pilot sample.
#' @param pi1_2 The sampling probabilities for the main sample (`dt1_2`).
#' @param method The ATE estimation method.
#' @return A numeric vector containing the combined ATE estimate and its variance.
get_tau_hat2 <- function(dt1_2, dt2, pi1_2, method){
  dt1_2$pi = pi1_2;
  N1_2 = nrow(dt1_2);n_m = nrow(dt2); N = n_m + N1_2
  # Simulate sampling for the main phase.
  dt1_2$V = rbinom(N1_2, 1, dt1_2$pi)
  dt1_2= filter(dt1_2, dt1_2$V==1)
  # Estimate ATE from the main sample (weighted).
  esATE_w(cbind(dt1_2$X1, dt1_2$U1), dt1_2$Y, dt1_2$Z, w = 1/dt1_2$pi, method, N) ->fit1_2
  # Estimate ATE from the pilot sample (unweighted).
  esATE(cbind(dt2$X1, dt2$U1), dt2$Y, dt2$Z, method) -> fit2
  # Combine the two estimates using inverse-variance weighting.
  (fit1_2[1]/fit1_2[2] + fit2$est/fit2$ve)/(1/fit1_2[2] + 1/fit2$ve) -> tau_hat
  (fit1_2[2]*fit2$ve)/(fit1_2[2]+fit2$ve) -> ve_hat
  cb_2 <- c(tau_hat, ve_hat)
  return(c(cb_2))
}


#' Main Simulation Function
#'
#' Runs a single replication of the simulation experiment, comparing all
#' specified sampling designs.
#'
#' @param N Total population size.
#' @param n_m Pilot sample size.
#' @param n Main sample size.
#' @param K Number of strata to create.
#' @param method ATE estimation method.
#' @param alpha Smoothing parameter for probability calculation.
#' @return A named numeric vector with ATE estimates and variances for all methods.
testf <- function(N, n_m, n, K, method, alpha){
  # Generate the full population data.
  Generatedt(N) -> dt1 -> dt
  
  #######
  # Create fixed strata based on outcome quantiles.
  te <- dt1 %>% dplyr::select(id, Y)
  te[order(te$Y), ] -> te
  R_O = rep(1:K, each = ceiling(N/K)) %>% as.factor()
  te$R_O =  R_O[1:N]
  dt1$R_O = te[order(te$id), "R_O"]
  
  ######
  # Calculate true and error-prone influence functions for the whole population.
  x = cbind(dt1$X1, dt1$U1); y = dt1$Y; A = dt1$Z
  
  # The true ATE and true influence functions (using observed X1 and unobserved U1).
  esATE(x, y, A, method) -> fit1
  dt1$phi = fit1$infl
  # Estimate the error-prone influence functions (using only observed X1).
  esATE(dt1$X1, y, A, method) -> fit1_c
  dt1$phi_m = fit1_c$infl
  
  # Calculate absolute IF.
  phi_bar = mean(dt1$phi)
  phi_bar_m = mean(dt1$phi_m)
  dt1$phi_a = abs(dt1$phi-phi_bar)
  dt1$phi_m_a = abs(dt1$phi_m-phi_bar_m)
  
  # Create strata based on error-prone influence functions.
  dt1 <- dt1[order(dt1$phi_m_a), ];
  R_M = rep(1:K, each = ceiling(N/K)) %>% as.factor()
  dt1$R_M = R_M[1:N]
  dt1 <- dt1[order(dt1$phi_a), ];
  
  #####
  # Define sampling probabilities for non-adaptive designs.
  # u: Uniform (Simple Random Sampling: SimRan)
  dt1$pi_u <- rep((n+n_m)/N, N)
  
  # b: Oracle (using true influence functions)
  getpi(dt1$phi_a, (n+n_m), alpha)*(n+n_m) -> dt1$pi_b
  
  ############
  # TWO-PHASE SAMPLING
  # Draw the pilot sample (middle phase).
  midphase_ind <- sample(1:N, size = n_m, replace = F)
  #### The pilot data
  dt2 <- dt1[midphase_ind,]
  #### The remaining data for the main sample
  dt1_2 <- dt1[-midphase_ind,]
  N1_2 = nrow(dt1_2)
  
  # --- Use pilot data (dt2) to inform sampling of main data (dt1_2) ---
  x2 = cbind(dt2$X1, dt2$U1); y2 = dt2$Y; A2 = dt2$Z
  # Estimate influence functions using the pilot data.
  esATE(x2, y2, A2, method) -> fit2
  dt2$phi_hat <- fit2$infl
  phi_hat_bar <- mean(dt2$phi_hat)
  dt2 <- dt2[order(dt2$phi_hat), ]
  # Create "optimal" strata based on pilot-estimated influence functions.
  R_L = rep(1:K, each = ceiling(n_m/K)) %>% as.factor()
  dt2$R_L = R_L[1:n_m]
  
  # Calculate intra-stratum variance of estimated influence functions.
  te <- c()
  for(k in 1:K){
    te[k] <- mean((dt2[dt2$R_L==k,"phi_hat"]-phi_hat_bar)^2)
  }
  te2 <- te
  # Calculate optimal stratum sampling probabilities based on these variances.
  getpi_stra(sqrt(te2), N1_2, n, K, p = rep(1/K, K), alpha) -> pi_m
  dt2$pi_m = NA
  for(k in 1:K){
    dt2[dt2$R_L==k, "pi_m"] = pi_m[k]
  }
  
  
  ##############
  # o: Outcome-dependent group (FixStrat)
  # Calculate allocation based on outcome strata, using pilot data variances.
  te <- c()
  for(k in 1:K){
    te[k] <- mean((dt2[dt2$R_O==k,"phi_hat"]-phi_hat_bar)^2)
  }
  te2 <- te
  getpi_stra(sqrt(te2), N1_2, n, K, p = rep(1/K, K), alpha) -> pi_O
  # Assign these probabilities to the main sample frame.
  pr_O <- data.frame(R_O = 1:K, p_O = pi_O)
  merge(dt1_2, pr_O, by="R_O") -> dt1_2_O
  merge(dt2, pr_O, by="R_O") -> dt2_O
  
  #######
  # m1: Stratified by error-prone IFs (m1)
  # Calculate allocation based on error-prone IF strata, using pilot data variances.
  te <- c()
  for(k in 1:K){
    te[k] <- mean((dt2[dt2$R_M==k,"phi_hat"]-phi_hat_bar)^2)
  }
  te2 <- te
  getpi_stra(sqrt(te2), N1_2, n, K, p = rep(1/K, K), alpha) -> p_M
  pr_M <- data.frame(R_M = 1:K, p_M = p_M)
  merge(dt1_2, pr_M, by="R_M") -> dt1_2_M1
  
  #######
  # m2: Stratified by error-prone IFs (m2)
  # Uses error-prone IF strata, but with probabilities `pi_m` derived from
  # the "optimal" pilot data strata (R_L).
  pr_M <- data.frame(R_M = 1:K, p_M = pi_m)
  merge(dt1_2, pr_M, by="R_M") -> dt1_2_M2
  
  ###########
  # z: Stratified using Neural Network
  # Train a neural net to predict the "optimal" stratum (R_L) from covariates.
  digit.x = dt2[, c("X1", "Y", "Z")]
  digit.y = dt2[, "R_L"]
  digit.ml <- train(x=digit.x, y=digit.y,
                    method="nnet",
                    tuneGrid=expand.grid(
                      .size=10,      # hidden neurons
                      .decay=0.1     # weight decay
                    ),
                    trControl=trainControl(method="none"),
                    MaxNWts=10000,   # maximum number of weights
                    maxit=500,       # maximum number of iterations
                    verbose = FALSE) # Suppress verbose output
  # Predict strata for the main sample frame and assign probabilities.
  dt1_2$R_Z <- predict(digit.ml, newdata = dt1_2)
  pr_Z <- data.frame(R_Z = 1:K, p_Z = pi_m)
  merge(dt1_2, pr_Z, by="R_Z") -> dt1_2_Z
  
  ##########
  # zm: Stratified using Neural Network with error-prone IF
  # Same as 'z', but include the error-prone IF as a predictor.
  digit.x = dt2[, c("X1", "Y", "Z", "phi_m")]
  digit.y = dt2[, "R_L"]
  digit.zm <- train(x=digit.x, y=digit.y,
                    method="nnet",
                    tuneGrid=expand.grid(
                      .size=10,      # hidden neurons
                      .decay=0.1     # weight decay
                    ),
                    trControl=trainControl(method="none"),
                    MaxNWts=10000,   # maximum number of weights
                    maxit=500,       # maximum number of iterations
                    verbose = FALSE) # Suppress verbose output
  dt1_2$R_ZM <- predict(digit.zm, newdata = dt1_2)
  pr_ZM <- data.frame(R_ZM = 1:K, p_ZM = pi_m)
  merge(dt1_2, pr_ZM, by="R_ZM") -> dt1_2_ZM
  
  
  ######################
  # f: Stratified using Random Forest
  # Train a random forest to predict the "optimal" stratum (R_L).
  fit_f <- randomForest(R_L~X1+Y+Z , data = dt2, proximity = T)
  dt1_2$R_F <- predict(fit_f, newdata = dt1_2)
  pr_F <- data.frame(R_F = 1:K, p_F = pi_m)
  merge(dt1_2, pr_F, by="R_F") -> dt1_2_F
  
  ######################
  # fm: Stratified using Random Forest with error-prone IF
  # Same as 'f', but include the error-prone IF as a predictor.
  fit_f <- randomForest(R_L~X1+Y+Z+phi_m, data = dt2, proximity = T)
  dt1_2$R_F <- predict(fit_f, newdata = dt1_2)
  pr_F <- data.frame(R_F = 1:K, p_FM = pi_m)
  merge(dt1_2, pr_F, by="R_F") -> dt1_2_FM
  
  ########################
  # Theoretical variances for comparison purposes.
  sum((dt1_2$phi_a)^2/(n/(N-n_m)))/N/N
  sum((dt1_2_M1$phi_a)^2/dt1_2_M1$p_M)/N/N
  sum((dt1_2_M2$phi_a)^2/dt1_2_M2$p_M)/N/N
  sum((dt1_2_Z$phi_a)^2/dt1_2_Z$p_Z)/N/N
  sum((dt1_2_ZM$phi_a)^2/dt1_2_ZM$p_ZM)/N/N
  sum((dt1_2_F$phi_a)^2/dt1_2_F$p_F)/N/N
  sum((dt1_2_FM$phi_a)^2/dt1_2_FM$p_FM)/N/N
  sum((dt1_2_O$phi_a)^2/dt1_2_O$p_O)/N/N
  
  ########
  # --- ESTIMATE ATE FOR EACH DESIGN ---
  # Single-phase designs
  get_tau_hat(dt1, dt1$pi_u, method) -> fit_u
  get_tau_hat(dt1, dt1$pi_b, method) -> fit_b
  
  # Two-phase designs
  get_tau_hat2(dt1_2_M1, dt2, dt1_2_M1$p_M, method) -> fit_M1
  get_tau_hat2(dt1_2_M2, dt2, dt1_2_M2$p_M, method) -> fit_M2
  get_tau_hat2(dt1_2_Z, dt2, dt1_2_Z$p_Z, method) -> fit_Z
  get_tau_hat2(dt1_2_ZM, dt2, dt1_2_ZM$p_ZM, method) -> fit_ZM
  get_tau_hat2(dt1_2_F, dt2, dt1_2_F$p_F, method) -> fit_F
  get_tau_hat2(dt1_2_FM, dt2, dt1_2_FM$p_FM, method) -> fit_FM
  get_tau_hat2(dt1_2_O, dt2, dt1_2_O$p_O, method) -> fit_O
  
  # Return all results as a single vector.
  return(c(fit_u, fit_b, fit_M1, fit_M2, fit_Z, fit_ZM, fit_F, fit_FM, fit_O))
}


#' Summarize Simulation Results
#'
#' Calculates bias, variance, and mean squared error (MSE) from the output
#' of multiple simulation replications.
#'
#' @param re A matrix or data frame where each row is a replication and each
#'   column corresponds to an estimator's result.
#' @return A matrix with "bias", "var", and "mse" for each method.
result_td <- function(re){
  # The true ATE value is hardcoded as 0.50 for the bias calculation.
  re = cbind(apply(re, 2, mean, na.rm = T)-0.50,
             apply(re, 2, var, na.rm = T),
             (apply(re, 2, mean, na.rm = T)-0.50)^2 +apply(re, 2, var, na.rm = T) )
  colnames(re) <- c("bias", "var", "mse")
  return(re)
}

# Example call to the main simulation function (commented out).
testf(N = 10000, 
      n_m = 500,
      n = 1000,
      K = 2,
      method = "ipw",
      alpha = 0.05)
