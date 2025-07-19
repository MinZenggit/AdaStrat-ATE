#############
# Simulation: main results in Table 1 and Figure 3

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



rm(list = ls())
library(dplyr)
library(ggplot2)
library(randomForest)
library(MASS)
library(foreach)
require(doSNOW)
library(caret)
name = "05-28-1test"
if(! dir.exists(paste0("result/",name))){
  dir.create(paste0("result/",name))
}

if(! dir.exists(paste0("result/",name))){
  dir.create(paste0("result/",name))
}

method_set <- c("aipw", "ipw", "reg")
N_set = c(10000)
n_m_set = c(500, 1000)
n_set = c(1000, 1200, 1500, 2000, 4000)
K_set = c(2, 4, 6, 8, 10)
alpha_set = c(0.05, 0.1, 0.2, 0.3)

expand.grid(method_set, N_set, n_m_set, n_set,K_set, alpha_set)->settings
colnames(settings) <- c("method", "N", "n_m", "n", "K", "alpha")

settings[order(settings$method, settings$N, 
               settings$n_m, settings$n, 
               settings$K, settings$alpha), ] -> settings
rownames(settings) <- c(1:nrow(settings))
cores = 4
trials =800
for(j in 1:nrow(settings)){
  method = settings$method[j]
  N = settings$N[j]
  n_m = settings$n_m[j]
  n = settings$n[j]
  K = settings$K[j]
  alpha = settings$alpha[j]
  cat(" Method=", method, ", N=", N,", n_m=", n_m,", n=", n, ", K=", K,  ", alpha=", alpha,"\n")
  starttime <- Sys.time()
  cls <- makeSOCKcluster(cores)
  registerDoSNOW(cls)
  pb <- txtProgressBar(max=trials, style=3, char = "*",)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
  result <- 
    foreach(i=c(1: trials),
            .combine = rbind, 
            .options.snow=opts,
            .errorhandling = "pass",
            .packages = c("dplyr","randomForest",
                          "nnet","caret",
                          "MASS"))%dopar%{
                            source("codes_simulation/somefunction17.R")
                            return(testf(N, 
                                         n_m,
                                         n,
                                         K,
                                         method,
                                         alpha))}
  stopCluster(cls)
  endtime <- Sys.time()
  cat( "\nRunning time: ", endtime-starttime, "\n")
  path = paste0("result/", name, "/", paste(N, n_m, n, K, alpha, method, sep = "_"), ".rdata")
  save(result, file = path)
}
# 
# library(stringr)
# isalpha <- function(str){
#   str <- ifelse(grepl("[a-zA-Z]", str), NA, str) %>% unlist %>% as.numeric()
#   return(str)
# }
# result_td <- function(re){
#   re = cbind((apply(re, 2, mean, na.rm = T)-0.50)/0.50 *100,
#              apply(re, 2, var, na.rm = T),
#              (apply(re, 2, mean, na.rm = T)-0.50)^2 +apply(re, 2, var, na.rm = T) )
#   colnames(re) <- c("bias", "var", "mse")
#   return(re)
# }
# 
# isalpha <- function(str){
#   str <- ifelse(grepl("[a-zA-Z]", str), NA, str) %>% unlist %>% as.numeric()
#   return(str)
# }
# r_td <- c()
# for(j in 1:nrow(settings)){
#   method = settings$method[j]
#   N = settings$N[j]
#   n_m = settings$n_m[j]
#   n = settings$n[j]
#   K = settings$K[j]
#   alpha = settings$alpha[j]
#   rrr <- c()
#   path = paste0("result/", name, "/", paste(N, n_m, n, K, alpha, method, sep = "_"), ".rdata")
#   # print(path)
#   load(file = path)
#   apply(result, 1, isalpha)%>%t() -> result
#   result_td(result[,  c(1,3,5,7,9,11, 13, 15, 17)]) %>% as.data.frame() -> rr
#   cbind(rr,
#         var_pot=colMeans(result[,c(1,3,5,7,9,11, 13, 15, 17)+1], na.rm = T),
#         method = rep(method, 9),
#         N = rep(N, 9),
#         n_m = rep(n_m,9),
#         n = rep(n, 9),
#         K = rep(K, 9),
#         alpha = rep(alpha, 9),
#         des = c("u", "b", "m1", "m2","z", "zm","f", "fm", "o"))->rrr
#   colnames(rrr) <- c("Bias", "Var", "MSE",
#                      "Var_est",
#                      "method", "N", "n_m", "n", "K", "alpha","des")
#   rbind(r_td, rrr) -> r_td
# }


