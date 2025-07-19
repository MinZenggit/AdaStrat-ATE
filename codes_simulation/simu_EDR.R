##################
# Simulations: EDR methods,
# Appendix A.3 (Figure S1, Table S1)
rm(list = ls())
library(dplyr)
library(ggplot2)
library(randomForest)
library(MASS)
library(nnet)
library(foreach)
require(doSNOW)
library(caret)

source("somefunction17.R")
source("EDR.R")
source("somefunctions_0601.R")

name = "06-02-1test"
if(! dir.exists(paste0("result/",name))){
  dir.create(paste0("result/",name))
}

N_set = c(20000)
n_m_set = c(500)
n_set = c(1000, 2000)
K_set = c(4, 8)
alpha_set = c(0.3)

expand.grid(N_set, n_m_set, n_set,K_set, alpha_set)->settings
colnames(settings) <- c("N", "n_m", "n", "K", "alpha")

settings[order(settings$N, settings$n_m, settings$n, 
               settings$K, settings$alpha), ] -> settings
rownames(settings) <- c(1:nrow(settings))


cores = 40; trials = 800

for(j in 1:nrow(settings)){
  method = settings$method[j]
  N = settings$N[j]
  n_m = settings$n_m[j]
  n = settings$n[j]
  K = settings$K[j]
  alpha = settings$alpha[j]
  cat("\n N=", N,", n_m=", n_m,", n=", n, ", K=", K,  ", alpha=", alpha,"\n")
  
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
                            source("somefunction17.R")
                            source("EDR.R")
                            source("somefunctions_0601.R")
                            return(testf(N, n_m, n, alpha, K))}
  apply(result, 1, isalpha)%>%t() -> result
  stopCluster(cls)
  endtime <- Sys.time()
  print("\n")
  print(apply(result[, 1:4], 2, var, na.rm = T))
  print(colMeans(result[, 5:8], na.rm = T))
  cat( "\nRunning time: ", endtime-starttime, "\n")
  path = paste0("result/", name, "/", paste(N, n_m, n, K, alpha, sep = "_"), ".rdata")
  save(result, file = path)
}

