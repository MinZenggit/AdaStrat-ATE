#######
# Evaluate the efficiency gain estimation
# Appendix A.4 (Table S2, Figure S2)
rm(list = ls())
library(dplyr)
library(ggplot2)
library(randomForest)
library(MASS)
library(foreach)
require(doSNOW)
library(caret)

name = "05-28-2test"
if(! dir.exists(paste0("result/",name))){
  dir.create(paste0("result/",name))
}

if(! dir.exists(paste0("result/",name))){
  dir.create(paste0("result/",name))
}

method_set <- c("aipw", "ipw", "reg")
N_set = c(10000)
n_m_set = c(500, 1000, 2000)
n_set = c(2000)
K_set = c(2, 4, 6, 8, 10)

expand.grid(method_set, N_set, n_m_set, n_set,K_set)->settings
colnames(settings) <- c("method", "N", "n_m", "n", "K")

settings[order(settings$method, settings$N, 
               settings$n_m, settings$n, 
               settings$K), ] -> settings
rownames(settings) <- c(1:nrow(settings))
cores = 40
trials = 800
for(j in 1:nrow(settings)){
  method = settings$method[j]
  N = settings$N[j]
  n_m = settings$n_m[j]
  n = settings$n[j]
  K = settings$K[j]
  cat(" Method=", as.character(method), ", N=", N,", n_m=", n_m,", n=", n, ", K=", K, "\n")
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
            .packages = c("dplyr",
                          "nnet","caret",
                          "MASS"))%dopar%{
                            source("somefunciton_rest.R")
                            return(testf(N, n_m, n, K, method))}
  stopCluster(cls)
  endtime <- Sys.time()
  cat( "\nRunning time: ", endtime-starttime, "\n")
  path = paste0("result/", name, "/", paste(N, n_m, n, K, method, sep = "_"), ".rdata")
  save(result, file = path)
}


result_td_r <- function(x){
  mean(result[, 6]) -> r_hat
  mean(result[, 5]) -> r
  var(result[, 1])/var(result[, 3]) -> r_em
  return(c(r, r_hat, r_em))
}

settings$r = settings$r_hat = settings$r_em = 0
result_td_r(result)
for(j in 1:nrow(settings)){
  method = settings$method[j]
  N = settings$N[j]
  n_m = settings$n_m[j]
  n = settings$n[j]
  K = settings$K[j]
  path = paste0("result/", name, "/", paste(N, n_m, n, K, method, sep = "_"), ".rdata")
  load(file = path)
  result_td_r(result) -> rr
  settings$r[j] = rr[1]
  settings$r_hat[j] = rr[2]
  settings$r_em[j] = rr[3]
}

