#############
# Test if the Propensity score model is wrong
# Appendix Table S6.
rm(list = ls())
library(dplyr)
library(ggplot2)
library(randomForest)
library(MASS)
library(foreach)
require(doSNOW)
library(caret)
source("somefunction_0602.R")
name = "06-02-2test"
if(! dir.exists(paste0("result/",name))){
  dir.create(paste0("result/",name))
}

if(! dir.exists(paste0("result/",name))){
  dir.create(paste0("result/",name))
}

method_set <- c("ipw", "aipw", "reg")
N_set = c(10000)
n_m_set = c(500)
n_set = c(1000, 2000)
K_set = c(4, 8)
alpha_set = c(0.3)

expand.grid(method_set, N_set, n_m_set, n_set,K_set, alpha_set)->settings
colnames(settings) <- c("method", "N", "n_m", "n", "K", "alpha")

settings[order(settings$method, settings$N, 
               settings$n_m, settings$n, 
               settings$K, settings$alpha), ] -> settings
rownames(settings) <- c(1:nrow(settings))
cores = 40
trials =800
for(j in 1:nrow(settings)){
  method = settings$method[j]
  N = settings$N[j]
  n_m = settings$n_m[j]
  n = settings$n[j]
  K = settings$K[j]
  alpha = settings$alpha[j]
  cat(" Method=", as.character(method), ", N=", N,", n_m=", n_m,", n=", n, ", K=", K,  ", alpha=", alpha,"\n")
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
                            source("somefunction_0602.R")
                            return(testf(N, 
                                         n_m,
                                         n,
                                         K,
                                         method,
                                         alpha))}
  apply(result, 1, isalpha)%>%t() -> result
  stopCluster(cls)
  endtime <- Sys.time()
  print("\n")
  print(apply(result[, c(1, 3, 5, 7, 9, 11, 13, 15, 17)], 2, var, na.rm = T))
  print(colMeans(result[, c(1, 3, 5, 7, 9, 11, 13, 15, 17) + 1], na.rm = T))
  cat( "\nRunning time: ", endtime-starttime, "\n")
  path = paste0("result/", name, "/", paste(N, n_m, n, K, alpha, method, sep = "_"), ".rdata")
  save(result, file = path)
}
# 
library(stringr)
isalpha <- function(str){
  str <- ifelse(grepl("[a-zA-Z]", str), NA, str) %>% unlist %>% as.numeric()
  return(str)
}
result_td <- function(re){
  re = cbind((apply(re, 2, mean, na.rm = T)-0.50)/0.50 *100,
             apply(re, 2, var, na.rm = T),
             (apply(re, 2, mean, na.rm = T)-0.50)^2 +apply(re, 2, var, na.rm = T) )
  colnames(re) <- c("bias", "var", "mse")
  return(re)
}
# 
# 
# 
# data.frame(result_td(result[, c(1,3,5,7,9,11, 13, 15, 17)]),
#            var_estimated = colMeans(result[, c(1,3,5,7,9,11, 13, 15, 17)+1])) -> rr
# rownames(rr) <- c("u", "b", "m1", "m2","z", "zm","f", "fm", "o")
# rr
#
if(! dir.exists(paste0("result_figure/",name))){
  dir.create(paste0("result_figure/",name))
}
r_td <- c()
for(j in 1:nrow(settings)){
  method = settings$method[j]
  N = settings$N[j]
  n_m = settings$n_m[j]
  n = settings$n[j]
  K = settings$K[j]
  alpha = settings$alpha[j]
  rrr <- c()
  path = paste0("result/", name, "/", paste(N, n_m, n, K, alpha, method, sep = "_"), ".rdata")
  # print(path)
  load(file = path)
  apply(result, 1, isalpha)%>%t() -> result
  result_td(result[,  c(1,3,5,7,9,11, 13, 15, 17)]) %>% as.data.frame() -> rr
  cbind(rr,
        var_pot=colMeans(result[,c(1,3,5,7,9,11, 13, 15, 17)+1], na.rm = T),
        method = rep(method, 9),
        N = rep(N, 9),
        n_m = rep(n_m,9),
        n = rep(n, 9),
        K = rep(K, 9),
        alpha = rep(alpha, 9),
        des = c("u", "b", "m1", "m2","z", "zm","f", "fm", "o"))->rrr
  colnames(rrr) <- c("Bias", "Var", "MSE",
                     "Var_est",
                     "method", "N", "n_m", "n", "K", "alpha","des")
  rbind(r_td, rrr) -> r_td
}
 
# ####################################
# # rm(p)
# alpha = 0.1
# for(method_i in "aipw"){
#   for(n_m_i in n_m_set){
#     for(K_i in K_set){
#       pdf(file = paste0("result_figure/", name, "/","n_m=",n_m_i, "_", "K=", K_i,"_", method_i, ".pdf"), 
#           onefile = F)
#       r_td %>% filter(method == method_i,n_m == n_m_i, K==K_i, alpha == 0.3,
#                       des %in%c("u", "b", "o", "fm")) %>%
#         ggplot(aes(x = n, y = MSE, color = des, shape = des)) + 
#         geom_line()+
#         geom_point(size = 5) +  # Adjust the size of the points as needed
#         theme_minimal() -> p
#       plot(p)
#       dev.off()
#       rm(p)
#     }
#   }
# }
# 
# 
# 
