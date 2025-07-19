##################
# Simulations: EDR methods,
# Appendix A.3 (Figure S1, Table S1)


rm(list = ls())
library(dplyr)
library(ggplot2)
library(randomForest)
library(MASS)
library(foreach)
require(doSNOW)
library(caret)
name = "06-02-1test"
if(! dir.exists(paste0("result/",name))){
  dir.create(paste0("result/",name))
}

N_set = c(10000)
n_m_set = c(500)
n_set = c(1000, 2000)
K_set = c(4, 8)
alpha_set = c(0.3)

expand.grid(N_set, n_m_set, n_set,K_set, alpha_set)->settings
colnames(settings) <- c("N", "n_m", "n", "K", "alpha")

settings[order(settings$N, settings$n_m, settings$n, 
               settings$K, settings$alpha), ] -> settings
rownames(settings) <- c(1:nrow(settings))

isalpha <- function(str){
  str <- ifelse(grepl("[a-zA-Z]", str), NA, str) %>% unlist %>% as.numeric()
  return(str)
}
result_td <- function(tau_est, tau_ve, tau_true){
  ci_up = tau_est + sqrt(tau_ve)*1.96
  ci_low = tau_est - sqrt(tau_ve)*1.96
  bias = mean(tau_est, na.rm = T) -tau_true
  rbias = bias/tau_true
  var_empirical = var(tau_est, na.rm = T)
  var_estimated = median(tau_ve, na.rm = T)
  mse = var_empirical+bias^2
  CP <- mean((tau_true <= ci_up) & (tau_true >= ci_low), na.rm = T)*100
  CIlength = mean(ci_up - ci_low, na.rm = T)
  return(c(bias = bias,
           var_empirical = var_empirical,
           mse = mse,
           var_estimated = var_estimated,
           rbias = rbias*100,
           CP = CP,
           CIlength = CIlength))
}

for(j in 1:nrow(settings)){
  N = settings$N[j]
  n_m = settings$n_m[j]
  n = settings$n[j]
  K = settings$K[j]
  alpha = settings$alpha[j]
  rrr <- c()
  path = paste0("result/", name, "/", paste(N, n_m, n, K, alpha, sep = "_"), ".rdata")
  # print(path)
  load(file = path)
  apply(result, 2, isalpha) -> result2
  rbind(result_td(result2[, 1], result2[, 5], tau_true = -0.61),
        result_td(result2[, 2], result2[, 6], tau_true = -0.61),
        result_td(result2[, 3], result2[, 7], tau_true = -0.61),
        result_td(result2[, 4], result2[, 8], tau_true = -0.61)) -> re
  rownames(re) <- c( "u","z", "f", "o")
  save(result, file = paste0("result_table/",name, "/", "n_m=",
                             n_m, "_n=",n,  "_K=", K,"_alpha=", alpha,".rdata"))
  re %>% write.csv(file = paste0("result_table/",name, "/", "n_m=",
                                 n_m, "_n=",n,  "_K=", K,"_alpha=", alpha,".csv"))
}

n_m_set = c(500)
n_set = c(1000, 2000)
K_set = c(4, 8)
alpha_set = c(0.3)
for(n_m in n_m_set){
  for(n in n_set){
    for(K in K_set){
      for(alpha in alpha_set){
        read.csv(file = paste0("result_table/", name, "/", "n_m=",
                               n_m, "_n=",n,  "_K=", K,"_alpha=", alpha,".csv")) %>% as.data.frame()  -> re
        re[c(1,4,3), ]%>%
          dplyr::select(c("rbias", "var_empirical", "CP", "CIlength") )%>%
          round(3) %>%
          write.csv(file = paste0("result_table/", name, "/cleanre/clean_n_m=",
                                  n_m, "/a=", alpha, "/", "n=",n,  "_K=", K, ".csv"))
      }
    }
  }
}





