#####################
# Combine the two sample estimators
rm(list = ls())
source("codes_simulation/somefunction16.R")
source("codes_ate_influence_functions/esATE.R")
source("codes_ate_influence_functions/esATE_w2.R")
library(dplyr)
library(randomForest)
library(nnet)
library(caret)
library(MASS)
library(foreach)
require(doSNOW)
set.seed(1)
name = "0905_test"
######
method_set <- c("aipw", "reg")
n_m_set = c(1000)
n_set = c(3000, 2000)
K_set = c(4, 8, 12)
alpha_set = c(0.3)
cores = 40
trials = 400

# n_m = 1000; n = 2000; method = "aipw"; alpha = 0.3; K = 8
for(method in method_set){
  for(n_m in n_m_set){
    for(n in n_set){
      for(K in K_set){
        for(alpha in alpha_set){
          cat(" Method=", method, ", n_m=", n_m,", n=", n, ", K=", K, ", alpha=",alpha, "\n")
          load(file = "dt.rdata")
          dt1 <- dt
          names(dt1) <- c("f.eid", "Ethnicity", "Sex", "Age", "Education", "PA_time", 
                          "BMI", "Alcohol", "Cigarette", "Highcholesterol", "Hypertension", 
                          "Diabete", "Fruit", "Vitamin", "Redmeat",
                          "Calcium", "Milk", 
                          "Type_PA_warking", "Type_PA_other", "Type_PA_strenuous", "Type_PA_LightDIY", 
                          "hExercise", "Type_PA_None", "blood_glucose", "BMD", 
                          "gpc_1", "gpc_2", "gpc_3", "gpc_4", "gpc_5", "gpc_6", "gpc_7", 
                          "gpc_8", "gpc_9", "gpc_10", "gpc_11", "gpc_12", "gpc_13", "gpc_14", 
                          "gpc_15", "gpc_16", "gpc_17", "gpc_18", "gpc_19", "gpc_20", "gpc_21", 
                          "gpc_22", "gpc_23", "gpc_24", "gpc_25", "gpc_26", "gpc_27", "gpc_28", 
                          "gpc_29", "gpc_30", "gpc_31", "gpc_32", "gpc_33", "gpc_34", "gpc_35", 
                          "gpc_36", "gpc_37", "gpc_38", "gpc_39", "gpc_40")
          Xvar_c <-c("Ethnicity", "Sex", "Age", "Education", "PA_time", "BMI",
                     "Alcohol", "Cigarette", "Highcholesterol", "Hypertension", 
                     "Diabete", "Calcium", "Milk")
          dt1$id = dt1$f.eid
          Xvar_e <- c("gpc_1", "gpc_2", "gpc_3", "gpc_4", 
                      "gpc_5", "gpc_6", "gpc_7", "gpc_8", "gpc_9", "gpc_10")
          Xvar <- c(Xvar_c, Xvar_e)
          # dt1[, Xvar] <- dt1[, Xvar] %>% apply(2, scale)
          Yvar <- "BMD"
          Avar <- "hExercise"
          X <- dt1[, Xvar] 
          X_c <- dt1[, Xvar_c]
          Y <- dt1[, Yvar] 
          A <- dt1[, Avar]
          N = nrow(dt1)
          #Influence function based group
          fit1 <- esATE(X, Y, A, "ipw")
          save(fit1, file = "0905_test/ipw.rdata")
          fit1 <- esATE(X, Y, A, "aipw")
          save(fit1, file = "0905_test/aipw.rdata")
          fit1 <- esATE(X, Y, A, "reg")
          save(fit1, file = "0905_test/reg.rdata")

          fit1_c <- esATE(X_c, Y, A, "ipw")
          save(fit1_c, file = "0905_test/ipw_c.rdata")
          fit1_c <- esATE(X_c, Y, A, "aipw")
          save(fit1_c, file = "0905_test/aipw_c.rdata")
          fit1_c <- esATE(X_c, Y, A, "reg")
          save(fit1_c, file = "0905_test/reg_c.rdata")
          if(method == "ipw"){
            load(file = "0905_test/ipw.rdata")
            load(file = "0905_test/ipw_c.rdata")
          }
          if(method == "aipw"){
            load(file = "0905_test/aipw.rdata")
            load(file = "0905_test/aipw_c.rdata")
          }
          if(method == "reg"){
            load(file = "0905_test/reg.rdata")
            load(file = "0905_test/reg_c.rdata")
          }
          cat(fit1_c$est/fit1$est)
          # cat(c(fit1$est - sqrt(fit1$ve)*1.96, fit1$est + sqrt(fit1$ve)*1.96))
          # cat(c(fit1_c$est - sqrt(fit1_c$ve)*1.96, fit1_c$est + sqrt(fit1_c$ve)*1.96))

          dt1$phi = fit1$infl
          dt1$phi_m = fit1_c$infl
          
          phi_bar = mean(dt1$phi)
          phi_bar_m = mean(dt1$phi_m)
          
          dt1$phi_a = abs(dt1$phi-phi_bar)
          dt1$phi_m_a = abs(dt1$phi_m-phi_bar_m)
          
          dt1 <- dt1[order(dt1$phi_m_a), ];
          R_M = rep(1:K, each = ceiling(N/K)) %>% as.factor()
          dt1$R_M = R_M[1:N]
          
          dt1 <- dt1[order(dt1$phi_a), ];
          R_T = rep(1:K, each = ceiling(N/K)) %>% as.factor()
          dt1$R_T = R_T[1:N]
          #####
          # uniform
          dt1$pi_u <- rep((n+n_m)/N, N)
          #####
          # best optimal oracle
          getpi(dt1$phi_a, (n+n_m), alpha)*(n+n_m) -> dt1$pi_b
          
          
          test_ukb <- function(dt1){
            #################
            # the pilot data/midphase data
            midphase_ind <- sample(1:N, size = n_m, replace = F)
            #### The interndata  $\mathcal P$
            dt2 <- dt1[midphase_ind,]
            #### The other first phase data 
            dt1_2 <- dt1[-midphase_ind,]
            N1_2 = nrow(dt1_2)
            
            x2 <- dt2[, Xvar]
            y2 <- dt2[, Yvar]
            A2 <- dt2[, Avar]
            esATE(x2, y2, A2, method) -> fit2
            dt2$phi_hat <- fit2$infl
            phi_hat_bar <- mean(dt2$phi_hat)
            dt2$phi_hat_a <- abs(dt2$phi_hat - phi_hat_bar)
            dt2 <- dt2[order(dt2$phi_hat_a), ]
            ###################
            # The clustering step: R_L
            # the orcale group result
            dt2 <- dt2[order(dt2$phi_hat_a), ]
            R_L = rep(1:K, each = ceiling(n_m/K)) %>% as.factor()
            dt2$R_L = R_L[1:n_m]
            te <- c()
            for(k in 1:K){
              te[k] <- mean((dt2[dt2$R_L==k,"phi_hat_a"]-phi_hat_bar)^2)
            }
            te2 <- te
            getpi_stra(sqrt(te2), N1_2, n, K, p = rep(1/K, K), alpha) -> pi_m
            ##############
            # The group optimal design
            rep(pi_m, each = ceiling(N1_2/K)) -> pi_1_0
            pi_1_0[1:N1_2] -> dt1_2$pi_1_0
            
            ####################
            #Lu's method
            ###################
            #Lu's group, K==8, without phi_m
            # (dt1$Type_PA_HeavyDIY +ifelse(dt1$edu>median(dt1$edu), 1, 0)*10+
            #    ifelse(dt1$age>median(dt1$age), 1, 0)*100) %>%
            #   strtoi(base = 2) ->R_O
            # dt1$R_O = R_O+1
            cut(dt1$BMD, breaks = quantile(dt1$BMD, probs = seq(0, 1, length.out = K + 1)),
                labels = 1:K) -> R_O
            dt1$R_O = R_O
            table(dt1$R_O)
            merge(dt1[, c("id", "R_O")], dt2) -> dt2
            merge(dt1[, c("id", "R_O")], dt1_2) -> dt1_2
            dt2 <- dt2[order(dt2$phi_hat_a), ]
            K_O = length(table(dt2$R_O))
            te <- c()
            for(k in 1:K_O){
              te[k] <- mean((dt2[dt2$R_O==k,"phi_hat"]-phi_hat_bar)^2)
            }
            te2 <- te
            getpi_stra(sqrt(te2), N1_2, n, K_O, p = rep(1/K_O, K_O), alpha = 0) -> pi_O
            pr_O <- data.frame(R_O = 1:K_O, p_O = pi_O)
            merge(dt1_2, pr_O, by="R_O", all = T) -> dt1_2_O
            ####################
            #Lu's group, K=8 with phi_m
            cut(dt1$phi_m, breaks = quantile(dt1$phi_m, probs = seq(0, 1, length.out = K + 1)),
                labels = 1:K) -> R_OM
            dt1$R_OM = R_OM
            table(dt1$R_OM)
            merge(dt1[, c("id", "R_OM")], dt2) -> dt2
            merge(dt1[, c("id", "R_OM")], dt1_2) -> dt1_2
            dt2 <- dt2[order(dt2$phi_hat_a), ]
            K_OM = length(table(dt2$R_OM))
            te <- c()
            for(k in 1:K_OM){
              te[k] <- mean((dt2[dt2$R_OM==k,"phi_hat"]-phi_hat_bar)^2)
            }
            te2 <- te
            getpi_stra(sqrt(te2), N1_2, n, K_OM, p = rep(1/K_OM, K_OM), alpha = 0) -> pi_OM
            pr_OM <- data.frame(R_OM = 1:K_OM, p_OM = pi_OM)
            merge(dt1_2, pr_OM, by="R_OM", all = T) -> dt1_2_OM
            
            ########################
            pr_M <- data.frame(R_M = 1:K, p_M = pi_m)
            merge(dt1_2, pr_M, by="R_M") -> dt1_2_M
            
            ######################
            ##The random foest method for classification, without "phi_m"
            formula = as.formula(paste0("R_L~", paste(c(Xvar_c, Yvar, Avar), collapse = "+")))
            fit_F <- randomForest(formula , data = dt2, proximity = T)
            dt1_2$R_F <- predict(fit_F, newdata = dt1_2)
            pr_F <- data.frame(R_F = 1:K, p_F = pi_m)
            merge(dt1_2, pr_F, by="R_F") -> dt1_2_F
            
            ######################
            ##The random foest method for classification
            formula = as.formula(paste0("R_L~", paste(c(Xvar_c, Yvar, Avar, "phi_m"), collapse = "+")))
            fit_FM <- randomForest(formula , data = dt2, proximity = T)
            dt1_2$R_FM <- predict(fit_FM, newdata = dt1_2)
            pr_FM <- data.frame(R_FM = 1:K, p_FM = pi_m)
            merge(dt1_2, pr_FM, by="R_FM") -> dt1_2_FM
            
            #######################
            ##Tree method for classification
            formula = as.formula(paste0("R_L~", paste(c(Xvar_c, Yvar, Avar), collapse = "+")))
            fit_Tr <- rpart(formula, data = dt2, cp = 0.005)
            dt1_2$R_Tr <- predict(fit_Tr, newdata = dt1_2, type = "class")
            pr_Tr <- data.frame(R_Tr = 1:K, p_Tr = pi_m)
            merge(dt1_2, pr_Tr, by="R_Tr") -> dt1_2_Tr
            
            #######################
            ##Tree method for classification
            formula = as.formula(paste0("R_L~", paste(c(Xvar_c, Yvar, Avar, "phi_m"), collapse = "+")))
            fit_TrM <- rpart(formula, data = dt2, cp = 0.005)
            dt1_2$R_TrM <- predict(fit_TrM, newdata = dt1_2, type = "class")
            pr_TrM <- data.frame(R_TrM = 1:K, p_TrM = pi_m)
            merge(dt1_2, pr_TrM, by="R_TrM") -> dt1_2_TrM
            
            #############
            sum((dt1_2$phi_a)^2/(n/N1_2))/N/N
            sum((dt1_2_M$phi_a)^2/dt1_2_M$p_M)/N/N
            sum((dt1_2_F$phi_a)^2/dt1_2_F$p_F)/N/N
            sum((dt1_2_FM$phi_a)^2/dt1_2_FM$p_FM)/N/N
            sum((dt1_2_Tr$phi_a)^2/dt1_2_Tr$p_Tr)/N/N
            sum((dt1_2_TrM$phi_a)^2/dt1_2_TrM$p_TrM)/N/N
            sum((dt1_2_O$phi_a)^2/dt1_2_O$p_O, na.rm = T)/N/N
            sum((dt1_2_OM$phi_a)^2/dt1_2_OM$p_OM, na.rm = T)/N/N
            
            get_tau_hat(dt1, dt1$pi_u, method) -> fit_u
            get_tau_hat(dt1, dt1$pi_b, method) -> fit_b
            get_tau_hat2(dt1_2_M, dt2, dt1_2_M$p_M, method) -> fit_M
            get_tau_hat2(dt1_2_F, dt2, dt1_2_F$p_F, method) -> fit_F
            get_tau_hat2(dt1_2_FM, dt2, dt1_2_FM$p_FM, method) -> fit_FM
            get_tau_hat2(dt1_2_Tr, dt2, dt1_2_Tr$p_Tr, method) -> fit_Tr
            get_tau_hat2(dt1_2_TrM, dt2, dt1_2_TrM$p_TrM, method) -> fit_TrM
            get_tau_hat2(dt1_2_O, dt2, dt1_2_O$p_O, method) -> fit_O
            get_tau_hat2(dt1_2_OM, dt2, dt1_2_OM$p_OM, method) -> fit_OM
            return(c(fit_u, fit_b, fit_M, fit_F, fit_FM, fit_Tr, fit_TrM, fit_O, fit_OM))
          }
          get_tau_hat <- function(dt, pi, method){
            dt$pi = pi
            N = nrow(dt)
            dt$V = rbinom(N, 1, dt$pi)
            dt_Z= filter(dt, dt$V==1)
            X <- dt_Z[, Xvar]
            Y <- dt_Z[, Yvar]
            A <- dt_Z[, Avar]
            esATE_w(X, Y, A, w = 1/dt_Z$pi, method, N) -> tau_hat
            return(tau_hat)
          }
          get_tau_hat2 <- function(dt1_2, dt2, pi1_2, method){
            dt1_2$pi = pi1_2;
            N1_2 = nrow(dt1_2);n_m = nrow(dt2); N = n_m + N1_2
            dt1_2$V = rbinom(N1_2, 1, dt1_2$pi)
            dt1_2= filter(dt1_2, dt1_2$V==1)
            esATE_w(dt1_2[, Xvar], dt1_2[, Yvar], dt1_2[, Avar], w = 1/dt1_2$pi, method, N) ->fit1_2
            esATE(dt2[, Xvar], dt2[, Yvar], dt2[, Avar], method) -> fit2
            (fit1_2[1]/fit1_2[2] + fit2$est/fit2$ve)/(1/fit1_2[2] + 1/fit2$ve) -> tau_hat
            (fit1_2[2]*fit2$ve)/(fit1_2[2]+fit2$ve) -> ve_hat
            cb_2 <- c(tau_hat, ve_hat)
            return(c(cb_2))
          }
          # test_ukb(dt1)
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
                                  "nnet", "caret","randomForest",
                                  "MASS", "glmnet", "rpart"))%dopar%{
                                    return(test_ukb(dt1))}
          stopCluster(cls)
          
          endtime <- Sys.time()
          cat( "\nRunning time: ", endtime-starttime, "\n")
          
          
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
          isalpha <- function(str){
            str <- ifelse(grepl("[a-d|f-z|A-Z]", str), NA, str) %>% unlist %>% as.numeric()
            return(str)
          }
          apply(result, 2, isalpha) -> result2
          apply(result2, 2, function(x){x[abs(x)>1] = NA;return(x)}) -> result2
          rbind(result_td(result2[, 1], result2[, 2], tau_true = fit1$est),
                result_td(result2[, 3], result2[, 4], tau_true = fit1$est),
                result_td(result2[, 5], result2[, 6], tau_true = fit1$est),
                result_td(result2[, 7], result2[, 8], tau_true = fit1$est),
                result_td(result2[, 9], result2[, 10], tau_true = fit1$est),
                result_td(result2[, 11], result2[, 12], tau_true = fit1$est),
                result_td(result2[, 13], result2[, 14], tau_true = fit1$est),
                result_td(result2[, 15], result2[, 16], tau_true = fit1$est),
                result_td(result2[, 17], result2[, 18], tau_true = fit1$est)) -> re
          rownames(re) <- c( "u","b", "m", "f", "fm", "tr", "trm", "o", "om")
          if(! dir.exists(paste0("result_UKB/",name))){
            dir.create(paste0("result_UKB/",name))
          }
          save(result, file = paste0("result_UKB/",name, "/", method, "_n_m=",
                                     n_m, "_n=",n,  "_K=", K,"_alpha=", alpha,".rdata"))
          re %>% write.csv(file = paste0("result_UKB/",name, "/", method, "_n_m=",
                                         n_m, "_n=",n,  "_K=", K,"_alpha=", alpha,".csv"))
        }
      }
    }
  }
}



