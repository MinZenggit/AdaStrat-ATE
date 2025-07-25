#####################
# Appendix A.5, Table S3
rm(list = ls())
source("somefunction16.R")
source("esATE.R")
source("esATE_w3.R")
library(dplyr)
library(randomForest)
library(nnet)
library(caret)
library(MASS)
library(foreach)
library(mice)
require(doSNOW)
set.seed(2023)
######
# Influence function based group
# fit1 <- esATE(X, Y, A, "ipw")
# save(fit1, file = "ipw.rdata")
# fit1 <- esATE(X, Y, A, "aipw")
# save(fit1, file = "aipw.rdata")
# fit1 <- esATE(X, Y, A, "reg")
# save(fit1, file = "reg.rdata")
# 
# fit1_c <- esATE(X_c, Y, A, "ipw")
# save(fit1_c, file = "ipw_c.rdata")
# fit1_c <- esATE(X_c, Y, A, "aipw")
# save(fit1_c, file = "aipw_c.rdata")
# fit1_c <- esATE(X_c, Y, A, "reg")
# save(fit1_c, file = "reg_c.rdata")

method_set <- c("aipw")
n_m_set = c(2000)
n_set = c(3000)
K_set = c(8, 12, 16)
alpha_set = c(0.3)
cores = 4
trials = 4
name = "impute_test_moregenes"
if(! dir.exists(paste0("result_UKB/",name))){
  dir.create(paste0("result_UKB/",name))
}

n_m = 1000; n = 2000; method = "aipw"; alpha = 0.2; K = 8
for(method in method_set){
  for(n_m in n_m_set){
    for(n in n_set){
      for(K in K_set){
        for(alpha in alpha_set){
          cat(" Method=", method, ", n_m=", n_m,", n=", n, ", K=", K, ", alpha=",alpha, "\n")
          load(file = "dt.rdata")
          dt1 <- dt 
          Xvar_c <-c("Ethnicity", "sex", "age", "edu", "Summins_activity", 
                     "BMI", "Alcohol", "Smoking_status", 
                     "highcholesterol", "hypertension", "Diabete",
                     "vitamin", "milk_new")
          dt1$id = dt1$f.eid
          Xvar_e <- c("gpc_1", "gpc_2", "gpc_3", "gpc_4", 
                      "gpc_5", "gpc_6", "gpc_7", "gpc_8", "gpc_9", "gpc_10",
                      "gpc_11", "gpc_12", "gpc_13", "gpc_14", "gpc_15", "gpc_16",
                      "gpc_17", "gpc_18", "gpc_19", "gpc_20", "gpc_21", "gpc_22",
                      "gpc_23", "gpc_24", "gpc_25")
          Xvar <- c(Xvar_c, Xvar_e)
          dt1[, Xvar] <- dt1[, Xvar] %>% apply(2, scale)
          Yvar <- "bone_density"
          Avar <- "Type_PA_HeavyDIY"
          X <- dt1[, Xvar] 
          X_c <- dt1[, Xvar_c]
          Y <- dt1[, Yvar] 
          A <- dt1[, Avar]
          N = nrow(dt1)
          ###Influence function based group
          # fit1 <- esATE(X, Y, A, "ipw")
          # save(fit1, file = paste0("result_UKB/", name, "/", "ipw.rdata"))
          # fit1 <- esATE(X, Y, A, "aipw")
          # save(fit1, file = paste0("result_UKB/", name, "/", "aipw.rdata"))
          # fit1 <- esATE(X, Y, A, "reg")
          # save(fit1, file = paste0("result_UKB/", name, "/", "reg.rdata"))
          # 
          # fit1_c <- esATE(X_c, Y, A, "ipw")
          # save(fit1_c, file = paste0("result_UKB/", name, "/", "ipw_c.rdata"))
          # fit1_c <- esATE(X_c, Y, A, "aipw")
          # save(fit1_c, file = paste0("result_UKB/", name, "/", "aipw_c.rdata"))
          # fit1_c <- esATE(X_c, Y, A, "reg")
          # save(fit1_c, file = paste0("result_UKB/", name, "/", "reg_c.rdata"))
          if(method == "ipw"){
            load(file = paste0("result_UKB/", name, "/", "ipw.rdata"))
            load(file = paste0("result_UKB/", name, "/", "ipw_c.rdata"))
          }
          if(method == "aipw"){
            load(file = paste0("result_UKB/", name, "/", "aipw.rdata"))
            load(file = paste0("result_UKB/", name, "/", "aipw_c.rdata"))
          }
          if(method == "reg"){
            load(file = paste0("result_UKB/", name, "/", "reg.rdata"))
            load(file = paste0("result_UKB/", name, "/", "reg_c.rdata"))
          }
          
          cat(c(fit1$est - sqrt(fit1$ve)*1.96, fit1$est + sqrt(fit1$ve)*1.96))
          cat(c(fit1_c$est - sqrt(fit1_c$ve)*1.96, fit1_c$est + sqrt(fit1_c$ve)*1.96))
          
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
          # best optimal
          getpi(dt1$phi_a, (n+n_m), alpha)*(n+n_m) -> dt1$pi_b
          
          
          test_ukb <- function(dt1){
            #################
            # the middle phase data
            midphase_ind <- sample(1:N, size = n_m, replace = F)
            #### The interndata
            dt2 <- dt1[midphase_ind,]
            #### The otherdata
            dt1_2 <- dt1[-midphase_ind,]
            fid_mis <- dt1$f.eid[midphase_ind]
            fid_misn <- dt1$f.eid[-midphase_ind]
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
            #Lu's group, K==8
            # (dt1$Type_PA_HeavyDIY +ifelse(dt1$edu>median(dt1$edu), 1, 0)*10+
            #    ifelse(dt1$age>median(dt1$age), 1, 0)*100) %>%
            #   strtoi(base = 2) ->R_O
            # dt1$R_O = R_O+1
            cut(dt1$bone_density, breaks = quantile(dt1$bone_density, probs = seq(0, 1, length.out = K + 1)),
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
            
            
            ######################
            ##The nnet method for classification
            digit.x = dt2[, c(Xvar_c, Yvar, Avar)]
            digit.y = dt2[, "R_L"]
            digit.ml <- train(x=digit.x, y=digit.y,
                              method="nnet",
                              tuneGrid=expand.grid(
                                # 10个隐藏神经元
                                .size=10,
                                # 衰变率
                                .decay=0.1
                              ),
                              trControl=trainControl(method="none"),
                              # 最大权重数量
                              MaxNWts=10000,
                              # 最大迭代次数
                              maxit=1000,
                              verbose = FALSE)
            dt1_2$R_Z <- predict(digit.ml, newdata = dt1_2)
            pr_Z <- data.frame(R_Z = 1:K, p_Z = pi_m)
            merge(dt1_2, pr_Z, by="R_Z") -> dt1_2_Z
            
            ########################
            pr_M <- data.frame(R_M = 1:K, p_M = pi_m)
            merge(dt1_2, pr_M, by="R_M") -> dt1_2_M
            
            ######################
            ##The nnet method for classification
            starttime <- Sys.time()
            digit.x = dt2[, c(Xvar_c, Yvar, Avar, "phi_m")]
            digit.y = dt2[, "R_L"]
            digit.zm <- train(x=digit.x, y=digit.y,
                              method="nnet",
                              tuneGrid=expand.grid(
                                # 10个隐藏神经元
                                .size=10,
                                # 衰变率
                                .decay=0.1
                              ),
                              trControl=trainControl(method="none"),
                              # 最大权重数量
                              MaxNWts=10000,
                              # 最大迭代次数
                              maxit=1000,
                              verbose = FALSE)
            dt1_2$R_ZM <- predict(digit.zm, newdata = dt1_2)
            pr_ZM <- data.frame(R_ZM = 1:K, p_ZM = pi_m)
            merge(dt1_2, pr_ZM, by="R_ZM") -> dt1_2_ZM
            endtime <- Sys.time()
            time_zm <- endtime - starttime
            
            ######################
            ##The random foest method for classification
            starttime <- Sys.time()
            formula = as.formula(paste0("R_L~", paste(c(Xvar_c, Yvar, Avar, "phi_m"), collapse = "+")))
            fit_f <- randomForest(formula , data = dt2, proximity = T)
            dt1_2$R_F <- predict(fit_f, newdata = dt1_2)
            pr_F <- data.frame(R_F = 1:K, p_F = pi_m)
            merge(dt1_2, pr_F, by="R_F") -> dt1_2_F
            endtime <- Sys.time()
            time_f <- endtime - starttime
            
            ##################
            # impute X_e methods
            starttime <- Sys.time()
            dt1_2[, Xvar_e] <- NA
            rbind(dt1_2[,c("id", Xvar, Yvar, Avar)], dt2[,c("id", Xvar, Yvar, Avar)]) -> dt11
            imputed_data <- mice(dt11[, c("id", Xvar, Yvar, Avar)], 
                                 m = 1,     
                                 maxit = 10,  
                                 printFlag = FALSE)
            dt1_imp <- complete(imputed_data, 1)
            esATE(dt1_imp[, Xvar], dt1_imp[, Yvar],
                  dt1_imp[, Avar], method) -> fit11
            dt1_imp$phi_impute <- fit11$infl
            merge(dt1_imp[, c("id", "phi_impute")], dt1, by = "id") -> dt11
            dt11 <- dt11[order(abs(dt11$phi_impute)), ];
            R_imp = rep(1:K, each = ceiling(N/K)) %>% as.factor()
            dt11$R_imp = R_imp[1:N]
            te <- c()
            for(k in 1:K){
              aa <- dt11 %>% filter(id %in% fid_mis, R_imp == k)
              te[k] <- mean(aa$phi_impute^2)
            }
            te2 <- te
            getpi_stra(sqrt(te2), N1_2, n, K, p = rep(1/K, K), alpha = 0) -> pi_imp
            pr_imp <- data.frame(R_imp = 1:K, p_imp = pi_imp)
            merge(dt11, pr_imp, by="R_imp") -> dt11_imp
            dt11_imp %>% filter(id %in% fid_misn) -> dt1_2_imp
            endtime <- Sys.time()
            time_imp <- endtime - starttime
            #############
            sum((dt1_2$phi_a)^2/(n/N1_2))/N/N
            sum((dt1_2_Z$phi_a)^2/dt1_2_Z$p_Z)/N/N
            sum((dt1_2_M$phi_a)^2/dt1_2_M$p_M)/N/N
            sum((dt1_2_ZM$phi_a)^2/dt1_2_ZM$p_ZM)/N/N
            sum((dt1_2_F$phi_a)^2/dt1_2_F$p_F)/N/N
            sum((dt1_2_imp$phi_a)^2/dt1_2_imp$p_imp)/N/N
            sum((dt1_2_O$phi_a)^2/dt1_2_O$p_O, na.rm = T)/N/N
            
            get_tau_hat(dt1, dt1$pi_u, method) -> fit_u
            get_tau_hat(dt1, dt1$pi_b, method) -> fit_b
            get_tau_hat2(dt1_2_Z, dt2, dt1_2_Z$p_Z, method) -> fit_Z
            get_tau_hat2(dt1_2_M, dt2, dt1_2_M$p_M, method) -> fit_M
            get_tau_hat2(dt1_2_ZM, dt2, dt1_2_ZM$p_ZM, method) -> fit_ZM
            get_tau_hat2(dt1_2_F, dt2, dt1_2_F$p_F, method) -> fit_F
            get_tau_hat2(dt1_2_imp, dt2, dt1_2_imp$p_imp, method) -> fit_imp
            get_tau_hat2(dt1_2_O, dt2, dt1_2_O$p_O, method) -> fit_O
            return(list(est <- c(fit_u, fit_b, fit_Z, fit_M, fit_ZM, fit_F, fit_imp, fit_O),
                        time <- as.numeric(c(time_zm, time_f, time_imp))))
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
                                  "MASS", "glmnet", "mice"))%dopar%{
                                    re1 <- test_ukb(dt1)
                                    return(c(re1[[1]], re1[[2]]))}
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
          times <- c( mean(result2[, 17], na.rm = T),  
                      mean(result2[, 18], na.rm = T), 
                      mean(result2[, 19], na.rm = T))
          print(times)
          apply(result2, 2, function(x){x[abs(x)>1] = NA;return(x)}) -> result2
          rbind(result_td(result2[, 1], result2[, 2], tau_true = fit1$est),
                result_td(result2[, 3], result2[, 4], tau_true = fit1$est),
                result_td(result2[, 5], result2[, 6], tau_true = fit1$est),
                result_td(result2[, 7], result2[, 8], tau_true = fit1$est),
                result_td(result2[, 9], result2[, 10], tau_true = fit1$est),
                result_td(result2[, 11], result2[, 12], tau_true = fit1$est),
                result_td(result2[, 13], result2[, 14], tau_true = fit1$est),
                result_td(result2[, 15], result2[, 16], tau_true = fit1$est)) -> re
          rownames(re) <- c( "u","b", "z", "m", "zm", "f", "imp", "o")
          if(! dir.exists(paste0("result_UKB/",name))){
            dir.create(paste0("result_UKB/",name))
          }
          save(result, file = paste0("result_UKB/",name, "/", method, "_n_m=",
                                     n_m, "_n=",n,  "_K=", K,"_alpha=", alpha,".rdata"))
          re %>% write.csv(file = paste0("result_UKB/",name, "/", method, "_n_m=",
                                         n_m, "_n=",n,  "_K=", K,"_alpha=", alpha,".csv"))
          write.csv(times, file = paste0("result_UKB/",name, "/", method, "_n_m=",
                                         n_m, "_n=",n,  "_K=", K,"_alpha=", alpha,"_time.csv"))
        }
      }
    }
  }
}



