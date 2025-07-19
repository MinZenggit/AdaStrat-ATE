##################
# If the propensity score model is wrong
source("esATE.R")
source("esATE_w2.R")
library(stringr)
isalpha <- function(str){
  str <- ifelse(grepl("[a-zA-Z]", str), NA, str) %>% unlist %>% as.numeric()
  return(str)
}
getpi <- function(phi_a, sn, alpha = 0.05){
  ## phi_a is sorted abs(phi)
  pi_a = phi_a/sum(phi_a)
  if(max(pi_a)<=(1/sn)){
    pi_a = (1-alpha)*pi_a + alpha * mean(pi_a)
    return(pi_a)
  }else{
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
    pi = pmin(H, phi_a)
    pi = pi/sum(pi)
    pi <- ifelse(pi>(1/sn), 1/sn, pi)
    pi = (1-alpha)*pi + alpha * mean(pi)
    return(pi)
  }
}
getpi_stra <- function(t, N, n, K, p, alpha = 0.05){ 
  ## phi_a is sorted abs(phi)
  pi_a = t/sum(t*p)*n/N
  if(max(pi_a)<=(1)){
    pi_a = (1-alpha)*pi_a + alpha * mean(pi_a)
    return(pi_a)
  }else{
    H <- re <- c()
    for(g in 1:(n*K/N-1)){
      H[g] = sum(t[1:(K-g)]*p[1:(K-g)])/(n/N-sum(p[(K-g+1): K]))
    }
    for(g in 1:(n*K/N-1)){
      (H[g] <= t[K-g+1]& H[g] >= t[K-g]) -> re[g]
    }
    g = which(re==T)[1]
    H[g] -> H
    pi = t/H
    pi <- ifelse(pi>1, 1, pi)
    pi = (1-alpha)*pi + alpha * mean(pi)
    return(pi)
  }
}

## 
# the true propensity score was not linear!!!, wrong propensity score model
Generatedt <- function(N){
  id <- 1:N
  X1 = runif(N, 0, 2)
  U1 = 0.5 + 0.5*X1 - 2*sin(X1) + 2*sign(sin(5*X1)) + rnorm(N, sd = 2)
  Y0 = - X1 - U1 + rnorm(N, sd = 2)
  Y1 = - X1 + 4*U1 + rnorm(N, sd = 2)
  # non linear propensity score model!!!!
  e = plogis(1 - 0.5*X1-0.5*U1 + U1*X1 - sin(X1) + sin(U1))
  Z <- rbinom(N, 1, e)
  Y = Z*Y1 + (1-Z)*Y0
  # ### additive form measurement
  Ym_a = 2*Y + X1 + rnorm(N) 
  # ### multiple form measurement
  Ym_m = Y * rnorm(N) 
  ### nonlinear form measurement
  Ym_s = plogis(Y/5 + rnorm(N, 0, 0.1))
  dt1 <- cbind(id, X1, U1, Y0, Y1, e, Z, Y, Ym_s, Ym_a, Ym_m) %>% as.data.frame()
  return(dt1)
}

# This function pick up the second phase data from first phase data
# by the design based probabilities and give an estimator of ATE.
get_tau_hat <- function(dt, pi, method){
  dt$pi = pi
  N = nrow(dt)
  dt$V = rbinom(N, 1, dt$pi)
  dt_Z= filter(dt, dt$V==1)
  esATE_w(cbind(dt_Z$X1, dt_Z$U1), dt_Z$Y, dt_Z$Z, w = 1/dt_Z$pi, method, N) -> tau_hat
  return(tau_hat)
}
#####
## Combine the two sub phase data
# get_tau_hat2 <- function(dt1_2, dt2, pi1_2, pi2, method){
#   dt1_2$pi = pi1_2; dt2$pi2 = pi2
#   N1_2 = nrow(dt1_2);n_m = nrow(dt2); N = n_m + N1_2
#   dt1_2$V = rbinom(N1_2, 1, dt1_2$pi)
#   dt1_2= filter(dt1_2, dt1_2$V==1)
#   dt1_2$pi*(1-n_m/N) + (n_m)/N -> dt1_2$pi_ad
#   dt2$pi_ad = dt2$pi2*(1-n_m/N)+ (n_m)/N
#   dtt <- rbind(dt1_2[, c("X1", "U1", "Y", "Z", "pi_ad")],
#                dt2[, c("X1", "U1", "Y", "Z", "pi_ad")])
#   esATE_w(cbind(dtt$X1, dtt$U1), dtt$Y, dtt$Z, w = 1/dtt$pi_ad, method, N) -> tau_hat
#   return(tau_hat)
# }
get_tau_hat2 <- function(dt1_2, dt2, pi1_2, method){
  dt1_2$pi = pi1_2;
  N1_2 = nrow(dt1_2);n_m = nrow(dt2); N = n_m + N1_2
  dt1_2$V = rbinom(N1_2, 1, dt1_2$pi)
  dt1_2= filter(dt1_2, dt1_2$V==1)
  esATE_w(cbind(dt1_2$X1, dt1_2$U1), dt1_2$Y, dt1_2$Z, w = 1/dt1_2$pi, method, N) ->fit1_2
  esATE(cbind(dt2$X1, dt2$U1), dt2$Y, dt2$Z, method) -> fit2
  (fit1_2[1]/fit1_2[2] + fit2$est/fit2$ve)/(1/fit1_2[2] + 1/fit2$ve) -> tau_hat
  (fit1_2[2]*fit2$ve)/(fit1_2[2]+fit2$ve) -> ve_hat
  cb_2 <- c(tau_hat, ve_hat)
  return(c(cb_2))
}

testf <- function(N, n_m, n, K, method, alpha){
  Generatedt(N) -> dt1 -> dt
  #######
  # Outcome dependent group
  te <- dt1 %>% dplyr::select(id, Y)
  te[order(te$Y), ] -> te
  R_O = rep(1:K, each = ceiling(N/K)) %>% as.factor()
  te$R_O =  R_O[1:N]
  dt1$R_O = te[order(te$id), "R_O"]
  ######
  # Influence function based group
  x = cbind(dt1$X1, dt1$U1); y = dt1$Y; A = dt1$Z
  esATE(x, y, A, method) -> fit1
  esATE(dt1$X1, y, A, method) -> fit1_c
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
  
  ############
  # the middle phase data
  midphase_ind <- sample(1:N, size = n_m, replace = F)
  #### The interndata
  dt2 <- dt1[midphase_ind,]
  #### The otherdata
  dt1_2 <- dt1[-midphase_ind,]
  N1_2 = nrow(dt1_2)
  
  x2 = cbind(dt2$X1, dt2$U1); y2 = dt2$Y; A2 = dt2$Z
  esATE(x2, y2, A2, method) -> fit2
  dt2$phi_hat <- fit2$infl
  phi_hat_bar <- mean(dt2$phi_hat)
  dt2$phi_hat_a <- abs(dt2$phi_hat - phi_hat_bar)
  dt2 <- dt2[order(dt2$phi_hat_a), ]
  R_L = rep(1:K, each = ceiling(n_m/K)) %>% as.factor()
  dt2$R_L = R_L[1:n_m]
  te <- c()
  for(k in 1:K){
    te[k] <- mean((dt2[dt2$R_L==k,"phi_hat"]-phi_hat_bar)^2)
  }
  te2 <- te
  getpi_stra(sqrt(te2), N1_2, n, K, p = rep(1/K, K), alpha) -> pi_m
  dt2$pi_m = NA
  for(k in 1:K){
    dt2[dt2$R_L==k, "pi_m"] = pi_m[k]
  }
  
  ##############
  # The group optimal design
  rep(pi_m, each = ceiling(N1_2/K)) -> pi_1_0
  pi_1_0[1:N1_2] -> dt1_2$pi_1_0
  
  
  ##############
  #Outcome dependent group 
  te <- c()
  for(k in 1:K){
    te[k] <- mean((dt2[dt2$R_O==k,"phi_hat"]-phi_hat_bar)^2)
  }
  te2 <- te
  getpi_stra(sqrt(te2), N1_2, n, K, p = rep(1/K, K), alpha) -> pi_O
  
  pr_O <- data.frame(R_O = 1:K, p_O = pi_O)
  merge(dt1_2, pr_O, by="R_O") -> dt1_2_O
  merge(dt2, pr_O, by="R_O") -> dt2_O
  
  #######
  # M1 model
  te <- c()
  for(k in 1:K){
    te[k] <- mean((dt2[dt2$R_M==k,"phi_hat"]-phi_hat_bar)^2)
  }
  te2 <- te
  getpi_stra(sqrt(te2), N1_2, n, K, p = rep(1/K, K), alpha) -> p_M
  pr_M <- data.frame(R_M = 1:K, p_M = p_M)
  merge(dt1_2, pr_M, by="R_M") -> dt1_2_M1
  
  #######
  # M2 model
  pr_M <- data.frame(R_M = 1:K, p_M = pi_m)
  merge(dt1_2, pr_M, by="R_M") -> dt1_2_M2
  ###########
  # Z model
  digit.x = dt2[, c("X1", "Y", "Z")]
  digit.y = dt2[, "R_L"]
  digit.ml <- train(x=digit.x, y=digit.y,
                    method="nnet",
                    tuneGrid=expand.grid(
                      # 5个隐藏神经元
                      .size=10,
                      # 衰变率
                      .decay=0.1
                    ),
                    trControl=trainControl(method="none"),
                    # 最大权重数量
                    MaxNWts=10000,
                    # 最大迭代次数
                    maxit=500)
  dt1_2$R_Z <- predict(digit.ml, newdata = dt1_2)
  pr_Z <- data.frame(R_Z = 1:K, p_Z = pi_m)
  merge(dt1_2, pr_Z, by="R_Z") -> dt1_2_Z
  
  ##########
  # ZM model
  digit.x = dt2[, c("X1", "Y", "Z", "phi_m")]
  digit.y = dt2[, "R_L"]
  digit.zm <- train(x=digit.x, y=digit.y,
                    method="nnet",
                    tuneGrid=expand.grid(
                      # 5个隐藏神经元
                      .size=10,
                      # 衰变率
                      .decay=0.1
                    ),
                    trControl=trainControl(method="none"),
                    # 最大权重数量
                    MaxNWts=10000,
                    # 最大迭代次数
                    maxit=500)
  dt1_2$R_ZM <- predict(digit.zm, newdata = dt1_2)
  pr_ZM <- data.frame(R_ZM = 1:K, p_ZM = pi_m)
  merge(dt1_2, pr_ZM, by="R_ZM") -> dt1_2_ZM
  
  
  ######################
  ##The random foest method for classification: F method
  fit_f <- randomForest(R_L~X1+Y+Z , data = dt2, proximity = T)
  dt1_2$R_F <- predict(fit_f, newdata = dt1_2)
  pr_F <- data.frame(R_F = 1:K, p_F = pi_m)
  merge(dt1_2, pr_F, by="R_F") -> dt1_2_F
  
  ######################
  ##The random foest method for classification: FM method
  fit_f <- randomForest(R_L~X1+Y+Z+phi_m, data = dt2, proximity = T)
  dt1_2$R_F <- predict(fit_f, newdata = dt1_2)
  pr_F <- data.frame(R_F = 1:K, p_FM = pi_m)
  merge(dt1_2, pr_F, by="R_F") -> dt1_2_FM
  
  # ######
  sum((dt1_2$phi_a)^2/(n/(N-n_m)))/N/N 
  sum((dt1_2_M1$phi_a)^2/dt1_2_M1$p_M)/N/N 
  sum((dt1_2_M2$phi_a)^2/dt1_2_M2$p_M)/N/N
  sum((dt1_2_Z$phi_a)^2/dt1_2_Z$p_Z)/N/N 
  sum((dt1_2_ZM$phi_a)^2/dt1_2_ZM$p_ZM)/N/N 
  sum((dt1_2_F$phi_a)^2/dt1_2_F$p_F)/N/N
  sum((dt1_2_FM$phi_a)^2/dt1_2_FM$p_FM)/N/N
  sum((dt1_2_O$phi_a)^2/dt1_2_O$p_O)/N/N
  ########
  get_tau_hat(dt1, dt1$pi_u, method) -> fit_u
  
  get_tau_hat(dt1, dt1$pi_b, method) -> fit_b
  
  get_tau_hat2(dt1_2_M1, dt2, dt1_2_M1$p_M, method) -> fit_M1
  
  get_tau_hat2(dt1_2_M2, dt2, dt1_2_M2$p_M, method) -> fit_M2
  
  get_tau_hat2(dt1_2_Z, dt2, dt1_2_Z$p_Z, method) -> fit_Z
  
  get_tau_hat2(dt1_2_ZM, dt2, dt1_2_ZM$p_ZM, method) -> fit_ZM
  
  get_tau_hat2(dt1_2_F, dt2, dt1_2_F$p_F, method) -> fit_F
  
  get_tau_hat2(dt1_2_FM, dt2, dt1_2_FM$p_FM, method) -> fit_FM
  
  get_tau_hat2(dt1_2_O, dt2, dt1_2_O$p_O, method) -> fit_O
  
  
  return(c(fit_u, fit_b, fit_M1, fit_M2, fit_Z, fit_ZM, fit_F, fit_FM, fit_O))
}

result_td <- function(re){
  re = cbind(apply(re, 2, mean, na.rm = T)-0.50,
             apply(re, 2, var, na.rm = T),
             (apply(re, 2, mean, na.rm = T)-0.50)^2 +apply(re, 2, var, na.rm = T) )
  colnames(re) <- c("bias", "var", "mse")
  return(re)
}
# testf(N, n_m, n, K, method, alpha)

