################################
# To evaluate the estimation of efficiency gain.

source("esATE.R")
source("esATE_w2.R")
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

Generatedt <- function(N){
  id <- 1:N
  X1 = runif(N, 0, 2)
  U1 = 0.5 + 0.5*X1 - 2*sin(X1) + 2*sign(sin(5*X1)) + rnorm(N, sd = 2)
  Y0 = - X1 - U1 + rnorm(N, sd = 2)
  Y1 = - X1 + 4*U1 + rnorm(N, sd = 2)
  e = plogis(1 - 0.5*X1-0.5*U1)
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

N = 10000; n_m = 1000; n = 1000; K = 10; method = "aipw"
testf <- function(N, n_m, n, K, method){
  Generatedt(N) -> dt1
  #############
  # Influence function based group
  x = cbind(dt1$X1, dt1$U1); y = dt1$Y; A = dt1$Z
  esATE(x, y, A, method) -> fit1
  dt1$phi = fit1$infl
  phi_bar = mean(dt1$phi)
  dt1$phi_a = abs(dt1$phi-phi_bar)
  dt1 <- dt1[order(dt1$phi_a), ]
  R_T = rep(1:K, each = ceiling(N/K)) %>% as.factor()
  dt1$R_T = R_T[1:N]
  #####
  # uniform
  dt1$pi_u <- rep(n/N, N)
  
  ######
  dt2 <- Generatedt(n_m)
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
    te[k] <- mean((dt2[dt2$R_L==k,"phi_hat_a"]-phi_hat_bar)^2)
  }
  te2 <- te
  getpi_stra(sqrt(te2), N, n, K, p = rep(1/K, K)) -> pi_m
########################################
##用交叉验证计算 efficiency gain
  ind = sample(1:n_m, size = n_m, replace = F)
  cv_f = 5
  thvar_z_hat <- thvar_f_hat<-thvar_u_hat <- c()
  for(cv_i in 1:cv_f){
    cat("fold ",cv_i, ", " )
    n_m/cv_f -> cv_n
    test_ind <- ind[(cv_n*cv_i-cv_n+1):(cv_n*cv_i)]
    dt2_train <- dt2[-test_ind, ]
    dt2_test <- dt2[test_ind, ]
    digit.x = dt2_train[, c("X1", "Y", "Z")]
    digit.y = dt2_train[, "R_L"]
    (digit.ml.test <- train(x=digit.x, y=digit.y,
                            method="nnet",
                            tuneGrid=expand.grid(
                              # 10个隐藏神经元
                              .size=10,
                              # 衰变率
                              .decay=0.1
                            ),
                            trControl=trainControl(method="none", verboseIter = T),
                            # 最大权重数量
                            MaxNWts=10000,
                            # 最大迭代次数
                            maxit=500)) %>% capture.output() -> output
    dt2_test$R_test <- predict(digit.ml.test, newdata = dt2_test)
    pr_test <- data.frame(R_test = 1:K, p_test = pi_m)
    merge(dt2_test, pr_test, by="R_test")  -> dt2_test_new
    mean((dt2_test_new$phi_hat- phi_hat_bar)^2/dt2_test_new$p_test)/N -> thvar_z_hat[cv_i]
    
    # fit_f <- randomForest(R_L~X1+Y+Z, data = dt2_train, proximity = T)
    # dt2_test$R_test <- predict(fit_f, newdata = dt2_test)
    # pr_test <- data.frame(R_test = 1:K, p_test = pi_m)
    # merge(dt2_test, pr_test, by="R_test") -> dt2_test_new
    # mean((dt2_test_new$phi_hat)^2/dt2_test_new$p_test)/N -> thvar_f_hat[cv_i]
    
    mean((dt2_test_new$phi_hat)^2/(n/N))/N -> thvar_u_hat[cv_i]
  }

  # (thvar_u_hat/thvar_f_hat) %>% mean() -> r_f_hat
  (thvar_u_hat/thvar_z_hat) %>% mean() -> r_z_hat
  
  
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
  dt1$R_Z <- predict(digit.ml, newdata = dt1)
  pr_Z <- data.frame(R_Z = 1:K, p_Z = pi_m)
  merge(dt1, pr_Z, by="R_Z") -> dt1_Z
  
  # fit_f <- randomForest(R_L~X1+Y+Z , data = dt2, proximity = T)
  # dt1$R_F <- predict(fit_f, newdata = dt1)
  # pr_F <- data.frame(R_F = 1:K, p_F = pi_m)
  # merge(dt1, pr_F, by="R_F") -> dt1_F
  
  #######
  sum((dt1$phi_a)^2/dt1$pi_u)/N/N -> thvar_u
  sum((dt1_Z$phi_a)^2/dt1_Z$p_Z)/N/N -> thvar_Z
  # sum((dt1_F$phi_a)^2/dt1_F$p_F)/N/N -> thvar_F
  r_Z = thvar_u/thvar_Z; 
  # r_F = thvar_u/thvar_F
  ########
  get_tau_hat(dt1, dt1$pi_u, method) -> tau_hat_u
  get_tau_hat(dt1_Z, dt1_Z$p_Z, method)-> tau_hat_Z
  # get_tau_hat(dt1_F, dt1_F$p_F, method)-> tau_hat_F
  
  # return(c(tau_hat_u, tau_hat_Z, tau_hat_F, 
  #          r_Z, r_F, r_z_hat, r_f_hat))
  
  return(c(tau_hat_u, tau_hat_Z, 
           r_Z,  r_z_hat))
}

# testf(N = 10000, n_m = 1000, n = 1000, K = 10, method = "aipw")

