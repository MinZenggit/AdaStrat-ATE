library(stringr)
isalpha <- function(str){
  str <- ifelse(grepl("[a-zA-Z]", str), NA, str) %>% unlist %>% as.numeric()
  return(str)
}
get_tau_edr <- function(dt1_2, pi1_2, R1_2, fit2){
  Y <- dt1_2$Y
  W <- 1/pi1_2; R = R1_2; A <- dt1_2$Z;
  V <- dt1_2$V <- rbinom(length(pi1_2), 1, pi1_2)
  
  dtva <- dt1_2 %>% filter(V==1)
  y <- dtva$Y; x <- dtva[,c("X1", "U1")] %>% as.matrix
  w <- 1/pi1_2[V==1]; r <- R[V==1]; a <- dtva$Z; p = unique(pi1_2)
  
  EDR_w(W, R, V, p, w, r, a, v, x, y) -> fit1_2
  
  (fit1_2$est/fit1_2$ve + fit2$est/fit2$ve)/(1/fit1_2$ve + 1/fit2$ve) -> est_hat
  (fit1_2$ve*fit2$ve)/(fit1_2$ve+fit2$ve) -> ve_hat
  cb_2 <- list(est = est_hat, ve = ve_hat)
  return(cb_2)
}
testf <- function(N, n_m, n, alpha, K){
  Generatedt(N) -> dt1 -> dtfull
  #######
  # Outcome dependent group
  te <- dt1 %>% dplyr::select(id, Y)
  te[order(te$Y), ] -> te
  R_O = rep(1:K, each = ceiling(N/K)) %>% as.factor()
  te$R_O =  R_O[1:N]
  dt1$R_O = te[order(te$id), "R_O"]
  ######
  
  ############
  # the middle phase data
  p = n_m/N;
  midphase_ind <- sample(1:N, size = n_m, replace = F)
  #### The interndata
  dt1$V = 0; dt1$V[midphase_ind] = 1
  #### The otherdata
  dt1_2 <- dt1[-midphase_ind,]
  N1_2 = nrow(dt1_2)
  
  Y <- dt1$Y; 
  W <- rep(1/p, N); R <- rep(1, N); A <- dt1$Z; V <- dt1$V
  
  # first consider the Y(A=1)
  dt2 <- dt1 %>% filter(V==1); y <- dt2$Y ;x <- dt2[,c("X1", "U1")] %>% as.matrix
  w <- rep(1/p, n_m); r <- rep(1, n_m); a <- dt2$Z; 
  
  EDR_w(W, R, V, p, w, r, a, v, x, y) -> fit2
  dt2$phi_hat <- fit2$inf.h1t
  phi_hat_bar <- mean(dt2$phi_hat)
  dt2 <- dt2[order(dt2$phi_hat), ]
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
  ######################
  ##The random foest method for classification: F method
  fit_f <- randomForest(R_L~X1+Y+Z , data = dt2, proximity = T)
  dt1_2$R_F <- predict(fit_f, newdata = dt1_2)
  pr_F <- data.frame(R_F = 1:K, p_F = pi_m)
  merge(dt1_2, pr_F, by="R_F") -> dt1_2_F
  ########
  get_tau_edr(dt1_2, pi1_2 = rep(n/(N-n_m), N-n_m), R1_2 = rep(1, N-n_m), fit2) -> fit_u
  get_tau_edr(dt1_2_Z, pi1_2 = dt1_2_Z$p_Z, R1_2 = dt1_2_Z$R_Z, fit2) -> fit_Z
  get_tau_edr(dt1_2_F, pi1_2 = dt1_2_F$p_F, R1_2 = dt1_2_F$R_F, fit2) -> fit_F
  get_tau_edr(dt1_2_O, pi1_2 = dt1_2_O$p_O, R1_2 = dt1_2_O$R_O, fit2) -> fit_O
  return(c(fit_u$est, fit_Z$est, fit_F$est, fit_O$est,
           fit_u$ve, fit_Z$ve, fit_F$ve, fit_O$ve))
}
testf(N, n_m, n, alpha, K)
