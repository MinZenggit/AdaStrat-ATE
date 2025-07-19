EDR_w <- function(W, R, V, p, w, r, a, v, x, y) {
  # mu_1
  N <- length(W)
  lm.out1 <- lm(y[which(a == 1)] ~ x[which(a == 1), ],
                weights = w[which(a == 1)])
  eta_hat <- lm.out1$coefficients
  mu1_hat <- cbind(1, x) %*% lm.out1$coefficients
  
  # pi
  glm.out <- glm(a ~ x, family = "binomial", weights = w)
  alpha_hat <- glm.out$coefficients
  pshat <- glm.out$fitted.values
  
  Exp_psi <- function(k, mu1x) {
    aa <- w * (a / pshat * (y - mu1x) - (a - pshat) / pshat * (mu1_hat - mu1x))
    mean(aa[r == k]) * p[k]
  }
  
  Un <- function(mu1) {
    unique_R <- unique(R)
    K <- length(unique_R)
    # 初始化 Exp_psi_values 向量
    Exp_psi_values <- numeric(K)
    for (i in 1:K) {
      k <- unique_R[i]
      idx <- which(r == k)
      aa <- w[idx] * (a[idx] / pshat[idx] * (y[idx] - mu1) - 
                        (a[idx] - pshat[idx]) / pshat[idx] * (mu1_hat[idx] - mu1))
      Exp_psi_values[i] <- mean(aa) * p[k]
    }
    Exp_psi_vec <- Exp_psi_values[match(R, unique_R)]
    # 计算 Un
    sum(a / pshat * (y - mu1) * w) - 
      sum(((a - pshat) / pshat * (mu1_hat - mu1)) * w) - 
      sum((V * W - 1) * Exp_psi_vec)
  }
  uniroot(Un, interval = c(-3, 3))$root -> mu1_est -> EST
  
  mu1_hat <- as.vector(mu1_hat)
  H <- c(mean(w * (mu1_hat - y) * a * (1 - pshat) / pshat),
         apply(w * (mu1_hat - y) * a * (1 - pshat) / pshat * x, 2, mean)) / mean(w)
  
  # Sigma_aa = E(S_a')
  dimx <- ncol(x)
  n <- nrow(x)
  E <- matrix(0, 1 + dimx, 1 + dimx)
  forE <- matrix(1, 1 + dimx, n)
  forE[2:(1 + dimx), ] <- t(x)
  for (jj in 1:(1 + dimx)) {
    for (kk in jj:(1 + dimx)) {
      E[jj, kk] <- E[kk, jj] <- mean(pshat * (1 - pshat) * forE[jj, ] * forE[kk, ] * w) / mean(w)
    }
  }
  (a - pshat) * H %*% ginv(E) %*% forE -> aaa
  
  ## H_t EDR
  B <- apply((1 - a / pshat) * cbind(1, x), 2, weighted.mean, w = w)
  B1mat <- cbind(a, x * a)
  B1 <- matrix(0, 1 + dimx, 1 + dimx)
  for (jj in 1:(1 + dimx)) {
    for (kk in jj:(1 + dimx)) {
      B1[jj, kk] <- B1[kk, jj] <- weighted.mean(B1mat[, jj] * B1mat[, kk], w)
    }
  }
  b1 <- B %*% ginv(B1)
  reg.resid1 <- reg.resid0 <- rep(0, n)
  reg.resid1[a == 1] <- (lm.out1$resid)
  apply(B1mat * matrix(b1, n, 1 + dimx, byrow = TRUE), 1, sum) * a * reg.resid1 -> bbb
  
  # psi_t
  (a / pshat * (y - mu1_est) - (a - pshat) / pshat * (mu1_hat - mu1_est)) -> ccc
  # Epsi_t
  dd <- w * ccc
  sapply(1:length(p), function(x) { mean(dd[r == x]) * p[x] }) -> e_psi
  sapply(r, function(x) { e_psi[x] }) -> ddd
  h1t <- as.numeric(ccc - ddd - aaa - bbb)
  H1t <- rep(0, N)
  H1t[V == 1] <- h1t
  H2t <- sapply(R, function(x) { e_psi[x] })
  
  INF.total <- V * W * H1t + H2t
  
  mean(INF.total^2) / N -> VE
  return(list(est = EST,
              ve = VE,
              inf.total = INF.total,
              inf.h1t = h1t,
              inf.H2t = H2t))
}
# EDR_w<-function(W, R, V, p, w, r, a, v, x, y){
#   # mu_1
#   N = length(W)
#   lm.out1<-lm(y[which(a==1)]~x[which(a==1),],
#               weights = w[which(a==1)])
#   eta_hat <- lm.out1$coefficients
#   mu1_hat<-cbind(1, x)%*%lm.out1$coefficients
# 
#   # pi
#   glm.out<-glm(a ~ x,family="binomial", weights = w)
#   alpha_hat <- glm.out$coefficients
#   pshat<-glm.out$fitted.values
# 
#   Exp_psi <- function(k, mu1x){
#     aa <- w * (a/pshat*(y-mu1x) - (a-pshat)/ pshat*(mu1_hat - mu1x))
#     mean(aa[r == k])*p[k]
#   }
# 
#   Un <- function(mu1){
#     sum(a/pshat*(y-mu1)*w) - sum(((a-pshat)/ pshat*(mu1_hat - mu1))*w)-
#       sum((V*W - 1) * unlist(lapply(R, Exp_psi, mu1x = mu1)))
#   }
#   uniroot(Un, interval =  c(-3, 3))$root -> mu1_est -> EST
# 
#   mu1_hat <- as.vector(mu1_hat)
#   H<-c(mean(w*(mu1_hat - y)*a*(1-pshat)/pshat),
#        apply(w*(mu1_hat - y)*a*(1-pshat)/pshat*x, 2 , mean))/mean(w)
# 
#   # Sigma_aa =  E(S_a')
#   dimx <- ncol(x); n = nrow(x)
#   E<-matrix(0,1+dimx,1+dimx)
#   forE<-matrix(1,1+dimx,n)
#   forE[2:(1+dimx),]<-t(x)
#   for(jj in 1:(1+dimx)){
#     for(kk in jj:(1+dimx)){
#       E[jj,kk]<-E[kk,jj]<-mean(pshat*(1-pshat)*forE[jj,]*forE[kk,]*w)/mean(w)
#     }
#   }
#   (a-pshat)*H%*%ginv(E)%*%forE -> aaa
#   ## H_t EDR
#   B<-apply((1 - a/pshat) * cbind(1,x),2, weighted.mean, w = w)
#   B1mat<-cbind(a, x*a)
#   B1<-matrix(0,1+dimx,1+dimx)
#   for(jj in 1:(1+dimx)){
#     for(kk in jj:(1+dimx)){
#       B1[jj,kk]<-B1[kk,jj]<-weighted.mean(B1mat[,jj]*B1mat[,kk], w)
#     }
#   }
#   b1<-B%*%ginv(B1)
#   reg.resid1<-reg.resid0<-rep(0,n)
#   reg.resid1[a==1]<-(lm.out1$resid)
#   apply(B1mat*matrix(b1,n,1+dimx,byrow=TRUE),1,sum)*a*reg.resid1 -> bbb
# 
#   # psi_t
#   (a/pshat*(y-mu1_est) - (a-pshat)/ pshat*(mu1_hat - mu1_est)) -> ccc
#   # Epsi_t
#   dd <- w * ccc
#   sapply(1:length(p), function(x){mean(dd[r == x])*p[x]}) -> e_psi
#   sapply(r, function(x){e_psi[x]}) -> ddd
#   h1t <- as.numeric(ccc - ddd -aaa - bbb)
#   H1t <- rep(0, N)
#   H1t[V==1] = h1t
#   H2t <- sapply(R, function(x){e_psi[x]})
# 
#   INF.total = V*W * H1t + H2t
# 
#   mean(INF.total^2)/N -> VE
#   return(list(est = EST,
#            ve = VE,
#            inf.total = INF.total,
#            inf.h1t = h1t,
#            inf.H2t = H2t))
# }

