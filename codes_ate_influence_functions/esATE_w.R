# estimation ate, only the point estimation

# Estimate Average Treatment Effect with weight (esATE_w)

# x the matrix of confounders fully observed  (n x dim(x))
# y the vector of outcome (n x 1)
# A the vector of binary treatment (n x 1)
# w is the vector of weitht (n x 1)
# method the estimation method for ATE
# select 1 method from c("reg","ipw","aipw")
#
# "reg": a linear regression imputation estimator of the ATE
#
# "ipw": the inverse probability of treatment weighting estimator of the ATE, where the propensity score follows a logistic regression model
#
# "aipw": the augmented inverse probability of treatment weighting estimator (Lunceford Davidian, 2004) of the ATE
#
#
# est: estimate of the ATE

esATE_w<-function(x, y, A, w, method){
  ## check argument
  method<-method[which(method%in%c("reg","ipw","aipw")==1)]
  if(length(method)==0){return(cat("Error: specify a valid argument for method"))}
  
  n<-length(y)#sample size 
  x<-as.matrix(x)
  if(method=="reg"){
    lm.y<- y
    lm.x<- x
    lm.w<- w
    dimx<-dim(lm.x)[2]
    lm.out1<-lm(lm.y[which(A==1)]~lm.x[which(A==1),],
                weights = lm.w[which(A==1)])
    reg.mu1<-cbind(1,lm.x)%*%lm.out1$coefficients
    lm.out0<-lm(lm.y[which(A==0)]~lm.x[which(A==0),],
                weights = lm.w[which(A==0)])
    reg.mu0<-cbind(1,lm.x)%*%lm.out0$coefficients
    reg<-mean(w*reg.mu1-w*reg.mu0)/mean(w)
    EST <- reg
  }
  if(method=="ipw"){
    Rx<-x
    Rx<-as.matrix(Rx)
    dimx<-dim(Rx)[2]
    glm.out<-glm(A~Rx,family="binomial", weights = w)
    pshat<-glm.out$fitted.values
    psw<-mean(w*y*A/pshat-w*y*(1-A)/(1-pshat))/mean(w)
    EST<-psw
  }
  if(method=="aipw"){
    lm.y<- y
    lm.x<- x
    lm.w<- w
    lm.out1<-lm(lm.y[which(A==1)]~lm.x[which(A==1),],
                weights = lm.w[which(A==1)])
    mu1<-cbind(1,lm.x)%*%lm.out1$coefficients
    lm.out0<-lm(lm.y[which(A==0)]~lm.x[which(A==0),],
                weights = lm.w[which(A==0)])
    mu0<-cbind(1,lm.x)%*%lm.out0$coefficients
    Rx<-x
    glm.out<-glm(A~Rx,family="binomial", weights = w)
    pshat<-glm.out$fitted.values
    psi.aug<-y*A/pshat-y*(1-A)/(1-pshat)+( mu1*(1-A/pshat)-mu0*(1-(1-A)/(1-pshat)))
    aipw<-mean(w*psi.aug)/mean(w)
    EST<-aipw
  }
  return(EST)
}

