#----------------------------------------------------------------------------------#
# Package: stocml                                                                  #
# stocml.svm(): The user interface for svm()                                       #
# Author: Xingguo Li                                                               #
# Email: <xingguo.leo@gmail.com>                                                   #
# Date: Aug 11th, 2015                                                             #
# Version: 0.1.0                                                                   #
#----------------------------------------------------------------------------------#

stocml.svm <- function(xx, 
                       yy, 
                       lambda = NULL,
                       nlambda = NULL,
                       lambda.min.ratio = NULL,
                       alg = "sdca",
                       prec = 1e-4,
                       max.ite = 1e4,
                       verbose = TRUE)
{
  n = nrow(xx)
  d = ncol(xx)
  if(verbose)
    cat("Linear classification \n")
  if(n==0 || d==0) {
    cat("No data input.\n")
    return(NULL)
  }
  if(alg!="sdca"){
    cat(" Wrong \"alg\" input. \n \"alg\" should be \"sdca\".\n", 
        alg,"does not exist. \n")
    return(NULL)
  }
  maxdf = max(n,d)
  
  if(!is.null(lambda)) nlambda = length(lambda)
  if(is.null(lambda)){
    if(is.null(nlambda))
      nlambda = 5
    if(is.null(lambda.min.ratio)){
      lambda.min.ratio = 0.25
    }
    lambda.max = max(abs(crossprod(xx,yy/n)))
    lambda.min = lambda.min.ratio*lambda.max
    lambda = exp(seq(log(lambda.max), log(lambda.min), length = nlambda))
    rm(lambda.max,lambda.min,lambda.min.ratio)
    gc()
  }
  begt=Sys.time()
  if (alg=="sdca")
    out = svm.sdca(yy, xx, lambda, nlambda, n, d, max.ite, prec, verbose)
  runt=Sys.time()-begt
  
  est = list()
  intcpt=matrix(0,nrow=1,ncol=nlambda)
  beta1=matrix(0,nrow=d,ncol=nlambda)
  for(k in 1:nlambda){
    beta1[,k] = out$beta[[k]]
  }
  
  est$alp = out$alp
  est$runt = out$runt
  est$beta = beta1
  #est$intercept = intcpt
  est$Y = yy
  est$X = xx
  est$lambda = lambda
  est$nlambda = nlambda
  #est$df = df
  est$alg = alg
  est$ite =out$ite
  est$verbose = verbose
  est$runtime = runt
  class(est) = "svm"
  return(est)
}

print.svm <- function(x, ...)
{  
  cat("\n SVM options summary: \n")
  cat(x$nlambda, " lambdas used:\n")
  print(signif(x$lambda,digits=3))
  cat("Method =", x$method, "\n")
  cat("Alg =", x$alg, "\n")
  cat("Degree of freedom:",min(x$df),"----->",max(x$df),"\n")
  if(units.difftime(x$runtime)=="secs") unit="secs"
  if(units.difftime(x$runtime)=="mins") unit="mins"
  if(units.difftime(x$runtime)=="hours") unit="hours"
  cat("Runtime:",x$runtime," ",unit,"\n")
}

plot.svm <- function(x, ...)
{
  matplot(x$lambda, t(x$beta), type="l", main="Regularization Path",
          xlab="Regularization Parameter", ylab="Coefficient")
}

coef.svm <- function(object, lambda.idx = c(1:3), beta.idx = c(1:3), ...)
{
  lambda.n = length(lambda.idx)
  beta.n = length(beta.idx)
  cat("\n Values of estimated coefficients: \n")
  cat(" index     ")
  for(i in 1:lambda.n){
    cat("",formatC(lambda.idx[i],digits=5,width=10),"")
  }
  cat("\n")
  cat(" lambda    ")
  for(i in 1:lambda.n){
    cat("",formatC(object$lambda[lambda.idx[i]],digits=4,width=10),"")
  }
  cat("\n")
  for(i in 1:beta.n){
    cat(" beta",formatC(beta.idx[i],digits=5,width=-5))
    for(j in 1:lambda.n){
      cat("",formatC(object$beta[beta.idx[i],lambda.idx[j]],digits=4,width=10),"")
    }
    cat("\n")
  }
}

predict.svm <- function(object, newdata, lambda.idx = c(1:3), Y.pred.idx = c(1:5), ...)
{
  pred.n = nrow(newdata)
  lambda.n = length(lambda.idx)
  Y.pred.n = length(Y.pred.idx)
  Y.pred = sign(newdata%*%object$beta[,lambda.idx])
  cat("\n Values of predicted responses: \n")
  cat("   index   ")
  for(i in 1:lambda.n){
    cat("",formatC(lambda.idx[i],digits=5,width=10),"")
  }
  cat("\n")
  cat("   lambda  ")
  for(i in 1:lambda.n){
    cat("",formatC(object$lambda[lambda.idx[i]],digits=4,width=10),"")
  }
  cat("\n")
  for(i in 1:Y.pred.n){
    cat("    Y",formatC(Y.pred.idx[i],digits=5,width=-5))
    for(j in 1:lambda.n){
      cat("",formatC(Y.pred[Y.pred.idx[i],j],digits=4,width=10),"")
    }
    cat("\n")
  }
  return(Y.pred)
}