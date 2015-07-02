#----------------------------------------------------------------------------------#
# Package: stocml                                                                  #
# stocml.ridge(): The user interface for ridge()                                   #
# Author: Xingguo Li                                                               #
# Email: <xingguo.leo@gmail.com>                                                   #
# Date: Jun 29th, 2015                                                             #
# Version: 0.1.0                                                                   #
#----------------------------------------------------------------------------------#

stocml.ridge <- function(X, 
                         Y, 
                         lambda = NULL,
                         nlambda = NULL,
                         lambda.min.ratio = NULL,
                         alg = "sdca",
                         design.sd = TRUE,
                         res.sd = FALSE,
                         prec = 1e-4,
                         max.ite = 1e4,
                         verbose = TRUE)
{
  n = nrow(X)
  d = ncol(X)
  if(verbose)
    cat("Linear regression. \n")
  if(n==0 || d==0) {
    cat("No data input.\n")
    return(NULL)
  }
  if(alg!="sdca" && alg!="spdc"){
    cat(" Wrong \"alg\" input. \n \"alg\" should be one of \"sdca\" and \"spdc\".\n", 
        alg,"does not exist. \n")
    return(NULL)
  }
  maxdf = max(n,d)
  if(design.sd){
    xm=matrix(rep(colMeans(X),n),nrow=n,ncol=d,byrow=T)
    x1=X-xm
    xinvc.vec=1/sqrt(colSums(x1^2)/(n-1))
    xx=x1%*%diag(xinvc.vec)
    ym=mean(Y)
    y1=Y-ym
    if(res.sd == TRUE){
      sdy=sqrt(sum(y1^2)/(n-1))
      yy=y1/sdy
    }else{
      sdy = 1
      yy = y1
    }
  }else{
    xinvc.vec = rep(1,d)
    sdy = 1
    xx = X
    yy = Y
  }
  
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
    #out = ridge.sdca(yy, xx, lambda, nlambda, n, d, max.ite, prec, verbose)
    out = ridge.sdca.call(yy, xx, lambda, nlambda, n, d, max.ite, prec, verbose)
  runt=Sys.time()-begt
  
  est = list()
  intcpt=matrix(0,nrow=1,ncol=nlambda)
  beta1=matrix(0,nrow=d,ncol=nlambda)
  if(design.sd){
    for(k in 1:nlambda){
      tmp.beta = out$beta[[k]]
      beta1[,k]=xinvc.vec*tmp.beta*sdy
      intcpt[k] = ym-as.numeric(xm[1,]%*%beta1[,k])+out$intcpt[k]*sdy
    }
  }else{
    for(k in 1:nlambda){
      beta1[,k]=out$beta[[k]]
      intcpt[k] = out$intcpt[k]
    }
  }
  
  #est$obj = out$obj
  est$runt = out$runt
  est$beta = beta1
  est$intercept = intcpt
  est$Y = Y
  est$X = X
  est$lambda = lambda
  est$nlambda = nlambda
  #est$df = df
  est$alg = alg
  est$ite =out$ite
  est$verbose = verbose
  est$runtime = runt
  class(est) = "ridge"
  return(est)
}

print.ridge <- function(x, ...)
{  
  cat("\n ridge options summary: \n")
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

plot.ridge <- function(x, ...)
{
  matplot(x$lambda, t(x$beta), type="l", main="Regularization Path",
          xlab="Regularization Parameter", ylab="Coefficient")
}

coef.ridge <- function(object, lambda.idx = c(1:3), beta.idx = c(1:3), ...)
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
  cat(" intercept ")
  for(i in 1:lambda.n){
    cat("",formatC(object$intercept[i],digits=4,width=10),"")
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

predict.ridge <- function(object, newdata, lambda.idx = c(1:3), Y.pred.idx = c(1:5), ...)
{
  pred.n = nrow(newdata)
  lambda.n = length(lambda.idx)
  Y.pred.n = length(Y.pred.idx)
  intcpt = matrix(rep(object$intercept[,lambda.idx],pred.n),nrow=pred.n,
                  ncol=lambda.n,byrow=T)
  Y.pred = newdata%*%object$beta[,lambda.idx] + intcpt
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