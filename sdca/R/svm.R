#----------------------------------------------------------------------------------#
# Package: stocml                                                                  #
# ridge(): SVM                                                                     #
# Author: Xingguo Li                                                               #
# Email: <xingguo.leo@gmail.com>                                                   #
# Date: Aug 12th, 2015                                                             #
# Version: 0.1.0                                                                   #
#----------------------------------------------------------------------------------#

svm.sdca <- function(Y, X, lambda, nlambda, n, d, max.ite, prec,verbose)
{
  if(verbose==TRUE){
    cat("Stochastic dual coordinate ascent\n")
  }
  beta = matrix(0,nrow=d,ncol=nlambda)
  beta.intcpt = rep(0,nlambda)
  alp = matrix(0,nrow=n,ncol=nlambda)
  ite.lamb = rep(0,nlambda)
  ite.in = rep(0,nlambda)
  runt = matrix(0,1,nlambda)
  #cat("X0",X[1],X[2],X[3],X[4],X[5])
  X = diag(c(Y))%*%X
  #cat("X1",X[1],X[2],X[3],X[4],X[5])
  X_rn = rowSums(X^2)
  X_maxrn = max(sqrt(X_rn))
  #cat(X_rn[1],X_rn[2],X_rn[3],X_rn[4],X_rn[5])
  str=.C("stocml_svm_sdca_idx", as.double(Y), as.double(X), as.double(X_maxrn), 
         as.double(beta), as.double(beta.intcpt), as.integer(n), as.integer(d), 
         as.integer(ite.lamb), as.integer(ite.in), 
         as.double(runt), as.double(lambda), as.integer(nlambda), 
         as.integer(max.ite), as.double(prec), as.double(X_rn),as.double(alp),
         PACKAGE="stocml")
  # idx = sample(c(1:n), n*10, replace = TRUE)
  # idx = idx-1
#   str=.C("stocml_svm_sdca", as.double(Y), as.double(X), as.double(X_maxrn), 
#          as.double(beta), as.double(beta.intcpt), as.integer(n), as.integer(d), 
#          as.integer(ite.lamb), as.integer(ite.in), 
#          as.double(runt), as.double(lambda), as.integer(nlambda), 
#          as.integer(max.ite), as.double(prec), as.integer(idx), PACKAGE="stocml")
  beta.list = vector("list", nlambda)
  for(i in 1:nlambda){
    beta.i = unlist(str[4])[((i-1)*d+1):(i*d)]
    beta.list[[i]] = beta.i
  }
  beta.intcpt = unlist(str[5])
  ite.lamb = unlist(str[8])
  ite.in = unlist(str[9])
  ite = list()
  ite[[1]] = ite.lamb
  ite[[2]] = ite.in
  runt = matrix(unlist(str[10]),ncol=nlambda,byrow = FALSE)
  return(list(beta=beta.list, intcpt = beta.intcpt, ite=ite, runt = runt, alp=alp))
}
