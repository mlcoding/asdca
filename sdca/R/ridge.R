#----------------------------------------------------------------------------------#
# Package: stocml                                                                  #
# ridge(): ridge                                                                   #
# Author: Xingguo Li                                                               #
# Email: <xingguo.leo@gmail.com>                                                   #
# Date: Jun 29th, 2015                                                             #
# Version: 0.1.0                                                                   #
#----------------------------------------------------------------------------------#

ridge.sdca <- function(Y, X, lambda, nlambda, n, d, max.ite, prec,verbose)
{
  if(verbose==TRUE){
    cat("Stochastic dual coordinate ascent\n")
  }
  beta = matrix(0,nrow=d,ncol=nlambda)
  beta.intcpt = rep(0,nlambda)
  ite.lamb = rep(0,nlambda)
  ite.in = rep(0,nlambda)
  runt = matrix(0,1,nlambda)
  X_maxrn = max(sqrt(rowSums(X^2)))
  idx = sample(c(1:n), n*10, replace = TRUE)
  idx = idx-1
#   str=.C("stocml_ridge_sdca_idx", as.double(Y), as.double(X), as.double(X_maxrn), 
#          as.double(beta), as.double(beta.intcpt), as.integer(n), as.integer(d), 
#          as.integer(ite.lamb), as.integer(ite.in), 
#          as.double(runt), as.double(lambda), as.integer(nlambda), 
#          as.integer(max.ite), as.double(prec), PACKAGE="stocml")
  str=.C("stocml_ridge_sdca", as.double(Y), as.double(X), as.double(X_maxrn), 
         as.double(beta), as.double(beta.intcpt), as.integer(n), as.integer(d), 
         as.integer(ite.lamb), as.integer(ite.in), 
         as.double(runt), as.double(lambda), as.integer(nlambda), 
         as.integer(max.ite), as.double(prec), as.integer(idx), PACKAGE="stocml")
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
  return(list(beta=beta.list, intcpt = beta.intcpt, ite=ite, runt = runt))
}
