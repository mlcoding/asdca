#----------------------------------------------------------------------------------#
# Package: stocml                                                                  #
# stocml(): The user interface for stocml()                                        #
# Author: Xingguo Li                                                               #
# Email: <xingguo.leo@gmail.com>                                                   #
# Date: Jun 29th, 2015                                                             #
# Version: 0.1.0                                                                   #
#----------------------------------------------------------------------------------#

stocml <- function(X, 
                   Y, 
                   lambda = NULL,
                   nlambda = NULL,
                   lambda.min.ratio = NULL,
                   family = "gaussian",
                   alg = "sdca",
                   standardize = FALSE,
                   design.sd = TRUE,
                   res.sd = FALSE,
                   prec = 1e-4,
                   max.ite = 1e4,
                   verbose = TRUE)
{
  if(family=="gaussian"){
    out = stocml.ridge(X = X, Y = Y, lambda = lambda, nlambda = nlambda, lambda.min.ratio = lambda.min.ratio,
                       alg = alg, design.sd = design.sd, res.sd = res.sd, prec = prec, max.ite = max.ite, verbose = verbose)
  }
  out$family = family
  return(out)
}
