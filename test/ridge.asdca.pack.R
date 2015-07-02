rm(list=ls())
library(stocml)
library(glmnet)

n = 1000
d = 20
prec1 = 1e-5
max.ite1 = 1e4
ite1.vec = rep(0,max.ite1)
prec2 = 1e-3
max.ite2 = n*5
verbose = TRUE
corX.flag = 0
# set.seed(101)
cor.X = 0.1
S0 = matrix(cor.X,d,d) + (1-cor.X)*diag(d)

R0 = chol(S0)
if(corX.flag==1){
  X = mvrnorm(n,rep(0,d),S0)/sqrt(n-1)
}else{
  X = scale(matrix(rnorm(n*d),n,d)%*%R0)/sqrt(n-1)*sqrt(n)
}
w0 = matrix(rnorm(d, mean = 0, sd = 1)*2,nrow=d,ncol=1)
Y = X%*%w0 + 1*rnorm(n)/1e1
nlambda = 10
lambda.max = max(abs(crossprod(X,Y)))/n*1.0001
lambda = exp(seq(log(lambda.max), log(lambda.max*0.05), length = nlambda))
out = stocml(X,Y,lambda,alg = "sdca",design.sd=TRUE,prec = prec1,max.ite = max.ite1)
w.hat = out$beta
ite = out$ite
intcpt.hat = out$intercept
out$runt
