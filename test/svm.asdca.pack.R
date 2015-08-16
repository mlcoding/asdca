rm(list=ls())
library(stocml)
library(glmnet)

n = 1000
d = 50
prec1 = 1e-5
prec2 = 1e-3
max.ite1 = 1e3
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
Y = sign(X%*%w0)
nlambda = 10
lambda.max = max(abs(crossprod(X,Y)))/n*1.0001
lambdas = exp(seq(log(lambda.max), log(lambda.max*0.05), length = nlambda))/10
out = stocml(X0,Y,lambda=lambdas,family="svm",alg = "sdca",prec = prec1,max.ite = max.ite1)
w.hat = out$beta
out$ite
out$runt
dif.vec1 = rep(0,nlambda)
dif.vec2 = rep(0,nlambda)
for(i in 1:nlambda){
  dif.vec1[i] = sum(abs(Y - sign(X0%*%w.hat[,i])))/2
  dif.vec2[i] = n-sum((diag(c(Y))%*%X0%*%w.hat[,i])>1)
}
# w0
plot(out)
dif.vec1
dif.vec2
