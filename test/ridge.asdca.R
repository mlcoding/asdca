rm(list=ls())
library(picasso)
library(glmnet)

n = 100
d = 5
prec1 = 1e-5
max.ite1 = 1e3
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
lambda = 0.1
R = max(max(sqrt(rowSums(X^2))),sqrt(11*n*lambda))
kappa = R^2/n-lambda
mu = lambda/2
rho = mu+kappa
eta = sqrt(mu/rho)
beta = (1-eta)/(1+eta)
lambdak = lambda+kappa
alp = matrix(0,nrow=n,ncol=1)
yw = w1 = w2 = matrix(0,nrow=d,ncol=1)
gap1 = 1
gap1.vec = rep(0,max.ite1)
gap2.mat = matrix(0,nrow=max.ite1,ncol=max.ite2)
ite1 = 0
while(gap1>prec1 && ite1<max.ite1){
  ite1 = ite1+1
  ite2 = 0
  gap2 = 1
  z = kappa*yw/lambdak
  y.tild = Y-X%*%z
  v = crossprod(X,alp)/(lambdak*n)
  idx = sample(c(1:n), max.ite2, replace = TRUE)
  while(gap2>prec2 && ite2<max.ite2){
    ite2 = ite2+1
    c_idx = idx[ite2]
    delta_alp = -(alp[c_idx]+crossprod(X[c_idx,],v)[1]-y.tild[c_idx])/
      (1+1*crossprod(X[c_idx,])[1]/(lambdak*n))
    alp[c_idx] = alp[c_idx]+delta_alp
    v = v + X[c_idx,]*delta_alp/(lambdak*n)
    w2 = v+z
    gap2 = (sum((X%*%w2-Y)^2)+sum((alp+Y)^2)-sum(Y^2))/(2*n)+lambdak*crossprod(w2,v)[1]
    gap2.mat[ite1,ite2] = gap2
    sprintf("ite=%d,gap=%.4f",ite2,gap2)
  }
  ite1.vec[ite1] = ite2
  yw = w2+beta*(w2-w1)
  gap1 = norm(w2-w1)
  w1 = w2
  gap1.vec[ite1] = gap1
  sprintf("ite=%d,gap=%.4f",ite1,gap1)
}
norm(w1-w0)
plot(gap1.vec,type="l",log="y")
gap1.vec[max.ite1]
