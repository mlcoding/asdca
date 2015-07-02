rm(list=ls())
library(picasso)
library(glmnet)

prec = 1e-5
max.ite = 1e4
n = 100
d = 5
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
Y = X%*%w0 + rnorm(n)/1e0
nlambda = 10
lambda.max = max(abs(crossprod(X,Y)))/n*1.0001
lambda = exp(seq(log(lambda.max), log(lambda.max*0.05), length = nlambda))
lambda = 1e-2
R = max(sqrt(rowSums(X^2)))
tau = sqrt(1/(n*lambda))/(2*R)
sigma = sqrt(n*lambda)/(2*R)
theta = 1-1/(n+R*sqrt(n/lambda))
alp1 = alp2 = matrix(rnorm(n),nrow=n,ncol=1)
w1 = w2 = w.bar = matrix(0,nrow=d,ncol=1)
u1 = u2 = crossprod(X,alp1)/n
idx = sample(c(1:n), max.ite, replace = TRUE)
gap = 1
gap.vec = rep(0,max.ite)
ite = 0
while(gap>prec && ite<max.ite){
  ite = ite+1
  c_idx = idx[ite]
  alp2[c_idx] = (alp1[c_idx]/sigma+crossprod(X[c_idx,],w.bar)[1]-Y[c_idx])/
    (1+1/sigma)
  w2 = (w1/tau-(u1+(alp2[c_idx]-alp1[c_idx])*X[c_idx,]))/(2*lambda+1/tau)
  u2 = u1+(alp2[c_idx]-alp1[c_idx])*X[c_idx,]/n
  w.bar = w2 + theta*(w2-w1)
  gap.w = norm(w2-w1)
  gap.u = norm(u2-u1)
  gap.alp = norm(alp2-alp1)
  w1 = w2
  u1 = u2 
  alp1 = alp2
  gap = max(gap.w,gap.u,gap.alp)
  gap.vec[ite] = gap
  sprintf("ite=%d,gap.w=%.4f,gap.u=%.4f,gap.alp=%.4f",
          ite,gap.w,gap.u,gap.alp)
}
norm(w1-w0)
plot(gap.vec,type="l",log="y")
gap.vec[max.ite]
w.hat = solve(crossprod(X),crossprod(X,Y))
norm(w0-w.hat)
