rm(list=ls())
library(picasso)
library(glmnet)

prec = 1e-4
max.ite = 1e4
n = 100
d = 5
verbose = TRUE
corX.flag = 0
# set.seed(101)
cor.X = 0.0
S0 = matrix(cor.X,d,d) + (1-cor.X)*diag(d)

R = chol(S0)
if(corX.flag==1){
  X = mvrnorm(n,rep(0,d),S0)/sqrt(n-1)
}else{
  X = scale(matrix(rnorm(n*d),n,d)%*%R)/sqrt(n-1)*sqrt(n)
}
w0 = matrix(rnorm(d, mean = 0, sd = 1)*2,nrow=d,ncol=1)
Y = X%*%w0 + 1*rnorm(n)/1e10
nlambda = 10
lambda.max = max(abs(crossprod(X,Y)))/n*1.0001
lambda = exp(seq(log(lambda.max), log(lambda.max*0.05), length = nlambda))
lambda = 1e-2
alp = matrix(rnorm(n)*0,nrow=n,ncol=1)
v = crossprod(X,alp)/(lambda*n)
idx = sample(c(1:n), max.ite, replace = TRUE)
gap = 1
gap.vec = rep(0,max.ite)
dif.vec = rep(0,max.ite)
ite = 0
X_rown = rowSums(X^2)
while(gap>prec && ite<max.ite){
  ite = ite+1
  c_idx = idx[ite]
  delta_alp = -(alp[c_idx]+crossprod(X[c_idx,],v)[1]-Y[c_idx])/
    (1+1*crossprod(X[c_idx,])[1]/(lambda*n))
  alp[c_idx] = alp[c_idx]+delta_alp
  v = v + X[c_idx,]*delta_alp/(lambda*n)
  w = v
  # gap = (sum((X%*%w-Y)^2)+sum((alp+Y)^2)-sum(Y^2))/(2*n)+lambda*crossprod(w,v)[1]
  gap = (sum((X%*%w-Y)^2)+sum((alp+Y)^2)-sum(Y^2)+crossprod(X_rown,alp^2)[1]/(lambda*n))/(2*n)+lambda*crossprod(w)[1]/2
  gap.vec[ite] = gap
  dif.vec[ite] = norm(w-w0)
  # cat("ite=",ite,", gap=",gap,"\n")
}
par(mfrow=c(2,1))
plot(gap.vec,type="l",log="y")
plot(dif.vec,type="l",log="y")
gap.vec[max.ite]
dif.vec[max.ite]
ridge.est = solve(crossprod(X)+lambda*diag(d),crossprod(X,Y))
norm(ridge.est-w0)
