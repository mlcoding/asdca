rm(list=ls())
library(picasso)
library(stocml)
library(glmnet)

smooth.svm = function(x, gamma){
  phi = rep(0,length(x))
  idx1 = x<=1-gamma
  phi[idx1] = 1-x[idx1]-gamma/2
  idx2 = x<1&x>1-gamma
  phi[idx2] = (1-x[idx2])^2/(2*gamma)
  phi
}

n = 1000
d = 50
prec1 = 1e-5
prec2 = 1e-3
max.ite2 = n*5
verbose = TRUE
corX.flag = 0
# set.seed(101)
cor.X = 0.1
S0 = matrix(cor.X,d,d) + (1-cor.X)*diag(d)

R0 = chol(S0)
if(corX.flag==1){
  X0 = mvrnorm(n,rep(0,d),S0)/sqrt(n-1)
}else{
  X0 = scale(matrix(rnorm(n*d),n,d)%*%R0)/sqrt(n-1)*sqrt(n)
}
w0 = matrix(rnorm(d, mean = 0, sd = 1)*2,nrow=d,ncol=1)
Y = sign(X0%*%w0)
nlambda = 10
lambda.max = max(abs(crossprod(X0,Y)))/n*1.0001
lambdas = exp(seq(log(lambda.max), log(lambda.max*0.05), length = nlambda))/10
X = diag(c(Y))%*%X0
X_rown = rowSums(X^2)
max.ite1 = 1e3
ite1.vec = matrix(0,nrow=nlambda,ncol=max.ite1)
dif.beta1 = rep(0,nlambda)
dif.beta2 = rep(0,nlambda)
act.alp = rep(0,nlambda)
runt.r = rep(0,nlambda)

yw = w1 = w2 = matrix(0,nrow=d,ncol=1)
for(i in 1:nlambda){
  alp = matrix(0,nrow=n,ncol=1)
  lambda = lambdas[i]
  gamma.svm = prec2*1
  R = max(max(sqrt(X_rown)),sqrt(11*n*lambda*1))
  kappa = R^2/(1*n)-lambda
  mu = lambda/2
  rho = mu+kappa
  eta = sqrt(mu/rho)
  beta = (1-eta)/(1+eta)
  lambdak = lambda+kappa
  gap1 = 1
  gap1.vec = rep(0,max.ite1)
  alp.act = rep(0,max.ite1)
  alp.up = rep(0,max.ite1)
  alp.actset = matrix(0,nrow=max.ite1,ncol=max.ite2)
  gap2.mat = matrix(0,nrow=max.ite1,ncol=max.ite2)
  dif.vec = rep(0,max.ite1)
  ite1 = 0
  # w2 = crossprod(X,alp)/(lambdak*n)
  epst = ((sum(smooth.svm(X%*%w2,gamma.svm))-sum(alp)+sum(alp^2)*gamma.svm/2)/n+lambda*crossprod(w2)[1])*(1+eta^(-2))
  begt = Sys.time()
  max.ite11 = ceiling(1+2*log(epst/prec1)/eta)
  #cat("lambda=",lambda,"gamma.svm=",gamma.svm,"R=",R,"kappa=",kappa,"mu=",mu,
  #    "rho=",rho,"eta=",eta,"beta=",beta,"epst=",epst,"max.ite11=",max.ite11,"\n")
  while(gap1>prec1 && ite1<max.ite11){
    ite1 = ite1+1
    ite2 = 0
    gap2 = 1
    z = kappa*yw/lambdak
    # w2 = z + crossprod(X,alp)/(lambdak*n)
    # idx = sample(c(1:n), n, replace = FALSE)
    idx = c(1:n)
    # idx = sample(c(1:n), max.ite2, replace = TRUE)
    prec2 = epst*eta/(2*(1+eta^(-2)))
    update=0
    idx_cnt = 0
    #cat(z[1],z[2],z[3],z[4],z[5],prec2)
    while(gap2>prec2 && ite2<max.ite2){
      ite2 = ite2+1
      # c_idx = idx[ite2]
      idx_cnt = ite2%%n
      if(idx_cnt==0) idx_cnt=n
      c_idx = idx[idx_cnt]
      delta_alp = max(-alp[c_idx],min(1-alp[c_idx],(1-crossprod(X[c_idx,],w2)[1]-gamma.svm*alp[c_idx])
                                      /(crossprod(X[c_idx,])[1]/(lambda*n)+gamma.svm)))
      if(delta_alp!=0){
        update = update+1
        alp.actset[ite1,ite2] = c_idx
      }
      alp[c_idx] = alp[c_idx]+delta_alp
      w2 = w2 + X[c_idx,]*delta_alp/(lambda*n)
      #cat("beta2",w2[1],w2[2],w2[3],w2[4],"\n")
      gap2 = (sum(smooth.svm(X%*%w2,gamma.svm))-sum(alp)+sum(alp^2)*gamma.svm/2)/n+lambda*crossprod(w2,w2-z)[1]
      gap2.mat[ite1,ite2] = gap2
      Xb=X%*%w2
      smo=smooth.svm(X%*%w2,gamma.svm)
      #Xb[1:5]
      #smo[1:5]
      #cat(c_idx,delta_alp,gap2)
    }
    alp.up[ite1] = update
    alp.act[ite1] = sum(alp!=0)
    ite1.vec[i,ite1] = ite2
    epst = (1-eta/2)^(ite1-1)*epst
    gap1 = (1+rho/mu)*epst+(rho*kappa*sum((w2-yw)^2))/(2*mu)
    yw = w2+beta*(w2-w1)
    # gap1 = norm(w2-w1)
    w1 = w2
    gap1.vec[ite1] = gap1
    dif.vec[ite1] = n-sum((X%*%w2)>1)
    # cat("i=",i,"ite1=",ite1,"ite2=",ite2,", gap1=",gap1,", gap2=",gap2,", prec2=",prec2,"\n")
  }
  runt.r[i] = Sys.time()-begt
  par(mfrow=c(2,2))
  plot(gap1.vec[1:ite1],type="l",log="y")
  plot(dif.vec[1:ite1],type="l",log="y")
  plot(alp.act[1:ite1],type="l",log="y")
  plot(alp.up[1:ite1],type="l",log="y")
  # plot(alp.actset[1:ite1],type="l",log="y")
  alp.up[1:ite1]
  sum(ite1.vec)
  dif.beta1[i] = dif.vec[ite1]
  dif.beta2[i] = sum(abs(Y - sign(X0%*%w2)))/2
  act.alp[i] = alp.act[ite1]
}
rowSums(ite1.vec)
dif.beta1
dif.beta2
act.alp
