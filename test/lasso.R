rm(list=ls())
library(picasso)
library(glmnet)

prec = 1e-4
max.ite = 1e4
n = 100
d = 50
verbose = TRUE
corX.flag = 0
# set.seed(101)
cor.X = 0.5
S0 = matrix(cor.X,d,d) + (1-cor.X)*diag(d)

R = chol(S0)
if(corX.flag==1){
  X = mvrnorm(n,rep(0,d),S0)/sqrt(n-1)
}else{
  X = scale(matrix(rnorm(n*d),n,d)%*%R)/sqrt(n-1)*sqrt(n)
}
w = c(2,rep(0,3),0,3,-2,0,0,0,1,0,rep(0,d-13),-1.5)
Y = X%*%w + rnorm(n) + 2
nlambda = 10
lambda.max = max(abs(crossprod(X,Y)))/n*1.0001
lambda = exp(seq(log(lambda.max), log(lambda.max*0.05), length = nlambda))
lambda2 = 0.1
gr.n =50
gr = list()
gr.size = rep(0,gr.n)
gr.d = d/gr.n
for(i in 1:gr.n){
  gr[[i]] = c((gr.d*(i-1)+1):(gr.d*i))
  gr.size[i] = length(gr[[i]])
}

lambda = NULL
nlambda = 20
lambda.min.ratio = 0.05
family = "gaussian"
method="lasso"
res.sd = FALSE
design.sd = FALSE
prec = 1e-4
max.ite = 1e4
alg = "stoc"
truncation=0.05
max.act.in=2
verbose = TRUE
lambda_max = max(abs(crossprod(X,Y)))/n
lambda_min = 0.05*lambda_max
lambda = exp(seq(log(lambda_max),log(lambda_min),length=20))

out.l1 = picasso(X, Y, nlambda = nlambda,lambda.min.ratio = lambda.min.ratio, 
                 family = family, method = "l1", alg=alg, design.sd=design.sd, 
                 max.act.in=max.act.in, truncation=truncation)
beta.l1 = out.l1$beta
ite.l1 = out.l1$ite
intcpt.l1 = out.l1$intercept
# out.l1$runt
out.l1$df
out.l1$lambda

out.gr = picasso(X, Y, lambda = lambda,lambda.min.ratio = lambda.min.ratio, 
                 family = family, method = "group", gr.d=2, alg=alg, design.sd=design.sd, 
                 max.act.in=max.act.in, truncation=truncation)
beta.gr = out.gr$beta
ite.gr = out.gr$ite
intcpt.gr = out.gr$intercept
# out.gr$runt
out.gr$df
gr.gr = out.gr$gr
grn.gr = out.gr$gr.n
grsize.gr = out.gr$gr.size

out.gr2 = picasso(X, Y, lambda = lambda,lambda.min.ratio = lambda.min.ratio, 
                 family = family, method = "group.mcp", gr.d=2, alg=alg, design.sd=design.sd, 
                 max.act.in=max.act.in, truncation=truncation)
beta.gr2 = out.gr2$beta
ite.gr2 = out.gr2$ite
intcpt.gr2 = out.gr2$intercept
# out.gr$runt
out.gr2$df
gr.gr2 = out.gr2$gr
grn.gr2 = out.gr2$gr.n
grsize.gr2 = out.gr2$gr.size

out.gr3 = picasso(X, Y, lambda = lambda,lambda.min.ratio = lambda.min.ratio, 
                 family = family, method = "group.scad", gr.d=2, alg=alg, design.sd=design.sd, 
                 max.act.in=max.act.in, truncation=truncation)
beta.gr3 = out.gr3$beta
ite.gr3 = out.gr3$ite
intcpt.gr3 = out.gr3$intercept
# out.gr$runt
out.gr3$df
gr.gr3 = out.gr3$gr
grn.gr3 = out.gr3$gr.n
grsize.gr3 = out.gr3$gr.size

out.scad = picasso(X, Y, nlambda = nlambda,lambda.min.ratio = lambda.min.ratio, 
                   family = family, method = "mcp",alg = alg,gamma=3, design.sd=design.sd, 
                   max.act.in=max.act.in, truncation=truncation)
beta.scad = out.scad$beta
ite.scad = out.scad$ite
intcpt.scad = out.scad$intercept
# out.scad$runt
out.scad$df

out.mcp = picasso(X, Y, nlambda = nlambda,lambda.min.ratio = lambda.min.ratio, 
                  family = family, method = "scad",alg = alg,gamma=3, design.sd=design.sd, 
                  max.act.in=max.act.in, truncation=truncation)
beta.mcp = out.mcp$beta
ite.mcp = out.mcp$ite
intcpt.mcp = out.mcp$intercept
# out.mcp$runt
out.mcp$df

graphics.off()
par(mfrow=c(3,2))
plot(out.l1)
plot(out.gr)
plot(out.mcp)
plot(out.gr2)
plot(out.scad)
plot(out.gr2)
# 
# begt = Sys.time()
# out.glm = glmnet(X, Y, family="gaussian",nlambda = nlambda,lambda.min.ratio = lambda.min.ratio)
# beta.glm = as.matrix(out.glm$beta)
# Sys.time()-begt
# out.glm$df
