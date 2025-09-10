install.packages("selectiveInference")
library(selectiveInference)

library(glmnet)
library(selectiveInference)

set.seed(43)
n = 200
p = 1010
sigma = 1
x = matrix(rnorm(n*p),n,p)
x = scale(x,TRUE,TRUE)
beta = c(0.3,0.3,rep(0,p-2))
y = x%*%beta + sigma*rnorm(n)
# first run glmnet
gfit = glmnet(x,y,standardize=FALSE)

# extract coef for a given lambda; note the 1/n factor!
# (and we don't save the intercept term)
lambda =10
beta = coef(gfit, x=x, y=y, s=lambda/n, exact=TRUE)[-1]
# compute fixed lambda p-values and selection intervals
out = fixedLassoInf(x,y,beta,lambda,sigma=sigma)
out
## as above, but use lar function instead to get initial 
## lasso fit (should get same results)
lfit = lar(x,y,normalize=FALSE)
beta = coef(lfit, s=lambda, mode="lambda")
out2 = fixedLassoInf(x, y, beta, lambda, sigma=sigma)
out2




set.seed(43)
n = 50
p = 10
sigma = 1
x = matrix(rnorm(n*p),n,p)
x=scale(x,TRUE,TRUE)
beta = c(3,2,rep(0,p-2))
y = x%*%beta + sigma*rnorm(n)
pf=c(rep(1,7),rep(.1,3))  #define penalty factors
pf=p*pf/sum(pf)   # penalty factors should be rescaled so they sum to p
xs=scale(x,FALSE,pf) #scale cols of x by penalty factors
# first run glmnet
gfit = glmnet(xs, y, standardize=FALSE)
# extract coef for a given lambda; note the 1/n factor!
# (and we don't save the intercept term)
lambda = .8
beta_hat = coef(gfit, x=xs, y=y, s=lambda/n, exact=TRUE)[-1]
# compute fixed lambda p-values and selection intervals
out = fixedLassoInf(xs,y,beta_hat,lambda,sigma=sigma)
#rescale conf points to undo the penalty factor
out$ci=t(scale(t(out$ci),FALSE,pf[out$vars]))
out
