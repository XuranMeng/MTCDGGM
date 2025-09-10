
l




#true null
dScore(Y, X, 4,family = 'gaussian',refit = T)$pval
?dScore
#binomial
Y  <-  sign(X%*%bt + rlogis(n)) 
#Y <- (Y+1)/2
#true alternative
dScore(Y, X, 1,family = 'binomial')$pval
#true null
dScore(Y, X, 5,family = 'binomial')$pval
dScore(Y, X, 5,family = 'binomial')

#poisson
Y  <-  rpois(n,exp(X%*%bt))
#true alternative
dScore(Y, X, 1,family = 'poisson')$pval
#true null
dScore(Y, X, 40,family = 'poisson',refit = T,nfolds = 10)$pval
