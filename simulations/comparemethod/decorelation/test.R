
library(mvtnorm)
dScore <- function(Y, X, coi, beta = NULL,theta0=0,maxit = 100000, nfolds = 5, lambda = NULL,
                   refit = TRUE, refit_ratio = 0.05){
  
  # Lots of error messages
  if(!is.numeric(Y)){
    stop('Y must be a numeric vector')
  }
  
  if(!is.numeric(X)){
    stop('X must be a numeric vector')
  }
  
  if(!is.numeric(coi) || (coi%%1!=0)){
    stop('Coeffcient of interest must be integer')
  }
  
  Y <-  as.matrix(Y)
  
  if(ncol(Y) != 1){
    stop('Y must be one dimensional')
  }
  
  if(nrow(Y) != nrow(X)){
    stop('Dimensions of X and Y do not match')
  }
  
  d = ncol(X) # Get dimensions
  n = nrow(X)
  
  if(ncol(X) < coi){
    stop('Coefficient of interest exceed dimensions')
  }
  
  
  if (is.null(beta) == FALSE){
    if(length(beta) != ncol(X)){
      stop('Dimensions of coefficient vector and data matrix do not match')
    }
  } else{
      cv.fit <- cv.glmnet(X,Y,maxit = maxit,nfolds = nfolds, family = 'gaussian', lambda = lambda)
      tmp <- which(cv.fit$glmnet.fit$lambda == cv.fit$lambda.min)
      beta <- cv.fit$glmnet.fit$beta[,tmp]
      if (sum(beta!=0) == 0){
        tmp = which(cv.fit$nzero >= 1)[1]
        beta <- cv.fit$glmnet.fit$beta[,tmp]
      }
      if (refit == TRUE){
        idx_nonzero = which(beta != 0)
        if(length(idx_nonzero) > floor(d*refit_ratio)){
          #warning(paste('hard-thresholding applied to select top',floor(d*refit_ratio),'variables',sep = ' '))
          tmp = which(cv.fit$nzero >= floor(d*refit_ratio))[1]
          beta <- cv.fit$glmnet.fit$beta[,tmp]
          idx_nonzero = which(beta != 0)
        }
        
        beta_nonzero <- glm(Y~X[,idx_nonzero]-1,family = 'gaussian')$coefficients
        beta[idx_nonzero] <- beta_nonzero
      }
    # beta <- fit_lasso_wrapper(X,Y,maxit,nfolds,lambda,refit,refit_ratio)
      
    } 
  
  RRR=test_func(beta,Y,X,coi,theta0,nfolds)
  U <- RRR[1]
  adjust0=RRR[2]
  pval <- 2*pnorm(-abs(U))
  
  results <- list(fitted_value = beta,
                  test_stats = U,
                  pval = as.vector(pval),
                  adjust0=adjust0)
  
  return(results)
}

  
test_func <- function(beta,Y,X,coi,theta0,nfolds){
  ahat <- beta[coi] 
  ghat <- beta[-coi]  
  Hs <- t(X)%*%X 
  n = nrow(X)
  la <- X[,coi] # Z_i in paper. Where Q_i = (Z_i,X_i) is all covariate data
  lg <- X[, -coi] # X_i in paper (note).
  
  sig2 = mean((Y - X%*%beta)^2)
  
  what <- fit_lasso_wrapper(lg,la,refit = F)
  
  res.n <- as.vector(Y-la*theta0-lg%*%ghat)
  S <- mean(res.n*la) - t(what)%*%(colMeans(res.n*lg))
  
  var <- (Hs[coi,coi] - t(what)%*%Hs[-coi,coi])/n
  
  U <- sqrt(n)*S/sqrt(var*sig2)
  return(c(U,S*sqrt(n)))
}

fit_lasso_wrapper <- function(X,Y,maxit = 10000,nfolds = 5,lambda = NULL,refit = FALSE,refit_ratio = 1){
  cv.fit <- cv.glmnet(X,Y,maxit = maxit,nfolds = nfolds, family = 'gaussian', lambda = lambda)
  tmp <- which(cv.fit$glmnet.fit$lambda == cv.fit$lambda.min)
  beta <- cv.fit$glmnet.fit$beta[,tmp]
  if (sum(beta!=0) == 0){
    tmp = which(cv.fit$nzero >= 1)[1]
    beta <- cv.fit$glmnet.fit$beta[,tmp]
  }
  if (refit == TRUE){
    idx_nonzero = which(beta != 0)
    if(length(idx_nonzero) > floor(d*refit_ratio)){
      tmp = which(cv.fit$nzero >= floor(d*refit_ratio))[1]
      beta <- cv.fit$glmnet.fit$beta[,tmp]
      idx_nonzero = which(beta != 0)
    }
    
    beta_nonzero <- glm(Y~X[,idx_nonzero]-1,family = 'gaussian')$coefficients
    beta[idx_nonzero] <- beta_nonzero
  }
  return(beta)
}




# remotes::install_github("huijiefeng/ScoreTest")



idx2block <- function(beta_idx,p=70,q=120) {
  block_size <- (p - 1) * (q + 1)  
  j <- ceiling(beta_idx / block_size)  
  idx_j <- ((beta_idx - 1) %% block_size) + 1 
  h <- ceiling(idx_j / (p - 1))              
  idx_k <- ((idx_j - 1) %% (p - 1)) + 1       
  k_list <- 1:p
  k_list <- k_list[-j]  
  k <- k_list[idx_k]    
  return(list(j = j, k = k, h = h))
}

block2idx <- function(j, k, h, p=70, q=120) {
  k_list <- 1:p
  k_list <- k_list[-j]
  idx_k <- which(k_list == k)  # 1~(p-1)
  idx_j <- (h - 1) * (p - 1) + idx_k  # 1~(p-1)*(q+1)
  block_size <- (p - 1) * (q + 1)
  beta_idx <- (j - 1) * block_size + idx_j
  return(beta_idx)
}



get_beta_j <- function(beta, j, p=70, q=120) {
  block_size <- (p - 1) * (q + 1)
  start_idx <- (j - 1) * block_size + 1
  end_idx   <- j * block_size
  beta_j <- beta[start_idx:end_idx]
  return(beta_j)
}

get_betaj_h <- function(beta_j, h, p=70, q=120) {
  block_size <- p-1
  start_idx <- (h - 1) * block_size + 1
  end_idx   <- h* block_size
  betaj_h <- beta_j[start_idx:end_idx]
  return(betaj_h)
}


make_X_y_nodewise <- function(Z, U, j) {
  # Z: n x p
  # U: n x (q+1)
  # j: response index
  n <- nrow(Z)
  p <- ncol(Z)
  q1 <- ncol(U)  
  idx_k <- setdiff(1:p, j)
  Z_minus_j <- Z[, idx_k, drop=FALSE]
  
  X <- NULL
  for (h in 1:q1) {
    X <- cbind(X, Z_minus_j * U[, h])
  }
  y <- Z[, j]
  colnames(X)=NULL
  return(list(X = X, y = y))
}



args <- commandArgs(trailingOnly = TRUE)
n <- as.numeric(args[1])

p=70
q=120
nu=0.3

library(doParallel)
cat("Starting parallel...")
cl <- makeCluster(15)
registerDoParallel(cl)


cat("Begin working...")
results<-foreach (rep0 = 1:200,.combine = 'rbind',.packages = c("glmnet","mvtnorm")) %dopar%{
  cat(rep0)


  path = paste0("datarepo/n=", n, "p=", p, "q=", q, "nu=", nu)

  Z0<-as.matrix(read.csv(paste0(path, "/rep",rep0,"/","Z0.csv"),head=FALSE))
  iU<-as.matrix(read.csv(paste0(path, "/rep",rep0,"/","iU.csv"),head=FALSE))

  beta<-as.matrix(read.csv(paste0(path, "/beta.csv"),head=FALSE))




  list_bias=c()
  list_count=c()
  list_rej=c()
  list_len=c()


  j=1
  Data=make_X_y_nodewise(Z0,iU,j)
  x=Data$X
  y=Data$y
  Result=dScore(y,x,coi=70,theta0=0,refit=TRUE)
  idxlist <- which(Result$fitted_value != 0)
  fitted_value=Result$fitted_value[Result$fitted_value!=0]
  pvalue=c()
  test_stat=c()
  for (idx in idxlist){
  Result0=dScore(y,x,coi=idx,theta0=0,refit=TRUE)
  pvalue=c(pvalue,Result0$pval)
  test_stat=c(test_stat,Result0$test_stats)
  }
  compare_deco=cbind(idxlist,fitted_value,pvalue,test_stat)

  write.csv(compare_deco,paste0(path,"/rep",rep0,"/",'compare_deco.csv'),row.names = FALSE)

  NULL
}

stopCluster(cl)


write.csv(result_bias,paste0(path,'/','postselection_bias.csv'),row.names = FALSE)

write.csv(result_len,paste0(path,'/','postselection_len.csv'),row.names = FALSE)

write.csv(result_count,paste0(path,'/','postselection_count.csv'),row.names = FALSE)

write.csv(result_rej,paste0(path,'/','postselection_rej.csv'),row.names = FALSE)
