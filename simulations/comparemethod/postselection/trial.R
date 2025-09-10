library(glmnet)
library(selectiveInference)
args <- commandArgs(trailingOnly = TRUE)
n <- as.numeric(args[1])
p=70
q=120
nu=0.3
lambda =80

cat(n)
#Select a relatively good lambda for different n

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


result_bias=c()
result_len=c()
result_count=c()
result_rej=c()
for (rep0 in 1:200){
cat(rep0)
if (n==800){
lambda=80
}
else if (n==400) {
   lambda=60
}
else {
   lambda=40
}

path = paste0("datarepo/n=", n, "p=", p, "q=", q, "nu=", nu)

Z0<-as.matrix(read.csv(paste0(path, "/rep",rep0,"/","Z0.csv"),head=FALSE))
iU<-as.matrix(read.csv(paste0(path, "/rep",rep0,"/","iU.csv"),head=FALSE))

beta<-as.matrix(read.csv(paste0(path, "/beta.csv"),head=FALSE))




list_bias=c()
list_count=c()
list_rej=c()
list_len=c()


for (j in 1:2){
Data=make_X_y_nodewise(Z0,iU,j)
x=Data$X
y=Data$y

gfit = glmnet(x,y,standardize=FALSE)

betapost = coef(gfit, x=x, y=y, s=lambda/n, exact=TRUE)[-1]
# compute fixed lambda p-values and selection intervals
out = fixedLassoInf(x,y,betapost,lambda,sigma=1,alpha=0.05)

if (70%in%out$vars){
coef1=out$coef0[out$vars==70]
CI1_low=out$ci[out$vars==70][1]
CI1_high=out$ci[out$vars==70][2]
}
else {
   coef1=NA
   CI1_low=NA
   CI1_high=NA
}
if (139 %in% out$vars){
coef2=out$coef0[out$vars==139]
CI2_low=out$ci[out$vars==139][1]
CI2_high=out$ci[out$vars==139][2]
}
else {
   coef2=NA
   CI2_low=NA
   CI2_high=NA
}
list_bias=c(list_bias,coef1+nu,coef2+nu)
list_len=c(list_len,CI1_high-CI1_low,CI2_high-CI2_low)
list_count=c(list_count,(-nu>CI1_low)&(-nu<CI1_high),(-nu>CI2_low)&(-nu<CI2_high))
list_rej=c(list_rej,(0>CI1_low)&(0<CI1_high),(0>CI2_low)&(0<CI2_high))
}
result_bias=rbind(result_bias,list_bias)
result_len=rbind(result_len,list_len)
result_count=rbind(result_count,list_count)
result_rej=rbind(result_rej,list_rej)

}


write.csv(result_bias,paste0(path,'/','postselection_bias.csv'),row.names = FALSE)

write.csv(result_len,paste0(path,'/','postselection_len.csv'),row.names = FALSE)

write.csv(result_count,paste0(path,'/','postselection_count.csv'),row.names = FALSE)

write.csv(result_rej,paste0(path,'/','postselection_rej.csv'),row.names = FALSE)





