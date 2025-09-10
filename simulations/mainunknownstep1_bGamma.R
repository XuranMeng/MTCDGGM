library(glmnet)
library(doParallel)
library(foreach)
Treptition=200

startrep = 1

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript mainunknownstep1_bGamma.R <n> [p] [q] [nu]")
}
n  <- as.numeric(args[1])
p = 70
q = 120


nu=0.3
dim=(p-1)*(q+1)



startrep=1

cl <- makeCluster(15)
registerDoParallel(cl)
results <- foreach(reptition = startrep:Treptition, .combine = "list", .packages = "glmnet") %dopar% {

path <- paste0("datarepounknown/n=", n,
               "p=", p,
               "q=", q,
               "nu=", nu,
               "/rep", reptition, "/")

Z0 <- as.matrix(read.csv(paste0(path, "Z0.csv"), header = FALSE))

iU <- as.matrix(read.csv(paste0(path, "iU.csv"), header = FALSE))

U <- iU[, -1]
A <- U
hatbGamma=NULL
for (j in 1:p){
    y=Z0[,j]

cvfit <- cv.glmnet(A, y, alpha = 1, intercept = FALSE)
beta_hat <- as.vector(coef(cvfit, s = "lambda.min"))[-1]
threshold <- cvfit$lambda.min 
beta_hat[abs(beta_hat) < threshold] <- 0
hatbGamma=cbind(hatbGamma,beta_hat)
}
write.table(hatbGamma,
            file = paste0(path, "hatbGamma.csv"),
            row.names = FALSE,
            col.names = FALSE,
            sep = ",")
return(NULL)
}
