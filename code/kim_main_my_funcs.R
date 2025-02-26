#####################################################################################  
#        R function for fast OSQP-based kernel SBW with Nystrom approximation       #                
#        by Kwhangho Kim (kkim@hcp.med.harvard.edu)                                 #
#        https://arxiv.org/abs/2311.00568                                           #
#####################################################################################

library(Matrix)
library(MASS)
# library(matrixStats)
library(osqp)
# library(kbal) # optional
library(RSpectra)
library(RcppParallel)
library(Rcpp)
# library(RcppArmadillo)

source('code/helper_funcs.R')
source('code/osqp-kernel-sbw.R')


# --------------------------------------------------------------------#
#                    Example (Hainmueller, 2012)                      #
# --------------------------------------------------------------------#
# made smalelr to deal with memory
n=1e+4
sig.123 <- diag(c(2,1,1))
sig.123[1,2] <- 1; sig.123[1,3] <- -1; sig.123[2,3] <- -0.5;
sig.123 <- forceSymmetric(sig.123)
beta_coef <- c(1,2,-2,-1,-0.5,1)

X.123 <- as.matrix(mvrnorm(n, mu = rep(0,3), Sigma = sig.123))
X.4 <- runif(n,-3,3)
X.5 <- rchisq(n,1)
X.6 <- rbinom(n,1,0.5)
X <- cbind(X.123, X.4, X.5, X.6)
A <- ifelse(X %*% matrix(beta_coef, ncol = 1) + rnorm(n,0,30) > 0,1,0)
Y <- (X.123[,1] + X.123[,2] + X.5)^2 + rnorm(n,0,1)

# ptm <- proc.time()
res <- osqp_kernel_sbw(X,A,Y,
                       delta.v=1e-4,
                       kernel.approximation = FALSE,
                       c = 100)
res[[1]]$y_hat
# with rank-restricted Nystrom approximation:
res <- osqp_kernel_sbw(X,A,Y,
                       delta.v=1e-4,
                       dim.reduction=TRUE,
                       c = 100, l=75, s=50)
mean(Y[A==1]) - res[[1]]$y_hat
# et <- proc.time() - ptm


# --------------------------------------------------------------------#
#                    more details about each step                     #
# --------------------------------------------------------------------#

# Step 1) Construct the approximated kernel basis
# + the sketch size c shouldn't be too large; typically recommend using 50~250 
#   (c=100 works well enough in the above example)
# + When using many covariates, consider the rank-restricted Nystrom approximation 
#   by setting dim.reduction=TRUE
#   with the regularization (l) and target rank (s) parameters
ptm <- proc.time()
X_ <- kernel.basis(X,A,Y,
                   kernel.approximation=TRUE,
                   c = 100)
# X_ <- kernel.basis(X,A,Y,
#                    kernel.approximation=TRUE,
#                    dim.reduction=TRUE,
#                    c = 100, l=75, s=50)
et <- proc.time() - ptm


# Step 2) Compute the sbw via osqp
# + Use smaller tolerance levels for higher accuracy
# + Polishing is an additional algorithm step where OSQP tries to compute a high-accuracy solution
high.acc.setting <- osqpSettings(alpha = 1.5, verbose = FALSE, 
                                 warm_start = TRUE, # use warm start when you iterate over multiple deltas
                                 polish=TRUE,# solution polishing  
                                 eps_abs = 5e-5, # use more strict absolute tolerance if needed
                                 eps_rel = 5e-5 # use more strict relative tolerance if needed
) 

# Example for multiple deltas
delta.v <- seq(0.0005, 0.01, by=0.0005)
res <- osqp_kernel_sbw(X,A,Y, 
                       delta.v=delta.v, X_=X_,
                       osqp.setting=high.acc.setting)
length(res)

