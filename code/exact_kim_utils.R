##################################################################################################
# Utility functions for OSQP-based kernel SBW with Nystrom approximation
##################################################################################################

##################################################################################################
# Build the full (n x n) Gaussian kernel matrix
# borrowing from the kbal R implementation by Chad Hazlett:
# "https://github.com/chadhazlett/KBAL/tree/master/R"

makeK <- function(allx, useasbases=NULL, b=NULL, linkernel = FALSE, scale = TRUE){
  N=nrow(allx)
  # If no "useasbasis" given, assume all observations are to be used.
  if(is.null(useasbases)) {useasbases = rep(1, N)}
  
  #default b is set to 2ncol to match kbal for now
  if (is.null(b)){ b=2*ncol(allx) }
  
  if(scale) {
    allx = scale(allx)
  } 
  bases = allx[useasbases==1, ]
  
  if(linkernel == TRUE) {
    K = allx
  } else {
    if(sum(useasbases) == N) {
      #symmetric K, build faster using
      K = kernel_parallel(X = allx, b = b)
    } else {
      #updated to build only triangle and mirror (4x faster)
      K = kernel_parallel_2(X = allx, Y = bases, b = b)
      #K = kernel_parallel_old(X = allx, Y = bases, b = b)
      
    }
    #old non parallel
    #new_gauss_kern(newx = allx, oldx = bases, b = b)
  }
  return(K)
}

kernel_parallel <- function(X, b) {
  .Call('_kbal_kernel_parallel', PACKAGE = 'kbal', X, b)
}

kernel_parallel_2 <- function(X, Y, b) {
  .Call('_kbal_kernel_parallel_2', PACKAGE = 'kbal', X, Y, b)
}

kernel_parallel_old <- function(X, Y, b) {
  .Call('_kbal_kernel_parallel_old', PACKAGE = 'kbal', X, Y, b)
}


##################################################################################################
# Construct the Gaussian kernel bases
# It can be done either based on
#    1) the full kernel matrix (n x n) (Hazlett, 2020)
# or 2) the reduced kernel matrix (n x l), l << n, based on the Nystrom approximation (Wang, 2019)

kernel.basis <- function(X,A,Y, 
                         kernel.approximation=TRUE,
                         dim.reduction=FALSE,
                         c = NULL, l=NULL, s=NULL, gamma=NULL, U.mat = NULL) {
  
  if (kernel.approximation) {
    if (is.null(c)) {
      # use some heuristics for choice of c 
      c = ifelse(n<1e+5, 100, 250) 
    }
    if (is.null(l)) {
      # set arbitrarily similar to Wang, 2019  
      l=round(c/2)
      s=round(l/2)
    }
    if (is.null(gamma)) {
      # heuristics by Hazlett, 2020  
      gamma = 1/(2*ncol(X)) 
    }
    
    set_c = sample(x = 1:n, size = c, replace = FALSE)
    set_c = sort(set_c)
    X <- scale(X)
    # --------------------------------------------------------
    # # slow version:  
    # f <- function(x,y) exp(-gamma*sum((x-y)^2))
    # C <- outer( 
    #   1:n, set_c, 
    #   Vectorize( function(i,j) f(X[i,], X[j,]) ) 
    # )
    # --------------------------------------------------------
    # C <- RBF_kernel_C(X, c, set_c)
    C <- RBF_kernel_C_parallel(X, c, set_c)
    # --------------------------------------------------------
    W <- C[set_c,]
    SVD_W <- svds(A=W, k=l)
    if (dim.reduction) {
      R <- C %*% SVD_W$u %*% diag(1/sqrt(SVD_W$d))
      SVD_R <- svds(A=R, k=s) 
      X_ <- R %*% SVD_R$v # B
    } else {
      X_ <- C %*% SVD_W$u # %*% diag(1/SVD_W$d)
      # lambda_ <- 1/SVD_W$d
    }
  } else{
    # O(n^2) spoce-, O(n^3) time-complexity
    gram.mat <- makeK(X)
    res = eigs_sym(gram.mat, c, which = "LM")  
    X_ <- if(is.null(U.mat)) {
      res$vectors %*% diag(1/sqrt(res$values))
    } else{
      U.mat %*% diag(1/sqrt(res$values))
    }
  }
  return(X_)
}

##################################################################################################
# Construct the SBW basis based on power series
# by using up to the k-th moment and k-th order interactions 
# (e.g., if K=2 and interactions=TRUE, then we use the 1st and 2nd moments, 
#  as well as all the pairwise products of the observed covariates)

power.basis <- function(X,A,Y,
                        K=2, interactions=FALSE) {
  for (k in 1:K){
    # add the k-th moment
    X.k <- X^k
    colnames(X.k) <- paste(colnames(X),".",k,sep = "")
    X_ <- cbind(X_, X.k)
  }
  
  if (interactions==TRUE & K >=2) {
    for (k in 2:K){
      # add the kth order interactions
      indx <- combn(colnames(X),k)
      int <- as.data.frame(do.call(cbind,
                                   lapply(split(indx, col(indx)), 
                                          function(x) rowProds(as.matrix(X[,x])))
      ))
      colnames(int) <- apply(indx, 2, function(x) paste(x,collapse="."))
      X_ <- cbind(X_, int)
    }
  }
  return(X_)
}

##################################################################################################
# The following function implements the OSQP-based kernel SBW 

osqp_kernel_sbw <- function(X,A,Y,
                            delta.v=0.005, 
                            X_=NULL,
                            osqp.setting=NULL,
                            basis="kernel", kernel.approximation=TRUE,
                            c = NULL, l=NULL, gamma=NULL, U.mat = NULL,
                            dim.reduction=FALSE, s=NULL,
                            K=2, interactions=FALSE) {
  
  #------------------------------------------------------------------------------------------
  # @X: covariates (n x d matrix)
  # @A: treatment (n x 1 vector)
  # @Y: outcome (n x 1 vector)
  # @X_: the basis matrix where (approximate) mean balancing between the treated and control groups is to be achieved
  #      (ncol(X_) = n)
  # @osqp.setting: osqp settings
  # @delta.v: tolerance level (either a single number or a numeric vector);
  #           the numeric vector is used for multiple iterations, each with different tolerance level
  # @basis: support kernel- or power series- (moment-) based constructions
  # + kernel basis (basis=='kernel')
  #     @kernel.approximation: if TRUE, use the Nystrom kernel approximation method proposed in \insertRef{wang2019scalable}
  #                            if FALSE, use the full gram matrix to compute eigenvectors as in \insertRef{hazlett2018}
  #     @c: sketch size for the standard Nystrom approximation. 
  #     @l: regularization parameter l < c. 
  #     @s: target rank for s < l the rank-restricted Nyström approximation
  #     @gamma: RBK parameter; set default value as 1/(2*dim(X))
  #     @U.mat: externel input for nxc eigenvector matrix; only used when kernel.approximation==FALSE
  # + power series basis (basis=='power')
  #     @K: use up to the K-th moment
  #     @interactions: add up to the K-th order interaction terms if TRUE
  #------------------------------------------------------------------------------------------
  
  
  res.list <- list()
  if (is.null(X_)) {
    if (basis=="kernel"){
      X_ <- kernel.basis(X,A,Y, 
                         kernel.approximation=kernel.approximation, 
                         c=c, l=l, gamma=gamma, U.mat=U.mat,
                         s=s, dim.reduction=dim.reduction)
    } 
    
    else if (basis=="power") {
      X_ <- power.basis(X,A,Y, 
                        K=K, interactions=interactions)
    } 
    
    else {
      stop("not available yet")
    }
  }
  
  nX <- ncol(X_)
  n1 <- sum(A); n0 <- sum(1-A)
  Xt <- X_[A==1,]; Xc <- X_[A==0,]
  Yt <- Y[A==1]; Yc <- Y[A==0]
  P.mat <- as(.symDiagonal(n=n0, x=1.), "dgCMatrix")
  q.vec <- rep(-1./n0,n0)
  A.mat <- Matrix(rbind(rep(1.,n0), 
                        P.mat,
                        t(Xc)), sparse = TRUE)
  if (is.null(osqp.setting)) {
    # use default one
    settings <- osqpSettings(alpha = 1.5, verbose = FALSE)  
  } else {
    settings <- osqp.setting
  }
  
  for (j in 1:length(delta.v)) {
    l.vec <- c(1., rep(0.,n0), 
               colMeans(Xt) - delta.v[j] * rep(1,nX))
    u.vec <- c(1., rep(1.,n0), 
               colMeans(Xt) + delta.v[j] * rep(1,nX))
    if (j==1) {
      model <- osqp(P.mat, q.vec, A.mat, l.vec, u.vec, settings)
    } else {
      model$Update(l = l.vec, u = u.vec)
    }
    res <- model$Solve()
    if (res$info$status != "solved") {
      warning(res$info$status)
    }
    res.list[[j]] <- list(w=res$x, y_hat=sum(res$x * Yc), t=res$info$solve_time)
  }
  return(res.list)
}   