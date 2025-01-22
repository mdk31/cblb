rm(list = ls())

library(reticulate)
use_virtualenv("r-reticulate", required = TRUE)

np <- import("numpy")
source_python("/Users/macbookair/Documents/Projects/cblb/Archive/gp_simu_gate.py")

#################    ##################################    #################
#################    ##################################    #################
#
# Data generation function
#
#################    ##################################    #################
#################    ##################################    #################


data_generation <- function(n,correct){
  
  #Data
  X1 <- rnorm(n,0,1)
  X2 <- rnorm(n,0,1)
  prt <- 1/(1+exp( -(-0.5 + X1 + X2) ))
  
  
  A <- rbinom(n,1,prt)
  
  truet <- effect_t <- 0.5
  
  
  Y0  <- X1 + X2 + rnorm(n,0,1)
  Y1  <- Y0 + effect_t 
  Y   <- Y0*(1-A) + Y1*(A)
  
  dta <- data.frame(Y,X1,X2,A)
  
  wrong2 <- log(abs(dta$X2))
  wrong1 <- (dta$X2)/(exp(dta$X1))
  
  
  Z1 <- correct*dta$X1 + (1-correct)*(wrong1)
  Z2 <- correct*dta$X2 + (1-correct)*(wrong2)
  Z <- cbind(Z1,Z2)
  colnames(Z) <- c("Z1","Z2")
  dta <- data.frame(dta,Z,Y1,Y0)
  
  
  return(dta)
  
}#end function


#################    ##################################    #################
#################    ##################################    #################
#
# Naive estimator - biased
#
#################    ##################################    #################
#################    ##################################    #################


ATE_est_means <- function(data){
  
  n <- nrow(data)
  
  #-------------Among controls
  e0 <- sum((1-data$A)*data$Y)/sum((1-data$A)) #E[Y_0]
  
  h0_0 <- (1-data$A) * as.numeric(t(data$Y -  e0) )
  
  one <- diag(n)
  g11 <- sum(1-data$A)/n
  g11inv <- 1/g11
  
  phi_mu0_hat <- g11inv * h0_0
  # B <- t(h0_0) %*% h0_0
  # A <- g11inv
  # temp0 <- glm(Y~1,data=data[data$A==0,])
  # diag(A%*%B%*%t(A))/n^2; diag(sandwich(temp0))
  # same as 
  # sum((phi_mu0_hat - mean(phi_mu0_hat))^2)/n^2; var(phi_mu0_hat)*((n-1)/n^2)
  # !!!! Andrew this is the same a your derivation sum(phi_mu0_hat^2)/n^2 + mean(phi_mu0_hat)^2 !!!!!
  
  se0 <- sqrt(var(phi_mu0_hat)/n) * sqrt((n-1)/n)
  
  #-------------Among treated
  e1 <- sum((data$A)*data$Y)/sum((data$A))
  
  h1_0 <- (data$A) * as.numeric(t(data$Y -  e1) )
  
  one <- diag(n)
  g11 <- sum(data$A)/n
  g11inv <- 1/g11
  
  phi_mu1_hat <- g11inv * h1_0
  # B <- t(h1_0) %*% h1_0
  # A <- g11inv
  # temp1 <- glm(Y~1,data=data[data$A==1,])
  # diag(A%*%B%*%t(A))/n^2; diag(sandwich(temp1))
  # same as 
  # sum((phi_mu1_hat - mean(phi_mu1_hat))^2)/n^2; var(phi_mu1_hat)*((n-1)/n^2)
  # !!!! Andrew this is the same a your derivation sum(phi_mu1_hat^2)/n^2 + mean(phi_mu1_hat)^2 !!!!!
  
  se1 <- sqrt(var(phi_mu1_hat)/n) * sqrt((n-1)/n)
  
  
  #-------------ATE
  
  ate <- e1 - e0 
  se_ate <- sqrt(var(phi_mu1_hat-phi_mu0_hat)/n) * sqrt((n-1)/n)
  
  return(list(mu0 = e0,
              se0 = se0,
              mu1 = e1,
              se1 = se1,
              ate_hat = ate,
              se_ate_hat = se_ate))
  
}


#################    ##################################    #################
#################    ##################################    #################
#
# Outcome regression - just for comparison
#
#################    ##################################    #################
#################    ##################################    #################


ATE_est_OR <- function(formula_m0,formula_m1,data){
  
  nn <- nrow(data)
  
  #-------------Among controls
  
  #regress among controls
  m0_only_cont <- glm(formula_m0,data=data[data$A==0,])
  
  Y <- data$Y
  X0 <- model.matrix(formula_m0, data = data,drop = FALSE)
  
  xbeta <- X0 %*% coef(m0_only_cont) # same as predict(m0_only_cont, newdata = data)
  e0 <- mean(xbeta) 
  
  h0_0 <- X0 *  (1-data$A) * as.numeric(t(Y -  xbeta) )
  h0_1 <- xbeta - e0
  
  g11 <- ( t(X0) %*% diag(1-data$A) %*% X0) / nn
  g11inv <- solve(g11)
  
  g21 <- apply(X0,2,sum)/nn
  
  g22 <- nn
  g22inv <- 1/nn
  
  phi_beta_hat <- g11inv %*% t(h0_0)
  # B <- t(h0_0) %*% h0_0
  # A <- g11inv
  # diag(A%*%B%*%t(A))/nn^2; diag(sandwich(m0_only_cont))
  # same as 
  # diag(phi_beta_hat %*% t(phi_beta_hat)/nn^2)
  
  phi_mu0 <- (xbeta - e0) - (t(phi_beta_hat)%*%g21)
  se0 <- sqrt(var(phi_mu0)/nn) * sqrt((nn-1)/nn)
  
  #-------------Among treated
  
  #regress among treated
  m1_only_cont <- glm(formula_m1,data=data[data$A==1,])
  
  Y <- data$Y
  X1 <- model.matrix(formula_m1, data = data,drop = FALSE)
  
  xbeta <- X1 %*% coef(m1_only_cont) # same as predict(m0_only_cont, newdata = data)
  e1 <- mean(xbeta) 
  
  h1_0 <- X1 *  (data$A) * as.numeric(t(Y -  xbeta) )
  h1_1 <- xbeta - e1
  
  g11 <- ( t(X1) %*% diag(data$A) %*% X1) / nn
  g11inv <- solve(g11)
  
  g21 <- apply(X1,2,sum)/nn
  
  g22 <- nn
  g22inv <- 1/nn
  
  phi_beta_hat <- g11inv %*% t(h1_0)
  # B <- t(h1_0) %*% h1_0
  # A <- g11inv
  # diag(A%*%B%*%t(A))/nn^2; diag(sandwich(m1_only_cont))
  # same as 
  # diag(phi_beta_hat %*% t(phi_beta_hat)/nn^2)
  
  phi_mu1 <- (xbeta - e1) - (t(phi_beta_hat)%*%g21)
  se1 <- sqrt(var(phi_mu1)/nn) * sqrt((nn-1)/nn)
  
  
  
  #-------------ATE
  
  ate <- e1 - e0 
  se_ate <- sqrt(var(phi_mu1-phi_mu0)/nn) * sqrt((nn-1)/nn)
  
  return(list(mu0 = e0,
              se0 = se0,
              mu1 = e1,
              se1 = se1,
              ate_hat = ate,
              se_ate_hat = se_ate))
  
}


#################    ##################################    #################
#################    ##################################    #################
#
# Kernel weights - https://davidahirshberg.bitbucket.io/static/amle-slides.pdf
# https://arxiv.org/pdf/2110.14831 
#################    ##################################    #################
#################    ##################################    #################

kernel_weights <- function(data,degree1,degree2,k1,k2,operator,penal){
  
  intervention <- data$A
  outcome <- data$Y
  confounders <- data.frame(data$X1,data$X2)
  
  
  #######################################################
  #######################################################
  #tune the hyperparametes
  #######################################################
  #######################################################
  
  
  
  t1 <-as.integer(intervention)
  t0 <-as.integer((1-intervention))
  
  
  y <- outcome
  Xtemp <- data.frame(confounders,y,intervention)
  
  
  Xsample <- Xtemp
  X1t <- Xsample[which(Xsample$intervention==1),]
  X0t <- Xsample[which(Xsample$intervention==0),]
  
  y1 <- y[which(Xsample$intervention==1)]
  y0 <- y[which(Xsample$intervention==0)]
  
  n1 <- length(which(Xsample$intervention==1))
  n0 <- length(which(Xsample$intervention==0))
  
  
  mX0 <- as.matrix(X0t[,1:(dim(X0t)[2]-2)])
  mX1 <- as.matrix(X1t[,1:(dim(X1t)[2]-2)])
  mY0 <- as.matrix(y0)
  mY1 <- as.matrix(y1)
  
  pyX0   <- np_array( np$array(mX0), dtype="float")
  pyX1   <- np_array( np$array(mX1), dtype="float")
  pyY0   <- np$array( mY0, dtype="float" )
  pyY1   <- np$array( mY1, dtype="float" )
  
  #source_python("rcode/gp_simu_gate.py")
  #source_python("gp_simu_gate.py")
  
  res.optim2_0    <- tryCatch(gp(pyX0,pyY0,
                                 degree1=degree1,
                                 degree2=degree2,
                                 k1=k1,
                                 k2=k2,
                                 operator=operator),
                              error=function(e) NULL)
  
  
  res.optim2_1    <- tryCatch(gp(pyX1,pyY1,
                                 degree1=degree1,
                                 degree2=degree2,
                                 k1=k1,
                                 k2=k2,
                                 operator=operator),
                              error=function(e) NULL)
  
  
  
  
  #compute K
  system_time_matrices <- system.time({
    matrix_eva <- as.matrix( confounders )
    # evaluation matrix
    res.optim2_1$par <- exp(res.optim2_1$gpr$kernel_$theta)
    res.optim2_0$par <- exp(res.optim2_0$gpr$kernel_$theta)
    
    K1 <- res.optim2_1$gpr$kernel_(matrix_eva)
    K0 <- res.optim2_0$gpr$kernel_(matrix_eva)
    # Gram matrices
    
    p1 <- as.numeric(res.optim2_1$gpr$predict(matrix_eva))
    p0 <- as.numeric(res.optim2_0$gpr$predict(matrix_eva))
    
    print("Building matrices")
    V <- rep(1/n,n)
    
    #Quadratic part
    I1 <- diag( t1 )
    I0 <- diag( t0 )
    I1KI1 <- I1%*%K1%*%I1 #maybe this can be improved.
    I0KI0 <- I0%*%K0%*%I0
    
    #Linear part
    # en^T K1 I1 + en^T K0 I0
    #K1 I1, K0 I0
    KI1 <- diag(t1)%*%K1
    KI0 <- diag(t0)%*%K0
    
    # en^T K1 I1, en^T K0 I0
    VKI1 <- t(V)%*%KI1
    VKI0 <- t(V)%*%KI0
    
    tol <- 1e-08
    
    sigma1 <- res.optim2_1$par[3]^2
    sigma0 <- res.optim2_0$par[3]^2
    
    Sigma <- sigma1*diag(t1) + sigma0*diag(t0)
    
    #Update Q
    Q <- (1/n^2)*( I1KI1 + I0KI0 + penal*Sigma )
    
    #Update c
    #Gurobi:
    c <- -2*(1/n^2)*(VKI1 + VKI0)
    
    rm(list = c("VKI1","VKI0"))
    
  })#end system.time
  
  #######################################################
  #######################################################
  #Solving QP
  #######################################################
  #######################################################
  
  print("Solving QP")
  
  model <- list()
  model$A          <- matrix(c( t1/n ,t0/n), nrow=2, byrow=T)
  model$rhs        <- c(1,1)
  model$modelsense <- "min"
  model$Q          <- Q
  model$obj        <- c
  model$sense      <- c("=")
  model$lb <- rep(tol,n)
  model$vtypes <- "C"
  Dmat <- Q  # Symmetric positive-definite matrix for the quadratic term
  dvec <- c  # Linear coefficients
  Amat <- t(matrix(c(t1/n, t0/n), nrow = 2, byrow = TRUE))
  # Amat <- rbind(t1/n, t0/n)  # Coefficients for the equality constraints
  bvec <- c(1, 1)  # Right-hand side values for the equality constraints
  meq <- 2  # N

  params <- list(Presolve=2,OutputFlag=0,QCPDual=0)
  
  # system_time_gurobi <- system.time(res <- tryCatch(gurobi(model,params),error=function(e) NULL))
  res <- quadprog::solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec, meq = meq)

  # same as AIPW just using the approximated weights instead of 1/p(A|X)
  # res$solution
  phi0 <- (1-data$A)*res$solution*(data$Y - p0) + p0
  phi1 <- (data$A)*res$solution*(data$Y - p1) + p1
  
  assertthat::assert_that(length(phi1) != 0)
  
  mp1 <- mean(phi1)
  mp0 <- mean(phi0)
  
  ate_akw <-mp1 - mp0
  se_ate_akw <- sqrt(var(phi1-phi0)/nrow(data))
  
  return(list(weights=res$solution,
              mp1=mp1,
              mp0=mp0,
              ate_akw=ate_akw,
              se_ate_akw=se_ate_akw))
  
}



#################    ##################################    #################
#################    ##################################    #################
#
# Simulation
#
#################    ##################################    #################
#################    ##################################    #################

n <- 1000
ate_akw <- se_ate_akw <- ate_kw <- ate_n <- se_ate_n <- se_ate_or <- ate_or <- NULL

#kernel weights parameters
degree1 <- 1
degree2 <- 1
k1 <- "poly"
k2 <- "poly"
operator <- "single" # here i'm considering a simple linear kernel just for comparison/true model

penal <- log(2) # keep it between 0.1 and log(k) where k is the number of features

m0 <- m1 <- Y ~ X1 + X2

itera <- 1
for(i in 1:itera){
  print(i)
  data <- data_generation(n,1)
  kw <- kernel_weights(data,degree1,degree2,k1,k2,operator,penal)
  
  #ate_kw[i] <- lm(Y  ~ A, data = data, weights=kw$x)$coeff[2]
  mu0 <- sum((1-data$A)*kw$weights*data$Y)/sum((1-data$A)*kw$weights)
  mu1 <- sum((data$A)*kw$weights*data$Y)/sum((data$A)*kw$weights)
  ate_kw[i] <- mu1 - mu0
  ate_akw[i] <- kw$ate_akw
  se_ate_akw[i] <- kw$se_ate_akw
  temp_n <- ATE_est_means(data)
  ate_n[i] <- temp_n$ate_hat
  se_ate_n[i] <- temp_n$se_ate_hat
  temp_or <- ATE_est_OR(m1,m0,data)
  ate_or[i] <- temp_or$ate_hat
  se_ate_or[i] <- temp_or$se_ate_hat
  
}

mean(ate_akw); sd(ate_akw); mean(se_ate_akw)
mean(ate_kw); sd(ate_kw)
mean(ate_n); sd(ate_n); mean(se_ate_n) #naive estimator - this is biased
mean(ate_or); sd(ate_or); mean(se_ate_or) #this is just for comparison



# > mean(ate_akw); sd(ate_akw); mean(se_ate_akw)
# [1] 0.5033198
# [1] 0.07457336
# [1] 0.07650964
# > mean(ate_kw); sd(ate_kw)
# [1] 0.5082344
# [1] 0.07454306
# > mean(ate_n); sd(ate_n); mean(se_ate_n) #this is biased
# [1] 1.964262
# [1] 0.1050023
# [1] 0.101024
# > mean(ate_or); sd(ate_or); mean(se_ate_or)
# [1] 0.5018955
# [1] 0.07427242
# [1] 0.07599452

