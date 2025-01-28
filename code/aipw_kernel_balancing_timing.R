library(data.table)
library(pbapply)
library(reticulate)

kernel_weights <- function(data,degree1,degree2,k1,k2,operator,penal, bootstrap_size=length(data)){
  
  intervention <- data$Tr
  outcome <- data$y
  confounders <- data.frame(data$X1,data$X2)
  n <- nrow(data)
  
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
  matrix_eva <- as.matrix( confounders )
  # evaluation matrix
  res.optim2_1$par <- exp(res.optim2_1$gpr$kernel_$theta)
  res.optim2_0$par <- exp(res.optim2_0$gpr$kernel_$theta)
  
  K1 <- res.optim2_1$gpr$kernel_(matrix_eva)
  K0 <- res.optim2_0$gpr$kernel_(matrix_eva)
  # Gram matrices
  
  p1 <- as.numeric(res.optim2_1$gpr$predict(matrix_eva))
  p0 <- as.numeric(res.optim2_0$gpr$predict(matrix_eva))
  
  # print("Building matrices")
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

  #######################################################
  #######################################################
  #Solving QP
  #######################################################
  #######################################################
  
  # print("Solving QP")
  
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
  bvec <- c(1, 1)  # Right-hand side values for the equality constraints
  meq <- 2  # N
  
  params <- list(Presolve=2,OutputFlag=0,QCPDual=0)

  res <- quadprog::solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec, meq = meq)

  phi0 <- (1-data$Tr)*res$solution*(data$y - p0) + p0
  phi1 <- (data$Tr)*res$solution*(data$y - p1) + p1
  
  return(list(phi1 = phi1, phi0 = phi0))
  
  # assertthat::assert_that(length(phi1) != 0)
  # 
  # mp1 <- mean(phi1)
  # mp0 <- mean(phi0)
  # 
  # ate_akw <- mp1 - mp0
  # se_ate_akw <- sqrt(var(phi1-phi0)/nrow(data))
  # 
  # return(res)
  
}


use_virtualenv("r-reticulate", required = TRUE)

np <- import("numpy")
source_python("/Users/macbookair/Documents/Projects/cblb/Archive/gp_simu_gate.py")

source('code/helper_funcs.R')

te <- 0.8
sigma <- 1
replications <- 25
degree1 <- 1
degree2 <- 1
k1 <- "poly"
k2 <- "poly"
operator <- "single" # here i'm considering a simple linear kernel just for comparison/true model
penal <- log(2) # keep it between 0.1 and log(k) where k is the number of features

base_nm <- 'aipw_kernel_balancing_timing'
image_path <- 'images'
dat_path <- 'data'

# Create the temporary directory
temp_dir <- file.path(dat_path, paste0(base_nm, '_tmp'))
if(!file.exists(temp_dir)){
  dir.create(temp_dir, recursive = TRUE)
}

img_tmp_dir <- file.path(image_path, paste0(base_nm, '_tmp'))
if(!file.exists(img_tmp_dir)){
  dir.create(img_tmp_dir, recursive = TRUE)
}

hyper_grid <- as.data.table(expand.grid(n = c(1000, 10000),
                                        B = c(100),
                                        prop_form = c('correct'),
                                        out_form = c('correct')))

seq_row <- seq_len(nrow(hyper_grid))

# CORRECT PROCEDURE----
if(file.exists(file.path(temp_dir, 'bootstrap.rds'))){
  cblb <- readRDS(file.path(temp_dir, 'bootstrap.rds'))
} else{
  cblb <- lapply(seq_row, function(i){
    grid_val <- hyper_grid[i]
    n <- grid_val$n
    prop_form <- grid_val$prop_form
    out_form <- grid_val$out_form
    B <- grid_val$B
    if(prop_form == 'correct'){
      prop_formula <- c('X1', 'X2')
    } else{
      prop_formula <- c('Z1', 'Z2')
    }
    if(out_form == 'correct'){
      out_formula <- c('Tr', 'X1', 'X2')
    } else{
      out_formula <- c('Tr', 'Z1', 'Z2')
    }
    subsets <- grid_val$subsets
    out <- pblapply(seq_len(replications), function(rp){
      set.seed(rp)
      dat <- kangschafer3(n = n, te = te, sigma = sigma, beta_overlap = 0.5)
      idx <- seq_len(n)
      time <- system.time({
        output <- kernel_weights(dat,degree1,degree2,k1,k2,operator,penal)
        phi1 <- output$phi1
        phi0 <- output$phi0
        mp1 <- mean(phi1)
        mp0 <- mean(phi0)
        ate_akw <- mp1 - mp0
        se_ate_akw <- sqrt(var(phi1-phi0)/nrow(dat))
        M <- rmultinom(n = B, size = n, prob = rep(1, n))
        
        blb_reps <- sapply(seq_len(B), function(bt){
          boot_phi1 <- M[, bt]*phi1
          boot_phi0 <- M[, bt]*phi0
          mean(boot_phi1) - mean(boot_phi0)
        })
        
        perc_ci <- boot:::perc.ci(blb_reps)
      })

      data.table(time_elapsed = time['elapsed'])
      }, cl = 1)

    out <- rbindlist(out)
    out[, `:=`(n = n,
               prop_form = prop_form,
               out_form = out_form,
               B = B)]
    out
  })
  cblb <- rbindlist(cblb)
  saveRDS(cblb, file.path(temp_dir, 'bootstrap.rds'))
}

# cblb[, .(mean(lower_ci <= te & upper_ci >= te)), by = c('n', 'type')]



# DISJOINT----
hyper_grid <- rbindlist(list(data.table(n = 10000, B = 100, prop_form = 'correct', out_form = 'correct', subsets = 10, gamma = 0.75),
                             data.table(n = 10000, B = 100, prop_form = 'correct', out_form = 'correct', subsets = 15, gamma = 0.70597),
                             data.table(n = 10000, B = 100, prop_form = 'correct', out_form = 'correct', subsets = 20, gamma = 0.67474),
                             data.table(n = 1000, B = 100, prop_form = 'correct', out_form = 'correct', subsets = 4, gamma = 0.79)))

seq_row <- seq_len(nrow(hyper_grid))

if(file.exists(file.path(temp_dir, 'disjoint_subsets.rds'))){
  cblb <- readRDS(file.path(temp_dir, 'disjoint_subsets.rds'))
} else{
  cblb <- lapply(seq_row, function(i){
    grid_val <- hyper_grid[i]
    n <- grid_val$n
    gamma <- grid_val$gamma
    b <- floor(n^gamma)
    part_idx <- seq_len(b)
    prop_form <- grid_val$prop_form
    out_form <- grid_val$out_form
    B <- grid_val$B
    if(prop_form == 'correct'){
      prop_formula <- c('X1', 'X2')
    } else{
      prop_formula <- c('Z1', 'Z2')
    }
    if(out_form == 'correct'){
      out_formula <- c('Tr', 'X1', 'X2')
    } else{
      out_formula <- c('Tr', 'Z1', 'Z2')
    }
    subsets <- grid_val$subsets
    
    out <- pblapply(seq_len(replications), function(rp){
      set.seed(rp)
      dat <- kangschafer3(n = n, te = te, sigma = sigma, beta_overlap = 0.5)
      time <- system.time({
        partitions <- make_partition(n = n, subsets = subsets, b = b, disjoint = TRUE)
        
        blb_out <- lapply(partitions, function(i){
          tmp_dat <- dat[i]
          output <- kernel_weights(tmp_dat,degree1,degree2,k1,k2,operator,penal)
          phi1 <- output$phi1
          phi0 <- output$phi0
          M <- rmultinom(n = B, size = n, prob = rep(1, b))
          
          blb_reps <- sapply(seq_len(B), function(bt){
            boot_phi1 <- M[, bt]*phi1
            boot_phi0 <- M[, bt]*phi0
            sum(boot_phi1)/n - sum(boot_phi0)/n
          })
          
          
          perc_ci <- boot:::perc.ci(blb_reps)
          return(data.table(lower_ci = perc_ci[4],
                            upper_ci = perc_ci[5],
                            estim = mean(blb_reps),
                            se = sd(blb_reps)))
        })
        
        
        blb_out <- rbindlist(blb_out)
        blb_out <- blb_out[, .(lower_ci = mean(lower_ci),
                               upper_ci = mean(upper_ci),
                               estim = mean(estim),
                               se = mean(se))]
        
      })
      data.table(time_elapsed = time['elapsed'])
    }, cl = 1)
    
    out <- rbindlist(out)
    out[, `:=`(n = n,
               gamma = gamma,
               subsets = subsets,
               prop_form = prop_form,
               out_form = out_form,
               B = B)]
    out
  })
  cblb <- rbindlist(cblb)
  saveRDS(cblb, file.path(temp_dir, 'disjoint_subsets.rds'))
}

# cblb[, .(mean(lower_ci <= te & upper_ci >= te)), by = c('n', 'gamma', 'subsets', 'out_form', 'prop_form', 'B')]



