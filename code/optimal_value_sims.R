library(data.table)
library(pbapply)
library(kernlab)

source('code/helper_funcs.R')

aol_dgp <- function(n){
  X <- lapply(1:5, function(i){
    dt <- data.frame(runif(n, -1, 1))
    names(dt) <- paste0('x', i)
    dt
  })
  X <- do.call(cbind, X)
  A <- ifelse(rbinom(n, size = 1, prob = 0.5) == 1, 1, -1)
  mu_y <- 0.5 + as.matrix(X) %*% c(0.5, 0.8, 0.3, -0.5, 0.7) + A*(0.2 - 0.6*X$x1 - 0.8*X$x2)
  y <- rnorm(n, mean = mu_y[, 1])
  X$y <- y
  X$A <- A
  return(X)
}

huber_hinge <- function(u, delta = 1) {
  ifelse(u >= 1, 0, ifelse(u >= -1, (1 - u)^2 / 4, -u))
}

train_aol <- function(dat){
  mu_y <- lm(y ~ x1 + x2 + x3 + x4 + x5 + A + A:x1 + A:x2, data = dat)
  
  pred_dat <- copy(dat)
  pred_dat$A <- -1
  mu_y0 <- predict(mu_y, pred_dat)
  
  pred_dat <- copy(dat)
  pred_dat$A <- 1
  mu_y1 <- predict(mu_y, pred_dat)
  
  g_tilde <- dat$y - (0.5 * mu_y0 + 0.5 * mu_y1)
  
  return(g_tilde)
}

# Objective function for optimization
aol_loss <- function(params, X, A, r_tilde, K_matrix, lambda) {
  n <- length(A)
  v <- params[1:n]  # Extract v coefficients
  b <- params[n+1]   # Extract intercept
  
  f_x <- K_matrix %*% v + b

  hinge_loss <- huber_hinge(A * sign(r_tilde) * f_x)
  
  loss_term <- mean(abs(r_tilde) / 0.5 * hinge_loss)
  reg_term <- (lambda / 2) * sum(v %*% K_matrix %*% v)
  
  return(loss_term + reg_term)
}


optimal_val <- 1.001
sigma <- 1
replications <- 1000
K <- 10

base_nm <- 'optimal_outcome'
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

hyper_grid <- as.data.table(expand.grid(n = c(1000),
                                        B = c(100)))
seq_row <- seq_len(nrow(hyper_grid))

# NAIVE NO REFIT PROCEDURE----
if(file.exists(file.path(temp_dir, 'naive_norefit.rds'))){
  cblb <- readRDS(file.path(temp_dir, 'naive_norefit.rds'))
} else{
  cblb <- lapply(seq_row, function(i){
    grid_val <- hyper_grid[i]
    n <- grid_val$n
    prop_form <- grid_val$prop_form
    out_form <- grid_val$out_form
    B <- grid_val$B
    
    subsets <- grid_val$subsets
    out <- pblapply(seq_len(replications), function(rp){
      set.seed(rp)
      dat <- aol_dgp(n = n)
      idx <- seq_len(n)
      lambda <- 0.01
      initial_params <- c(rep(0, n), 0)  # Initial v and b
      
      X <- dat[, c('x1', 'x2', 'x3', 'x4', 'x5')]
      A <- dat$A
      K_matrix <- kernelMatrix(vanilladot(), as.matrix(X))
      r_tilde <- train_aol(dat)

      opt_result <- optim(
        par = initial_params,
        fn = aol_loss,
        X = X,
        A = A,
        r_tilde = r_tilde,
        K_matrix = K_matrix,
        lambda = lambda,
        method = "L-BFGS-B"
      )
      
      # Extract optimized parameters
      v_opt <- opt_result$par[1:n]
      b_opt <- opt_result$par[n+1]
      decision_boundary <- K_matrix %*% v_opt + b_opt
      estim_opt_regime <- ifelse(decision_boundary > 0, 1, -1)

      # full_estim <- sum(dat$y/0.5*(dat$A == estim_opt_regime))/n
      M <- rmultinom(n = B, size = n, prob = rep(1, n))
      
      blb_reps <- sapply(seq_len(B), function(bt){
        sum(M[, bt]*dat$y/0.5*(dat$A == estim_opt_regime))/n
      })

      perc_ci <- boot:::perc.ci(blb_reps)
      return(data.table(lower_ci = c(perc_ci[4]),
                        upper_ci = c(perc_ci[5]),
                        type = c('Percentile'),
                        estim = c(mean(blb_reps)),
                        se = c(sd(blb_reps))))
      }, cl = 5)

    out <- rbindlist(out)
    out[, `:=`(n = n,
               B = B)]
    out
  })
  cblb <- rbindlist(cblb)
  saveRDS(cblb, file.path(temp_dir, 'naive_norefit.rds'))
}
cblb[, .(mean(lower_ci <= optimal_val & upper_ci >= optimal_val)), by = c('n', 'type')]

# DISJOINT----
hyper_grid <- data.table(n = c(1000),
                         subsets = c(5),
                         B = 100)
hyper_grid[, `:=`(gamma = calculate_gamma(n, subsets))]
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
    subsets <- grid_val$subsets
    
    out <- pblapply(seq_len(replications), function(rp){
      set.seed(rp)
      dat <- aol_dgp(n = n)
      idx <- seq_len(b)
      lambda <- 0.01
      partitions <- make_partition(n = n, subsets = subsets, b = b, disjoint = TRUE)
      
      blb_out <- lapply(partitions, function(i){
        tmp_dat <- dat[i, ]
        initial_params <- c(rep(0, b), 0)  # Initial v and b
      
        X <- tmp_dat[, c('x1', 'x2', 'x3', 'x4', 'x5')]
        A <- tmp_dat$A
        K_matrix <- kernelMatrix(vanilladot(), as.matrix(X))
        r_tilde <- train_aol(tmp_dat)
        
        opt_result <- optim(
          par = initial_params,
          fn = aol_loss,
          X = X,
          A = A,
          r_tilde = r_tilde,
          K_matrix = K_matrix,
          lambda = lambda,
          method = "L-BFGS-B"
        )
        
        # Extract optimized parameters
        v_opt <- opt_result$par[1:b]
        b_opt <- opt_result$par[b+1]
        decision_boundary <- K_matrix %*% v_opt + b_opt
        estim_opt_regime <- ifelse(decision_boundary > 0, 1, -1)
        
        M <- rmultinom(n = B, size = n, prob = rep(1, b))
        blb_reps <- sapply(seq_len(B), function(bt){
          sum(M[, bt]*tmp_dat$y/0.5*(tmp_dat$A == estim_opt_regime))/n
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
      blb_out
    }, cl = 5)
    
    out <- rbindlist(out)
    out[, `:=`(n = n,
               gamma = gamma,
               subsets = subsets,
               B = B)]
    out
  })
  cblb <- rbindlist(cblb)
  saveRDS(cblb, file.path(temp_dir, 'disjoint_subsets.rds'))
}

cblb[, .(mean(lower_ci <= optimal_val & upper_ci >= optimal_val)), by = c('n', 'gamma', 'subsets','B')]
