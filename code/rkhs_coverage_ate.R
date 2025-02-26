library(data.table)
library(pbapply)
library(ggplot2)
source('code/helper_funcs.R')

compute_gcv <- function(Y, K, lambda_values) {
  n <- nrow(K)
  gcv_scores <- numeric(length(lambda_values))
  for (i in 1:length(lambda_values)) {
    lambda <- lambda_values[i]
    
    K_lambda <- K + n*lambda*diag(n)
    inv_K_lambda <- solve(K_lambda)
    
    S_lambda <- K %*% inv_K_lambda
    trace_S <- sum(diag(S_lambda))
    
    residuals <- Y %*% S_lambda - Y

    gcv_scores[i] <- mean((residuals/(1-trace_S/n))^2)
  }
  
  return(gcv_scores)
}

# Try running with larger n BLB
# Multiple smaller lambdas BLB
# Look at kernels for full n = 100 and subset of size b = 100 from n = 1000
# Try with linear kernel (vanilladot), no intercept (if possible) (Do this first)

d_vals <- seq(-0.5, 0.5, length.out = 10)
true_dr <- 1.2*d_vals

replications <- 1000
r <- 100

base_nm <- 'rkhs_coverage_ate'
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

# FULL SIMULATIONS----
if(file.exists(file.path(temp_dir, 'full_coverage.rds'))){
  full <- readRDS(file.path(temp_dir, 'full_coverage.rds'))
} else{
  full <- lapply(seq_row, function(i){
    grid_val <- hyper_grid[i]
    n <- grid_val$n
    B <- grid_val$B
    idx <- seq_len(n)

    out <- pblapply(seq_len(replications), function(rp){
      set.seed(rp)
      dat <- rkhs_sim(n = n)
      rbf <- kernlab::polydot()
      
      K_DD <- kernlab::kernelMatrix(rbf, as.matrix(dat[, c('D')]))
      K_XX <- kernlab::kernelMatrix(rbf, as.matrix(dat[, c('X1', 'X2', 'V')]))
      K_Dd <- kernlab::kernelMatrix(rbf, as.matrix(dat[, c('D')]), as.matrix(d_vals))

      lambda_values <- 10^seq(-4, 4, length.out = 20)
      gcv_gamma <- compute_gcv(dat$Y, K_DD * K_XX, lambda_values)
      lambda_1 <- lambda_values[which.min(gcv_gamma)]
      # 
      first_inverse <- solve(K_DD *K_XX + n * lambda_1 * diag(n))
      first_term <- t(dat$Y) %*% first_inverse

      ate_terms <- lapply(seq_len(nrow(dat)), function(id){
        K_Xx <- kernlab::kernelMatrix(rbf, as.matrix(dat[, c('X1', 'X2', 'V')]),
                                      y = as.matrix(dat[rep(id, length(d_vals)), c('X1', 'X2', 'V')]))
        first_term %*% (K_Dd*K_Xx)
      })

      ate_terms <- do.call(rbind, ate_terms)
      M <- rmultinom(B, size = n, prob = rep(1, n))
      boot_reps <- lapply(1:B, function(bt){
        (M[, bt] %*% ate_terms)/n
      })

      boot_reps <- do.call(rbind, boot_reps)
      perc_ci <- apply(boot_reps, 2, boot:::perc.ci)
      return(data.table(
        lower_ci = perc_ci[4,],
        upper_ci = perc_ci[5,],
        tr_vals = d_vals,
        true_dr = true_dr
      ))
      
    }, cl = 5) 
    out <- rbindlist(out)
    out[, `:=`(n = n,
               B = B)]
    out
  })
  full <- rbindlist(full)
  saveRDS(full, file.path(temp_dir, 'full_coverage.rds'))
}
full[, .(coverage = mean(lower_ci <= true_dr & upper_ci >= true_dr)), by = c('n', 'B', 'tr_vals', 'true_dr')]

# FULL SIMULATIONS REFIT----
if(file.exists(file.path(temp_dir, 'full_coverage_refit.rds'))){
  full <- readRDS(file.path(temp_dir, 'full_coverage_refit.rds'))
} else{
  full <- lapply(seq_row, function(i){
    grid_val <- hyper_grid[i]
    n <- grid_val$n
    B <- grid_val$B
    idx <- seq_len(n)
    
    out <- pblapply(seq_len(replications), function(rp){
      set.seed(rp)
      dat <- rkhs_sim(n = n)
      rbf <- kernlab::polydot()
      
      K_DD <- kernlab::kernelMatrix(rbf, as.matrix(dat[, c('D')]))
      K_XX <- kernlab::kernelMatrix(rbf, as.matrix(dat[, c('X1', 'X2', 'V')]))
      K_Dd <- kernlab::kernelMatrix(rbf, as.matrix(dat[, c('D')]), as.matrix(d_vals))
      
      lambda_values <- 10^seq(-4, 4, length.out = 20)
      gcv_gamma <- compute_gcv(dat$Y, K_DD * K_XX, lambda_values)
      lambda_1 <- lambda_values[which.min(gcv_gamma)]

      boot_reps <- lapply(1:B, function(bt){
        boot_dat <- dat[sample(idx, n, replace = TRUE), ]
        K_DD <- kernlab::kernelMatrix(rbf, as.matrix(boot_dat[, c('D')]))
        K_XX <- kernlab::kernelMatrix(rbf, as.matrix(boot_dat[, c('X1', 'X2', 'V')]))
        K_Dd <- kernlab::kernelMatrix(rbf, as.matrix(boot_dat[, c('D')]), as.matrix(d_vals))
        
        first_inverse <- solve(K_DD *K_XX + n * lambda_1 * diag(n))
        first_term <- t(boot_dat$Y) %*% first_inverse
        
        ate_terms <- lapply(seq_len(nrow(boot_dat)), function(id){
          K_Xx <- kernlab::kernelMatrix(rbf, as.matrix(boot_dat[, c('X1', 'X2', 'V')]),
                                        y = as.matrix(boot_dat[rep(id, length(d_vals)), c('X1', 'X2', 'V')]))
          first_term %*% (K_Dd*K_Xx)
        })
        ate_terms <- do.call(rbind, ate_terms)
        colMeans(ate_terms)
      })
      
      boot_reps <- do.call(rbind, boot_reps)
      perc_ci <- apply(boot_reps, 2, boot:::perc.ci)
      return(data.table(
        lower_ci = perc_ci[4,],
        upper_ci = perc_ci[5,],
        tr_vals = d_vals,
        true_dr = true_dr
      ))
      
    }, cl = 5) 
    out <- rbindlist(out)
    out[, `:=`(n = n,
               B = B)]
    out
  })
  full <- rbindlist(full)
  saveRDS(full, file.path(temp_dir, 'full_coverage_refit.rds'))
}
full[, .(coverage = mean(lower_ci <= true_dr & upper_ci >= true_dr)), by = c('n', 'B', 'tr_vals', 'true_dr')]

# CBLB SCALE BASELINE----
hyper_grid <- data.table(n = c(1000, 10000),
                         subsets = c(5, 20),
                         B = 100)
hyper_grid[, `:=`(gamma = calculate_gamma(n, subsets))]
seq_row <- seq_len(nrow(hyper_grid))

if(file.exists(file.path(temp_dir, 'blb_coverage_baseline.rds'))){
  cblb <- readRDS(file.path(temp_dir, 'blb_coverage_baseline.rds'))
} else{
  cblb <- lapply(seq_row, function(i){
    grid_val <- hyper_grid[i]
    n <- grid_val$n
    gamma <- grid_val$gamma
    b <- floor(n^gamma)
    part_idx <- seq_len(b)
    B <- grid_val$B
    subsets <- grid_val$subsets
    
    out <- pblapply(seq_len(replications), function(rp){
      set.seed(rp)
      dat <- rkhs_sim(n = n)
      rbf <- kernlab::polydot()
      partitions <- make_partition(n = n, subsets = subsets, b = b, disjoint = TRUE)
      blb_out <- lapply(partitions, function(i){
        tmp_dat <- dat[i, ]
        
        K_DD <- kernlab::kernelMatrix(rbf, as.matrix(tmp_dat[, c('D')]))
        K_XX <- kernlab::kernelMatrix(rbf, as.matrix(tmp_dat[, c('X1', 'X2', 'V')]))
        K_Dd <- kernlab::kernelMatrix(rbf, as.matrix(tmp_dat[, c('D')]), as.matrix(d_vals))
        
        lambda_values <- 10^seq(-4, 4, length.out = 20)
        gcv_gamma <- compute_gcv(tmp_dat$Y, K_DD * K_XX, lambda_values)
        lambda_1 <- lambda_values[which.min(gcv_gamma)]

        first_inverse <- solve(K_DD *K_XX + b *lambda_1 * diag(b))
        first_term <- t(tmp_dat$Y) %*% first_inverse
        ate_terms <- lapply(seq_len(nrow(tmp_dat)), function(id){
          K_Xx <- kernlab::kernelMatrix(rbf, as.matrix(tmp_dat[, c('X1', 'X2', 'V')]),
                                        y = as.matrix(tmp_dat[rep(id, length(d_vals)), c('X1', 'X2', 'V')]))
          
          first_term %*% (K_Dd*K_Xx)
        })
        ate_terms <- do.call(rbind, ate_terms)
        M <- rmultinom(B, size = n, prob = rep(1, b))
        blb_reps <- lapply(1:B, function(bt){
          (M[, bt] %*% ate_terms)/n
        })
        blb_reps <- do.call(rbind, blb_reps)
        
        perc_ci <- apply(blb_reps, 2, boot:::perc.ci)
        return(data.table(lower_ci = perc_ci[4, ],
                          upper_ci = perc_ci[5, ],
                          tr_vals = d_vals,
                          true_dr = true_dr))
      })
      blb_out <- rbindlist(blb_out)
      blb_out <- blb_out[, .(lower_ci = mean(lower_ci),
                             upper_ci = mean(upper_ci)), by = c('true_dr', 'tr_vals')]    

    }, cl = 5)
    
    out <- rbindlist(out)
    out[, `:=`(n = n,
               gamma = gamma,
               subsets = subsets,
               B = B)]
    out
  })
  cblb <- rbindlist(cblb)
  saveRDS(cblb, file.path(temp_dir, 'blb_coverage_baseline.rds'))
}

cblb[, .(coverage = mean(lower_ci <= true_dr & upper_ci >= true_dr)), by = c('n', 'B', 'tr_vals', 'subsets', 'gamma')]


