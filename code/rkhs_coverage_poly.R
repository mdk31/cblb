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
v <- 0
true_dr <- 1.2*d_vals

full_lambda_1 <- 0.0004291934
full_lambda_2 <- 0.001757511

replications <- 1000
r <- 100

base_nm <- 'rkhs_coverage'
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

hyper_grid <- as.data.table(expand.grid(n = c(100),
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
      rbf <- kernlab::rbfdot()
      
      K_DD <- kernlab::kernelMatrix(rbf, as.matrix(dat[, c('D')]))
      K_VV <- kernlab::kernelMatrix(rbf, as.matrix(dat[, c('V')]))
      browser()
      K_XX <- kernlab::kernelMatrix(rbf, as.matrix(dat[, c('X1', 'X2')]))
      K_Dd <-
        kernlab::kernelMatrix(rbf, as.matrix(dat[, c('D')]), as.matrix(d_vals))
      K_Vv <-
        kernlab::kernelMatrix(rbf, as.matrix(dat[, c('V')]), as.matrix(rep(v, length(d_vals))))
      
      lambda_values <- 10^seq(-3, 3, length.out = 10)
      gcv_gamma <- compute_gcv(dat$Y, K_DD * K_VV * K_XX, lambda_values)
      lambda_1 <- lambda_values[which.min(gcv_gamma)]
      
      gcv_mu_x <- compute_gcv(K_VV, K_XX, lambda_values)  
      lambda_2 <- lambda_values[which.min(gcv_mu_x)]

      first_inverse <- solve(K_DD * K_VV * K_XX + n * 0.0001 * diag(n))
      second_inverse <- solve(K_VV + n * 0.0001 * diag(n))
      
      M <- rmultinom(B, size = n, prob = rep(1, n))
      boot_reps <- lapply(1:B, function(bt) {
        (M[, bt] * dat$Y) %*% first_inverse %*% (K_Dd * K_Vv * (K_XX %*% second_inverse %*% K_Vv))
      })
      
      boot_reps <- do.call(rbind, boot_reps)
      perc_ci <- apply(boot_reps, 2, boot:::perc.ci)
      return(data.table(
        lower_ci = perc_ci[4,],
        upper_ci = perc_ci[5,],
        tr_vals = d_vals,
        true_dr = true_dr
      ))
      
    }, cl = 1) 
    out <- rbindlist(out)
    out[, `:=`(n = n,
               B = B)]
    out
  })
  full <- rbindlist(full)
  saveRDS(full, file.path(temp_dir, 'full_coverage.rds'))
}
full[, .(coverage = mean(lower_ci <= true_dr & upper_ci >= true_dr)), by = c('n', 'B', 'tr_vals', 'true_dr')]

# CBLB SCALE BASELINE----
hyper_grid <- as.data.table(expand.grid(n = c(1000),
                                        subsets = c(10),
                                        B = c(100)))
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
      rbf <- kernlab::rbfdot()
      partitions <- make_partition(n = n, subsets = subsets, b = b, disjoint = TRUE)
      blb_out <- lapply(partitions, function(i){
        tmp_dat <- dat[i, ]
        
        K_DD <- kernlab::kernelMatrix(rbf, as.matrix(tmp_dat[, c('D')]))
        K_VV <- kernlab::kernelMatrix(rbf, as.matrix(tmp_dat[, c('V')]))
        K_XX <- kernlab::kernelMatrix(rbf, as.matrix(tmp_dat[, c('X1', 'X2')]))
        K_Dd <- kernlab::kernelMatrix(rbf, as.matrix(tmp_dat[, c('D')]), as.matrix(d_vals))
        K_Vv <- kernlab::kernelMatrix(rbf, as.matrix(tmp_dat[, c('V')]), as.matrix(rep(v, length(d_vals))))
        
        # lambda_values <- 10^seq(-3, 3, length.out = 10)
        # gcv_gamma <- compute_gcv(tmp_dat$Y, K_DD*K_VV*K_XX, lambda_values)
        # lambda_1 <- lambda_values[which.min(gcv_gamma)]
        # 
        # gcv_mu_x <- compute_gcv(K_VV, K_XX, lambda_values)
        # lambda_2 <- lambda_values[which.min(gcv_mu_x)]
        first_inverse <- solve(K_DD * K_VV * K_XX + b * 0.00001 * diag(b))
        second_inverse <- solve(K_VV + b * 0.00001 * diag(b))
        
        M <- rmultinom(B, size = n, prob = rep(1, b))
        blb_reps <- lapply(1:B, function(bt){
          (M[, bt] * tmp_dat$Y) %*% first_inverse %*% (K_Dd * K_Vv * (K_XX %*% second_inverse %*% K_Vv))
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
      
    }, cl = 1)
    
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


# CBLB REWEIGHTING----
hyper_grid <- as.data.table(expand.grid(n = c(1000),
                                        subsets = c(10),
                                        B = c(100)))
hyper_grid[, `:=`(gamma = calculate_gamma(n, subsets))]
seq_row <- seq_len(nrow(hyper_grid))

if(file.exists(file.path(temp_dir, 'blb_coverage.rds'))){
  cblb <- readRDS(file.path(temp_dir, 'blb_coverage.rds'))
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
      rbf <- kernlab::rbfdot()
      partitions <- make_partition(n = n, subsets = subsets, b = b, disjoint = TRUE)
      blb_out <- lapply(partitions, function(i){
        tmp_dat <- dat[i, ]
        
        K_DD <- kernlab::kernelMatrix(rbf, as.matrix(tmp_dat[, c('D')]))
        K_VV <- kernlab::kernelMatrix(rbf, as.matrix(tmp_dat[, c('V')]))
        K_XX <- kernlab::kernelMatrix(rbf, as.matrix(tmp_dat[, c('X1', 'X2')]))
        K_Dd <- kernlab::kernelMatrix(rbf, as.matrix(tmp_dat[, c('D')]), as.matrix(d_vals))
        K_Vv <- kernlab::kernelMatrix(rbf, as.matrix(tmp_dat[, c('V')]), as.matrix(rep(v, length(d_vals))))
        
        lambda_values <- 10^seq(-3, 3, length.out = 10)
        gcv_gamma <- compute_gcv(tmp_dat$Y, K_DD*K_VV*K_XX, lambda_values)
        lambda_1 <- lambda_values[which.min(gcv_gamma)]
        
        gcv_mu_x <- compute_gcv(K_VV, K_XX, lambda_values)
        lambda_2 <- lambda_values[which.min(gcv_mu_x)]
        
        M <- rmultinom(B, size = n, prob = rep(1, b))
        blb_reps <- lapply(1:B, function(bt){
          W <- M[, bt]
          W_sqrt <- sqrt(W)
          W_diag <- diag(W_sqrt)
          Y_w <- W*tmp_dat$Y
          K_DD_w <- W_diag %*% K_DD %*% W_diag
          K_VV_w <- W_diag %*% K_VV %*% W_diag
          K_XX_w <- W_diag %*% K_XX %*% W_diag
          # K_Dd_w <- W_diag %*% K_Dd %*% W_diag
          # K_Vv_w <- W_diag %*% K_Vv %*% W_diag
          
          first_inverse_w <- solve(K_DD_w*K_VV_w*K_XX_w + n*lambda_1*diag(b))
          second_inverse_w <- solve(K_VV_w + n*lambda_2*diag(b))
          t(Y_w) %*% (first_inverse_w) %*% (K_Dd*K_Vv*(K_XX_w %*% second_inverse_w) %*% K_Vv)
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
  saveRDS(cblb, file.path(temp_dir, 'blb_coverage.rds'))
}

cblb[, .(coverage = mean(lower_ci <= true_dr & upper_ci >= true_dr)), by = c('n', 'B', 'tr_vals', 'subsets', 'gamma')]


# CBLB REWEIGHTING SCALE----
# hyper_grid <- as.data.table(expand.grid(n = c(1000, 50000),
#                                         subsets = c(1),
#                                         B = c(100)))
hyper_grid <- data.table(n = c(1000, 1000),
                         subsets = c(5, 10),
                         B = 100)
hyper_grid[, `:=`(gamma = calculate_gamma(n, subsets))]
seq_row <- seq_len(nrow(hyper_grid))

if(file.exists(file.path(temp_dir, 'blb_coverage_rescale.rds'))){
  cblb <- readRDS(file.path(temp_dir, 'blb_coverage_rescale.rds'))
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
      rbf <- kernlab::rbfdot()
      partitions <- make_partition(n = n, subsets = subsets, b = b, disjoint = TRUE)
      blb_out <- lapply(partitions, function(i){
        tmp_dat <- dat[i, ]
        
        K_DD <- kernlab::kernelMatrix(rbf, as.matrix(tmp_dat[, c('D')]))
        K_VV <- kernlab::kernelMatrix(rbf, as.matrix(tmp_dat[, c('V')]))
        K_XX <- kernlab::kernelMatrix(rbf, as.matrix(tmp_dat[, c('X1', 'X2')]))
        K_Dd <- kernlab::kernelMatrix(rbf, as.matrix(tmp_dat[, c('D')]), as.matrix(d_vals))
        K_Vv <- kernlab::kernelMatrix(rbf, as.matrix(tmp_dat[, c('V')]), as.matrix(rep(v, length(d_vals))))
        
        lambda_values <- 10^seq(-3, 3, length.out = 10)
        gcv_gamma <- compute_gcv(tmp_dat$Y, K_DD*K_VV*K_XX, lambda_values)
        lambda_1 <- lambda_values[which.min(gcv_gamma)]
        
        gcv_mu_x <- compute_gcv(K_VV, K_XX, lambda_values)
        lambda_2 <- lambda_values[which.min(gcv_mu_x)]
        
        M <- rmultinom(B, size = n, prob = rep(1, b))
        blb_reps <- lapply(1:B, function(bt){
          W <- M[, bt]
          W_sqrt <- sqrt(W)
          W_diag <- diag(W_sqrt)
          Y_w <- W*tmp_dat$Y
          K_DD_w <- W_diag %*% K_DD %*% W_diag
          K_VV_w <- W_diag %*% K_VV %*% W_diag
          K_XX_w <- W_diag %*% K_XX %*% W_diag
          # K_Dd_w <- apply(K_Dd, 2, function(x) x*W)
          # K_Vv_w <- apply(K_Dd, 2, function(x) x*W)

          first_inverse_w <- solve(K_DD_w*K_VV_w*K_XX_w + n*(b/n)*lambda_1*diag(b))
          second_inverse_w <- solve(K_VV_w + n*(b/n)*lambda_2*diag(b))
          t(Y_w) %*% (first_inverse_w) %*% (K_Dd*K_Vv*(K_XX_w %*% second_inverse_w) %*% K_Vv)
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
  saveRDS(cblb, file.path(temp_dir, 'blb_coverage_rescale.rds'))
}

cblb[, .(coverage = mean(lower_ci <= true_dr & upper_ci >= true_dr)), by = c('n', 'B', 'tr_vals', 'subsets', 'gamma')]



# CBLB REWEIGHTING FULL LAMBDA----
hyper_grid <- as.data.table(expand.grid(n = c(1000),
                                        subsets = c(10),
                                        B = c(100)))
hyper_grid[, `:=`(gamma = calculate_gamma(n, subsets))]
seq_row <- seq_len(nrow(hyper_grid))

if(file.exists(file.path(temp_dir, 'blb_coverage_fulllambda.rds'))){
  cblb <- readRDS(file.path(temp_dir, 'blb_coverage_fulllambda.rds'))
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
      rbf <- kernlab::rbfdot()
      partitions <- make_partition(n = n, subsets = subsets, b = b, disjoint = TRUE)
      blb_out <- lapply(partitions, function(i){
        tmp_dat <- dat[i, ]
        
        K_DD <- kernlab::kernelMatrix(rbf, as.matrix(tmp_dat[, c('D')]))
        K_VV <- kernlab::kernelMatrix(rbf, as.matrix(tmp_dat[, c('V')]))
        K_XX <- kernlab::kernelMatrix(rbf, as.matrix(tmp_dat[, c('X1', 'X2')]))
        K_Dd <- kernlab::kernelMatrix(rbf, as.matrix(tmp_dat[, c('D')]), as.matrix(d_vals))
        K_Vv <- kernlab::kernelMatrix(rbf, as.matrix(tmp_dat[, c('V')]), as.matrix(rep(v, length(d_vals))))
        
        # lambda_values <- 10^seq(-3, 3, length.out = 10)
        # gcv_gamma <- compute_gcv(tmp_dat$Y, K_DD*K_VV*K_XX, lambda_values)
        # lambda_1 <- lambda_values[which.min(gcv_gamma)]
        # 
        # gcv_mu_x <- compute_gcv(K_VV, K_XX, lambda_values)
        # lambda_2 <- lambda_values[which.min(gcv_mu_x)]
        
        M <- rmultinom(B, size = n, prob = rep(1, b))
        blb_reps <- lapply(1:B, function(bt){
          W <- M[, bt]
          W_sqrt <- sqrt(W)
          W_diag <- diag(W_sqrt)
          Y_w <- W*tmp_dat$Y
          K_DD_w <- W_diag %*% K_DD %*% W_diag
          K_VV_w <- W_diag %*% K_VV %*% W_diag
          K_XX_w <- W_diag %*% K_XX %*% W_diag
          # K_Dd_w <- apply(K_Dd, 2, function(x) x*W)
          # K_Vv_w <- apply(K_Dd, 2, function(x) x*W)
          
          first_inverse_w <- solve(K_DD_w*K_VV_w*K_XX_w + n*full_lambda_1*diag(b))
          second_inverse_w <- solve(K_VV_w + n*full_lambda_2*diag(b))
          t(Y_w) %*% (first_inverse_w) %*% (K_Dd*K_Vv*(K_XX_w %*% second_inverse_w) %*% K_Vv)
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
  saveRDS(cblb, file.path(temp_dir, 'blb_coverage_fulllambda.rds'))
}

cblb[, .(coverage = mean(lower_ci <= true_dr & upper_ci >= true_dr)), by = c('n', 'B', 'tr_vals', 'subsets', 'gamma')]



