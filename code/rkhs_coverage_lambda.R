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
# Try with linear kernel (vanilladot), no intercept (if possible)

d_vals <- seq(-0.5, 0.5, length.out = 10)
v <- 0
true_dr <- 1.2*d_vals

replications <- 1000
r <- 100

base_nm <- 'rkhs_coverage_lambda'
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


# CBLB SCALE BASELINE----
hyper_grid <- as.data.table(expand.grid(n = c(1000),
                                        subsets = c(10),
                                        lambda = 10^(-5:-10),
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
    lambda <- grid_val$lambda
    
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

        first_inverse <- solve(K_DD * K_VV * K_XX + b * lambda * diag(b))
        second_inverse <- solve(K_VV + b * lambda * diag(b))
        
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
      
    }, cl = 5)
    
    out <- rbindlist(out)
    out[, `:=`(n = n,
               gamma = gamma,
               subsets = subsets,
               lambda = lambda,
               B = B)]
    out
  })
  cblb <- rbindlist(cblb)
  saveRDS(cblb, file.path(temp_dir, 'blb_coverage_baseline.rds'))
}

cblb[, .(coverage = mean(lower_ci <= true_dr & upper_ci >= true_dr)), by = c('n', 'B', 'tr_vals', 'subsets', 
                                                                             'gamma', 'lambda')]


# CBLB SCALE LAMBDA VARIATION----
hyper_grid <- as.data.table(expand.grid(n = c(1000),
                                        subsets = c(10),
                                        lambda = 10^seq(0, -10, length.out=20),
                                        scaling = c('n', 'b'),
                                        B = c(100)))
hyper_grid[, `:=`(gamma = calculate_gamma(n, subsets))]
seq_row <- seq_len(nrow(hyper_grid))

if(file.exists(file.path(temp_dir, 'blb_coverage_lambda_variation.rds'))){
  cblb <- readRDS(file.path(temp_dir, 'blb_coverage_lambda_variation.rds'))
} else{
  cblb <- lapply(seq_row, function(i){
    grid_val <- hyper_grid[i]
    n <- grid_val$n
    gamma <- grid_val$gamma
    b <- floor(n^gamma)
    part_idx <- seq_len(b)
    B <- grid_val$B
    subsets <- grid_val$subsets
    lambda <- grid_val$lambda
    scaling <- grid_val$scaling
    
    if(scaling == 'n'){
      scaler <- n
    } else{
      scaler <- b
    }
    
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
        
        first_inverse <- solve(K_DD * K_VV * K_XX + scaler * lambda * diag(b))
        second_inverse <- solve(K_VV + scaler * lambda * diag(b))
        
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
      
    }, cl = 5)
    
    out <- rbindlist(out)
    out[, `:=`(n = n,
               gamma = gamma,
               subsets = subsets,
               lambda = lambda,
               scaling = scaling,
               B = B)]
    out
  })
  cblb <- rbindlist(cblb)
  saveRDS(cblb, file.path(temp_dir, 'blb_coverage_lambda_variation.rds'))
}

View(cblb[, .(coverage = mean(lower_ci <= true_dr & upper_ci >= true_dr)), by = c('n', 'B', 'tr_vals', 'subsets', 
                                                                             'gamma', 'lambda', 'scaling')])


