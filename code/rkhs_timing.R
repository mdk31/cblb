library(data.table)
library(pbapply)
library(ggplot2)
source('code/helper_funcs.R')

compute_gcv <- function(Y, K, lambda_values) {
  n <- nrow(K)
  gcv_scores <- numeric(length(lambda_values))
  for (i in 1:length(lambda_values)) {
    lambda <- lambda_values[i]
    
    K_lambda <- K + lambda * diag(n)
    
    inv_K_lambda <- solve(K_lambda)
    
    S_lambda <- K %*% inv_K_lambda
    trace_S <- sum(diag(S_lambda))
    
    residuals <- (diag(n) - S_lambda) %*% Y
    
    gcv_scores[i] <- sum(residuals^2) / (trace_S / n)^2
  }
  
  return(gcv_scores)
}

# Generate evaluation points
# Get true effects
# Do full method (see where bottleneck is in large n)
# See how well this does in terms of bias, RMSE, etc.
# From there, think how BLB can improve this

d_vals <- seq(-1.1, 2, length.out = 10)
v <- 0
true_dr <- 1.2*d_vals

replications <- 25
r <- 100

base_nm <- 'rkhs_timing'
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
if(file.exists(file.path(temp_dir, 'full_timing.rds'))){
  cblb <- readRDS(file.path(temp_dir, 'full_timing.rds'))
} else{
  cblb <- lapply(seq_row, function(i){
    grid_val <- hyper_grid[i]
    n <- grid_val$n
    B <- grid_val$B
    idx <- seq_len(n)

    out <- pblapply(seq_len(replications), function(rp){
      set.seed(rp)
      dat <- rkhs_sim(n = n)
      time <- system.time({
        rbf <- kernlab::rbfdot()

        K_DD <- kernlab::kernelMatrix(rbf, as.matrix(dat[, c('D')]))
        K_VV <- kernlab::kernelMatrix(rbf, as.matrix(dat[, c('V')]))
        K_XX <- kernlab::kernelMatrix(rbf, as.matrix(dat[, c('X1', 'X2')]))
        K_Dd <- kernlab::kernelMatrix(rbf, as.matrix(dat[, c('D')]), as.matrix(d_vals))
        K_Vv <- kernlab::kernelMatrix(rbf, as.matrix(dat[, c('V')]), as.matrix(rep(v, length(d_vals))))
        
        lambda_values <- 10^seq(-3, 3, length.out = 10)
        gcv_gamma <- compute_gcv(dat$Y, K_DD*K_VV*K_XX, lambda_values)
        lambda_1 <- lambda_values[which.min(gcv_gamma)]
        
        gcv_mu_x <- compute_gcv(K_VV, K_XX, lambda_values)  # Assuming X[,1] as a proxy for X
        lambda_2 <- lambda_values[which.min(gcv_mu_x)]
        
        first_inverse <- solve(K_DD*K_VV*K_XX + n*lambda_1*diag(n))
        second_inverse <- solve(K_VV + n*lambda_2*diag(n))
        
        M <- rmultinom(B, size = n, prob = rep(1, n))
        boot_reps <- lapply(1:B, function(bt){
          (M[, bt]*dat$Y) %*% first_inverse %*% (K_Dd*K_Vv*(K_XX %*% second_inverse %*% K_Vv))
        })
        boot_reps <- do.call(rbind, boot_reps)
        perc_ci <- apply(boot_reps, 2, boot:::perc.ci)
      })
      data.table(time_elapsed = time['elapsed'])
    }, cl = 1) 
    out <- rbindlist(out)
    out[, `:=`(n = n,
               B = B)]
    out
  })
  cblb <- rbindlist(cblb)
  saveRDS(cblb, file.path(temp_dir, 'full_timing.rds'))
}
print(cblb)

# CBLB TIMING----
hyper_grid <- as.data.table(expand.grid(n = c(1000),
                                        subsets = c(5),
                                        B = c(100)))
hyper_grid[, `:=`(gamma = calculate_gamma(n, subsets))]
seq_row <- seq_len(nrow(hyper_grid))

if(file.exists(file.path(temp_dir, 'blb_timing.rds'))){
  cblb <- readRDS(file.path(temp_dir, 'blb_timing.rds'))
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
      time <- system.time({
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
          
          gcv_mu_x <- compute_gcv(K_VV, K_XX, lambda_values)  # Assuming X[,1] as a proxy for X
          lambda_2 <- lambda_values[which.min(gcv_mu_x)]
          
          first_inverse <- solve(K_DD*K_VV*K_XX + b*lambda_1*diag(b))
          second_inverse <- solve(K_VV + b*lambda_2*diag(b))
          
          M <- rmultinom(B, size = n, prob = rep(1, b))
          blb_reps <- lapply(1:B, function(bt){
            (M[, bt]*tmp_dat$Y) %*% first_inverse %*% (K_Dd*K_Vv*(K_XX %*% second_inverse %*% K_Vv))
          })
          blb_reps <- do.call(rbind, blb_reps)

          perc_ci <- apply(blb_reps, 2, boot:::perc.ci)
          return(data.table(lower_ci = perc_ci[4, ],
                            upper_ci = perc_ci[5, ],
                            tr_vals = d_vals))
        })
        blb_out <- rbindlist(blb_out)
        blb_out <- blb_out[, .(lower_ci = mean(lower_ci),
                               upper_ci = mean(upper_ci)), by = 'tr_vals']     
        
      })
      data.table(time_elapsed = time['elapsed'])
    }, cl = 1)
    
    out <- rbindlist(out)
    out[, `:=`(n = n,
               gamma = gamma,
               subsets = subsets,
               B = B)]
    out
  })
  cblb <- rbindlist(cblb)
  saveRDS(cblb, file.path(temp_dir, 'blb_timing.rds'))
}

full <- readRDS(file.path(temp_dir, 'full_timing.rds'))
cblb <- readRDS(file.path(temp_dir, 'blb_timing.rds'))

cblb[, `:=`(type = paste0('BLB, subsets = ', subsets, ', gamma = ', gamma),
            gamma = NULL,
            subsets = NULL)]
full[, `:=`(type = 'Full')]
timing <- rbindlist(list(full, cblb))

for(i in c(1000, 50000)){
  filt <- timing[n == i]
  p <- ggplot(filt, aes(x = type, y = time_elapsed)) +
    geom_boxplot() +
    ggtitle(paste0('n = ', i)) +
    ylab('Time Elapsed (seconds)') +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  ggsave(paste('n_', i, '.png'), plot = p, width = 6, height = 4, dpi = 300)
  print(p)
}
