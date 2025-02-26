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
cl <- 5

replications <- 1
r <- 100

base_nm <- 'rkhs_coverage_comparison'
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
                                        B = c(100)))
hyper_grid[, `:=`(gamma = calculate_gamma(n, subsets))]
seq_row <- seq_len(nrow(hyper_grid))

if(file.exists(file.path(temp_dir, 'comparison.rds'))){
  cblb <- readRDS(file.path(temp_dir, 'comparison.rds'))
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
      K_DD_full <- kernlab::kernelMatrix(rbf, as.matrix(dat[, c('D')]))
      K_VV_full <- kernlab::kernelMatrix(rbf, as.matrix(dat[, c('V')]))
      K_XX_full <- kernlab::kernelMatrix(rbf, as.matrix(dat[, c('X1', 'X2')]))
      lambda_values <- 10^seq(-3, 3, length.out = 15)
      gcv_gamma <- compute_gcv(dat$Y, K_DD_full * K_VV_full * K_XX_full, lambda_values)
      lambda_1 <- lambda_values[which.min(gcv_gamma)]
      
      gcv_mu_x <- compute_gcv(K_VV_full, K_XX_full, lambda_values)  
      lambda_2 <- lambda_values[which.min(gcv_mu_x)]
      
      first_inverse <- solve(K_DD_full * K_VV_full * K_XX_full + n * lambda_1 * diag(n))
      second_inverse <- solve(K_VV_full + n * lambda_1 * diag(n))

      full_out <- purrr::map2(list(K_DD_full, K_VV_full, K_XX_full, first_inverse, second_inverse), 
                              list('DD', 'VV', 'XX', 'first_inverse', 'second_inverse'), ~{
        data.table(values = c(.x), matrix_type = .y)
      })
      full_out <- rbindlist(full_out)
      full_out[, `:=`(rep = rp, type = 'full')]
      partitions <- make_partition(n = n, subsets = subsets, b = b, disjoint = TRUE)
      blb_out <- lapply(seq_along(partitions), function(i){
        tmp_dat <- dat[partitions[[i]], ]
        
        K_DD <- kernlab::kernelMatrix(rbf, as.matrix(tmp_dat[, c('D')]))
        K_VV <- kernlab::kernelMatrix(rbf, as.matrix(tmp_dat[, c('V')]))
        K_XX <- kernlab::kernelMatrix(rbf, as.matrix(tmp_dat[, c('X1', 'X2')]))
        lambda_values <- 10^seq(-3, 3, length.out = 15)
        gcv_gamma <- compute_gcv(tmp_dat$Y, K_DD * K_VV * K_XX, lambda_values)
        lambda_1 <- lambda_values[which.min(gcv_gamma)]
        
        gcv_mu_x <- compute_gcv(K_VV, K_XX, lambda_values)  
        lambda_2 <- lambda_values[which.min(gcv_mu_x)]
        
        first_inverse <- solve(K_DD * K_VV * K_XX + b * lambda_1 * diag(b))
        second_inverse <- solve(K_VV + b * lambda_2 * diag(b))
        blb_out_in <- purrr::map2(list(K_DD, K_VV, K_XX, first_inverse, second_inverse), 
                                  list('DD', 'VV', 'XX', 'first_inverse', 'second_inverse'), ~{
          data.table(values = c(.x), matrix_type = .y)
        })
        blb_out_in <- rbindlist(blb_out_in)
        blb_out_in[, `:=`(subset = i, type = 'subset')]
        return(blb_out_in)
      })
      return(rbindlist(list(full_out, rbindlist(blb_out))))
    }, cl = cl)
    
    out <- rbindlist(out)
    out
  })
  cblb <- rbindlist(cblb)
  saveRDS(cblb, file.path(temp_dir, 'comparison.rds'))
}


reps <- unique(cblb$rep)
for(r in reps){
  filt <- cblb[(rep == r) | (rep == 1 & type == 'full')]
  p <- ggplot(filt, aes(x = values, fill = type)) +
    geom_density(alpha = 0.3) +
    scale_fill_manual(values = c("blue", "red")) +
    scale_color_manual(values = c("blue", "red")) +
    facet_wrap(~ matrix_type) + 
    theme_minimal() +
    labs(title = paste0("Overlaid Histogram by Category", ', replication ', r),
         x = "Value",
         y = "Count",
         fill = "Category")
  print(p)
}


