library(data.table)
library(pbapply)
library(KernSmooth)
source('code/helper_funcs.R')

# Generate evaluation points
# Get true effects
# Do full method (see where bottleneck is in large n)
# See how well this does in terms of bias, RMSE, etc.
# From there, think how BLB can improve this

a.vals <- seq(-3.3, 3.4, length.out = 100)
true_dr <- 0.75*a.vals

te <- 0.8
sigma <- 1
replications <- 25
r <- 100

base_nm <- 'dose_response_timing'
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

hyper_grid <- as.data.table(expand.grid(n = c(1000, 50000),
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
      dat <- continuous_treatment_sim(n = n, sigma = sigma, beta_overlap = 0.5)
      time <- system.time({
        model_output <- estimate_models(data = dat, a.vals = a.vals)
        
        boot_reps <- replicate(B, {
          boot_idx <- sample(idx, size = n, replace = TRUE)
          boot_dat <- dat[boot_idx, ]
          boot_pseudo <- model_output$pseudo.out[boot_idx]
          approx(locpoly(boot_dat$exposure, boot_pseudo, bandwidth = model_output$bandwidth), 
                 xout = a.vals)$y
          
        })
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
hyper_grid <- as.data.table(expand.grid(n = c(1000, 50000),
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
      dat <- continuous_treatment_sim(n = n, sigma = sigma, beta_overlap = 0.5)
      time <- system.time({
        partitions <- make_partition(n = n, subsets = subsets, b = b, disjoint = TRUE)
        blb_out <- lapply(partitions, function(i){
          tmp_dat <- dat[i, ]
          model_output <- estimate_models(data = tmp_dat, a.vals = a.vals)
          
          blb_reps <- replicate(B, {
            boot_idx <- sample(part_idx, size = n, replace = TRUE)
            boot_dat <- tmp_dat[boot_idx, ]
            boot_pseudo <- model_output$pseudo.out[boot_idx]
            approx(locpoly(boot_dat$exposure, boot_pseudo, bandwidth = model_output$bandwidth), 
                   xout = a.vals)$y
            
          })
          perc_ci <- apply(blb_reps, 2, boot:::perc.ci)
          return(data.table(lower_ci = perc_ci[4, ],
                            upper_ci = perc_ci[5, ],
                            tr_vals = a.vals))
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

for(i in c(1000, 50000))
