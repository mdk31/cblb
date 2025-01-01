library(data.table)
library(pbapply)
library(caret)

source('code/helper_funcs.R')

te <- 0.8
sigma <- 1
reps <- 5
r <- 100
max_subsets <- 100
n <- 100e3
gamma <- 0.7
b <- round(n^gamma)
cl <- 4

idx <- seq_len(n)

base_nm <- 'ci_convergence'
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

calculate_estimate <- function(dt, n = NULL, weights = NULL){
  if(is.null(n)){
    n <- nrow(dt)
  }
  m <- glm(y ~ Tr + X1 + X2, data = dt, weights = weights)
  g <- glm(Tr ~ X1 + X2, data = dt, family = 'binomial', weights = weights)
  
  pred_dat <- data.table::copy(dt)
  
  dt$prop_score <- predict(g, type = 'response')
  
  pred_dat[, `:=`(Tr = 1)]
  dt$m1 <- predict(m, newdata = pred_dat)
  
  pred_dat[, `:=`(Tr = 0)]
  dt$m0 <- predict(m, newdata = pred_dat)
  phi1 <- (dt$Tr/dt$prop_score)*(dt$y - dt$m1) + dt$m1
  phi0 <- (1 - dt$Tr)/(1 - dt$prop_score)*(dt$y - dt$m0) + dt$m0
  if(!is.null(weights)){
    phi1 <- weights*phi1
    phi0 <- weights*phi0
  }
  
  sum(phi1)/n - sum(phi0)/n

}

# TRUE CI
if(!file.exists(file.path(temp_dir, 'true_bootstrap.rds'))){
  bootstrap_cis <- pblapply(seq_len(reps), function(i){
  set.seed(i)
  dat <- kangschafer3(n = n, te = te, sigma = sigma, beta_overlap = 0.5)
  
  boots <- sapply(seq_len(5000), function(q){
    boot_dat <- dat[sample(idx, replace = TRUE)]
    return(calculate_estimate(boot_dat))
  })

  data.table(lower = boot:::perc.ci(boots)[4], upper = boot:::perc.ci(boots)[5],
             reps = i)
  
}, cl = cl)

bootstrap_cis <- rbindlist(bootstrap_cis)
saveRDS(bootstrap_cis, file.path(temp_dir, 'true_bootstrap.rds'))
}


# bc-CBLB SIMULATIONS----
if(file.exists(file.path(temp_dir, 'ci_convergence.rds'))){
  cblb <- readRDS(file.path(temp_dir, 'ci_convergence.rds'))
} else{
  cblb <- pblapply(seq_len(reps), function(rp){
      set.seed(rp)
      dat <- kangschafer3(n = n, te = te, sigma = sigma, beta_overlap = 0.5)
      tau_hat_full <- calculate_estimate(dat)
      partitions <- make_partition(n = n, subsets = max_subsets, b = b, disjoint = FALSE)
      
      blb_out <- lapply(partitions, function(i){
        tmp_dat <- dat[i]
        M <- rmultinom(n = r, size = n, prob = rep(1, b))
        
        blb_reps <- sapply(seq_len(r), function(bt){
          return(calculate_estimate(dt = tmp_dat, n = n, weights = M[, bt]))
        })
        
        bias <- mean(blb_reps) - tau_hat_full
        blb_bias_correct <- blb_reps - bias
        perc_ci <- boot:::perc.ci(blb_bias_correct)
        return(data.table(lower_ci = perc_ci[4],
                          upper_ci = perc_ci[5]))
      })

      
      blb_out <- rbindlist(blb_out)
      blb_out <- blb_out[, .(lower = cumsum(lower_ci)/seq_len(max_subsets),
                             upper = cumsum(upper_ci)/seq_len(max_subsets),
                             reps = rp,
                             subset_num = seq_len(max_subsets))]
      blb_out
    }, cl = cl)
  
  cblb <- rbindlist(cblb)
  saveRDS(cblb, file.path(temp_dir, 'ci_convergence.rds'))
}

# CLASSIC BLB
if(file.exists(file.path(temp_dir, 'ci_convergence_blb.rds'))){
  cblb <- readRDS(file.path(temp_dir, 'ci_convergence_blb.rds'))
} else{
  cblb <- pblapply(seq_len(reps), function(rp){
    set.seed(rp)
    dat <- kangschafer3(n = n, te = te, sigma = sigma, beta_overlap = 0.5)
    partitions <- make_partition(n = n, subsets = max_subsets, b = b, disjoint = FALSE)
    
    blb_out <- lapply(partitions, function(i){
      tmp_dat <- dat[i]
      M <- rmultinom(n = r, size = n, prob = rep(1, b))
      
      blb_reps <- sapply(seq_len(r), function(bt){
        return(calculate_estimate(dt = tmp_dat, n = n, weights = M[, bt]))
      })
      
      perc_ci <- boot:::perc.ci(blb_reps)
      return(data.table(lower_ci = perc_ci[4],
                        upper_ci = perc_ci[5]))
    })
    
    blb_out <- rbindlist(blb_out)
    blb_out <- blb_out[, .(lower = cumsum(lower_ci)/seq_len(max_subsets),
                           upper = cumsum(upper_ci)/seq_len(max_subsets),
                           reps = rp,
                           subset_num = seq_len(max_subsets))]
    blb_out
  }, cl = cl)
  
  cblb <- rbindlist(cblb)
  saveRDS(cblb, file.path(temp_dir, 'ci_convergence_blb.rds'))
}


truth <- readRDS(file.path(temp_dir, 'true_bootstrap.rds'))
bcblb <- readRDS(file.path(temp_dir, 'ci_convergence.rds'))
blb <- readRDS(file.path(temp_dir, 'ci_convergence_blb.rds'))

plot <- ggplot(bcblb, aes(x = subset_num)) +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  geom_hline(data = truth, aes(yintercept = lower), color = 'red') +
  geom_hline(data = truth, aes(yintercept = upper), color = 'red') +
  facet_wrap(~ reps) + 
  xlab('Number of Subsets') + 
  ylab('ATE') +
  ggtitle('bc-BLB Convergence') +
  ylim(c(0.725, 0.85)) +
  theme_minimal()

ggsave(file.path(image_path, 'bcblb_convergence.png'), plot = plot, width = 8, height = 6, dpi = 400)

plot <- ggplot(blb, aes(x = subset_num)) +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  geom_hline(data = truth, aes(yintercept = lower), color = 'red') +
  geom_hline(data = truth, aes(yintercept = upper), color = 'red') +
  facet_wrap(~ reps) +
  ggtitle('BLB Convergence') +
  ylab('ATE') +
  xlab('Number of Subsets') +
  ylim(c(0.725, 0.85)) +
  theme_minimal()

ggsave(file.path(image_path, 'blb_convergence.png'), plot = plot, width = 8, height = 6, dpi = 400)
