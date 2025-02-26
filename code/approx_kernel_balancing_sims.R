library(data.table) 
library(pbapply)
library(Rcpp)

source('code/helper_funcs.R')
source('code/osqp-kernel-sbw.R')
Rcpp::sourceCpp("code/RBF_kernel_C_parallel.cpp")

te <- 0.8
sigma <- 1
replications <- 1000
degree1 <- 1
degree2 <- 1
k1 <- "poly"
k2 <- "poly"
operator <- "single" # here i'm considering a simple linear kernel just for comparison/true model
penal <- log(2) # keep it between 0.1 and log(k) where k is the number of features

base_nm <- 'approx_kernel_balancing'
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

hyper_grid <- as.data.table(expand.grid(n = c(10000),
                                        B = c(100),
                                        prop_form = c('correct'),
                                        out_form = c('correct')))

seq_row <- seq_len(nrow(hyper_grid))

# PARAMETER BIAS----
if(file.exists(file.path(temp_dir, 'param_values.rds'))){
  cblb <- readRDS(file.path(temp_dir, 'param_values.rds'))
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
      output <- osqp_kernel_sbw(X = as.matrix(dat[, c('X1', 'X2')]),
                                A = dat$Tr,
                                Y = dat$y,
                                delta.v=1e-4,
                                kernel.approximation = TRUE,
                                c = 100)
      data.table(att = mean(dat$y[dat$Tr == 1]) - output[[1]]$y_hat)
    }, cl = 5)
    
    out <- rbindlist(out)
    out[, `:=`(n = n,
               prop_form = prop_form,
               out_form = out_form,
               B = B)]
    out
  })
  cblb <- rbindlist(cblb)
  saveRDS(cblb, file.path(temp_dir, 'param_values.rds'))
}

summary(cblb)

# CORRECT PROCEDURE W/ NYDSTROM APPROX----
if(file.exists(file.path(temp_dir, 'bootstrap_approx.rds'))){
  cblb <- readRDS(file.path(temp_dir, 'bootstrap_approx.rds'))
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
      n0 <- sum(1-dat$Tr)
      idx <- seq_len(n)
      
      blb_reps <- sapply(seq_len(B), function(bt){
        boot_dat <- dat[sample(idx, n, replace = TRUE)]
        output <- osqp_kernel_sbw(X = as.matrix(boot_dat[, c('X1', 'X2')]),
                                  A = boot_dat$Tr,
                                  Y = boot_dat$y,
                                  delta.v=1e-4,
                                  kernel.approximation = TRUE,
                                  c = 100)
        mean(dat$y[dat$Tr == 1]) - output[[1]]$y_hat
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
               prop_form = prop_form,
               out_form = out_form,
               B = B)]
    out
  })
  cblb <- rbindlist(cblb)
  saveRDS(cblb, file.path(temp_dir, 'bootstrap_approx.rds'))
}

cblb[, .(mean(lower_ci <= te & upper_ci >= te)), by = c('n', 'type')]

# FULL PROCEDURE W/ NYDSTROM APPROX----
if(file.exists(file.path(temp_dir, 'full_bootstrap_approx.rds'))){
  cblb <- readRDS(file.path(temp_dir, 'full_bootstrap_approx.rds'))
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
      n0 <- sum(1-dat$Tr)
      n1 <- n - n0
      idx <- seq_len(n)
      output <- osqp_kernel_sbw(X = as.matrix(dat[, c('X1', 'X2')]),
                                A = dat$Tr,
                                Y = dat$y,
                                delta.v=1e-4,
                                kernel.approximation = TRUE,
                                c = 100)
      M <- rmultinom(n = B, size = n, prob = rep(1, n))
      
      blb_reps <- sapply(seq_len(B), function(bt){
        weights <- M[, bt]
        treat <- split(weights, dat$Tr)
        sum(treat[['1']]*dat$y[dat$Tr == 1])/n1 - sum(treat[['0']]*output[[1]]$w*dat$y[dat$Tr == 0])
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
               prop_form = prop_form,
               out_form = out_form,
               B = B)]
    out
  })
  cblb <- rbindlist(cblb)
  saveRDS(cblb, file.path(temp_dir, 'full_bootstrap_approx.rds'))
}

cblb[, .(mean(lower_ci <= te & upper_ci >= te)), by = c('n', 'type')]


# DISJOINT----
hyper_grid <- rbindlist(list(data.table(n = 10000, B = 100, prop_form = 'correct', out_form = 'correct', subsets = 10, gamma = 0.75),
                             data.table(n = 10000, B = 100, prop_form = 'correct', out_form = 'correct', subsets = 15, gamma = 0.70597),
                             data.table(n = 10000, B = 100, prop_form = 'correct', out_form = 'correct', subsets = 20, gamma = 0.67474)))

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
      n1 <- sum(dat$Tr == 1)
      n0 <- n - n1
      partitions <- make_partition(n = n, subsets = subsets, b = b, disjoint = TRUE)
      
      blb_out <- lapply(partitions, function(i){
        tmp_dat <- dat[i]
        b1 <- sum(tmp_dat$Tr == 1)
        b0 <- b - b1
        output <- osqp_kernel_sbw(X = as.matrix(tmp_dat[, c('X1', 'X2')]),
                                  A = tmp_dat$Tr,
                                  Y = tmp_dat$y,
                                  delta.v=1e-4,
                                  kernel.approximation = TRUE,
                                  c = 100)
        
        # output <- osqp_kernel_sbwN0(X = as.matrix(tmp_dat[, c('X1', 'X2')]),
        #                           A = tmp_dat$Tr,
        #                           Y = tmp_dat$y,
        #                           delta.v=1e-4,
        #                           kernel.approximation = TRUE,
        #                           c = 100)

        M <- rmultinom(n = B, size = n, prob = rep(1, b))
        
        blb_reps <- sapply(seq_len(B), function(bt){
          weights <- M[, bt]
          treat <- split(weights, tmp_dat$Tr)
          sum(treat[['1']]*tmp_dat$y[tmp_dat$Tr == 1])/n1 - sum(treat[['0']]*output[[1]]$w*tmp_dat$y[tmp_dat$Tr == 0]*b0)/n0
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
               prop_form = prop_form,
               out_form = out_form,
               B = B)]
    out
  })
  cblb <- rbindlist(cblb)
  saveRDS(cblb, file.path(temp_dir, 'disjoint_subsets.rds'))
}

cblb[, .(mean(lower_ci <= te & upper_ci >= te)), by = c('n', 'gamma', 'subsets', 'out_form', 'prop_form', 'B')]




