library(data.table)
library(pbapply)
library(CDML)
library(xgboost)

source('code/helper_funcs.R')

te <- 0.8
sigma <- 1
replications <- 25
r <- 100

base_nm <- 'timing_cdml'
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


hyper_grid <- as.data.table(expand.grid(n = c(10000, 50000),
                                        gamma = c(0.7, 0.8),
                                        subsets = c(2, 5, 10),
                                        prop_form = c('correct', 'wrong'),
                                        out_form = c('correct', 'wrong'),
                                        cores = c(1)))

seq_row <- seq_len(nrow(hyper_grid))

# CBLB SIMULATIONS----
if(file.exists(file.path(temp_dir, 'timing.rds'))){
  cblb <- readRDS(file.path(temp_dir, 'timing.rds'))
} else{
  cblb <- lapply(seq_row, function(i){
    grid_val <- hyper_grid[i]
    n <- grid_val$n
    gamma <- grid_val$gamma
    b <- round(n^gamma)
    part_idx <- seq_len(b)
    prop_form <- grid_val$prop_form
    out_form <- grid_val$out_form
    cores <- grid_val$cores
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
    
    out <- lapply(seq_len(replications), function(rp){
      set.seed(rp)
      dat <- kangschafer3(n = n, te = te, sigma = sigma, beta_overlap = 0.5)
      time <- system.time({
        m_train <- as.matrix(dat[, out_formula, with = FALSE])
        g_train <- as.matrix(dat[, prop_formula, with = FALSE])
        # Full sample
        m <- xgboost(data = m_train, label = dat$y, verbose = 0, nrounds = 15,
                     params = list(objective = 'reg:squarederror'))
        g <- xgboost(data = g_train, label = dat$Tr, verbose = 0, nrounds = 15, 
                     params = list(objective = 'binary:logistic'))
        
        newdata <- data.frame(Tr = 1, X1 = dat$X1, X2 = dat$X2, Z1 = dat$Z1, Z2 = dat$Z2)
        newdata <- as.matrix(newdata[, out_formula])
        m1 <- predict(m, newdata = newdata)

        newdata <- data.frame(Tr = 0, X1 = dat$X1, X2 = dat$X2, Z1 = dat$Z1, Z2 = dat$Z2)
        newdata <- as.matrix(newdata[, out_formula])
        m0 <- predict(m, newdata = newdata)
        
        prop_score <- predict(g, newdata = g_train)

        tau_hat_full <- estimate_cdml(A = dat$Tr, Y = dat$y, mu1 = m1, mu0 = m0, pi0 = 1 - prop_score, 
                      pi1 = prop_score, functional = CDML:::functionals_info[['ATE']][['fun']],
                      representer = CDML:::functionals_info[['ATE']][['rep']])

        partitions <- make_partition(n = n, subsets = subsets, b = round(n^gamma), disjoint = FALSE)
        
        blb_out <- parallel::mclapply(partitions, function(i){
          tmp_dat <- dat[i]
          m_train <- as.matrix(tmp_dat[, out_formula, with = FALSE])
          g_train <- as.matrix(tmp_dat[, prop_formula, with = FALSE])
          M <- rmultinom(n = r, size = n, prob = rep(1, b))
          
          blb_reps <- sapply(seq_len(r), function(bt){
            # Full sample
            m <- xgboost(data = m_train, label = tmp_dat$y, verbose = 0, nrounds = 15,
                         params = list(objective = 'reg:squarederror'), weight = M[, bt])
            g <- xgboost(data = g_train, label = tmp_dat$Tr, verbose = 0, nrounds = 15, 
                         params = list(objective = 'binary:logistic'), weight = M[, bt])
            
            newdata <- data.frame(Tr = 1, X1 = tmp_dat$X1, X2 = tmp_dat$X2, Z1 = tmp_dat$Z1, Z2 = tmp_dat$Z2)
            newdata <- as.matrix(newdata[, out_formula])
            m1 <- predict(m, newdata = newdata)
            
            newdata <- data.frame(Tr = 0, X1 = tmp_dat$X1, X2 = tmp_dat$X2, Z1 = tmp_dat$Z1, Z2 = tmp_dat$Z2)
            newdata <- as.matrix(newdata[, out_formula])
            m0 <- predict(m, newdata = newdata)
            
            prop_score <- predict(g, newdata = g_train)
            estimate_cdml(A = tmp_dat$Tr, Y = tmp_dat$y, mu1 = m1, mu0 = m0, 
                          pi0 = 1 - prop_score, pi1 = prop_score, 
                          functional = CDML:::functionals_info[['ATE']][['fun']],
                          representer = CDML:::functionals_info[['ATE']][['rep']],
                          weights = M[, bt])
            
            
          })
          
          bias <- mean(blb_reps) - tau_hat_full
          blb_bias_correct <- blb_reps - bias
          blb_ci <- boot:::perc.ci(blb_bias_correct)
        }, mc.cores = cores)
        # Full sample
      })
      data.table(time_elapsed = time['elapsed'])
    })
    
    out <- rbindlist(out)
    out[, `:=`(n = n,
               gamma = gamma,
               subsets = subsets,
               prop_form = prop_form,
               out_form = out_form,
               cores = cores)]
    out
  })
  cblb <- rbindlist(cblb)
  saveRDS(cblb, file.path(temp_dir, 'timing.rds'))
}

times <- cblb[, .(med_time = median(time_elapsed)), by = c('out_form', 'prop_form', 'subsets', 'cores', 'n', 'gamma')]
setorder(times, -med_time)
times
