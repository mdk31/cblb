library(data.table)
library(caret)
library(pbapply)
library(xgboost)

source('code/helper_funcs.R')

te <- 0.8
sigma <- 1
replications <- 1000
r <- 100
K <- 10
n <- 1e6

base_nm <- 'verylargen_crossfit'
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

hyper_grid <- as.data.table(expand.grid(n = c(n),
                                        B = c(100),
                                        prop_form = c('correct'),
                                        out_form = c('correct')))

seq_row <- seq_len(nrow(hyper_grid))

# CORRECT PROCEDURE----
if(file.exists(file.path(temp_dir, 'bootstrap.rds'))){
  cblb <- readRDS(file.path(temp_dir, 'bootstrap.rds'))
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
      idx <- seq_len(nrow(dat))
      folds <- split(idx, sample(rep(1:K, length.out = length(idx))))
      crossfit <- lapply(folds, function(test_idx){
        train_idx <- setdiff(idx, test_idx)
        train_dat <- dat[train_idx]
        test_dat <- dat[-train_idx]
        
        m_train <- as.matrix(train_dat[, out_formula, with = FALSE])
        g_train <- as.matrix(train_dat[, prop_formula, with = FALSE])
        
        g_test <- as.matrix(test_dat[, prop_formula, with = FALSE])
        
        m <- xgboost(data = m_train, label = train_dat$y, verbose = 0, nrounds = 15,
                     params = list(objective = 'reg:squarederror'))
        g <- xgboost(data = g_train, label = train_dat$Tr, verbose = 0, nrounds = 15, 
                     params = list(objective = 'binary:logistic'))
        
        prop_score <- predict(g, newdata = g_test)
        
        newdata <- data.frame(Tr = 1, X1 = test_dat$X1, X2 = test_dat$X2, Z1 = test_dat$Z1, Z2 = test_dat$Z2)
        newdata <- as.matrix(newdata[, out_formula])
        m1 <- predict(m, newdata = newdata)
        
        newdata <- data.frame(Tr = 0, X1 = test_dat$X1, X2 = test_dat$X2, Z1 = test_dat$Z1, Z2 = test_dat$Z2)
        newdata <- as.matrix(newdata[, out_formula])
        m0 <- predict(m, newdata = newdata)
        
        test_dat$prop_score <- prop_score
        test_dat$m1 <- m1
        test_dat$m0 <- m0
        
        test_dat
      })
      crossfit <- rbindlist(crossfit)
      
      M <- rmultinom(n = B, size = n, prob = rep(1, n))
      
      blb_reps <- sapply(seq_len(B), function(bt){
        phi1 <- M[, bt]*((crossfit$Tr/crossfit$prop_score)*(crossfit$y - crossfit$m1) + crossfit$m1)
        phi0 <- M[, bt]*((1 - crossfit$Tr)/(1 - crossfit$prop_score)*(crossfit$y - crossfit$m0) + crossfit$m0)
        sum(phi1)/n - sum(phi0)/n
        
      })
        
      perc_ci <- boot:::perc.ci(blb_reps)
      return(data.table(lower_ci = perc_ci[4],
                        upper_ci = perc_ci[5],
                        estim = mean(blb_reps),
                        se = sd(blb_reps)))
      }, cl = 4)
    
    out <- rbindlist(out)
    out[, `:=`(n = n,
               prop_form = prop_form,
               out_form = out_form,
               B = B)]
    out
  })
  cblb <- rbindlist(cblb)
  saveRDS(cblb, file.path(temp_dir, 'bootstrap.rds'))
}

cblb[, .(mean(lower_ci <= te & upper_ci >= te)), by = c('n')]


# CROSSFIT----
hyper_grid <- rbindlist(list(data.table(n = n, B = 100, prop_form = 'correct', out_form = 'correct', subsets = 3, gamma = 0.91),
                        data.table(n = n, B = 100, prop_form = 'correct', out_form = 'correct', subsets = 15, gamma = 0.8),
                        data.table(n = n, B = 100, prop_form = 'correct', out_form = 'correct', subsets = 12, gamma = 0.82),
                        data.table(n = n, B = 100, prop_form = 'correct', out_form = 'correct', subsets = 63, gamma = 0.7)))

seq_row <- seq_len(nrow(hyper_grid))

if(file.exists(file.path(temp_dir, 'crossfit_disjoint_subsets.rds'))){
  cblb <- readRDS(file.path(temp_dir, 'crossfit_disjoint_subsets.rds'))
} else{
  cblb <- lapply(seq_row, function(i){
    grid_val <- hyper_grid[i]
    n <- grid_val$n
    gamma <- grid_val$gamma
    b <- round(n^gamma)
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
      partitions <- make_partition(n = n, subsets = subsets, b = round(n^gamma), disjoint = TRUE)
      
      blb_out <- lapply(partitions, function(i){
        tmp_dat <- dat[i]
        idx <- seq_len(nrow(tmp_dat))
        folds <- split(idx, sample(rep(1:K, length.out = length(idx))))
        crossfit <- lapply(folds, function(test_idx){
          train_idx <- setdiff(idx, test_idx)
          train_dat <- tmp_dat[train_idx]
          test_dat <- tmp_dat[-train_idx]
          
          m_train <- as.matrix(train_dat[, out_formula, with = FALSE])
          g_train <- as.matrix(train_dat[, prop_formula, with = FALSE])
          
          g_test <- as.matrix(test_dat[, prop_formula, with = FALSE])
          
          m <- xgboost(data = m_train, label = train_dat$y, verbose = 0, nrounds = 15,
                       params = list(objective = 'reg:squarederror'))
          g <- xgboost(data = g_train, label = train_dat$Tr, verbose = 0, nrounds = 15, 
                       params = list(objective = 'binary:logistic'))
          
          prop_score <- predict(g, newdata = g_test)
          
          newdata <- data.frame(Tr = 1, X1 = test_dat$X1, X2 = test_dat$X2, Z1 = test_dat$Z1, Z2 = test_dat$Z2)
          newdata <- as.matrix(newdata[, out_formula])
          m1 <- predict(m, newdata = newdata)
          
          newdata <- data.frame(Tr = 0, X1 = test_dat$X1, X2 = test_dat$X2, Z1 = test_dat$Z1, Z2 = test_dat$Z2)
          newdata <- as.matrix(newdata[, out_formula])
          m0 <- predict(m, newdata = newdata)
          
          test_dat$prop_score <- prop_score
          test_dat$m1 <- m1
          test_dat$m0 <- m0
          
          test_dat
        })
        crossfit <- rbindlist(crossfit)
        
        M <- rmultinom(n = B, size = n, prob = rep(1, b))
        
        blb_reps <- sapply(seq_len(B), function(bt){
          phi1 <- M[, bt]*((crossfit$Tr/crossfit$prop_score)*(crossfit$y - crossfit$m1) + crossfit$m1)
          phi0 <- M[, bt]*((1 - crossfit$Tr)/(1 - crossfit$prop_score)*(crossfit$y - crossfit$m0) + crossfit$m0)
          sum(phi1)/n - sum(phi0)/n
          
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
    }, cl = 4)
    
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
  saveRDS(cblb, file.path(temp_dir, 'crossfit_disjoint_subsets.rds'))
}

cblb[, .(mean(lower_ci <= te & upper_ci >= te)), by = c('n', 'gamma', 'subsets', 'out_form', 'prop_form', 'B')]



