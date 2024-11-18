library(data.table)
library(caret)

source('code/helper_funcs.R')

te <- 0.8
sigma <- 1
replications <- 1000
r <- 100

base_nm <- 'coverage_sims'
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
                                        B = c(25, 50, 100),
                                        prop_form = c('correct', 'wrong'),
                                        out_form = c('correct', 'wrong')))

seq_row <- seq_len(nrow(hyper_grid))


# CBLB SIMULATIONS----
if(file.exists(file.path(temp_dir, 'coverage.rds'))){
  cblb <- readRDS(file.path(temp_dir, 'coverage.rds'))
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
      # Full sample
      m_train <- as.matrix(dat[, out_formula, with = FALSE])
      g_train <- as.matrix(dat[, prop_formula, with = FALSE])

      m <- xgboost(data = m_train, label = dat$y, verbose = 0, nrounds = 15,
                   params = list(objective = 'reg:squarederror'))
      g <- xgboost(data = g_train, label = dat$Tr, verbose = 0, nrounds = 15, 
                   params = list(objective = 'binary:logistic'))
      
      prop_score <- predict(g, newdata = g_train)
      
      newdata <- data.frame(Tr = 1, X1 = dat$X1, X2 = dat$X2, Z1 = dat$Z1, Z2 = dat$Z2)
      newdata <- as.matrix(newdata[, out_formula])
      m1 <- predict(m, newdata = newdata)
      
      newdata <- data.frame(Tr = 0, X1 = dat$X1, X2 = dat$X2, Z1 = dat$Z1, Z2 = dat$Z2)
      newdata <- as.matrix(newdata[, out_formula])
      m0 <- predict(m, newdata = newdata)
      
      full_dat <- copy(dat)
      full_dat$prop_score <- prop_score
      full_dat$m1 <- m1
      full_dat$m0 <- m0
      
      phi1_full <- (full_dat$Tr/full_dat$prop_score)*(full_dat$y - full_dat$m1) + full_dat$m1
      phi0_full <- (1 - full_dat$Tr)/(1 - full_dat$prop_score)*(full_dat$y - full_dat$m0) + full_dat$m0
      tau_hat_full <- mean(phi1_full) - mean(phi0_full)
      
      partitions <- make_partition(n = n, subsets = subsets, b = round(n^gamma), disjoint = FALSE)
      
      blb_out <- lapply(partitions, function(i){
        tmp_dat <- dat[i]
        m_train <- as.matrix(tmp_dat[, out_formula, with = FALSE])
        g_train <- as.matrix(tmp_dat[, prop_formula, with = FALSE])
        
        M <- rmultinom(n = B, size = n, prob = rep(1, b))
        
        blb_reps <- sapply(seq_len(B), function(bt){
          m <- xgboost(data = m_train, label = tmp_dat$y, verbose = 0, nrounds = 15,
                       params = list(objective = 'reg:squarederror'), weight = M[, bt])
          g <- xgboost(data = g_train, label = tmp_dat$Tr, verbose = 0, nrounds = 15, 
                       params = list(objective = 'binary:logistic'), weight = M[, bt])
          
          prop_score <- predict(g, newdata = g_train)
          
          newdata <- data.frame(Tr = 1, X1 = tmp_dat$X1, X2 = tmp_dat$X2, Z1 = tmp_dat$Z1, Z2 = tmp_dat$Z2)
          newdata <- as.matrix(newdata[, out_formula])
          m1 <- predict(m, newdata = newdata)
          
          newdata <- data.frame(Tr = 0, X1 = tmp_dat$X1, X2 = tmp_dat$X2, Z1 = tmp_dat$Z1, Z2 = tmp_dat$Z2)
          newdata <- as.matrix(newdata[, out_formula])
          m0 <- predict(m, newdata = newdata)
          
          tmp_dat$prop_score <- prop_score
          tmp_dat$m1 <- m1
          tmp_dat$m0 <- m0
          
          phi1 <- M[, bt]*((tmp_dat$Tr/tmp_dat$prop_score)*(tmp_dat$y - tmp_dat$m1) + tmp_dat$m1)
          phi0 <- M[, bt]*((1 - tmp_dat$Tr)/(1 - tmp_dat$prop_score)*(tmp_dat$y - tmp_dat$m0) + tmp_dat$m0)
          sum(phi1)/n - sum(phi0)/n
        })
        
        bias <- mean(blb_reps) - tau_hat_full
        blb_bias_correct <- blb_reps - bias
        perc_ci <- boot:::perc.ci(blb_bias_correct)
        return(data.table(lower_ci = perc_ci[4],
                          upper_ci = perc_ci[5],
                          estim = mean(blb_bias_correct),
                          se = sd(blb_bias_correct)))
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
  saveRDS(cblb, file.path(temp_dir, 'coverage.rds'))
}

# COVERAGE----
zip <- copy(cblb)
zip[, `:=`(cent = abs(estim - te)/se)]
setnames(zip, old = c('prop_form', 'out_form'), new = c('Propensity Model', 'Outcome Model'))

zip[, `:=`(rank = as.integer(cut(cent, quantile(cent, probs = seq(0, 1, by = 0.01)), include.lowest = TRUE))), 
    by = c('n', 'gamma', 'subsets', 'Propensity Model', 'Outcome Model', 'B')]

zip[, `:=`(n = format(n, scientific = FALSE, big.mark = ','),
           covered = fifelse(lower_ci <= te & upper_ci >= te, 'Coverer', 'Non-coverer'))]
zip_labels <- zip[, .(perc_cover = round(mean(covered == 'Coverer'), 3)),
                  by = c('n', 'gamma', 'subsets', 'Propensity Model', 'Outcome Model', 'B')]

value_grid <- expand.grid(n = unique(zip$n),
                          gamma = unique(zip$gamma),
                          subsets = unique(zip$subsets), stringsAsFactors = FALSE)

ggdat <- copy(zip)
# Keep B = 100 for simulations

for(row_idx in seq_len(nrow(value_grid))){
  row <- value_grid[row_idx, ]
  n_val <- row$n
  gamma_val <- row$gamma
  subsets_val <- row$subsets
  nm <- paste0('zip_plot_n_', n_val, '_subset_', subsets_val, '_gamma_', gamma_val, '.pdf')
  title <- bquote(paste(s == .(subsets_val), ' and ', gamma == .(gamma_val), ' and ', n == .(n_val)))
  ggsub <- ggdat[n == n_val & subsets == subsets_val & gamma == gamma_val & B == 100]
  label_sub <- zip_labels[n == n_val & subsets == subsets_val & gamma == gamma_val & B == 100]
  p <- ggplot(ggsub, aes(y = rank)) +
    geom_segment(aes(x = lower_ci, y = rank, xend = upper_ci, yend = rank, color = covered)) +
    facet_grid(`Propensity Model` ~ `Outcome Model`, labeller = label_both) +
    geom_vline(aes(xintercept = te), color = 'yellow', size = 0.5, linetype = 'dashed') +
    ylab('Fractional Centile of |z|') +
    xlab('95% Confidence Intervals') +
    theme_bw() +
    scale_y_continuous(breaks = c(5, 50, 95)) +
    scale_color_discrete(name = "Coverage") +
    geom_text(x = 0.75, y = 50, aes(label = perc_cover), data = label_sub, size = 4) +
    ggtitle(title)
  
  print(p)
  ggsave(file.path(image_path, nm), height = 9, width = 7)
}

# B robustness checks
ggdat <- copy(zip)

for(row_idx in seq_len(nrow(value_grid))){
  row <- value_grid[row_idx, ]
  n_val <- row$n
  gamma_val <- row$gamma
  subsets_val <- row$subsets
  nm <- paste0('zip_plot_n_', n_val, '_subset_', subsets_val, '_gamma_', gamma_val, '_B_25.pdf')
  title <- bquote(paste(s == .(subsets_val), ' and ', gamma == .(gamma_val), ' and ', n == .(n_val)))
  ggsub <- ggdat[n == n_val & subsets == subsets_val & gamma == gamma_val & B == 25]
  label_sub <- zip_labels[n == n_val & subsets == subsets_val & gamma == gamma_val & B == 25]
  p <- ggplot(ggsub, aes(y = rank)) +
    geom_segment(aes(x = lower_ci, y = rank, xend = upper_ci, yend = rank, color = covered)) +
    facet_grid(`Propensity Model` ~ `Outcome Model`, labeller = label_both) +
    geom_vline(aes(xintercept = te), color = 'yellow', size = 0.5, linetype = 'dashed') +
    ylab('Fractional Centile of |z|') +
    xlab('95% Confidence Intervals') +
    theme_bw() +
    scale_y_continuous(breaks = c(5, 50, 95)) +
    scale_color_discrete(name = "Coverage") +
    geom_text(x = 0.75, y = 50, aes(label = perc_cover), data = label_sub, size = 4) +
    ggtitle(title)
  
  print(p)
  ggsave(file.path(image_path, nm), height = 9, width = 7)
}
