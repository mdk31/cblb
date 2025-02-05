library(data.table)
library(caret)
library(pbapply)
library(e1071)

source('code/helper_funcs.R')

te <- 0.8
sigma <- 1
replications <- 1000
K <- 10
r <- 100

base_nm <- 'svm_dml'
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
                                        B = c(100)))

seq_row <- seq_len(nrow(hyper_grid))


# FULL SIMULATIONS----
if(file.exists(file.path(temp_dir, 'full_coverage.rds'))){
  cblb <- readRDS(file.path(temp_dir, 'full_coverage.rds'))
} else{
  cblb <- lapply(seq_row, function(i){
    grid_val <- hyper_grid[i]
    n <- grid_val$n
    B <- grid_val$B
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
        
        m <- svm(y ~ Tr + X1 + X2, data = train_dat, kernel = 'linear', type = 'eps-regression')
        g <- svm(Tr ~ X1 + X2, data = train_dat, kernel = 'linear', type = 'C-classification', probability = TRUE)

        prop_score <- attr(predict(g, newdata = test_dat, probability = TRUE), 'probabilities')[, 1]
        
        newdata <- data.frame(Tr = 1, X1 = test_dat$X1, X2 = test_dat$X2, Z1 = test_dat$Z1, Z2 = test_dat$Z2)
        m1 <- predict(m, newdata = newdata)
        
        newdata <- data.frame(Tr = 0, X1 = test_dat$X1, X2 = test_dat$X2, Z1 = test_dat$Z1, Z2 = test_dat$Z2)
        m0 <- predict(m, newdata = newdata)
        
        test_dat$prop_score <- prop_score
        test_dat$m1 <- m1
        test_dat$m0 <- m0
        
        test_dat
      })
      crossfit <- rbindlist(crossfit)
      M <- rmultinom(n = B, size = n, prob = rep(1, n))
      
      boot_reps <- sapply(seq_len(B), function(bt){
        phi1 <- M[, bt]*((crossfit$Tr/crossfit$prop_score)*(crossfit$y - crossfit$m1) + crossfit$m1)
        phi0 <- M[, bt]*((1 - crossfit$Tr)/(1 - crossfit$prop_score)*(crossfit$y - crossfit$m0) + crossfit$m0)
        sum(phi1)/n - sum(phi0)/n
        })
      
      boot_ci <- boot:::perc.ci((boot_reps))
      blb_out <- data.table(lower_ci = boot_ci[4],
                             upper_ci = boot_ci[5],
                             estim = mean(boot_reps),
                             se = mean(boot_reps))
      blb_out
    }, cl = 4)

    out <- rbindlist(out)
    out[, `:=`(n = n,
               B = B)]
    out
  })
  cblb <- rbindlist(cblb)
  saveRDS(cblb, file.path(temp_dir, 'full_coverage.rds'))
}

cblb[, .(mean(lower_ci <= te & upper_ci >= te)), by = c('n')]



# CBLB SIMULATIONS----
hyper_grid <- as.data.table(expand.grid(n = c(10000),
                                        subsets = c(5, 10, 15),
                                        B = c(100)))
hyper_grid <- rbindlist(list(hyper_grid,
                             data.table(n = 50000, subsets = 50, B = 100)))
hyper_grid[, `:=`(gamma = calculate_gamma(n, subsets))]
seq_row <- seq_len(nrow(hyper_grid))

if(file.exists(file.path(temp_dir, 'coverage.rds'))){
  cblb <- readRDS(file.path(temp_dir, 'coverage.rds'))
} else{
  cblb <- lapply(seq_row, function(i){
    grid_val <- hyper_grid[i]
    n <- grid_val$n
    gamma <- grid_val$gamma
    b <- floor(n^gamma)
    idx <- seq_len(b)
    B <- grid_val$B
    subsets <- grid_val$subsets
    
    out <- pblapply(seq_len(replications), function(rp){
      set.seed(rp)
      dat <- kangschafer3(n = n, te = te, sigma = sigma, beta_overlap = 0.5)
      partitions <- make_partition(n = n, subsets = subsets, b = b, disjoint = TRUE)
      
      blb_out <- lapply(partitions, function(i){
        tmp_dat <- dat[i]
        folds <- split(idx, sample(rep(1:K, length.out = length(idx))))
        crossfit <- lapply(folds, function(test_idx){
          train_idx <- setdiff(idx, test_idx)
          train_dat <- tmp_dat[train_idx]
          test_dat <- tmp_dat[-train_idx]
          
          m <- svm(y ~ Tr + X1 + X2, data = train_dat, kernel = 'linear', type = 'eps-regression')
          g <- svm(Tr ~ X1 + X2, data = train_dat, kernel = 'linear', type = 'C-classification', probability = TRUE)
          prop_score <- attr(predict(g, newdata = test_dat, probability = TRUE), 'probabilities')[, 1]
          
          newdata <- data.frame(Tr = 1, X1 = test_dat$X1, X2 = test_dat$X2, Z1 = test_dat$Z1, Z2 = test_dat$Z2)
          m1 <- predict(m, newdata = newdata)
          
          newdata <- data.frame(Tr = 0, X1 = test_dat$X1, X2 = test_dat$X2, Z1 = test_dat$Z1, Z2 = test_dat$Z2)
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
               B = B)]
    out
  })
  cblb <- rbindlist(cblb)
  saveRDS(cblb, file.path(temp_dir, 'coverage.rds'))
}

cblb[, .(mean(lower_ci <= te & upper_ci >= te)), by = c('n', 'gamma', 'subsets')]

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
