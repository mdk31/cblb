library(data.table)
library(pbapply)
library(KernSmooth)
source('code/helper_funcs.R')

# Generate evaluation points
# Get true effects
# Do full method (see where bottleneck is in large n)
# See how well this does in terms of bias, RMSE, etc.
# From there, think how BLB can improve this

continuous_treatment_sim <- function(n, sigma = 1, beta_overlap = 0.5){
  X1 <- rnorm(n, 0, 1)
  X2 <- rnorm(n, 0, 1)
  logit_lambda <- 1/(1 + exp(-beta_overlap*(X1 + X2)))
  
  exposure <- rbeta(n, logit_lambda, 1 - logit_lambda)/20
  
  y  <- exposure*(0.3*X1 - 0.13^2*exposure) + 0.5*X1 + 0.7*X2 + rnorm(n, 0, sigma)
  out <- cbind(y, exposure, X1, X2)
  out <- as.data.frame(out)
  return(out)

}

estimate_models <- function(data, a.vals){
  
  # Construct data for predictions
  la.new <- rbind(cbind(data[, c("X1", "X2")], exposure = data$exposure),
                  cbind(data[rep(1:nrow(data), length(a.vals)), c("X1", "X2")], 
                        exposure  = rep(a.vals, rep(nrow(data), length(a.vals)))))
  
  l.new <- la.new[, -3]
  
  # Fit models
  pimod <- lm(exposure ~ X1 + X2, data = data)
  pimod.vals <- predict(pimod, newdata = l.new)
  sq.res <- (data$exposure - pimod.vals[1:nrow(data)])^2
  
  pi2mod <- lm(sq.res ~ X1 + X2, data = data)
  pi2mod.vals <- predict(pi2mod, newdata = l.new)
  
  mumod <- lm(y ~ exposure + X1 + X2, data = data)
  muhat.vals <- predict(mumod, newdata = la.new, type = "response")
  
  # Compute standardized residuals
  a.std <- (la.new$exposure - pimod.vals) / sqrt(pi2mod.vals)
  
  # Kernel density estimation for pihat 
  approx.fn <- function(x, y, z) { predict(smooth.spline(x, y), x = z)$y }
  pihat.vals <- approx.fn(density(a.std[1:nrow(data)])$x, density(a.std[1:nrow(data)])$y, a.std)
  pihat <- pihat.vals[1:nrow(data)]
  
  pihat.mat <- matrix(pihat.vals[-(1:nrow(data))], nrow = nrow(data), ncol = length(a.vals))
  varpihat <- approx.fn(a.vals, apply(pihat.mat, 2, mean), data$exposure)
  varpihat.mat <- matrix(rep(apply(pihat.mat, 2, mean), nrow(data)), byrow = TRUE, nrow = nrow(data))
  
  muhat <- muhat.vals[1:nrow(data)]
  muhat.mat <- matrix(muhat.vals[-(1:nrow(data))], nrow = nrow(data), ncol = length(a.vals))
  mhat <- approx.fn(a.vals, apply(muhat.mat, 2, mean), data$exposure)
  mhat.mat <- matrix(rep(apply(muhat.mat, 2, mean), nrow(data)), byrow = TRUE, nrow = nrow(data))
  
  # Compute pseudo outcome
  pseudo.out <- (data$y - muhat) / (pihat / varpihat) + mhat
  
  # Select bandwidth via cross-validation
  kern <- function(x) { dnorm(x) }
  w.fn <- function(bw) {
    w.avals <- NULL
    for (a.val in a.vals) {
      a.std <- (data$exposure - a.val) / bw
      kern.std <- kern(a.std) / bw
      w.avals <- c(w.avals, mean(a.std^2 * kern.std) * (kern(0) / bw) /
                     (mean(kern.std) * mean(a.std^2 * kern.std) - mean(a.std * kern.std)^2))
    }
    return(w.avals / nrow(data))
  }
  hatvals <- function(bw) { approx(a.vals, w.fn(bw), xout = data$exposure)$y }
  cts.eff <- function(out, bw) { 
    approx(locpoly(data$exposure, out, bandwidth = bw), xout = data$exposure)$y }
  
  h.opt <- optimize(function(h) {
    hats <- hatvals(h)
    mean(((pseudo.out - cts.eff(pseudo.out, bw = h)) / (1 - hats))^2)
  }, c(0.01, 50), tol = 0.01)$minimum
  
  # Final effect estimate
  # est <- approx(locpoly(data$exposure, pseudo.out, bandwidth = h.opt), xout = a.vals)$y
  return(list(tr_vals = a.vals, pseudo.out = pseudo.out, bandwidth = h.opt))
}

estimate_dose_response_est <- function(data){
  ### INPUT: l is an n*p matrix, a and y are vectors of length n
  ### l = matrix of covariates
  ### a = vector of treatment values
  ### y = vector of observed outcomes
  
  # Set up evaluation points & matrices for predictions
  a.min <- min(a); a.max <- max(a)
  a.vals <- seq(a.min, a.max, length.out = 100)
  la.new <- rbind(cbind(l, a), cbind(l[rep(1:n, length(a.vals)), ], a = rep(a.vals, rep(n, length(a.vals)))))
  l.new <- la.new[, -dim(la.new)[2]]
  
  # Fit linear models instead of SuperLearner
  pimod <- lm(a ~ ., data = as.data.frame(cbind(a, l)))
  pimod.vals <- predict(pimod, newdata = as.data.frame(l.new))
  sq.res <- (a - pimod.vals[1:n])^2
  
  pi2mod <- lm(sq.res ~ ., data = as.data.frame(cbind(sq.res, l)))
  pi2mod.vals <- predict(pi2mod, newdata = as.data.frame(l.new))
  
  mumod <- glm(y ~ ., family = binomial, data = as.data.frame(cbind(y, l, a)))
  muhat.vals <- predict(mumod, newdata = as.data.frame(la.new), type = "response")
  
  # Construct estimated pi/varpi and mu/m values
  approx.fn <- function(x, y, z) {
    predict(smooth.spline(x, y), x = z)$y
  }
  
  a.std <- (la.new$a - pimod.vals) / sqrt(pi2mod.vals)
  pihat.vals <- approx.fn(density(a.std[1:n])$x, density(a.std[1:n])$y, a.std)
  pihat <- pihat.vals[1:n]
  pihat.mat <- matrix(pihat.vals[-(1:n)], nrow = n, ncol = length(a.vals))
  varpihat <- approx.fn(a.vals, apply(pihat.mat, 2, mean), a)
  varpihat.mat <- matrix(rep(apply(pihat.mat, 2, mean), n), byrow = TRUE, nrow = n)
  
  muhat <- muhat.vals[1:n]
  muhat.mat <- matrix(muhat.vals[-(1:n)], nrow = n, ncol = length(a.vals))
  mhat <- approx.fn(a.vals, apply(muhat.mat, 2, mean), a)
  mhat.mat <- matrix(rep(apply(muhat.mat, 2, mean), n), byrow = TRUE, nrow = n)
  
  # Form adjusted/pseudo outcome xi
  pseudo.out <- (y - muhat) / (pihat / varpihat) + mhat
  
  # Leave-one-out cross-validation to select bandwidth
  library(KernSmooth)
  kern <- function(x) { dnorm(x) }
  
  w.fn <- function(bw) {
    w.avals <- NULL
    for (a.val in a.vals) {
      a.std <- (a - a.val) / bw
      kern.std <- kern(a.std) / bw
      w.avals <- c(w.avals, mean(a.std^2 * kern.std) * (kern(0) / bw) /
                     (mean(kern.std) * mean(a.std^2 * kern.std) - mean(a.std * kern.std)^2))
    }
    return(w.avals / n)
  }
  
  hatvals <- function(bw) { approx(a.vals, w.fn(bw), xout = a)$y }
  cts.eff <- function(out, bw) { approx(locpoly(a, out, bw), xout = a)$y }
  
  # Optimize bandwidth selection
  h.opt <- optimize(function(h) {
    hats <- hatvals(h)
    mean(((pseudo.out - cts.eff(pseudo.out, bw = h)) / (1 - hats))^2)
  }, c(0.01, 50), tol = 0.01)$minimum
  
  # Estimate effect curve with optimal bandwidth
  est <- approx(locpoly(a, pseudo.out, bandwidth = h.opt), xout = a.vals)$y
  
}

dat <- continuous_treatment_sim(n = 1000)
# estimate_models(dat)
# estimate_dose_response_est(dat)

te <- 0.8
sigma <- 1
replications <- 1000
r <- 100

base_nm <- 'dose_response'
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
                                        B = c(100)))

seq_row <- seq_len(nrow(hyper_grid))


# FULL SIMULATIONS----
if(file.exists(file.path(temp_dir, 'full_coverage.rds'))){
  cblb <- readRDS(file.path(temp_dir, 'full_coverage.rds'))
} else{
  cblb <- lapply(seq_row, function(i){
    grid_val <- hyper_grid[i]
    n <- grid_val$n
    part_idx <- seq_len(b)
    B <- grid_val$B
    subsets <- grid_val$subsets

    out <- pblapply(seq_len(replications), function(rp){
      set.seed(rp)
      dat <- kangschafer3(n = n, te = te, sigma = sigma, beta_overlap = 0.5)
      model_out <- estimate_models()
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
      blb_out <- blb_out[, .(lower_ci = boot_ci[4],
                             upper_ci = boot_ci[5],
                             estim = mean(boot_reps),
                             se = mean(boot_reps))]
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


# CBLB SIMULATIONS----
hyper_grid <- as.data.table(expand.grid(n = c(10000),
                                        subsets = c(10),
                                        B = c(100)))
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
    part_idx <- seq_len(b)
    B <- grid_val$B
    subsets <- grid_val$subsets
    
    out <- pblapply(seq_len(replications), function(rp){
      set.seed(rp)
      dat <- continuous_treatment_sim(n = n, sigma = sigma, beta_overlap = 0.5)
      a.min <- min(dat$exposure)
      a.max <- max(dat$exposure)
      a.vals <- seq(a.min, a.max, length.out = 100)
      true_dr <- a.vals*(0 - 0.13^2*a.vals)
      partitions <- make_partition(n = n, subsets = subsets, b = b, disjoint = FALSE)
      
      blb_out <- lapply(partitions, function(i){
        tmp_dat <- dat[i, ]
        idx <- seq_len(nrow(tmp_dat))
        model_output <- estimate_models(tmp_dat)
        blb_reps <- replicate(B, {
          boot_idx <- sample(idx, size = n, replace = TRUE)
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
      blb_out[, `:=`(true_dr = true_dr)]
      blb_out <- blb_out[, .(coverage = sum(lower_ci <= true_dr & upper_ci >= true_dr),
                  rep = rp,
                  total_vals = length(tr_vals))]
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

# cblb[, `:=`(xval = rep(1:100, 1000),
#             rp = rep(1:100, each = 1000))]
# cblb[, .(coverage = mean(lower_ci <= true_dr & upper_ci >= true_dr)), by = c('n', 'gamma', 'subsets', 'B')]


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
