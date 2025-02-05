
calculate_gamma <- function(n, subsets){
  soln <- 1 - log(subsets)/log(n)
  return(truncate_to_n(soln, 5))
}

causal_blb <- function(data, y_formula, prop_formula, y_method, prop_method, r = 100, 
                       subsets = 5, b = NULL, m_tune = NULL, prop_tune = NULL,
                       y_args = NULL, prop_args = NULL, cores = 1){
  
  if(data.table::is.data.table(data) == FALSE){
    data <- as.data.table(data)
  }
  n <- nrow(data)
  
  if(is.null(b)){
    b <- nrow(data)^0.7
  }
  W <- all.vars(prop_formula)[1]
  Y <- all.vars(y_formula)[1]
  
  data$Tr2 <- factor(data[[W]])
  levels(data$Tr2) <- make.names(levels(data$Tr2))
  prop_formula <- update.formula(as.formula(prop_formula), Tr2 ~ .)
  
  
  m <- do.call(caret::train, c(list(form = as.formula(y_formula),
                                    data = data,
                                    method = y_method,
                                    trControl = caret::trainControl(method = 'none'),
                                    tuneGrid = m_tune), y_args))
  
  g <- do.call(caret::train, c(list(form = as.formula(prop_formula),
                                    data = data,
                                    method = prop_method,
                                    trControl = caret::trainControl(method = 'none', classProbs = TRUE),
                                    tuneGrid = prop_tune), prop_args))
  

  # Calculate full sample
  # Full sample
  prop_score <- predict(g, type = 'prob')[, 1]
  pred_data <- data.table::copy(data)
  pred_data[[W]] <- 1
  m1 <- predict(m, pred_data)
  pred_data[[W]] <- 0
  m0 <- predict(m, pred_data)
  pred_data[[W]] <- data[[W]]
  
  phi1_full <- (pred_data[[W]]/prop_score)*(pred_data[[Y]] - m1) + m1
  phi0_full <- (1 - pred_data[[W]])/(1 - prop_score)*(pred_data[[Y]] - m0) + m0
  
  tau_hat_full <- mean(phi1_full) - mean(phi0_full)
  partitions <- make_partition(n = n, subsets = subsets, b = b, disjoint = FALSE)
  
  blb_out <- parallel::mclapply(partitions, function(i){
    tmp_dat <- data[i]
    M <- rmultinom(n = r, size = n, prob = rep(1, b))
    
    blb_reps <- sapply(seq_len(r), function(bt){
      m <- do.call(caret::train, c(list(form = as.formula(y_formula),
                                        data = tmp_dat,
                                        method = y_method,
                                        trControl = caret::trainControl(method = 'none'),
                                        tuneGrid = m_tune, weights = M[, bt]), y_args))
      
      g <- do.call(caret::train, c(list(form = as.formula(prop_formula),
                                        data = tmp_dat,
                                        method = prop_method,
                                        trControl = caret::trainControl(method = 'none', classProbs = TRUE),
                                        tuneGrid = prop_tune, weights = M[, bt]), prop_args))
      
      prop_score <- predict(g, type = 'prob')[, 1]
      pred_data <- data.table::copy(tmp_dat)
      pred_data[[W]] <- 1
      m1 <- predict(m, pred_data)
      pred_data[[W]] <- 0
      m0 <- predict(m, pred_data)
      pred_data[[W]] <- tmp_dat[[W]]
      
      phi1 <- M[, bt]*((pred_data[[W]]/prop_score)*(pred_data[[Y]] - m1) + m1)
      phi0 <- M[, bt]*((1 - pred_data[[W]])/(1 - prop_score)*(pred_data[[Y]] - m0) + m0)
      sum(phi1)/n - sum(phi0)/n
    })
    
    bias <- mean(blb_reps) - tau_hat_full
    blb_bias_correct <- blb_reps - bias
    perc_ci <- boot:::perc.ci(blb_bias_correct)
    return(data.table::data.table(lower_ci = perc_ci[4],
                                  upper_ci = perc_ci[5]))
  }, mc.cores = cores)
  
  
  blb_out <- data.table::rbindlist(blb_out)
  blb_out <- blb_out[, .(lower_ci = mean(lower_ci),
                         upper_ci = mean(upper_ci))]
  blb_out
}


continuous_treatment_sim <- function(n, sigma = 1, beta_overlap = 0.5){
  X1 <- rnorm(n, 0, 1)
  X2 <- rnorm(n, 0, 1)
  exposure <- X1 + X2 + rnorm(n, 0, 1)
  
  y  <- exposure*0.75 + 1.5*X1 + 1.125*X2 + rnorm(n, 0, sigma)
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
  
  pi2mod <- glm(sq.res ~ X1 + X2, data = data, family = gaussian(link = "log"))
  pi2mod.vals <-  predict(pi2mod, newdata = l.new, type = "response")
  
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
  return(list(pseudo.out = pseudo.out, bandwidth = h.opt))
}

inv_logit <- function(x){
  1/(1 + exp(-x))
}

kangschafer3 <- function(n, te, beta_overlap = 0.5, sigma) {
  X1 <- rnorm(n, 0, 1)
  X2 <- rnorm(n, 0, 1)
  prt <- 1/(1 + exp(-beta_overlap*(X1 + X2)))
  
  Z1 <- exp(X1/2)
  Z2 <- 0.25*(X1*X2)
  
  Tr <- rbinom(n, 1, prt)
  
  Y0  <- X1 + X2 + rnorm(n, 0, sigma)
  Y1  <- Y0 + te
  y   <- Y0*(1-Tr) + Y1*(Tr)
  out <- cbind(y, Tr, X1, X2, Z1, Z2, Y1, Y0)
  out <- data.table::as.data.table(out)
  return(out)
}


make_partition <- function(n, subsets, b, disjoint = FALSE){
  part_idx <- seq(1, n, by = 1)
  if(disjoint){
    # Generate disjoint sets
    # Permute indices
    partition <- sample(part_idx, n)
    totality <- b*subsets
    stopifnot(totality <= n)
    # Exclude rest of sample
    partition <- partition[1:totality]
    partition <- split(partition, f = rep(1:subsets, each = b))
  } else{
    partition <- replicate(subsets, {
      sample(part_idx, size = b, replace = FALSE)
    }, simplify = FALSE)
  }
  partition
}

truncate_to_n <- function(number, n) {
  factor <- 10^n
  truncated <- trunc(number * factor) / factor
  return(truncated)
}


