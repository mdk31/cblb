
causal_blb <- function(data, y_formula, prop_formula, y_method, prop_method, r = 100, 
                       subsets = 5, b = NULL, m_tune = NULL, prop_tune = NULL,
                       y_args = NULL, prop_args = NULL, cores = 1){
  
  if(data.table::is.data.table(data) == FALSE){
    data <- as.data.table(data)
  }
  
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
