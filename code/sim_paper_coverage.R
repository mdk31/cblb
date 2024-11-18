library(data.table)
library(pbapply)
library(xgboost)
library(caret)

source('code/helper_funcs.R')

te <- 0.8
sigma <- 1
replications <- 1000
r <- 100

base_nm <- 'sim_paper'
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

hyper_grid <- as.data.table(expand.grid(n = c(50000),
                                        gamma = c(0.8),
                                        subsets = c(5, 10),
                                        prop_form = c('correct', 'wrong'), stringsAsFactors = FALSE))

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
    if(prop_form == 'correct'){
      prop_formula <- c('X1', 'X2')
    } else{
      prop_formula <- c('Z1', 'Z2')
    }
    subsets <- grid_val$subsets

    out <- pblapply(seq_len(replications), function(rp){
      set.seed(rp)
      dat <- kangschafer3(n = n, te = te, sigma = sigma, beta_overlap = 0.5)
      n <- nrow(dat)
      n0 <- sum(dat$Tr == 0)
      n1 <- n - n0

      partitions <- make_partition(n = n, subsets = subsets, b = round(n^gamma), disjoint = FALSE)
      
      blb_out <- lapply(partitions, function(i){
        tmp_dat <- dat[i]
        g_train <- as.matrix(tmp_dat[, prop_formula, with = FALSE])
        # Make weights
        # g <- xgboost(data = g_train, label = tmp_dat$Tr, verbose = 0, nrounds = 15, 
        #              params = list(objective = 'binary:logistic'))
        prop_form <- paste0('Tr ~ ', paste(prop_formula, collapse = ' + '))
        g <- glm(as.formula(prop_form), family = 'binomial', data = tmp_dat)
        prop_score <- predict(g, type = 'response')
        wts <- 1/prop_score*(tmp_dat$Tr == 1) + (1/(1-prop_score))*(tmp_dat$Tr == 0)
        ns <- abs(tapply(wts, tmp_dat$Tr, sum))
        ns <- ns['1']*(tmp_dat$Tr == 1) + (tmp_dat$Tr == 0)*ns['0']
        wts <- wts/ns
        
        ws <- split(wts, tmp_dat$Tr)
        ys <- split(tmp_dat$y, tmp_dat$Tr)
        w0 <- ws[['0']]
        w1 <- ws[['1']]
        
        M1 <- rmultinom(r, n1, prob = w1)
        y1 <- colSums(ys[['1']]*M1)/n1
        M0 <- rmultinom(r, n0, prob = w0)
        y0 <- colSums(ys[['0']]*M0)/n0

        perc_ci <- boot:::perc.ci(y1 - y0)
        return(data.table(lower_ci = perc_ci[4],
                          upper_ci = perc_ci[5],
                          estim = mean(y1 - y0),
                          se = sd(y1 - y0)))
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
               prop_form = prop_form)]
    out
  })
  cblb <- rbindlist(cblb)
  saveRDS(cblb, file.path(temp_dir, 'coverage.rds'))
}

# COVERAGE----
zip <- copy(cblb)
zip[, `:=`(cent = abs(estim - te)/se)]
setnames(zip, old = c('prop_form'), new = c('Propensity Model'))

zip[, `:=`(rank = as.integer(cut(cent, quantile(cent, probs = seq(0, 1, by = 0.01)), include.lowest = TRUE))), 
    by = c('n', 'gamma', 'subsets', 'Propensity Model')]

zip[, `:=`(n = format(n, scientific = FALSE, big.mark = ','),
           covered = fifelse(lower_ci <= te & upper_ci >= te, 'Coverer', 'Non-coverer'))]
zip_labels <- zip[, .(perc_cover = round(mean(covered == 'Coverer'), 3)),
                  by = c('n', 'gamma', 'subsets', 'Propensity Model')]

value_grid <- expand.grid(n = unique(zip$n),
                          gamma = unique(zip$gamma),
                          subsets = unique(zip$subsets), stringsAsFactors = FALSE)

ggdat <- copy(zip)

n_val <- row$n
gamma_val <- row$gamma

nm <- paste0('SIM_zip_plot_n_', n_val, '_gamma_', gamma_val, '.pdf')
title <- bquote(paste(gamma == .(gamma_val), ' and ', n == .(n_val)))
ggsub <- ggdat[n == n_val & gamma == gamma_val]
label_sub <- zip_labels[n == n_val & gamma == gamma_val]
p <- ggplot(ggsub, aes(y = rank)) +
  geom_segment(aes(x = lower_ci, y = rank, xend = upper_ci, yend = rank, color = covered)) +
  facet_grid(~`Propensity Model` ~ subsets, labeller = label_both) +
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


for(row_idx in seq_len(nrow(value_grid))){

}


