library(data.table)
library(ggplot2)

source('code/helper_funcs.R')

te <- 0.8
sigma <- 1
replications <- 25
r <- 100

image_path <- 'images'
dat_path <- 'data'

base_nm <- 'timing_full'

temp_dir <- file.path(dat_path, paste0(base_nm, '_tmp'))
full <- as.data.table(readRDS(file.path(temp_dir, 'timing.rds')))
full[, `:=`(type = 'regular', cores = 1, subsets = 1, gamma = 1)]

base_nm <- 'timing_blb'
temp_dir <- file.path(dat_path, paste0(base_nm, '_tmp'))
blb <- as.data.table(readRDS(file.path(temp_dir, 'timing.rds')))
# blb[, `:=`(type = paste0('Cores: ', cores, '; Subsets = ', subsets, '; gamma = ', gamma))]
# blb[, `:=`(cores = NULL, subsets = NULL, gamma = NULL)]
blb[, `:=`(type = 'blb')]

out <- rbindlist(list(blb, full), use.names = TRUE)
out[, `:=`(cores = factor(cores))]

ns <- unique(out$n)
for(i in ns){
  
  tmp <- out[n == i & prop_form == 'correct' & out_form == 'correct']
  reg <- grepl('regular', tmp$type, ignore.case=TRUE)
  full <- mean(tmp[reg][['time_elapsed']])
  tmp <- tmp[!reg]
  
  p <- ggplot(tmp, aes(x = cores, y = time_elapsed)) +
    geom_boxplot() +
    facet_grid(subsets ~ gamma, labeller = label_both) +
    theme_minimal() +
    geom_hline(yintercept = full, color='red') +
    theme(axis.text.x = element_text(angle = 90)) +
    ylab('Timing (seconds)')
  
  print(p)
  
  ggsave(file.path(image_path, paste0('timing_comparison_n', i, '.pdf')), height = 9, width = 7)
}

