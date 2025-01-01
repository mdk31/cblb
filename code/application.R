library(arrow)
library(data.table)
library(ggplot2)
library(dplyr)
library(knitr)
library(kableExtra)

source('code/helper_funcs.R')

base_dir <- '/Users/macbookair/Documents/Projects/cblb/Bootstrap'

df <- read_parquet(file.path(base_dir, 'Data/environmental_simulated_data.parquet'))
dim(df)
colnames(df)

# Extract treatment variable
treatment <- df$pm25

# Create confounders dataframe excluding specific columns
confounders <- df[, !names(df) %in% c("zip", "year", "pm25", "sim_outcome")]

# Standardize numeric columns
confounders <- as.data.frame(lapply(confounders, function(x) {
  if(is.numeric(x)) {
    return((x - mean(x)) / sd(x))
  }
  return(x)
}))

# Define air quality effect function
air_quality_effect <- function(x) {
  -0.05 * pmax(10 - x, 0)^2 + 0.3 * x
}


# # Create plot using ggplot2
# ggplot(data.frame(x = x, y = 10 * (y/max(y) - 1)), aes(x = x, y = y)) +
#   geom_line() +
#   labs(
#     x = "PM 2.5 standard (lower mean higher air quality)",
#     y = "% Reduction in deaths",
#     title = "True causal effect of changes in air quality standards"
#   ) +
#   theme_minimal()

for (c in names(confounders)) {
  corr_treatment_c <- cor(treatment, confounders[[c]])
  cat(sprintf("Correlation between treatment and %s: %f\n", c, corr_treatment_c))
}


# Define mortality coefficients
mortality_coeffs <- c(
  mean_bmi = 0.02,
  smoke_rate = 0.02,
  hispanic = 0.02,
  pct_blk = 5,
  medhouseholdincome = -0.05,
  medianhousevalue = -0.04,
  poverty = 0.02,
  education = -5,
  popdensity = 0.0,
  pct_owner_occ = -0.005,
  summer_tmmx = -0.025,
  winter_tmmx = -0.025,
  summer_rmax = 0.025,
  winter_rmax = 0.025
)

confounders <- confounders[, names(confounders) %in% names(mortality_coeffs)]
which(!names(mortality_coeffs) %in% names(confounders))

# Calculate outcome
sim_outcome <- as.matrix(confounders) %*% mortality_coeffs
sim_outcome <- as.vector(sim_outcome) + air_quality_effect(treatment)
sim_outcome <- sim_outcome - quantile(sim_outcome, 0.01) + 10
sim_outcome[sim_outcome < 5] <- 5
hist(sim_outcome)

# # Combine treatment and confounders into a single dataframe for regression
# regression_data <- data.frame(
#   outcome = sim_outcome,
#   treatment = treatment,
#   confounders
# )
# 
# # Run linear regression
# model <- lm(outcome ~ ., data = regression_data)
# 
# # View summary of regression results
# summary(model)
# 
# # Get confidence interval for treatment effect
# treatment_ci <- confint(model)["treatment", ]
# print(paste("95% Confidence Interval:", 
#             round(treatment_ci[1], 4), "to", 
#             round(treatment_ci[2], 4)))



binary_treatment <- as.numeric(df$pm25 > 9)

# Combine treatment and confounders into a single dataframe for regression
binary_regression_data <- data.frame(
  outcome = sim_outcome,
  treatment = binary_treatment,
  confounders
)

# Run linear regression
binary_model <- lm(outcome ~ ., data = binary_regression_data)

# View summary of regression results
summary(binary_model)
cis <- confint(binary_model)['treatment', ]

## Causal-BLB Analysis

bcblb_out <- lapply(c(1, 5, 10), function(s){
  causal_blb(binary_regression_data, 
             y_method = 'glm', 
             prop_method = 'glm', 
             y_formula = outcome ~ treatment + mean_bmi + smoke_rate + hispanic + pct_blk + medhouseholdincome + medianhousevalue + poverty + education + popdensity + pct_owner_occ + summer_tmmx + winter_tmmx + summer_rmax + winter_rmax, 
             prop_formula = treatment ~ mean_bmi, 
             subsets = s, 
             r = 100,
             cores = 4)
})

bcblb_out <- rbindlist(bcblb_out)

bcblb_out$method <- 'Bootstrap'
bcblb_out$subsets <- c(1, 5, 10)
bcblb_out <- rbind(bcblb_out, data.frame(method = 'Binary model', lower_ci = cis[1], upper_ci = cis[2],
                                         subsets = NA))
saveRDS(bcblb_out, 'data/applications.rds')

kable(bcblb_out, format = "latex", booktabs = TRUE, align = "c", caption = "Sample Data Table",
      digits = 3) %>%
  kable_styling(latex_options = c("striped", "hold_position"), font_size = 10)
# 
# tf1 <- tempfile(fileext = ".parquet")
# write_parquet(data.frame(df), tf1)
