rm(list = ls())
cat("\014")

library(trelliscopejs)
library(cmdstanr)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

locf <- FALSE

file_path <- if_else(isTRUE(locf), 
                     "iv_2cmt_linear_covariates_time_varying/Generic/Data/iv_2cmt_prop_covariates_time_varying_generic_locf.csv",
                     "iv_2cmt_linear_covariates_time_varying/Generic/Data/iv_2cmt_prop_covariates_time_varying_generic_nocb.csv")

nonmem_data <- read_csv(
  file_path,
  na = ".") %>% 
  rename_all(tolower) %>% 
  rename(ID = "id",
         DV = "dv") %>% 
  mutate(cmt = 2,
         DV = if_else(is.na(DV), 5555555, DV),    # This value can be anything except NA. It'll be indexed away 
         bloq = if_else(is.na(bloq), -999, bloq)) # This value can be anything except NA. It'll be indexed away 

## Summary of BLOQ values
nonmem_data %>%
  filter(evid == 0) %>%
  group_by(ID) %>%
  summarize(lloq = unique(lloq),
            n_obs = n(),
            n_bloq = sum(bloq)) %>%
  filter(n_bloq > 0)

(p1 <- ggplot(nonmem_data %>%
                group_by(ID) %>%
                mutate(Dose = factor(max(amt, na.rm = TRUE))) %>%
                ungroup() %>%
                filter(mdv == 0)) +
    geom_line(mapping = aes(x = time, y = DV, group = ID, color = Dose)) +
    geom_point(mapping = aes(x = time, y = DV, group = ID, color = Dose)) +
    scale_color_discrete(name = "Dose (mg)") +
    scale_y_continuous(name = latex2exp::TeX("Drug Conc. $(\\mu g/mL)$"),
                       limits = c(NA, NA),
                       trans = "log10") + 
    scale_x_continuous(name = "Time (d)",
                       breaks = seq(0, max(nonmem_data$time), by = 14),
                       labels = seq(0, max(nonmem_data$time), by = 14),
                       limits = c(0, max(nonmem_data$time))) +
    theme_bw(18) +
    theme(axis.text = element_text(size = 14, face = "bold"),
          axis.title = element_text(size = 18, face = "bold"),
          axis.line = element_line(linewidth = 2),
          legend.position = "bottom"))

# p1 +
#   facet_wrap(~ID, scales = "free_y", labeller = label_both, ncol = 4)

n_subjects <- nonmem_data %>%  # number of individuals
  distinct(ID) %>%
  count() %>%
  deframe()

n_total <- nrow(nonmem_data)   # total number of records

i_obs <- nonmem_data %>%
  mutate(row_num = 1:n()) %>%
  filter(evid == 0) %>%
  select(row_num) %>%
  deframe()

n_obs <- length(i_obs)

subj_start <- nonmem_data %>%
  mutate(row_num = 1:n()) %>%
  group_by(ID) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(row_num) %>%
  deframe()

subj_end <- c(subj_start[-1] - 1, n_total)

wt <- nonmem_data %>% 
  pull(wt)

race_asian <- nonmem_data %>% 
  pull(race_asian)

egfr <- nonmem_data %>% 
  pull(egfr)

stan_data <- list(n_subjects = n_subjects,
                  n_total = n_total,
                  n_obs = n_obs,
                  i_obs = i_obs,
                  ID = nonmem_data$ID,
                  amt = nonmem_data$amt,
                  cmt = nonmem_data$cmt,
                  evid = nonmem_data$evid,
                  rate = nonmem_data$rate,
                  ii = nonmem_data$ii,
                  addl = nonmem_data$addl,
                  ss = nonmem_data$ss,
                  time = nonmem_data$time,
                  dv = nonmem_data$DV,
                  subj_start = subj_start,
                  subj_end = subj_end,
                  lloq = nonmem_data$lloq,
                  bloq = nonmem_data$bloq,
                  wt = wt,
                  race_asian = race_asian,
                  egfr = egfr,
                  location_tvcl = 0.25,
                  location_tvvc = 3,
                  location_tvq = 1,
                  location_tvvp = 4,
                  scale_tvcl = 1,
                  scale_tvvc = 1,
                  scale_tvq = 1,
                  scale_tvvp = 1,
                  scale_omega_cl = 0.4,
                  scale_omega_vc = 0.4,
                  scale_omega_q = 0.4,
                  scale_omega_vp = 0.4,
                  lkj_df_omega = 2,
                  scale_sigma_p = 0.5,
                  prior_only = 0,
                  solver = 1,
                  no_gq_predictions = 0)

model <- cmdstan_model(
  "iv_2cmt_linear_covariates_time_varying/Generic/Stan/Fit/iv_2cmt_prop_all_solvers_covariates_time_varying_generic.stan",
  cpp_options = list(stan_threads = TRUE))

fit_analytical <- model$sample(
  data = stan_data,
  seed = 11235,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = parallel::detectCores()/4,
  iter_warmup = 500,
  iter_sampling = 1000,
  adapt_delta = 0.8,
  refresh = 500,
  max_treedepth = 10,
  save_warmup = TRUE,
  output_dir = "iv_2cmt_linear_covariates_time_varying/Generic/Stan/Fits/Output",
  output_basename = "analytical_all_prop",
  init = function() list(TVCL = rlnorm(1, log(0.25), 0.3),
                         TVVC = rlnorm(1, log(3), 0.3),
                         TVQ = rlnorm(1, log(1), 0.3),
                         TVVP = rlnorm(1, log(4), 0.3),
                         theta_cl_wt = rnorm(1), 
                         theta_vc_wt = rnorm(1), 
                         theta_q_wt = rnorm(1), 
                         theta_vp_wt = rnorm(1), 
                         theta_vc_race_asian = rnorm(1), 
                         theta_cl_egfr = rnorm(1),
                         omega = rlnorm(4, log(0.3), 0.3),
                         sigma_p = rlnorm(1, log(0.2), 0.3)))

fit_analytical$save_object(
  "iv_2cmt_linear_covariates_time_varying/Generic/Stan/Fits/analytical_all_prop.rds")

fit_analytical$save_data_file(
  dir = "iv_2cmt_linear_covariates_time_varying/Generic/Stan/Fits/Stan_Data",
  basename = "analytical_all_prop", timestamp = FALSE, 
  random = FALSE)

stan_data$solver <- 2
stan_data$cmt <- stan_data$cmt - 1
fit_mat_exp <- model$sample(
  data = stan_data,
  seed = 11235,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = parallel::detectCores()/4,
  iter_warmup = 500,
  iter_sampling = 1000,
  adapt_delta = 0.8,
  refresh = 500,
  max_treedepth = 10,
  save_warmup = TRUE,
  output_dir = "iv_2cmt_linear_covariates_time_varying/Generic/Stan/Fits/Output",
  output_basename = "mat_exp_all_prop",
  init = function() list(TVCL = rlnorm(1, log(0.25), 0.3),
                         TVVC = rlnorm(1, log(3), 0.3),
                         TVQ = rlnorm(1, log(1), 0.3),
                         TVVP = rlnorm(1, log(4), 0.3),
                         theta_cl_wt = rnorm(1), 
                         theta_vc_wt = rnorm(1), 
                         theta_q_wt = rnorm(1), 
                         theta_vp_wt = rnorm(1), 
                         theta_vc_race_asian = rnorm(1), 
                         theta_cl_egfr = rnorm(1),
                         omega = rlnorm(4, log(0.3), 0.3),
                         sigma_p = rlnorm(1, log(0.2), 0.3)))

fit_mat_exp$save_object(
  "iv_2cmt_linear_covariates_time_varying/Generic/Stan/Fits/mat_exp_all_prop.rds")

fit_mat_exp$save_data_file(
  dir = "iv_2cmt_linear_covariates_time_varying/Generic/Stan/Fits/Stan_Data",
  basename = "mat_exp_all_prop", timestamp = FALSE, 
  random = FALSE)

stan_data$solver <- 3
fit_rk45 <- model$sample(
  data = stan_data,
  seed = 11235,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = parallel::detectCores()/4,
  iter_warmup = 500,
  iter_sampling = 1000,
  adapt_delta = 0.8,
  refresh = 500,
  max_treedepth = 10,
  save_warmup = TRUE,
  output_dir = "iv_2cmt_linear_covariates_time_varying/Generic/Stan/Fits/Output",
  output_basename = "rk45_all_prop",
  init = function() list(TVCL = rlnorm(1, log(0.25), 0.3),
                         TVVC = rlnorm(1, log(3), 0.3),
                         TVQ = rlnorm(1, log(1), 0.3),
                         TVVP = rlnorm(1, log(4), 0.3),
                         theta_cl_wt = rnorm(1), 
                         theta_vc_wt = rnorm(1), 
                         theta_q_wt = rnorm(1), 
                         theta_vp_wt = rnorm(1), 
                         theta_vc_race_asian = rnorm(1), 
                         theta_cl_egfr = rnorm(1),
                         omega = rlnorm(4, log(0.3), 0.3),
                         sigma_p = rlnorm(1, log(0.2), 0.3)))

fit_rk45$save_object(
  "iv_2cmt_linear_covariates_time_varying/Generic/Stan/Fits/rk45_all_prop.rds")

fit_rk45$save_data_file(
  dir = "iv_2cmt_linear_covariates_time_varying/Generic/Stan/Fits/Stan_Data",
  basename = "rk45_all_prop", timestamp = FALSE, 
  random = FALSE)

stan_data$solver <- 4
fit_bdf <- model$sample(
  data = stan_data,
  seed = 11235,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = parallel::detectCores()/4,
  iter_warmup = 500,
  iter_sampling = 1000,
  adapt_delta = 0.8,
  refresh = 500,
  max_treedepth = 10,
  save_warmup = TRUE,
  output_dir = "iv_2cmt_linear_covariates_time_varying/Generic/Stan/Fits/Output",
  output_basename = "bdf_all_prop",
  init = function() list(TVCL = rlnorm(1, log(0.25), 0.3),
                         TVVC = rlnorm(1, log(3), 0.3),
                         TVQ = rlnorm(1, log(1), 0.3),
                         TVVP = rlnorm(1, log(4), 0.3),
                         theta_cl_wt = rnorm(1), 
                         theta_vc_wt = rnorm(1), 
                         theta_q_wt = rnorm(1), 
                         theta_vp_wt = rnorm(1), 
                         theta_vc_race_asian = rnorm(1), 
                         theta_cl_egfr = rnorm(1),
                         omega = rlnorm(4, log(0.3), 0.3),
                         sigma_p = rlnorm(1, log(0.2), 0.3)))

fit_bdf$save_object(
  "iv_2cmt_linear_covariates_time_varying/Generic/Stan/Fits/bdf_all_prop.rds")

fit_bdf$save_data_file(
  dir = "iv_2cmt_linear_covariates_time_varying/Generic/Stan/Fits/Stan_Data",
  basename = "bdf_all_prop", timestamp = FALSE, 
  random = FALSE)

