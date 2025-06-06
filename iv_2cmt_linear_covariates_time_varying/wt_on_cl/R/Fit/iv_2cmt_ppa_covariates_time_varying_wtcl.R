rm(list = ls())
cat("\014")

library(trelliscopejs)
library(cmdstanr)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

nonmem_data <- read_csv(
  "iv_2cmt_linear_covariates_time_varying/wt_on_cl/Data/nonmem_IV.csv",
  na = ".") %>% 
  rename_all(tolower) %>% 
  rename(ID = "id",
         DV = "dv") %>% 
  mutate(cmt = 2, # using the analytical solution
         lloq = 0,
         DV = if_else(is.na(DV), 5555555, DV),    # This value can be anything except NA. It'll be indexed away 
         bloq = 0)                                # This value can be anything except NA. It'll be indexed away 

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
                mutate(ID = factor(ID)) %>%
                ungroup() %>%
                filter(mdv == 0)) +
    geom_line(mapping = aes(x = time, y = DV, group = ID, color = ID)) +
    geom_point(mapping = aes(x = time, y = DV, group = ID, color = ID)) +
    scale_color_discrete(guide = "none") +
    scale_y_continuous(name = latex2exp::TeX("Drug Conc. $(\\mu g/mL)$"),
                       limits = c(NA, NA),
                       trans = "log10") + 
    scale_x_continuous(name = "Time (units)",
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
                  location_tvcl = 1,
                  location_tvvc = 10,
                  location_tvq = 2,
                  location_tvvp = 10,
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
                  scale_sigma_a = 3,
                  lkj_df_sigma = 2,
                  prior_only = 0,
                  no_gq_predictions = 0)

model <- cmdstan_model(
  "iv_2cmt_linear_covariates_time_varying/wt_on_cl/Stan/Fit/iv_2cmt_ppa_covariates_time_varying_wtcl.stan",
  cpp_options = list(stan_threads = TRUE))

fit <- model$sample(
  data = stan_data,
  seed = 11235,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = parallel::detectCores()/4,
  iter_warmup = 500,
  iter_sampling = 1000,
  adapt_delta = 0.8,
  refresh = 100,
  max_treedepth = 10,
  save_warmup = TRUE,
  output_dir = "iv_2cmt_linear_covariates_time_varying/wt_on_cl/Stan/Fits/Output",
  output_basename = "andras",
  init = function() list(TVCL = rlnorm(1, log(0.25), 0.3),
                         TVVC = rlnorm(1, log(3), 0.3),
                         TVQ = rlnorm(1, log(1), 0.3),
                         TVVP = rlnorm(1, log(4), 0.3),
                         theta_cl_wt = rnorm(1), 
                         omega = rlnorm(4, log(0.3), 0.3),
                         sigma = rlnorm(2, log(0.2), 0.3)))

fit$save_object(
  str_c("iv_2cmt_linear_covariates_time_varying/wt_on_cl/Stan/Fits/", 
        "andras", 
        ".rds"))

fit$save_data_file(
  dir = "iv_2cmt_linear_covariates_time_varying/Generic/Stan/Fits/Stan_Data",
  basename = "andras", timestamp = FALSE, 
  random = FALSE)
