rm(list = ls())
cat("\014")

library(trelliscopejs)
library(cmdstanr)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

nonmem_data <- read_csv(
  "depot_1cmt_linear_covariates/Data/depot_1cmt_ppa_covariates.csv",
  na = ".") %>% 
  rename_all(tolower) %>% 
  rename(ID = "id",
         DV = "dv") %>% 
  mutate(DV = if_else(is.na(DV), 5555555, DV),    # This value can be anything except NA. It'll be indexed away 
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
    scale_y_continuous(name = latex2exp::TeX("$Drug\\;Conc.\\;(\\mu g/mL)$"),
                       limits = c(NA, NA),
                       trans = "log10") +
    scale_x_continuous(name = "Time (d)",
                       breaks = seq(0, 216, by = 24),
                       labels = seq(0, 216/24, by = 24/24),
                       limits = c(0, NA)) +
    theme_bw(18) +
    theme(axis.text = element_text(size = 14, face = "bold"),
          axis.title = element_text(size = 18, face = "bold"),
          axis.line = element_line(linewidth = 2),
          legend.position = "bottom"))

# p1 +
#   facet_wrap(~ID, scales = "free_y", labeller = label_both, ncol = 4)
# 
# p1 +
#   facet_trelliscope(~ID, scales = "free_y", ncol = 2, nrow = 2)


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
  group_by(ID) %>% 
  distinct(wt) %>% 
  ungroup() %>% 
  pull(wt)

cmppi <- nonmem_data %>% 
  group_by(ID) %>% 
  distinct(cmppi) %>% 
  ungroup() %>% 
  pull(cmppi)

egfr <- nonmem_data %>% 
  group_by(ID) %>% 
  distinct(egfr) %>% 
  ungroup() %>% 
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
                  cmppi = cmppi,
                  egfr = egfr,
                  location_tvcl = 1,
                  location_tvvc = 8,
                  location_tvka = 0.8,
                  scale_tvcl = 1,
                  scale_tvvc = 1,
                  scale_tvka = 1,
                  scale_omega_cl = 0.4,
                  scale_omega_vc = 0.4,
                  scale_omega_ka = 0.4,
                  lkj_df_omega = 2,
                  scale_sigma_p = 0.5,
                  scale_sigma_a = 0.5,
                  lkj_df_sigma = 2,
                  prior_only = 0)

model <- cmdstan_model(
  "depot_1cmt_linear_covariates/Stan/Fit/depot_1cmt_ppa_no_threading_covariates.stan")

fit <- model$sample(data = stan_data,
                    seed = 11235,
                    chains = 4,
                    parallel_chains = 4,
                    iter_warmup = 500,
                    iter_sampling = 1000,
                    adapt_delta = 0.8,
                    refresh = 500,
                    max_treedepth = 10,
                    init = function() list(TVCL = rlnorm(1, log(1), 0.3),
                                           TVVC = rlnorm(1, log(8), 0.3),
                                           TVKA = rlnorm(1, log(0.8), 0.3),
                                           theta_cl_wt = rnorm(1), 
                                           theta_vc_wt = rnorm(1), 
                                           theta_ka_cmppi = rnorm(1), 
                                           theta_cl_egfr = rnorm(1),
                                           omega = rlnorm(3, log(0.3), 0.3),
                                           sigma = rlnorm(2, log(0.2), 0.3)))

fit$save_object(
  "depot_1cmt_linear_covariates/Stan/Fits/depot_1cmt_ppa_no_threading_covariates.rds")