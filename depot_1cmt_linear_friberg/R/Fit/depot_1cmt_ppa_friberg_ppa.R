rm(list = ls())
cat("\014")

# library(trelliscopejs)
library(patchwork)
library(cmdstanr)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

nonmem_data <- read_csv("depot_1cmt_linear_friberg/Data/depot_1cmt_ppa_friberg_ppa.csv",
                        na = ".") %>% 
  rename_all(tolower) %>% 
  rename(ID = "id",
         DV = "dv") %>% 
  mutate(DV = if_else(is.na(DV), 5555555, DV),    # This value can be anything except NA. It'll be indexed away 
         bloq = if_else(is.na(bloq), -999, bloq)) # This value can be anything except NA. It'll be indexed away 

## Summary of BLOQ values
nonmem_data %>%
  filter(evid == 0) %>%
  group_by(ID, cmt) %>%
  summarize(lloq = unique(lloq),
            n_obs = n(),
            n_bloq = sum(bloq)) %>%
  filter(n_bloq > 0)

(p_pk <- ggplot(nonmem_data %>% 
                  group_by(ID) %>% 
                  mutate(Dose = factor(max(amt))) %>% 
                  ungroup() %>% 
                  filter(!is.na(DV), cmt == 2, bloq == 0) %>% 
                  mutate(cmt = factor(case_when(cmt == 2 ~ "PK",
                                                cmt == 4 ~ "PD",
                                                TRUE ~ "Dosing"),
                                      levels = c("PK", "PD", "Dosing")))) +
    geom_point(mapping = aes(x = time, y = DV, group = ID, color = Dose)) +
    geom_line(mapping = aes(x = time, y = DV, group = ID, color = Dose)) +
    theme_bw(18) +
    scale_y_continuous(name = latex2exp::TeX("Drug Conc. $(\\mu g/mL)$"),
                       trans = "log10") + 
    scale_x_continuous(name = "Time (h)",
                       breaks = seq(0, max(nonmem_data$time), by = 24),
                       labels = seq(0, max(nonmem_data$time), by = 24),
                       limits = c(0, max(nonmem_data$time))) +
    facet_wrap(~ cmt, scales = "free"))

(p_pd <- ggplot(nonmem_data %>% 
                  group_by(ID) %>% 
                  mutate(Dose = factor(max(amt))) %>% 
                  ungroup() %>% 
                  filter(!is.na(DV), cmt == 4, bloq == 0) %>% 
                  mutate(cmt = factor(case_when(cmt == 2 ~ "PK",
                                                cmt == 4 ~ "PD",
                                                TRUE ~ "Dosing"),
                                      levels = c("PK", "PD", "Dosing")))) +
    geom_point(mapping = aes(x = time, y = DV, group = ID, color = Dose)) +
    geom_line(mapping = aes(x = time, y = DV, group = ID, color = Dose)) +
    theme_bw(18) +
    scale_y_continuous(name = latex2exp::TeX("Neutrophils $(\\times 10^9/L)"),
                       trans = "identity") +
    scale_x_continuous(name = "Time (h)",
                       breaks = seq(0, max(nonmem_data$time), by = 24),
                       labels = seq(0, max(nonmem_data$time), by = 24),
                       limits = c(0, max(nonmem_data$time))) +
    facet_wrap(~ cmt, scales = "free"))

p_pk + 
  p_pd +
  plot_layout(guides = 'collect') &
  theme(legend.position = "bottom")

# p_pk +
#   facet_trelliscope(~ID, scales = "free_y", ncol = 2, nrow = 2)
# p_pd +
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
                  scale_sigma_a = 5,
                  lkj_df_sigma = 2,
                  location_tvmtt = 100,
                  location_tvcirc0 = 5,
                  location_tvgamma = 0.2,
                  location_tvalpha = 3e-4,
                  scale_tvmtt = 1,
                  scale_tvcirc0 = 1,
                  scale_tvgamma = 1,
                  scale_tvalpha = 1,
                  scale_omega_mtt = 0.4,
                  scale_omega_circ0 = 0.4,
                  scale_omega_gamma = 0.4,
                  scale_omega_alpha = 0.4,
                  lkj_df_omega_pd = 2,
                  scale_sigma_p_pd = 0.5,
                  scale_sigma_a_pd = 0.5,
                  lkj_df_sigma_pd = 2,
                  prior_only = 0,
                  no_gq_predictions = 0)

model <- cmdstan_model(
  "depot_1cmt_linear_friberg/Stan/Fit/depot_1cmt_ppa_friberg_ppa.stan",
  cpp_options = list(stan_threads = TRUE))

fit <- model$sample(
  data = stan_data,
  seed = 112358,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = parallel::detectCores()/4,
  iter_warmup = 500,
  iter_sampling = 1000,
  adapt_delta = 0.8,
  refresh = 100,
  max_treedepth = 10, 
  output_dir = "depot_1cmt_linear_friberg/Stan/Fits/Output",
  output_basename = "ppa_ppa",
  init = function() list(TVCL = rlnorm(1, log(0.6), 0.3),
                         TVVC = rlnorm(1, log(18), 0.3),
                         TVKA = rlnorm(1, log(1), 0.3),
                         omega = rlnorm(3, log(0.3), 0.3),
                         sigma = c(rlnorm(1, log(0.2), 0.3), 
                                   rlnorm(1, log(0.8), 0.3)),
                         TVMTT = rlnorm(1, log(120), 0.3),
                         TVCIRC0 = rlnorm(1, log(5), 0.3),
                         TVGAMMA = rlnorm(1, log(0.2), 0.3),
                         TVALPHA = rlnorm(1, log(3e-4), 0.3),
                         omega_pd = rlnorm(4, log(0.35), 0.3),
                         sigma_pd = c(rlnorm(1, log(0.2), 0.3), 
                                      rlnorm(1, log(0.8), 0.3))))

fit$save_object("depot_1cmt_linear_friberg/Stan/Fits/depot_1cmt_ppa_friberg_ppa.rds")
fit$save_data_file(dir = "depot_1cmt_linear_friberg/Stan/Fits/Stan_Data",
                   basename = "ppa_ppa", timestamp = FALSE, random = FALSE)
