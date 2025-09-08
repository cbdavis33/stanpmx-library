rm(list = ls())
cat("\014")

library(patchwork)
library(cmdstanr)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

nonmem_data <- read_csv("depot_1cmt_linear_ir3/Data/depot_1cmt_prop_ir3_prop.csv",
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
                                                cmt == 3 ~ "PD",
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
                  filter(!is.na(DV), cmt == 3, bloq == 0) %>% 
                  mutate(cmt = factor(case_when(cmt == 2 ~ "PK",
                                                cmt == 3 ~ "PD",
                                                TRUE ~ "Dosing"),
                                      levels = c("PK", "PD", "Dosing")))) +
    geom_point(mapping = aes(x = time, y = DV, group = ID, color = Dose)) +
    geom_line(mapping = aes(x = time, y = DV, group = ID, color = Dose)) +
    theme_bw(18) +
    scale_y_continuous(name = latex2exp::TeX("Response (unit)"),
                       trans = "log10") + 
    scale_x_continuous(name = "Time (h)",
                       breaks = seq(0, max(nonmem_data$time), by = 24),
                       labels = seq(0, max(nonmem_data$time), by = 24),
                       limits = c(0, max(nonmem_data$time))) +
    facet_wrap(~ cmt, scales = "free"))

p_pk + 
  p_pd +
  plot_layout(guides = 'collect') &
  theme(legend.position = "bottom")

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
                  location_tvcl = 0.6,
                  location_tvvc = 13,
                  location_tvka = 0.9,
                  scale_tvcl = 1,
                  scale_tvvc = 1,
                  scale_tvka = 1,
                  scale_omega_cl = 0.4,
                  scale_omega_vc = 0.4,
                  scale_omega_ka = 0.4,
                  lkj_df_omega_pk = 2,
                  scale_sigma_p_pk = 0.5,
                  location_tvkin = 5,
                  location_tvkout = 0.2,
                  location_tvsc50 = 10,
                  location_tvsmax = 5,
                  scale_tvkin = 1,
                  scale_tvkout = 1,
                  scale_tvsc50 = 1,
                  scale_tvsmax = 1,
                  scale_omega_kin = 0.4,
                  scale_omega_kout = 0.4,
                  scale_omega_sc50 = 0.4,
                  scale_omega_smax = 0.4,
                  lkj_df_omega_pd = 2,
                  scale_sigma_p_pd = 0.5,
                  prior_only = 0,
                  no_gq_predictions = 0,
                  coupled = 0,
                  solver = 1) # rk45 = 1, bdf = 2, adams = 3 (adams is non-coupled only)

model <- cmdstan_model(
  "depot_1cmt_linear_ir3/Stan/Fit/depot_1cmt_prop_ir3_prop_all_solvers.stan",
  cpp_options = list(stan_threads = TRUE))

fit_rk45 <- model$sample(
  data = stan_data,
  seed = 11235,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = parallel::detectCores()/4,
  iter_warmup = 500,
  iter_sampling = 100,
  adapt_delta = 0.8,
  refresh = 10,
  max_treedepth = 10,
  init = function() 
    with(stan_data,
         list(TVCL = rlnorm(1, log(location_tvcl), scale_tvcl),
              TVVC = rlnorm(1, log(location_tvvc), scale_tvvc),
              TVKA = rlnorm(1, log(location_tvka), scale_tvka),
              omega_pk = abs(rnorm(3, 0, c(scale_omega_cl,
                                        scale_omega_vc,
                                        scale_omega_ka))),
              sigma_p_pk = abs(rnorm(1, 0, scale_sigma_p_pk)),
              TVKIN = rlnorm(1, log(location_tvkin), scale_tvkin),
              TVKOUT = rlnorm(1, log(location_tvkout), scale_tvkout),
              TVSC50 = rlnorm(1, log(location_tvsc50), scale_tvsc50),
              TVSMAX = rlnorm(1, log(location_tvsmax), scale_tvsmax),
              omega_pd = abs(rnorm(4, 0, c(scale_omega_kin,
                                           scale_omega_kout,
                                           scale_omega_sc50,
                                           scale_omega_smax))),
              sigma_p_pd = abs(rnorm(1, 0, scale_sigma_p_pd)))))


fit_rk45$save_object("depot_1cmt_linear_ir3/Stan/Fits/depot_1cmt_prop_ir3_prop_rk45.rds")

stan_data$solver <- 3
fit_adams <- model$sample(
  data = stan_data,
  seed = 11235,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = parallel::detectCores()/4,
  iter_warmup = 500,
  iter_sampling = 100,
  adapt_delta = 0.8,
  refresh = 10,
  max_treedepth = 10,
  init = function() 
    with(stan_data,
         list(TVCL = rlnorm(1, log(location_tvcl), scale_tvcl),
              TVVC = rlnorm(1, log(location_tvvc), scale_tvvc),
              TVKA = rlnorm(1, log(location_tvka), scale_tvka),
              omega_pk = abs(rnorm(3, 0, c(scale_omega_cl,
                                           scale_omega_vc,
                                           scale_omega_ka))),
              sigma_p_pk = abs(rnorm(1, 0, scale_sigma_p_pk)),
              TVKIN = rlnorm(1, log(location_tvkin), scale_tvkin),
              TVKOUT = rlnorm(1, log(location_tvkout), scale_tvkout),
              TVSC50 = rlnorm(1, log(location_tvsc50), scale_tvsc50),
              TVSMAX = rlnorm(1, log(location_tvsmax), scale_tvsmax),
              omega_pd = abs(rnorm(4, 0, c(scale_omega_kin,
                                           scale_omega_kout,
                                           scale_omega_sc50,
                                           scale_omega_smax))),
              sigma_p_pd = abs(rnorm(1, 0, scale_sigma_p_pd)))))

fit_adams$save_object(
  "depot_1cmt_linear_ir3/Stan/Fits/depot_1cmt_prop_ir3_prop_adams.rds")

stan_data$solver <- 2
fit_bdf <- model$sample(
  data = stan_data,
  seed = 11235,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = parallel::detectCores()/4,
  iter_warmup = 500,
  iter_sampling = 100,
  adapt_delta = 0.8,
  refresh = 10,
  max_treedepth = 10,
  init = function() 
    with(stan_data,
         list(TVCL = rlnorm(1, log(location_tvcl), scale_tvcl),
              TVVC = rlnorm(1, log(location_tvvc), scale_tvvc),
              TVKA = rlnorm(1, log(location_tvka), scale_tvka),
              omega_pk = abs(rnorm(3, 0, c(scale_omega_cl,
                                           scale_omega_vc,
                                           scale_omega_ka))),
              sigma_p_pk = abs(rnorm(1, 0, scale_sigma_p_pk)),
              TVKIN = rlnorm(1, log(location_tvkin), scale_tvkin),
              TVKOUT = rlnorm(1, log(location_tvkout), scale_tvkout),
              TVSC50 = rlnorm(1, log(location_tvsc50), scale_tvsc50),
              TVSMAX = rlnorm(1, log(location_tvsmax), scale_tvsmax),
              omega_pd = abs(rnorm(4, 0, c(scale_omega_kin,
                                           scale_omega_kout,
                                           scale_omega_sc50,
                                           scale_omega_smax))),
              sigma_p_pd = abs(rnorm(1, 0, scale_sigma_p_pd)))))

fit_bdf$save_object(
  "depot_1cmt_linear_ir3/Stan/Fits/depot_1cmt_prop_ir3_prop_bdf.rds")


############################################################################

stan_data$coupled <- 1
stan_data$solver <- 1
fit_rk45_coupled <- model$sample(
  data = stan_data,
  seed = 11235,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = parallel::detectCores()/4,
  iter_warmup = 500,
  iter_sampling = 100,
  adapt_delta = 0.8,
  refresh = 10,
  max_treedepth = 10,
  init = function() 
    with(stan_data,
         list(TVCL = rlnorm(1, log(location_tvcl), scale_tvcl),
              TVVC = rlnorm(1, log(location_tvvc), scale_tvvc),
              TVKA = rlnorm(1, log(location_tvka), scale_tvka),
              omega_pk = abs(rnorm(3, 0, c(scale_omega_cl,
                                           scale_omega_vc,
                                           scale_omega_ka))),
              sigma_p_pk = abs(rnorm(1, 0, scale_sigma_p_pk)),
              TVKIN = rlnorm(1, log(location_tvkin), scale_tvkin),
              TVKOUT = rlnorm(1, log(location_tvkout), scale_tvkout),
              TVSC50 = rlnorm(1, log(location_tvsc50), scale_tvsc50),
              TVSMAX = rlnorm(1, log(location_tvsmax), scale_tvsmax),
              omega_pd = abs(rnorm(4, 0, c(scale_omega_kin,
                                           scale_omega_kout,
                                           scale_omega_sc50,
                                           scale_omega_smax))),
              sigma_p_pd = abs(rnorm(1, 0, scale_sigma_p_pd)))))


fit_rk45_coupled$save_object(
  "depot_1cmt_linear_ir3/Stan/Fits/depot_1cmt_prop_ir3_prop_rk45_coupled.rds")


stan_data$solver <- 2
fit_bdf_coupled <- model$sample(
  data = stan_data,
  seed = 11235,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = parallel::detectCores()/4,
  iter_warmup = 500,
  iter_sampling = 100,
  adapt_delta = 0.8,
  refresh = 10,
  max_treedepth = 10,
  init = function() 
    with(stan_data,
         list(TVCL = rlnorm(1, log(location_tvcl), scale_tvcl),
              TVVC = rlnorm(1, log(location_tvvc), scale_tvvc),
              TVKA = rlnorm(1, log(location_tvka), scale_tvka),
              omega_pk = abs(rnorm(3, 0, c(scale_omega_cl,
                                           scale_omega_vc,
                                           scale_omega_ka))),
              sigma_p_pk = abs(rnorm(1, 0, scale_sigma_p_pk)),
              TVKIN = rlnorm(1, log(location_tvkin), scale_tvkin),
              TVKOUT = rlnorm(1, log(location_tvkout), scale_tvkout),
              TVSC50 = rlnorm(1, log(location_tvsc50), scale_tvsc50),
              TVSMAX = rlnorm(1, log(location_tvsmax), scale_tvsmax),
              omega_pd = abs(rnorm(4, 0, c(scale_omega_kin,
                                           scale_omega_kout,
                                           scale_omega_sc50,
                                           scale_omega_smax))),
              sigma_p_pd = abs(rnorm(1, 0, scale_sigma_p_pd)))))

fit_bdf_coupled$save_object(
  "depot_1cmt_linear_ir3/Stan/Fits/depot_1cmt_prop_ir3_prop_bdf_coupled.rds")

