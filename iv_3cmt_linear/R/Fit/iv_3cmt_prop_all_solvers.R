rm(list = ls())
cat("\014")

library(trelliscopejs)
library(cmdstanr)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

solver <- 1 # mat exp = 1, rk45 = 2, bdf = 3

nonmem_data <- read_csv("iv_3cmt_linear/Data/iv_3cmt_prop.csv",
                        na = ".") %>% 
  rename_all(tolower) %>% 
  rename(ID = "id",
         DV = "dv") %>% 
  mutate(DV = if_else(is.na(DV), 5555555, DV),    # This value can be anything except NA. It'll be indexed away 
         bloq = if_else(is.na(bloq), -999, bloq), # This value can be anything except NA. It'll be indexed away 
         cmt = if_else(solver == 1, 2, 1))

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

p1 +
  facet_wrap(~ID, scales = "free_y", labeller = label_both, ncol = 4)

p1 +
  facet_trelliscope(~ID, scales = "free_y", ncol = 2, nrow = 2)


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
                  location_tvcl = 2,
                  location_tvvc = 5,
                  location_tvq1 = 2,
                  location_tvvp1 = 3.5,
                  location_tvq2 = 0.5,
                  location_tvvp2 = 4,
                  scale_tvcl = 1,
                  scale_tvvc = 1,
                  scale_tvq1 = 1,
                  scale_tvvp1 = 1,
                  scale_tvq2 = 1,
                  scale_tvvp2 = 1,
                  scale_omega_cl = 0.4,
                  scale_omega_vc = 0.4,
                  scale_omega_q1 = 0.4,
                  scale_omega_vp1 = 0.4,
                  scale_omega_q2 = 0.4,
                  scale_omega_vp2 = 0.4,
                  lkj_df_omega = 2,
                  scale_sigma_p = 0.5,
                  prior_only = 0,
                  solver = 1)

model <- cmdstan_model("iv_3cmt_linear/Stan/Fit/iv_3cmt_prop_all_solvers.stan",
                       cpp_options = list(stan_threads = TRUE))

fit_mat_exp <- model$sample(
  data = stan_data,
  seed = 11235,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = parallel::detectCores()/4,
  iter_warmup = 300,
  iter_sampling = 200,
  adapt_delta = 0.8,
  refresh = 50,
  max_treedepth = 10,
  init = function() list(TVCL = rlnorm(1, log(2), 0.3),
                         TVVC = rlnorm(1, log(5), 0.3),
                         TVQ1 = rlnorm(1, log(2), 0.3),
                         TVVP1 = rlnorm(1, log(3.5), 0.3),
                         TVQ2 = rlnorm(1, log(0.8), 0.3),
                         TVVP2 = rlnorm(1, log(4), 0.3),
                         omega = rlnorm(6, log(0.3), 0.3),
                         sigma_p = rlnorm(1, log(0.2), 0.3)))

fit_mat_exp$save_object("iv_3cmt_linear/Stan/Fits/iv_3cmt_prop_mat_exp.rds")


stan_data$solver <- 2
fit_rk45 <- model$sample(
  data = stan_data,
  seed = 11235,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = parallel::detectCores()/4,
  iter_warmup = 300,
  iter_sampling = 200,
  adapt_delta = 0.8,
  refresh = 50,
  max_treedepth = 10,
  init = function() list(TVCL = rlnorm(1, log(2), 0.3),
                         TVVC = rlnorm(1, log(5), 0.3),
                         TVQ1 = rlnorm(1, log(2), 0.3),
                         TVVP1 = rlnorm(1, log(3.5), 0.3),
                         TVQ2 = rlnorm(1, log(0.8), 0.3),
                         TVVP2 = rlnorm(1, log(4), 0.3),
                         omega = rlnorm(6, log(0.3), 0.3),
                         sigma_p = rlnorm(1, log(0.2), 0.3)))

fit_rk45$save_object("iv_3cmt_linear/Stan/Fits/iv_3cmt_prop_rk45.rds")

stan_data$solver <- 3
fit_bdf <- model$sample(
  data = stan_data,
  seed = 11235,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = parallel::detectCores()/4,
  iter_warmup = 300,
  iter_sampling = 200,
  adapt_delta = 0.8,
  refresh = 50,
  max_treedepth = 10,
  init = function() list(TVCL = rlnorm(1, log(2), 0.3),
                         TVVC = rlnorm(1, log(5), 0.3),
                         TVQ1 = rlnorm(1, log(2), 0.3),
                         TVVP1 = rlnorm(1, log(3.5), 0.3),
                         TVQ2 = rlnorm(1, log(0.8), 0.3),
                         TVVP2 = rlnorm(1, log(4), 0.3),
                         omega = rlnorm(6, log(0.3), 0.3),
                         sigma_p = rlnorm(1, log(0.2), 0.3)))

fit_bdf$save_object("iv_3cmt_linear/Stan/Fits/iv_3cmt_prop_bdf.rds")

