rm(list = ls())
cat("\014")

library(patchwork)
library(cmdstanr)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

nonmem_data <- read_csv("depot_1cmt_linear_ir2/Data/depot_1cmt_exp_ir2_exp.csv",
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

# (p_pk <- ggplot(nonmem_data %>% 
#                   group_by(ID) %>% 
#                   mutate(Dose = factor(max(amt))) %>% 
#                   ungroup() %>% 
#                   filter(!is.na(DV), cmt == 2, bloq == 0) %>% 
#                   mutate(cmt = factor(case_when(cmt == 2 ~ "PK",
#                                                 cmt == 3 ~ "PD",
#                                                 TRUE ~ "Dosing"),
#                                       levels = c("PK", "PD", "Dosing")))) +
#     geom_point(mapping = aes(x = time, y = DV, group = ID, color = Dose)) +
#     geom_line(mapping = aes(x = time, y = DV, group = ID, color = Dose)) +
#     theme_bw(18) +
#     scale_y_continuous(name = latex2exp::TeX("Drug Conc. $(\\mu g/mL)$"),
#                        trans = "log10") + 
#     scale_x_continuous(name = "Time (h)",
#                        breaks = seq(0, max(nonmem_data$time), by = 24),
#                        labels = seq(0, max(nonmem_data$time), by = 24),
#                        limits = c(0, max(nonmem_data$time))) +
#     facet_wrap(~ cmt, scales = "free"))
# 
# (p_pd <- ggplot(nonmem_data %>% 
#                   group_by(ID) %>% 
#                   mutate(Dose = factor(max(amt))) %>% 
#                   ungroup() %>% 
#                   filter(!is.na(DV), cmt == 3, bloq == 0) %>% 
#                   mutate(cmt = factor(case_when(cmt == 2 ~ "PK",
#                                                 cmt == 3 ~ "PD",
#                                                 TRUE ~ "Dosing"),
#                                       levels = c("PK", "PD", "Dosing")))) +
#     geom_point(mapping = aes(x = time, y = DV, group = ID, color = Dose)) +
#     geom_line(mapping = aes(x = time, y = DV, group = ID, color = Dose)) +
#     theme_bw(18) +
#     scale_y_continuous(name = latex2exp::TeX("Response (unit)"),
#                        trans = "log10") + 
#     scale_x_continuous(name = "Time (h)",
#                        breaks = seq(0, max(nonmem_data$time), by = 24),
#                        labels = seq(0, max(nonmem_data$time), by = 24),
#                        limits = c(0, max(nonmem_data$time))) +
#     facet_wrap(~ cmt, scales = "free"))
# 
# p_pk + 
#   p_pd +
#   plot_layout(guides = 'collect') &
#   theme(legend.position = "bottom")

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
                  scale_sigma_pk = 0.5,
                  location_tvkin = 5,
                  location_tvkout = 0.2,
                  location_tvic50 = 10,
                  scale_tvkin = 1,
                  scale_tvkout = 1,
                  scale_tvic50 = 1,
                  scale_omega_kin = 0.4,
                  scale_omega_kout = 0.4,
                  scale_omega_ic50 = 0.4,
                  lkj_df_omega_pd = 2,
                  scale_sigma_pd = 0.5,
                  prior_only = 1,
                  no_gq_predictions = 0)

model <- cmdstan_model(
  "depot_1cmt_linear_ir2/Stan/Fit/depot_1cmt_exp_ir2_exp.stan",
  cpp_options = list(stan_threads = TRUE))

prior <- model$sample(
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
  init = function() 
    with(stan_data,
         list(TVCL = rlnorm(1, log(location_tvcl), scale_tvcl),
              TVVC = rlnorm(1, log(location_tvvc), scale_tvvc),
              TVKA = rlnorm(1, log(location_tvka), scale_tvka),
              omega_pk = abs(rnorm(3, 0, c(scale_omega_cl,
                                           scale_omega_vc,
                                           scale_omega_ka))),
              sigma_pk = abs(rnorm(1, 0, scale_sigma_pk)),
              TVKIN = rlnorm(1, log(location_tvkin), scale_tvkin),
              TVKOUT = rlnorm(1, log(location_tvkout), scale_tvkout),
              TVIC50 = rlnorm(1, log(location_tvic50), scale_tvic50),
              omega_pd = abs(rnorm(3, 0, c(scale_omega_kin,
                                           scale_omega_kout,
                                           scale_omega_ic50))),
              sigma_pd = abs(rnorm(1, 0, scale_sigma_pd)))))

prior_df <- prior$draws(format = "draws_df") %>%  
  drop_na()

prior_preds_summary <- prior_df %>%
  spread_draws(c(ipred, epred, dv_ppc)[i]) %>% 
  mean_qi(ipred, epred, dv_ppc, .width = 0.95) %>%
  mutate(DV = nonmem_data$DV[nonmem_data$evid == 0][i],
         bloq = nonmem_data$bloq[nonmem_data$evid == 0][i],
         lloq = nonmem_data$lloq[nonmem_data$evid == 0][i],
         ID = nonmem_data$ID[nonmem_data$evid == 0][i],
         time = nonmem_data$time[nonmem_data$evid == 0][i],
         cmt = nonmem_data$cmt[nonmem_data$evid == 0][i]) %>% 
  rowwise() %>% 
  mutate(lower_for_bloq_sim = min(ipred.lower, lloq),
         upper_for_bloq_sim = min(ipred.upper, lloq),
         DV = case_when(bloq == 0 ~ DV, 
                        bloq == 1 ~ runif(1, lower_for_bloq_sim, 
                                          upper_for_bloq_sim),
                        .default = NA_real_)) %>% 
  ungroup() %>% 
  select(-contains("for_bloq_sim"))

tmp <- ggplot(prior_preds_summary) +
  ggforce::facet_wrap_paginate(~ID + cmt, labeller = label_both,
                               nrow = 3, ncol = 4,
                               page = 1, scales = "free_y")

for(i in 1:ggforce::n_pages(tmp)){
  print(ggplot() +
          geom_ribbon(data = prior_preds_summary,
                      mapping = aes(x = time, ymin = dv_ppc.lower,
                                    ymax = dv_ppc.upper, group = ID),
                      fill = "blue", alpha = 0.25, show.legend = FALSE) +
          geom_point(data = prior_preds_summary, # %>% filter(bloq == 0),
                     mapping = aes(x = time, y = DV, group = ID, color = factor(bloq)),
                     size = 2, show.legend = FALSE) +
          scale_color_manual(values = c("0" = "black", "1" = "green")) +
          scale_y_continuous(name = "Conc. or Response",
                             limits = c(NA, NA),
                             trans = "identity") +
          scale_x_continuous(name = "Time (h)",
                             breaks = seq(0, max(nonmem_data$time), by = 24),
                             labels = seq(0, max(nonmem_data$time), by = 24),
                             limits = c(0, max(nonmem_data$time))) +
          theme_bw() +
          theme(axis.text = element_text(size = 14, face = "bold"),
                axis.title = element_text(size = 18, face = "bold"),
                legend.position = "bottom") +
          ggforce::facet_wrap_paginate(~ ID + cmt, labeller = label_both,
                                       nrow = 3, ncol = 4,
                                       page = i, scales = "free"))
  
}

cols_pk <- nonmem_data %>% 
  filter(evid == 0, bloq == 0) %>% 
  mutate(row_num = row_number()) %>% 
  filter(cmt == 2) %>% 
  pull(row_num)

cols_pd <- nonmem_data %>% 
  filter(evid == 0, bloq == 0) %>% 
  mutate(row_num = row_number()) %>% 
  filter(cmt == 3) %>% 
  pull(row_num)

(ppc_dens_overlay_pk <- 
    ppc_dens_overlay(nonmem_data %>% 
                       filter(evid == 0, bloq == 0, cmt == 2) %>% 
                       pull(DV), 
                     prior_df %>%
                       select(contains("epred[")) %>% 
                       select(-(nonmem_data %>% 
                                  filter(evid == 0) %>% 
                                  mutate(i = row_number()) %>% 
                                  filter(bloq == 1) %>% 
                                  pull(i))) %>% 
                       select(all_of(cols_pk)) %>% 
                       slice_sample(n = 1000) %>% 
                       as.matrix()) +
    scale_x_continuous(trans = "log10",
                       limits = c(min(nonmem_data %>% 
                                        filter(evid == 0, bloq == 0, cmt == 2) %>% 
                                        pull(DV))/1000, NA)))

(ppc_dens_overlay_pd <- 
    ppc_dens_overlay(nonmem_data %>% 
                       filter(evid == 0, bloq == 0, cmt == 3) %>% 
                       pull(DV), 
                     prior_df %>%
                       select(contains("epred[")) %>% 
                       select(-(nonmem_data %>% 
                                  filter(evid == 0) %>% 
                                  mutate(i = row_number()) %>% 
                                  filter(bloq == 1) %>% 
                                  pull(i))) %>% 
                       select(all_of(cols_pd)) %>% 
                       slice_sample(n = 1000) %>% 
                       as.matrix()) +
    scale_x_continuous(trans = "log10",
                       limits = c(min(nonmem_data %>% 
                                        filter(evid == 0, bloq == 0, cmt == 3) %>% 
                                        pull(DV))/1000, NA)))
