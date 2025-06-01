rm(list = ls())
cat("\014")

library(bayesplot)
library(vpc)
library(cmdstanr)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

nonmem_data <- read_csv("iv_1cmt_linear/Data/iv_1cmt_prop.csv",
                        na = ".") %>% 
  rename_all(tolower) %>% 
  rename(ID = "id",
         DV = "dv") %>% 
  mutate(DV = if_else(is.na(DV), 5555555, DV),    # This value can be anything except NA. It'll be indexed away 
         bloq = if_else(is.na(bloq), -999, bloq), # This value can be anything except NA. It'll be indexed away 
         cmt = 1)

## Summary of BLOQ values
nonmem_data %>%
  filter(evid == 0) %>%
  group_by(ID) %>%
  summarize(lloq = unique(lloq),
            n_obs = n(),
            n_bloq = sum(bloq)) %>%
  filter(n_bloq > 0)

# (p1 <- ggplot(nonmem_data %>%
#                 group_by(ID) %>%
#                 mutate(Dose = factor(max(amt, na.rm = TRUE))) %>%
#                 ungroup() %>%
#                 filter(mdv == 0)) +
#     geom_line(mapping = aes(x = time, y = DV, group = ID, color = Dose)) +
#     geom_point(mapping = aes(x = time, y = DV, group = ID, color = Dose)) +
#     scale_color_discrete(name = "Dose (mg)") +
#     scale_y_continuous(name = latex2exp::TeX("$Drug\\;Conc.\\;(\\mu g/mL)$"),
#                        limits = c(NA, NA),
#                        trans = "log10") +
#     scale_x_continuous(name = "Time (d)",
#                        breaks = seq(0, 216, by = 24),
#                        labels = seq(0, 216/24, by = 24/24),
#                        limits = c(0, NA)) +
#     theme_bw(18) +
#     theme(axis.text = element_text(size = 14, face = "bold"),
#           axis.title = element_text(size = 18, face = "bold"),
#           axis.line = element_line(linewidth = 2),
#           legend.position = "bottom"))

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
                  location_tvcl = 0.5,
                  location_tvvc = 4,
                  scale_tvcl = 1,
                  scale_tvvc = 1,
                  scale_omega_cl = 0.4,
                  scale_omega_vc = 0.4,
                  lkj_df_omega = 2,
                  scale_sigma_p = 0.5,
                  prior_only = 1,
                  no_gq_predictions = 0) 

model <- cmdstan_model("iv_1cmt_linear/Stan/Fit/iv_1cmt_prop.stan",
                       cpp_options = list(stan_threads = TRUE))

prior <- model$sample(data = stan_data,
                      seed = 112358,
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
                                  omega = abs(rnorm(2, 0, c(scale_omega_cl, 
                                                            scale_omega_vc))),
                                  sigma_p = abs(rnorm(1, 0, scale_sigma_p)))))

prior_df <- prior$draws(format = "draws_df") %>%  
  drop_na()

prior_preds_summary <- prior_df %>%
  spread_draws(c(ipred, epred, dv_ppc)[i]) %>% 
  mean_qi(ipred, epred, dv_ppc, .width = 0.95) %>%
  mutate(DV = nonmem_data$DV[nonmem_data$evid == 0][i],
         bloq = nonmem_data$bloq[nonmem_data$evid == 0][i],
         lloq = nonmem_data$lloq[nonmem_data$evid == 0][i],
         ID = nonmem_data$ID[nonmem_data$evid == 0][i],
         time = nonmem_data$time[nonmem_data$evid == 0][i]) %>% 
  rowwise() %>% 
  mutate(lower_for_bloq_sim = min(ipred.lower, lloq),
         upper_for_bloq_sim = min(ipred.upper, lloq),
         DV = case_when(bloq == 0 ~ DV, 
                        bloq == 1 ~ runif(1, lower_for_bloq_sim, 
                                          upper_for_bloq_sim),
                        .default = NA_real_)) %>% 
  ungroup() %>% 
  select(-contains("for_bloq_sim"))

tmp <- ggplot(prior_preds_summary, aes(x = time, group = ID)) +
  ggforce::facet_wrap_paginate(~ID, labeller = label_both,
                               nrow = 2, ncol = 2,
                               page = 1, scales = "free")

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
          scale_y_continuous(name = latex2exp::TeX("$Drug\\;Conc.\\;(\\mu g/mL)$"),
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
          ggforce::facet_wrap_paginate(~ ID, labeller = label_both,
                                       nrow = 2, ncol = 2,
                                       page = i, scales = "free"))
  
}

ppc_dens_overlay(nonmem_data %>% 
                   filter(evid == 0, bloq == 0) %>% 
                   pull(DV), 
                 prior_df %>%
                   select(contains("epred[")) %>% 
                   select(-(nonmem_data %>% 
                              filter(evid == 0) %>% 
                              mutate(i = row_number()) %>% 
                              filter(bloq == 1) %>% 
                              pull(i))) %>% 
                   slice_sample(n = 1000) %>% 
                   as.matrix()) +
  scale_x_continuous(trans = "log10",
                     limits = c(min(nonmem_data %>% 
                                      filter(evid == 0, bloq == 0) %>% 
                                      pull(DV))/1000, NA)) 

## Prior VPCs

preds_col <- prior_df %>% 
  spread_draws(epred[i]) %>% 
  summarize(epred_mean = mean(epred)) %>% 
  ungroup() %>% 
  bind_cols(nonmem_data %>% 
              select(ID, time, evid) %>% 
              filter(evid == 0)) %>% 
  select(ID, time, epred_mean, -i)

sim <- prior_df %>%
  spread_draws(epred_stan[i], epred[i]) %>% 
  ungroup() %>% 
  mutate(ID = nonmem_data[nonmem_data$evid == 0, ]$ID[i],
         time = nonmem_data[nonmem_data$evid == 0, ]$time[i],
         mdv = nonmem_data[nonmem_data$evid == 0, ]$mdv[i],
         evid = nonmem_data[nonmem_data$evid == 0, ]$evid[i],
         bloq = nonmem_data[nonmem_data$evid == 0, ]$bloq[i]) %>% 
  # mutate(epred = if_else(bloq == 1, 0, epred)) %>% 
  select(ID, sim = ".draw", time, epred_stan, epred) %>% 
  arrange(ID, sim, time) %>% 
  left_join(preds_col, by = c("ID", "time"))

obs <- nonmem_data %>% 
  select(ID, time, dv = "DV", amt, evid, mdv, bloq, lloq) %>% 
  filter(evid == 0) %>% 
  # mutate(dv = if_else(bloq == 1, 0, dv)) %>%
  mutate(dv = if_else(bloq == 1, lloq, dv)) %>% 
  left_join(preds_col, by = c("ID", "time"))

(p_vpc <- vpc(sim = sim, obs = obs,
              sim_cols = list(idv = "time",
                              dv = "epred",
                              id = "ID",
                              pred = "epred_mean",
                              sim = "sim"),
              obs_cols = list(dv = "dv",
                              idv = "time",
                              id = "ID",
                              pred = "epred_mean"),
              show = list(obs_dv = TRUE, 
                          obs_ci = FALSE,
                          pi = FALSE,
                          pi_as_area = FALSE,
                          pi_ci = TRUE,
                          obs_median = TRUE,
                          sim_median = TRUE,
                          sim_median_ci = FALSE),
              log_y = TRUE,
              lloq = 1) +
    # scale_y_continuous(name = "Drug Conc. (ug/mL)",
    #                    trans = "identity") +
    scale_x_continuous(name = "Time (h)",
                       breaks = seq(0, 168, by = 24),
                       labels = seq(0, 168, by = 24)) +
    theme_bw())

(p_pcvpc <- vpc(sim = sim, obs = obs,
                sim_cols = list(idv = "time",
                                dv = "epred",
                                id = "ID",
                                pred = "epred_mean",
                                sim = "sim"),
                obs_cols = list(dv = "dv",
                                idv = "time",
                                id = "ID",
                                pred = "epred_mean"),
                pred_corr = TRUE,
                show = list(obs_dv = TRUE, 
                            obs_ci = FALSE,
                            pi = FALSE,
                            pi_as_area = FALSE,
                            pi_ci = TRUE,
                            obs_median = TRUE,
                            sim_median = TRUE,
                            sim_median_ci = FALSE),
                log_y = TRUE) +
    # scale_y_continuous(name = "Drug Conc. (ug/mL)",
    #                    trans = "log10") +
    scale_x_continuous(name = "Time (h)",
                       breaks = seq(0, 168, by = 24),
                       labels = seq(0, 168, by = 24)) +
    theme_bw())

