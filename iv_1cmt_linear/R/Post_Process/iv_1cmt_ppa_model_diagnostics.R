rm(list = ls())
cat("\014")

library(gganimate)
library(loo)
library(vpc)
library(npde)
library(cmdstanr)
library(tidybayes)
library(bayesplot)
library(patchwork)
library(posterior)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

nonmem_data <- read_csv("iv_1cmt_linear/Data/iv_1cmt_ppa.csv",
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

## Read in fit
fit <- read_rds("iv_1cmt_linear/Stan/Fits/iv_1cmt_ppa.rds")

draws_df <- fit$draws(format = "draws_df")

## DV vs. PRED/EPRED/IPRED/DV_ppc. All point estimates are posterior means. That
## can be changed to medians if you prefer
# PRED = f(TVs, x, eta = 0) - the "typical" subject. Not a fan of this
# EPRED_stan = f(TVs, x, eta = eta_new) - This is the "EPRED" in the Staniverse.
#   Simulates eta_new ~ multi_normal(0, Omega)
# EPRED = EPRED_stan + error - This is the "EPRED" from NONMEM
# IPRED = f(TVs, x, eta = eta_i) - Uses eta values for the observed subjects
# DV_ppc = IPRED + error
# Note: for BLOQ values in this section, I impute a value for DV only for 
# plotting purposes. It can easily be filtered out or done differently

post_preds_summary <- draws_df %>%
  spread_draws(c(pred, epred_stan, epred, ipred, dv_ppc)[i]) %>% 
  mean_qi(pred, epred_stan, epred, ipred, dv_ppc, .width = 0.95) %>%
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

# post_preds_summary %>%
#   filter(bloq == 0) %>%
#   mutate(in_epred_stan_interval = between(DV, epred_stan.lower, epred_stan.upper),
#          in_epred_interval = between(DV, epred.lower, epred.upper),
#          in_ipred_interval = between(DV, ipred.lower, ipred.upper),
#          in_ppc_interval = between(DV, dv_ppc.lower, dv_ppc.upper)) %>%
#   summarize(across(c(in_epred_stan_interval, in_epred_interval,
#                      in_ipred_interval, in_ppc_interval), mean))

p_dv_vs_pred <- ggplot(post_preds_summary, # %>% filter(bloq == 0), 
                       aes(x = pred, y = DV, color = factor(bloq))) +
  geom_point() +
  scale_color_manual(guide = "none",
                     values = c("0" = "black", "1" = "green")) +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, color = "blue", linewidth = 1.5) +
  geom_smooth(color = "red", se = FALSE, linewidth = 1.5) +
  xlab("Population Predictions (PRED)") +
  ylab("Observed Concentration") +
  theme(axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        axis.line = element_line(linewidth = 2))

p_dv_vs_epred <- ggplot(post_preds_summary, # %>% filter(bloq == 0),  
                        aes(x = epred, y = DV, color = factor(bloq))) +
  geom_point() +
  scale_color_manual(guide = "none",
                     values = c("0" = "black", "1" = "green")) +
  geom_errorbarh(aes(xmin = epred.lower, xmax = epred.upper)) +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, color = "blue", linewidth = 1.5) +
  geom_smooth(color = "red", se = FALSE, linewidth = 1.5) +
  xlab("Population Predictions (EPRED)") +
  ylab("Observed Concentration") +
  theme(axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        axis.line = element_line(linewidth = 2))

p_dv_vs_ipred <- ggplot(post_preds_summary, # %>% filter(bloq == 0),  
                        aes(x = pred, y = DV, color = factor(bloq))) +
  geom_point() +
  scale_color_manual(guide = "none",
                     values = c("0" = "black", "1" = "green")) +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, color = "blue", linewidth = 1.5) +
  geom_smooth(color = "red", se = FALSE, linewidth = 1.5) +
  xlab("Individual Predictions (IPRED)") +
  ylab("Observed Concentration") +
  theme(axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        axis.line = element_line(linewidth = 2))

p_dv_vs_ppc <- ggplot(post_preds_summary, # %>% filter(bloq == 0), 
                      aes(x = dv_ppc, y = DV, color = factor(bloq))) +
  geom_point() +
  scale_color_manual(guide = "none",
                     values = c("0" = "black", "1" = "green")) +
  geom_errorbarh(aes(xmin = dv_ppc.lower, xmax = dv_ppc.upper)) +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, color = "blue", linewidth = 1.5) +
  geom_smooth(color = "red", se = FALSE, linewidth = 1.5) +
  xlab("Individual Predictions (PPC)") +
  ylab("Observed Concentration") +
  theme(axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        axis.line = element_line(linewidth = 2))

p_dv_vs_pred +
  p_dv_vs_epred +
  p_dv_vs_ipred +
  p_dv_vs_ppc

(p_dv_vs_pred +
    scale_y_log10() +
    scale_x_log10()) +
  (p_dv_vs_epred +
     scale_y_log10() +
     scale_x_log10()) +
  (p_dv_vs_ipred +
     scale_y_log10() +
     scale_x_log10()) +
  (p_dv_vs_ppc +
     scale_y_log10() +
     scale_x_log10())


## Individual Plots

tmp <- ggplot(post_preds_summary) +
  ggforce::facet_wrap_paginate(~ID, labeller = label_both,
                               nrow = 3, ncol = 4,
                               page = 1, scales = "free_y")

for(i in 1:ggforce::n_pages(tmp)){
  print(ggplot() +
          geom_lineribbon(data = post_preds_summary,
                          mapping = aes(x = time, y = ipred, ymin = ipred.lower,
                                        ymax = ipred.upper, group = ID),
                          fill = "blue", color = "blue", linewidth = 2,
                          alpha = 0.5, show.legend = FALSE) +
          geom_lineribbon(data = post_preds_summary,
                          mapping = aes(x = time, ymin = dv_ppc.lower,
                                        ymax = dv_ppc.upper, group = ID),
                          fill = "blue", linewidth = 2, 
                          alpha = 0.25, show.legend = FALSE) +
          geom_lineribbon(data = post_preds_summary,
                          mapping = aes(x = time, y = epred, ymin = epred.lower,
                                        ymax = epred.upper, group = ID),
                          fill = "red", color = "red", linewidth = 2,
                          alpha = 0.25, show.legend = FALSE) +
          # geom_line(data = post_preds_summary,
          #                 mapping = aes(x = time, y = pred, group = ID),
          #                 linewidth = 2,
          #                 show.legend = FALSE) +
          geom_point(data = post_preds_summary, # %>% filter(bloq == 0),
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
                                       nrow = 3, ncol = 4,
                                       page = i, scales = "free"))
  
}

# Those plots look a little weird, since they only take into account observation
# timepoints. We can make predictions at a denser grid to show the multiple 
# dosing and make better-looking plots

end_inf_data <- nonmem_data %>% 
  filter(evid == 1) %>% 
  realize_addl() %>% 
  mutate(tinf = amt/rate,
         time = time + tinf) %>% 
  select(-tinf)

new_data_to_simulate <- nonmem_data %>% 
  # filter(evid == 0) %>%
  group_by(ID) %>% 
  slice(c(1, n())) %>% 
  # expand(time = seq(time[1], time[2], by = 0.5)) %>%
  expand(time = seq(time[1], time[2], by = 12/24)) %>%
  ungroup() %>% 
  bind_rows(end_inf_data) %>% 
  mutate(amt = 0,
         evid = 2,
         rate = 0,
         addl = 0,
         ii = 0,
         cmt = 1,
         mdv = 1, 
         ss = 0,
         lloq = 1,
         bloq = 0,
         DV = NA_real_,
         c = NA_character_) %>% 
  select(c, ID, time, everything())

observed_data_to_simulate <- nonmem_data %>% 
  mutate(evid = if_else(evid == 1, 1, 2))

new_data <- new_data_to_simulate %>% 
  bind_rows(observed_data_to_simulate) %>% 
  arrange(ID, time, evid) %>% 
  distinct(ID, time, .keep_all = TRUE) %>% 
  select(-DV, -lloq, -bloq)

n_subjects <- new_data %>%  # number of individuals
  distinct(ID) %>% 
  count() %>% 
  deframe()

n_time_new <- nrow(new_data) # total number of time points at which to predict

subj_start <- new_data %>% 
  mutate(row_num = 1:n()) %>% 
  group_by(ID) %>% 
  slice_head(n = 1) %>%
  ungroup() %>% 
  select(row_num) %>% 
  deframe()

subj_end <- c(subj_start[-1] - 1, n_time_new)  

stan_data <- list(n_subjects = n_subjects,
                  n_time_new = n_time_new,
                  time = new_data$time,
                  amt = new_data$amt,
                  cmt = new_data$cmt,
                  evid = new_data$evid,
                  rate = new_data$rate,
                  ii = new_data$ii,
                  addl = new_data$addl,
                  ss = new_data$ss,
                  subj_start = subj_start,
                  subj_end = subj_end,
                  t_1 = 70,
                  t_2 = 84,
                  want_auc_cmax = 0)

model <- cmdstan_model(
  "iv_1cmt_linear/Stan/Predict/iv_1cmt_ppa_predict_observed_subjects.stan")

# This object can get VERY big with the dense grid, so I sometimes thin a little
# if I'm just plotting. It won't make much difference. For actual inference, 
# don't thin. To do no thinning, do thin_draws(1)
preds <- model$generate_quantities(fit$draws() %>%
                                     thin_draws(10),
                                   data = stan_data,
                                   parallel_chains = 4,
                                   seed = 1234)

preds_df <- preds$draws(format = "draws_df")

# rm(list = setdiff(ls(), c("preds_df", "new_data", "nonmem_data")))
# gc()

post_preds <- preds_df %>%
  spread_draws(c(pred, epred_stan, epred, ipred, dv)[i]) 

post_preds_summary <- post_preds %>% 
  mean_qi(pred, epred_stan, epred, ipred, dv, .width = 0.95) %>%
  mutate(ID = new_data$ID[i],
         time = new_data$time[i]) %>% 
  left_join(nonmem_data %>% 
              filter(evid == 0) %>%
              select(ID, time, bloq, lloq, DV), 
            by = c("ID", "time")) %>%
  rowwise() %>% 
  # Again, simulate a DV for the BLOQ observations, purely for plotting purposes
  mutate(lower_for_bloq_sim = if_else(is.na(bloq), 0, min(ipred.lower, lloq)),
         upper_for_bloq_sim = if_else(is.na(bloq), 1, min(ipred.upper, lloq)),
         DV = case_when(bloq == 0 ~ DV, 
                        bloq == 1 ~ runif(1, lower_for_bloq_sim, 
                                          upper_for_bloq_sim),
                        .default = NA_real_)) %>% 
  ungroup() %>% 
  select(ID, time, everything(), -i, -contains("for_bloq_sim"))

post_preds <- post_preds %>%
  ungroup() %>% 
  mutate(ID = new_data$ID[i],
         time = new_data$time[i]) %>% 
  left_join(nonmem_data %>% 
              filter(evid == 0) %>%
              select(ID, time, DV), 
            by = c("ID", "time")) %>%
  select(ID, time, everything())

tmp <- ggplot(post_preds_summary, aes(x = time, group = ID)) +
  ggforce::facet_wrap_paginate(~ID, labeller = label_both,
                               nrow = 2, ncol = 2,
                               page = 1, scales = "free")

for(i in 1:ggforce::n_pages(tmp)){
  print(ggplot() +
          geom_lineribbon(data = post_preds_summary,
                          mapping = aes(x = time, y = epred, ymin = epred.lower,
                                        ymax = epred.upper, group = ID),
                          fill = "red", color = "red", linewidth = 2,
                          alpha = 0.25, show.legend = FALSE) +
          geom_lineribbon(data = post_preds_summary,
                          mapping = aes(x = time, ymin = dv.lower,
                                        ymax = dv.upper, group = ID),
                          fill = "blue", linewidth = 2, 
                          alpha = 0.25, show.legend = FALSE) +
          geom_lineribbon(data = post_preds_summary,
                          mapping = aes(x = time, y = ipred, ymin = ipred.lower,
                                        ymax = ipred.upper, group = ID),
                          fill = "blue", color = "blue", linewidth = 2,
                          alpha = 0.5, show.legend = FALSE) +
          geom_point(data = post_preds_summary, # %>% filter(bloq == 0),
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


for(i in 1:ggforce::n_pages(tmp)){
  print(ggplot(data = post_preds %>% 
                 sample_draws(ndraws = 50)) +
          geom_line(mapping = aes(x = time, y = epred_stan, group = .draw),
                    color = "red", linewidth = 0.5,
                    alpha = 0.15, show.legend = FALSE) +
          # geom_line(mapping = aes(x = time, y = dv, group = .draw),
          #           color = "green", linewidth = 0.5,
          #           alpha = 0.15, show.legend = FALSE) +
          geom_line(mapping = aes(x = time, y = ipred, group = .draw),
                    color = "blue", linewidth = 0.5,
                    alpha = 0.15, show.legend = FALSE) +
          geom_point(data = post_preds_summary, # %>% filter(bloq == 0),
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

# Check some observations to see how they're covered by their posterior
zxc <- draws_df %>%
  spread_draws(c(pred, epred_stan, epred, ipred, dv_ppc)[i]) %>% 
  mutate(DV = nonmem_data$DV[nonmem_data$evid == 0][i],
         bloq = nonmem_data$bloq[nonmem_data$evid == 0][i],
         lloq = nonmem_data$lloq[nonmem_data$evid == 0][i],
         ID = nonmem_data$ID[nonmem_data$evid == 0][i],
         time = nonmem_data$time[nonmem_data$evid == 0][i],
         DV = na_if(DV, 5555555)) %>% 
  left_join(post_preds_summary %>% 
              filter(bloq == 1) %>% 
              select(ID, time, DV),
            by = c("ID", "time")) %>% 
  mutate(DV = coalesce(DV.x, DV.y), .keep = "unused")

zxc %>% 
  sample_draws(ndraws = 16, draw = "i") %>%
  pivot_longer(cols = c(dv_ppc, epred)) %>%
  ggplot(aes(x = value, group = name)) +
  geom_histogram(aes(fill = name), alpha = 0.5, position = position_identity()) + 
  scale_fill_manual(name = "Posterior Draw Type",
                    values = c("dv_ppc" = "blue", "epred" = "red"),
                    breaks = c("dv_ppc", "epred"),
                    labels = c("PPC", "EPRED")) +
  geom_vline(aes(xintercept = DV, color = factor(bloq)),
             show.legend = FALSE) +
  scale_color_manual(values = c("0" = "black", "1" = "green")) +
  scale_x_continuous(name = "Posterior Draws") +
  theme_bw() +
  theme(legend.position = "bottom") +
  facet_wrap(~i, ncol = 4, scales = "free",
             labeller = labeller(i = function(x) str_c("Observation #", x))) 


## VPCs. Also see ../VPC/iv_1cmt_ppa_vpc.R for an alternate but equivalent
## way to do VPCs

preds_col <- draws_df %>% 
  spread_draws(epred[i]) %>% 
  summarize(epred_mean = mean(epred)) %>% 
  ungroup() %>% 
  bind_cols(nonmem_data %>% 
              filter(evid == 0) %>% 
              select(ID, time, evid)) %>% 
  filter(evid == 0) %>% 
  select(ID, time, epred_mean, -i)

sim <- draws_df %>%
  spread_draws(epred[i]) %>% 
  ungroup() %>% 
  mutate(ID = nonmem_data$ID[nonmem_data$evid == 0][i],
         time = nonmem_data$time[nonmem_data$evid == 0][i],
         mdv = nonmem_data$mdv[nonmem_data$evid == 0][i],
         evid = nonmem_data$evid[nonmem_data$evid == 0][i],
         bloq = nonmem_data$bloq[nonmem_data$evid == 0][i]) %>% 
  filter(evid == 0) %>%
  # mutate(epred = if_else(bloq == 1, 0, epred)) %>%
  select(ID, sim = ".draw", time, epred) %>% 
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
                              pred = "pred",
                              sim = "sim"),
              obs_cols = list(dv = "dv",
                              idv = "time",
                              id = "ID",
                              pred = "pred"),
              show = list(obs_dv = TRUE, 
                          obs_ci = TRUE,
                          pi = TRUE,
                          pi_as_area = FALSE,
                          pi_ci = TRUE,
                          obs_median = TRUE,
                          sim_median = TRUE,
                          sim_median_ci = TRUE),
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
                            obs_ci = TRUE,
                            pi = TRUE,
                            pi_as_area = FALSE,
                            pi_ci = TRUE,
                            obs_median = TRUE,
                            sim_median = TRUE,
                            sim_median_ci = TRUE),
                log_y = TRUE) +
    # scale_y_continuous(name = "Drug Conc. (ug/mL)",
    #                    trans = "log10") +
    scale_x_continuous(name = "Time (h)",
                       breaks = seq(0, 168, by = 24),
                       labels = seq(0, 168, by = 24)) +
    theme_bw())

p_vpc + 
  p_pcvpc


## NPDEs. Also see ../NPDE/iv_1cmt_ppa_npde.R for an alternate but 
## equivalent way to do NPDEs

data_obs <- nonmem_data %>% 
  filter(evid == 0) %>% 
  mutate(DV = if_else(bloq == 1, lloq, DV)) 

simdata <- draws_df %>%
  spread_draws(epred[i]) %>% 
  ungroup() %>% 
  mutate(ID = data_obs$ID[data_obs$evid == 0][i],
         time = data_obs$time[data_obs$evid == 0][i]) %>% 
  left_join(data_obs %>% 
              filter(evid == 0) %>%
              select(ID, time, DV), 
            by = c("ID", "time")) %>%
  select(ID, time, everything()) %>% 
  arrange(.draw, ID, time)

my_npde <- autonpde(data_obs %>% 
                      filter(evid == 0) %>% 
                      select(ID, amt, time, DV), 
                    simdata %>% 
                      select(ID, time, epred), 
                    1, 3, 4, boolsave = FALSE, cens.method = "cdf")

results_npde <- data_obs %>% 
  select(ID, lloq, bloq, time, DV) %>% 
  bind_cols(my_npde@results@res %>% 
              tibble() %>% 
              rename(ewres = ydobs) %>% 
              select(ypred, ycomp, ewres, npde))

# NPDE histogram
results_npde %>% 
  ggplot() + 
  geom_histogram(aes(x = npde, y = after_stat(density)), fill = "gray") +
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1),
                color = "red", size = 1) +
  scale_x_continuous(name = "NPDE") +
  theme_bw() +
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())

# NPDE QQ-plot
results_npde %>% 
  ggplot(aes(sample = npde)) + 
  geom_qq() +
  geom_abline() +
  theme_bw()

# NPDE vs. Population Prediction
results_npde %>% 
  ggplot(aes(x = ypred, y = npde)) +
  geom_point() +
  scale_x_continuous(name = "Population Prediction (EPRED)") +
  scale_y_continuous(name = "NPDE") +
  geom_hline(yintercept = 0, color = "blue", linewidth = 1.5) +
  geom_hline(aes(yintercept = qnorm(0.025)), 
             linetype = "dashed", color = "blue",
             linewidth = 1.25) +
  geom_hline(aes(yintercept = qnorm(0.975)), 
             linetype = "dashed", color = "blue",
             linewidth = 1.25) +
  geom_smooth(se = FALSE, color = "red", linewidth = 1.25) +
  theme_bw()

# NPDE vs. Time
results_npde %>% 
  ggplot(aes(x = time, y = npde)) +
  geom_point() +
  scale_x_continuous(name = "Time (h)") +
  scale_y_continuous(name = "NPDE") +
  geom_hline(yintercept = 0, color = "blue", linewidth = 1.5) +
  geom_hline(aes(yintercept = qnorm(0.025)), 
             linetype = "dashed", color = "blue",
             linewidth = 1.25) +
  geom_hline(aes(yintercept = qnorm(0.975)), 
             linetype = "dashed", color = "blue",
             linewidth = 1.25) +
  geom_smooth(se = FALSE, color = "red", linewidth = 1.25) +
  theme_bw()


## Residuals
### EWRES
# EWRES density plot
results_npde %>% 
  ggplot() + 
  geom_density(aes(x = ewres), fill = "gray") + 
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1),
                color = "red", size = 1) +
  scale_x_continuous(name = "EWRES") +
  theme_bw() +
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())

# EWRES QQ-plot
results_npde %>% 
  ggplot(aes(sample = ewres)) + 
  geom_qq() +
  geom_abline() +
  theme_bw()

# EWRES vs. Population Prediction
results_npde %>% 
  ggplot(aes(x = ypred, y = ewres)) +
  geom_point() +
  scale_x_continuous(name = "Population Prediction (EPRED)") +
  scale_y_continuous(name = "EWRES") +
  geom_hline(yintercept = 0, color = "blue", linewidth = 1.5) +
  geom_hline(aes(yintercept = qnorm(0.025)), 
             linetype = "dashed", color = "blue",
             linewidth = 1.25) +
  geom_hline(aes(yintercept = qnorm(0.975)), 
             linetype = "dashed", color = "blue",
             linewidth = 1.25) +
  geom_smooth(se = FALSE, color = "red", linewidth = 1.25) +
  theme_bw()

# EWRES vs. Time
results_npde %>% 
  ggplot(aes(x = time, y = ewres)) +
  geom_point() +
  scale_x_continuous(name = "Time (h)") +
  scale_y_continuous(name = "EWRES") +
  geom_hline(yintercept = 0, color = "blue", linewidth = 1.5) +
  geom_hline(aes(yintercept = qnorm(0.025)), 
             linetype = "dashed", color = "blue",
             linewidth = 1.25) +
  geom_hline(aes(yintercept = qnorm(0.975)), 
             linetype = "dashed", color = "blue",
             linewidth = 1.25) +
  geom_smooth(se = FALSE, color = "red", linewidth = 1.25) +
  theme_bw()


### IWRES
## Residuals and epsilon shrinkage = 1 - sd(IWRES) - I think there are two ways 
## to look at this, both of which have value:
##   1) by draw 
##   2) using posterior summaries

# Look at it by draw first
residuals <- draws_df %>%
  spread_draws(c(iwres, ipred)[i]) %>% 
  mutate(time = nonmem_data$time[nonmem_data$evid == 0][i],
         ID = nonmem_data$ID[nonmem_data$evid == 0][i],
         bloq = nonmem_data$bloq[nonmem_data$evid == 0][i])

(shrinkage_eps <- residuals %>% 	
    filter(bloq == 0) %>% 
    group_by(.draw) %>% 	
    summarize(sd_iwres = sd(iwres)) %>% 	
    ungroup() %>% 	
    mutate(shrinkage = 1 - sd_iwres) %>% 	
    mean_qi(shrinkage))

residuals %>% 
  filter(bloq == 0) %>%
  sample_draws(4) %>% 
  ggplot(aes(sample = iwres)) + 
  geom_qq() +
  geom_abline() +
  theme_bw() +
  facet_wrap(~.draw, labeller = label_both)

some_residuals_qq <- residuals %>% 
  filter(bloq == 0) %>%
  sample_draws(100) %>% 
  ggplot(aes(sample = iwres)) + 
  geom_qq() +
  geom_abline() +
  theme_bw() +
  transition_manual(.draw)

animate(some_residuals_qq, nframes = 100, width = 384, height = 384, res = 96, 
        dev = "png", type = "cairo")

residuals %>% 
  filter(bloq == 0) %>%
  sample_draws(4) %>%
  ggplot(aes(x = iwres)) + 
  geom_histogram(aes(y = after_stat(density))) +
  geom_density(color = "blue", linewidth = 1.5) +
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1),
                color = "red", size = 1) +
  theme_bw() +
  xlab("IWRES") +
  facet_wrap(~ .draw, labeller = label_both)

residuals %>% 
  filter(bloq == 0) %>%
  sample_draws(9) %>%
  mutate(qn_lower = qnorm(0.025),
         qn_upper = qnorm(0.975)) %>% 
  ggplot(aes(x = time, y = iwres)) + 
  geom_point() +
  geom_hline(yintercept = 0, color = "blue", linewidth = 1.5) +
  geom_hline(aes(yintercept = qn_lower), linetype = "dashed", color = "blue",
             linewidth = 1.25) +
  geom_hline(aes(yintercept = qn_upper), linetype = "dashed", color = "blue",
             linewidth = 1.25) +
  geom_smooth(se = FALSE, color = "red", linewidth = 1.25) +
  theme_bw() +
  ylab("IWRES") +
  xlab("Time (h)") +
  facet_wrap(~ .draw, labeller = label_both)

residuals %>% 
  filter(bloq == 0) %>%
  sample_draws(9) %>%
  ggplot(aes(x = ipred, y = iwres)) + 
  geom_point() +
  geom_hline(yintercept = 0, color = "blue", linewidth = 1.5) +
  geom_hline(aes(yintercept = qnorm(0.025)), linetype = "dashed", color = "blue",
             linewidth = 1.25) +
  geom_hline(aes(yintercept = qnorm(0.975)), linetype = "dashed", color = "blue",
             linewidth = 1.25) +
  geom_smooth(se = FALSE, color = "red", linewidth = 1.25) +
  theme_bw() +
  ylab("IWRES") +
  xlab("Individual Prediction") +
  facet_wrap(~ .draw, labeller = label_both)

some_residuals_scatter <- residuals %>% 
  filter(bloq == 0) %>%
  sample_draws(30) %>% 
  ggplot(aes(x = time, y = iwres)) + 
  geom_point() +
  geom_hline(yintercept = 0, color = "blue", linewidth = 1.5) +
  geom_hline(aes(yintercept = qnorm(0.025)), linetype = "dashed", color = "blue",
             linewidth = 1.25) +
  geom_hline(aes(yintercept = qnorm(0.975)), linetype = "dashed", color = "blue",
             linewidth = 1.25) +
  geom_smooth(se = FALSE, color = "red", linewidth = 1.25) +
  theme_bw() +
  transition_manual(.draw)

animate(some_residuals_scatter, nframes = 30, width = 384, height = 384, 
        res = 96, dev = "png", type = "cairo", fps = 5)

# |iwres| vs. IPRED
residuals %>% 
  filter(bloq == 0) %>%
  sample_draws(9) %>% 
  mutate(abs_iwres = abs(iwres)) %>% 
  ggplot(aes(x = ipred, y = abs_iwres)) +
  geom_point() + 
  geom_smooth(method = "loess", se = FALSE, color = "red", linewidth = 1.25) +
  theme_bw() +
  ylab(("|IWRES|")) +
  xlab("Individual Prediction") +
  facet_wrap(~ .draw, labeller = label_both)


# Now look at it using posterior point estimates (I use posterior means, but 
# posterior medians would be fine)

iwres_point <- residuals %>%
  filter(bloq == 0) %>%
  group_by(i, time) %>% 
  summarize(iwres_point = mean(iwres)) %>% 
  ungroup()

iwres_point %>% 
  summarize(shrinkage = 1 - sd(iwres_point))

iwres_point %>% 
  ggplot(aes(x = time, y = iwres_point)) +
  theme_bw() +
  geom_point() +
  geom_hline(yintercept = 0, color = "blue", linewidth = 1.5) +
  scale_y_continuous(name = "IWRES") +
  scale_x_continuous(name = "Time (h)",
                     limits = c(NA, NA))

iwres_point %>% 
  ggplot(aes(sample = iwres_point)) + 
  geom_qq() +
  geom_abline() +
  theme_bw()


## Eta-shrinkage
## I think there are three ways to look at this, all of which have value:
##   1) by draw - this is similar to "Enhanced Method for Diagnosing 
##        Pharmacometric Models: Random Sampling from Conditional Distributions"
##        from the Lixoft group. I maybe shouldn't call it shrinkage, but it has 
##        a similar idea but calculate 1 - sd(eta)/omega within each posterior 
##        draw and get a "distribution" of eta-shrinkage
##   2) using posterior summaries, like a posterior mean for eta-hat (analogous 
##        to EBEs) (1 - sd(individual posterior means)/(posterior mean of omega))
##   3) Using the method from "Bayesian Measures of Explained Variance and
##        Pooling in Multilevel (Hierarchical) Models" by Gelman and Pardoe 
##        except using standard deviations instead of variances

# A function to visualize the shrinkage
plot_shrinkage <- function(.variable, std_dev, ind_params, ...){
  
  dots <- list(...)
  
  x_label <- glue::glue("$\\eta_{[.variable]}$", .open = "[", .close = "]")
  
  p_1 <- ggplot() +
    geom_histogram(aes(x = ind_params$ind_params,
                       y = after_stat(density))) +
    stat_function(fun = dnorm, 
                  args = list(mean = 0, sd = std_dev),
                  color = "red", linewidth = 1) +
    scale_x_continuous(name = latex2exp::TeX(x_label),
                       limits = c(qnorm(c(0.001, 0.995), 
                                        0, std_dev))) +
    theme_bw()
  
  
  if(is.numeric(dots$.draw)){
    
    p_1 <- p_1 +
      labs(subtitle = str_c("draw = ", dots$.draw)) +
      theme(plot.subtitle = element_text(hjust = 0.5))
    
  }
  
  return(p_1)
  
}

# This function stacks multiple shrinkage plots (e.g. by draw)
plot_shrinkage_multiple <- function(draw){
  
  data_shrinkage_by_draw %>% 
    filter(.draw == draw) %>% 
    pmap(plot_shrinkage) %>% 
    wrap_plots(nrow = 1)
  
}

# Look at it by draw first
(shrinkage <- draws_df %>%
    gather_draws(`eta_.*`[ID], regex = TRUE) %>% 
    group_by(.draw, .variable) %>% 
    summarize(sd_eta = sd(.value)) %>% 
    ungroup() %>% 
    mutate(.variable = str_remove(.variable, "eta_")) %>% 
    left_join(draws_df %>%
                gather_draws(`omega_.*`, regex = TRUE) %>% 
                ungroup() %>% 
                arrange(.draw) %>% 
                rename(omega = ".value") %>% 
                mutate(.variable = str_remove(.variable, "omega_")),
              by = c(".draw", ".variable")) %>% 
    mutate(shrinkage = 1 - sd_eta/omega) %>% 
    select(.draw, .variable, shrinkage) %>% 
    group_by(.variable) %>% 
    mean_qi())

draws_for_shrinkage <- draws_df %>%
  sample_draws(5)

data_shrinkage_by_draw <- draws_for_shrinkage %>%
  gather_draws(`omega_.*`, regex = TRUE) %>%
  mutate(.variable = str_remove(.variable, "omega_") %>%
           toupper()) %>%
  rename(std_dev = .value) %>%
  ungroup() %>% 
  inner_join(draws_for_shrinkage %>%
               gather_draws(`eta_.*`[ID], regex = TRUE) %>%
               ungroup()  %>% 
               mutate(.variable = str_remove(.variable, "eta_") %>%
                        toupper()) %>%
               rename(ind_params = `.value`) %>%
               select(-ID) %>%
               group_by(.draw) %>%
               nest(ind_params = ind_params) %>%
               ungroup(),
             by = c(".draw", ".chain", ".iteration", ".variable")) %>%
  arrange(.draw, .variable) %>%
  select(-.chain, -.iteration)

map(data_shrinkage_by_draw %>% 
      distinct(.draw) %>% 
      pull(.draw), 
    .f = plot_shrinkage_multiple) %>% 
  wrap_plots(ncol = 1)


# Now look at it using posterior point estimates
draws_df %>%
  gather_draws(`eta_.*`[ID], regex = TRUE) %>% 
  summarize(estimate = mean(.value)) %>% 
  ungroup() %>% 
  group_by(.variable) %>% 
  summarize(std_dev = sd(estimate)) %>% 
  ungroup() %>% 
  mutate(.variable = str_remove(.variable, "eta_")) %>% 
  left_join(draws_df %>%
              gather_draws(`omega_.*`, regex = TRUE) %>%
              summarize(estimate = mean(.value)) %>% 
              ungroup() %>% 
              mutate(.variable = str_remove(.variable, "omega_")),
            by = ".variable") %>% 
  mutate(shrinkage = 1 - std_dev/estimate)

data_shrinkage_with_point_estimates <- draws_df %>% 
  gather_draws(`omega_.*`, regex = TRUE) %>%
  summarize(std_dev = mean(.value)) %>% 
  mutate(.variable = str_remove(.variable, "omega_") %>% 
           toupper()) %>% 
  inner_join(draws_df %>%
               gather_draws(`eta_.*`[ID], regex = TRUE) %>% 
               summarize(ind_params = mean(.value)) %>% 
               ungroup() %>% 
               mutate(.variable = str_remove(.variable, "eta_") %>%
                        toupper()) %>%
               select(-ID) %>% 
               nest(ind_params = ind_params), 
             by = ".variable")

pmap(data_shrinkage_with_point_estimates, 
     plot_shrinkage) %>% 
  wrap_plots(nrow = 1)

# Gelman and Pardoe's method - for the denominator, using the sd of the etas in 
# each draw and then calculating the mean of those sds. If you use the mean of 
# the omegas, it's exactly the same as the point estimate method above
draws_df %>%
  gather_draws(`eta_.*`[ID], regex = TRUE) %>% 
  summarize(estimate = mean(.value)) %>% 
  ungroup() %>% 
  group_by(.variable) %>% 
  summarize(numerator = sd(estimate)) %>% 
  ungroup() %>% 
  inner_join(draws_df %>%
               gather_draws(`eta_.*`[ID], regex = TRUE) %>% 
               group_by(.draw, .variable) %>% 
               summarize(sd_eta = sd(.value)) %>% 
               ungroup() %>% 
               group_by(.variable) %>% 
               summarize(denominator = mean(sd_eta)),
             by = ".variable") %>% 
  mutate(.variable = str_remove(.variable, "eta_"),
         shrinkage = 1 - numerator/denominator) %>% 
  select(.variable, shrinkage)


# And maybe it's just cool to look at individual posterior densities on top of 
# the density of the population parameter
blah <- draws_df %>%
  gather_draws(CL[ID], VC[ID], TVCL, TVVC) %>%
  ungroup() %>%
  mutate(across(c(ID, .variable), as.factor))

ggplot() +
  geom_density(data = blah %>% 
                 filter(is.na(ID)) %>% 
                 mutate(.variable = blah %>% 
                          filter(is.na(ID)) %>% 
                          pull(.variable) %>% 
                          str_remove("TV"),
                        ID = "Population"),
               mapping = aes(x = .value, group = .variable,
                             fill = .variable), color = "black",
               position = "identity") +
  stat_density(data = blah %>% 
                 filter(!is.na(ID)), 
               mapping = aes(x = .value, group = ID, color = ID),
               geom = "line", position = "identity") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~.variable, ncol = 1, scales = "free") +
  theme(legend.position = "bottom") +
  scale_fill_manual(name = "Population Parameter",
                    values = c("red", "blue", "green")) +
  guides(color = "none") +
  ggtitle("Individual Parameter Posterior Densities") 

## Check Leave-One-Out Cross-Validation - this can also be used for model 
## comparison, but here if a bunch of Pareto-k diagnostic values are high, then
## it's a sign the model isn't so good
fit_loo <- fit$loo()
fit_loo
plot(fit_loo, label_points = TRUE)
