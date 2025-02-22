rm(list = ls())
cat("\014")

library(trelliscopejs)
library(cmdstanr)
library(tidybayes)
library(posterior)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

locf <- TRUE

file_path <- if_else(isTRUE(locf), 
                     "iv_2cmt_linear_covariates_time_varying/Generic/Data/iv_2cmt_prop_covariates_time_varying_generic_locf.csv",
                     "iv_2cmt_linear_covariates_time_varying/Generic/Data/iv_2cmt_prop_covariates_time_varying_generic_nocb.csv")


fit <- if(isTRUE(locf)){
  read_rds("iv_2cmt_linear_covariates_time_varying/Generic/Stan/Fits/prop_fixed_allometric_locf.rds")
}else{
  read_rds("iv_2cmt_linear_covariates_time_varying/Generic/Stan/Fits/prop_fixed_allometric_nocb.rds")
}

stan_data <- jsonlite::read_json(
  str_c("iv_2cmt_linear_covariates_time_varying/Generic/Stan/Fits/Stan_Data/", 
        "prop_fixed_allometric_locf.json")) %>% 
  map(function(x) if(is.list(x)) as_vector(x) else x)

nonmem_data <- read_csv(
  file_path,
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

tmp <- nonmem_data %>% 
  # filter(evid == 0) %>%
  group_by(ID) %>% 
  slice(c(1, n())) %>% 
  mutate(time_sequence = list(seq(time[1], time[2], by = 2))) %>%
  slice(1) %>% 
  ungroup() %>% 
  select(ID, time_sequence) 

new_data_to_simulate <- nonmem_data %>% 
  inner_join(tmp, by = "ID") %>%
  group_by(ID) %>% 
  # mutate(time_ahead = lead(time, default = last(time))) %>% 
  mutate(time_ahead = lead(time, default = last(time) + 1)) %>% 
  rowwise() %>% 
  mutate(time_new = list(time_sequence[time_sequence >= time & time_sequence < time_ahead])) %>% 
  unnest(cols = c(time_new)) %>% 
  ungroup() %>% 
  select(-time_sequence, -time_ahead, -time) %>% 
  rename(time = time_new) %>% 
  filter(evid == 0) 

observed_data_to_simulate <- nonmem_data %>% 
  mutate(evid = if_else(evid == 1, 1, 2)) 

new_data <- new_data_to_simulate %>% 
  bind_rows(observed_data_to_simulate) %>% 
  arrange(ID, time, evid) %>% 
  distinct(ID, time, .keep_all = TRUE) %>% 
  select(-DV, -lloq, -bloq) %>% 
  mutate(cmt = 2) # analytical solution

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

wt <- new_data %>% 
  pull(wt)

race_asian <- new_data %>% 
  pull(race_asian)

egfr <- new_data %>% 
  pull(egfr)

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
                  wt = wt,
                  race_asian = race_asian,
                  egfr = egfr,
                  t_1 = 144,
                  t_2 = 168,
                  want_auc_cmax = 0,
                  theta_cl_wt = stan_data$theta_cl_wt,
                  theta_vc_wt = stan_data$theta_vc_wt,
                  theta_q_wt = stan_data$theta_q_wt,
                  theta_vp_wt = stan_data$theta_vp_wt)

model <- cmdstan_model(
  "iv_2cmt_linear_covariates_time_varying/Generic/Stan/Predict/iv_2cmt_prop_predict_observed_subjects_covariates_fixed_allometric_time_varying.stan")

preds <- model$generate_quantities(fit$draws() %>%
                                     thin_draws(100),
                                   data = stan_data,
                                   parallel_chains = 4,
                                   seed = 1234)

preds_df <- preds$draws(format = "draws_df")

rm(list = setdiff(ls(), c("preds_df", "new_data", "nonmem_data")))

post_preds_summary <- preds_df %>%
  spread_draws(c(pred, ipred, dv)[i]) %>%
  mean_qi(pred, ipred, dv) %>% 
  mutate(ID = new_data$ID[i],
         time = new_data$time[i],
         wt = new_data$wt[i],
         cmppi = new_data$cmppi[i],
         egfr = new_data$egfr[i]) %>% 
  left_join(nonmem_data %>% 
              filter(mdv == 0) %>%
              select(ID, time, DV), 
            by = c("ID", "time")) %>% 
  select(ID, time, everything(), -i)

ggplot(post_preds_summary, aes(x = time, group = ID)) +
  geom_ribbon(aes(ymin = dv.lower, ymax = dv.upper),
              fill = "blue", alpha = 0.25, show.legend = FALSE) +
  geom_ribbon(aes(ymin = ipred.lower, ymax = ipred.upper),
              fill = "blue", alpha = 0.5, show.legend = FALSE) +
  geom_line(aes(y = ipred), linetype = 1, size = 1.15) +
  geom_line(aes(y = pred), linetype = 2, size = 1.05) +
  geom_point(aes(y = DV), size = 2, show.legend = FALSE, 
             color = "red") +
  scale_y_continuous(name = latex2exp::TeX("Drug Conc. $(\\mu g/mL)$"),
                     trans = "log10",
                     limits = c(NA, NA)) +
  scale_x_continuous(name = "Time (d)",
                     breaks = seq(0, max(nonmem_data$time), by = 14),
                     labels = seq(0, max(nonmem_data$time), by = 14),
                     limits = c(0, max(nonmem_data$time))) +
  theme_bw() +
  theme(axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "bottom") +
  facet_trelliscope(~ ID, nrow = 3, ncol = 4, scales = "free")


tmp <- ggplot(post_preds_summary, aes(x = time, group = ID)) +
  geom_line(aes(y = ipred), linetype = 1, size = 1.15) +
  ggforce::facet_wrap_paginate(~ID, labeller = label_both,
                               nrow = 2, ncol = 2,
                               page = 1, scales = "free")

for(i in 1:ggforce::n_pages(tmp)){
  print(ggplot(post_preds_summary, aes(x = time, group = ID)) +
          geom_ribbon(aes(ymin = dv.lower, ymax = dv.upper),
                      fill = "blue", alpha = 0.25, show.legend = FALSE) +
          geom_ribbon(aes(ymin = ipred.lower, ymax = ipred.upper),
                      fill = "blue", alpha = 0.5, show.legend = FALSE) +
          geom_line(aes(y = ipred), linetype = 1, linewidth = 1.15) +
          geom_line(aes(y = pred), linetype = 2, linewidth = 1.05) +
          geom_point(aes(y = DV), size = 2, show.legend = FALSE, 
                     color = "red") +
          scale_y_log10(name = latex2exp::TeX("Drug Conc. $(\\mu g/mL)$"),
                        limits = c(NA, NA)) +
          scale_x_continuous(name = "Time (d)",
                             breaks = seq(0, max(nonmem_data$time), by = 14),
                             labels = seq(0, max(nonmem_data$time), by = 14),
                             limits = c(0, max(nonmem_data$time))) +
          theme_bw() +
          theme(axis.text = element_text(size = 14, face = "bold"),
                axis.title = element_text(size = 18, face = "bold"),
                legend.position = "bottom") +
          ggforce::facet_wrap_paginate(~ ID, labeller = label_both,
                                       nrow = 2, ncol = 2,
                                       page = i, scales = "free"))
  
}


## TODO: work from here down
# ## Individual estimates (posterior mean)
# 
# if(stan_data$want_cmax_auc){
#   est_ind <- preds_df %>%
#     spread_draws(c(auc_ss, c_max, t_max)[ID]) %>% 
#     mean_qi() %>% 
#     select(ID, auc_ss, c_max, t_max) %>% 
#     inner_join(post_preds_summary %>% 
#                  filter(time == 168) %>% 
#                  select(ID, c_trough = "ipred") %>% 
#                  distinct(),
#                by = "ID")
# }
# 
# 
# 
# est <- preds_df %>%
#   spread_draws(c(CL, VC, KA)[i]) %>% 
#   mean_qi() %>% 
#   select(i, CL, VC, KA)
# 
# 
# 
# # est_ind <- preds_df %>%
# #   spread_draws(c(CL, VC, KA, 
# #                  auc_ss, c_max, t_max, t_half)[ID]) %>% 
# #   mean_qi() %>% 
# #   select(ID, CL, VC, KA, 
# #          auc_ss, c_max, t_max, t_half) %>% 
# #   inner_join(post_preds_summary %>% 
# #                filter(time == 168) %>% 
# #                select(ID, c_trough = "ipred") %>% 
# #                distinct(),
# #              by = "ID")
# # 
# # est_ind
