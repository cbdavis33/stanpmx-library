rm(list = ls())
cat("\014")

library(cmdstanr)
library(tidybayes)
library(posterior)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

fit <- read_rds("depot_1cmt_linear_ir2/Stan/Fits/depot_1cmt_exp_ir2_exp.rds")

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

new_data_to_simulate_pk <- nonmem_data %>% 
  # filter(evid == 0) %>%
  group_by(ID) %>% 
  slice(c(1, n())) %>% 
  # expand(time = seq(time[1], time[2], by = 0.5)) %>%
  expand(time = seq(time[1], time[2], by = 1)) %>%
  ungroup() %>% 
  mutate(amt = 0,
         evid = 2,
         rate = 0,
         addl = 0,
         ii = 0,
         cmt = 2,
         mdv = 1, 
         ss = 0,
         lloq = 1,
         bloq = 0,
         DV = NA_real_,
         c = NA_character_) %>% 
  select(c, ID, time, everything())

new_data_to_simulate_pd <- nonmem_data %>% 
  # filter(evid == 0) %>%
  group_by(ID) %>% 
  slice(c(1, n())) %>% 
  # expand(time = seq(time[1], time[2], by = 0.5)) %>%
  expand(time = seq(time[1], time[2], by = 1)) %>%
  ungroup() %>% 
  mutate(amt = 0,
         evid = 2,
         rate = 0,
         addl = 0,
         ii = 0,
         cmt = 3,
         mdv = 1, 
         ss = 0,
         lloq = 1,
         bloq = 0,
         DV = NA_real_,
         c = NA_character_) %>% 
  select(c, ID, time, everything())

observed_data_to_simulate <- nonmem_data %>% 
  mutate(evid = if_else(evid == 1, 1, 2))

new_data <- new_data_to_simulate_pk %>% 
  bind_rows(new_data_to_simulate_pd, 
            observed_data_to_simulate) %>% 
  arrange(ID, time, evid, cmt) %>% 
  distinct(ID, time, cmt, .keep_all = TRUE) %>% 
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
                  t_1 = 144,
                  t_2 = 168,
                  want_auc_cmax = 1)

model <- cmdstan_model(
  "depot_1cmt_linear_ir2/Stan/Predict/depot_1cmt_exp_ir2_exp_predict_observed_subjects.stan")

preds <- model$generate_quantities(fit$draws() %>%
                                     thin_draws(10),
                                   data = stan_data,
                                   parallel_chains = 4,
                                   seed = 1234) 

preds_df <- preds$draws(format = "draws_df")

rm(list = setdiff(ls(), c("preds_df", "new_data", "nonmem_data", "stan_data")))
gc()

## This is cleaner, but can be really slow for bigger datasets, so I'll do them
## separately here
# post_preds <- preds_df %>%
#   spread_draws(c(epred_stan, epred, ipred, dv)[i])

post_epred <- preds_df %>%
  spread_draws(epred[i]) 

post_epred_stan <- preds_df %>%
  spread_draws(epred_stan[i]) 

post_ipred <- preds_df %>%
  spread_draws(ipred[i]) 

post_dv <- preds_df %>%
  spread_draws(dv[i]) 

post_preds <- post_epred %>% 
  inner_join(post_epred_stan) %>% 
  inner_join(post_ipred) %>% 
  inner_join(post_dv) 

post_preds_summary <- post_preds %>% 
  mean_qi(epred_stan, epred, ipred, dv, .width = 0.95) %>%
  mutate(ID = new_data$ID[i],
         time = new_data$time[i],
         cmt = new_data$cmt[i]) %>% 
  left_join(nonmem_data %>% 
              filter(evid == 0) %>%
              select(ID, time, cmt, bloq, lloq, DV), 
            by = c("ID", "time", "cmt")) %>%
  mutate(cmt = case_when(cmt == 2 ~ "PK",
                         cmt == 3 ~ "PD",
                         .default = NA_character_) %>% 
           factor(levels = c("PK", "PD"))) %>% 
  filter(!is.na(cmt)) %>% 
  rowwise() %>% 
  # Note: for BLOQ values, I impute a value for DV only for plotting purposes.
  # It can easily be filtered out or done differently
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
         time = new_data$time[i],
         cmt = new_data$cmt[i]) %>% 
  left_join(nonmem_data %>% 
              filter(evid == 0) %>%
              select(ID, time, cmt, DV), 
            by = c("ID", "time", "cmt"),
            relationship = "many-to-many") %>%
  select(ID, time, everything())

tmp <- ggplot(post_preds_summary) +
  ggforce::facet_wrap_paginate(~ID + cmt, labeller = label_both,
                               nrow = 3, ncol = 4,
                               page = 1, scales = "free_y")

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
          ggforce::facet_wrap_paginate(~ ID + cmt, labeller = label_both,
                                       nrow = 3, ncol = 4,
                                       page = i, scales = "free"))
  
  
}

for(i in 1:ggforce::n_pages(tmp)){
  print(ggplot(data = post_preds %>% 
                 sample_draws(ndraws = 50) %>% 
                 mutate(cmt = case_when(cmt == 2 ~ "PK",
                                        cmt == 3 ~ "PD",
                                        .default = NA_character_) %>% 
                          factor(levels = c("PK", "PD"))) %>% 
                 filter(!is.na(cmt))) +
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
          ggforce::facet_wrap_paginate(~ ID + cmt, labeller = label_both,
                                       nrow = 3, ncol = 4,
                                       page = i, scales = "free"))
  
}

## Individual estimates (posterior mean)
est_ind <- if(stan_data$want_auc_cmax){
  preds_df %>%
    spread_draws(c(CL, VC, KA, 
                   auc_ss, c_max, t_max, t_half,
                   KIN, KOUT, IC50, t_max_pd, r_max)[ID]) %>% 
    mean_qi() %>% 
    select(ID, CL, VC, KA, 
           auc_ss, c_max, t_max, t_half,
           KIN, KOUT, IC50, t_max_pd, r_max) %>% 
    inner_join(post_preds_summary %>% 
                 filter(time == 168, cmt == "PK") %>% 
                 select(ID, c_trough = "ipred") %>% 
                 distinct(),
               by = "ID")
}else{
  preds_df %>%
    spread_draws(c(CL, VC, KA, t_half,
                   KIN, KOUT, IC50)[ID]) %>% 
    mean_qi() %>% 
    select(ID, CL, VC, KA, t_half,
           KIN, KOUT, IC50) %>% 
    inner_join(post_preds_summary %>% 
                 filter(time == 168, cmt == "PK") %>% 
                 select(ID, c_trough = "ipred") %>% 
                 distinct(),
               by = "ID")
}

est_ind %>% 
  as_tibble()

