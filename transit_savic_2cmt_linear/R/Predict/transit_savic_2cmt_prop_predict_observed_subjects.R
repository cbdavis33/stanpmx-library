rm(list = ls())
cat("\014")

library(trelliscopejs)
library(cmdstanr)
library(tidybayes)
library(posterior)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

fit <- read_rds("transit_savic_2cmt_linear/Stan/Fits/transit_savic_2cmt_prop.rds")

nonmem_data <- read_csv("transit_savic_2cmt_linear/Data/transit_savic_2cmt_prop.csv",
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

new_data_to_simulate <- nonmem_data %>% 
  # filter(evid == 0) %>%
  group_by(ID) %>% 
  slice(c(1, n())) %>% 
  expand(time = seq(time[1], time[2], by = 0.5)) %>%
  # expand(time = seq(time[1], time[2], by = 2)) %>%
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

new_data_dose <- new_data %>% 
  filter(evid == 1)

n_dose <- new_data_dose %>% 
  nrow()

subj_start_dose <- new_data_dose %>% 
  mutate(row_num = 1:n()) %>% 
  group_by(ID) %>% 
  slice_head(n = 1) %>%
  ungroup() %>% 
  select(row_num) %>% 
  deframe()

subj_end_dose <- c(subj_start_dose[-1] - 1, n_dose) 

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
                  n_dose = n_dose,
                  dosetime = new_data_dose$time,
                  doseamt = new_data_dose$amt,
                  subj_start_dose = subj_start_dose,
                  subj_end_dose = subj_end_dose,
                  t_1 = 0,
                  t_2 = 24)

model <- cmdstan_model(
  "transit_savic_2cmt_linear/Stan/Predict/transit_savic_2cmt_prop_predict_observed_subjects.stan")

preds <- model$generate_quantities(fit,
                                   data = stan_data,
                                   parallel_chains = 4,
                                   seed = 1234) 

# preds <- model$generate_quantities(fit$draws() %>%
#                                      thin_draws(100),
#                                    data = stan_data,
#                                    parallel_chains = 4,
#                                    seed = 1234)

preds_df <- preds$draws(format = "draws_df")

rm(list = setdiff(ls(), c("preds_df", "new_data", "nonmem_data")))

post_preds_summary <- preds_df %>%
  spread_draws(pred[i], ipred[i], dv[i]) %>%
  mean_qi(pred, ipred, dv) %>% 
  mutate(ID = new_data$ID[i],
         time = new_data$time[i]) %>% 
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
  scale_x_continuous(name = "Time (h)",
                     breaks = seq(0, 216, by = 24),
                     labels = seq(0, 216, by = 24),
                     limits = c(0, NA)) +
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
          geom_line(aes(y = ipred), linetype = 1, size = 1.15) +
          geom_line(aes(y = pred), linetype = 2, size = 1.05) +
          geom_point(aes(y = DV), size = 2, show.legend = FALSE, 
                     color = "red") +
          scale_y_continuous(name = latex2exp::TeX("Drug Conc. $(\\mu g/mL)$"),
                             limits = c(NA, NA),
                             trans = "identity") +
          scale_x_continuous(name = "Time (h)",
                             breaks = seq(0, 168, by = 12),
                             labels = seq(0, 168, by = 12),
                             limits = c(0, 48)) +
          theme_bw() +
          theme(axis.text = element_text(size = 14, face = "bold"),
                axis.title = element_text(size = 18, face = "bold"),
                legend.position = "bottom") +
          ggforce::facet_wrap_paginate(~ ID, labeller = label_both,
                                       nrow = 2, ncol = 2,
                                       page = i, scales = "free"))
  
}

## Individual estimates (posterior mean)
est_ind <- preds_df %>%
  spread_draws(c(CL, VC, Q, VP, KA, NTR, MTT, 
               auc_ss, c_max, t_max, 
               t_half_alpha, t_half_terminal)[ID]) %>% 
  mean_qi() %>% 
  select(ID, CL, VC, Q, VP, KA, NTR, MTT, 
         auc_ss, c_max, t_max, starts_with("t_half")) %>% 
  inner_join(post_preds_summary %>% 
               filter(time == 48) %>% 
               select(ID, c_trough = "ipred") %>% 
               distinct(),
             by = "ID")

est_ind
