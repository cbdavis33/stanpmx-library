rm(list = ls())
cat("\014")

# library(trelliscopejs)
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
                  t_2 = 168)

model <- cmdstan_model(
  "depot_1cmt_linear_ir2/Stan/Predict/depot_1cmt_exp_ir2_exp_predict_observed_subjects.stan")

preds <- model$generate_quantities(fit %>%
                                     as_draws_array() %>%
                                     thin_draws(10),
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
         cmt = new_data$cmt[i]) %>% 
  left_join(nonmem_data %>% 
              filter(mdv == 0) %>%
              select(ID, time, cmt, DV), 
            by = c("ID", "time", "cmt")) %>% 
  mutate(cmt = case_when(cmt == 2 ~ "PK",
                         cmt == 3 ~ "PD",
                         TRUE ~ NA_character_) %>% 
           factor(levels = c("PK", "PD"))) %>% 
  select(ID, time, everything(), -i) %>% 
  filter(!is.na(cmt))

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
  facet_trelliscope(~ ID + cmt, nrow = 3, ncol = 4, scales = "free")


tmp <- ggplot(post_preds_summary) +
  ggforce::facet_wrap_paginate(~ID + cmt, labeller = label_both,
                               nrow = 3, ncol = 4,
                               page = 1, scales = "free_y")

for(i in 1:ggforce::n_pages(tmp)){
  print(ggplot() +
          geom_ribbon(data = post_preds_summary, 
                      mapping = aes(x = time, ymin = ipred.lower, 
                                    ymax = ipred.upper, group = ID),
                      fill = "blue", alpha = 0.5, show.legend = FALSE) +
          geom_ribbon(data = post_preds_summary, 
                      mapping = aes(x = time, ymin = dv.lower, 
                                    ymax = dv.upper, group = ID),
                      fill = "blue", alpha = 0.25, show.legend = FALSE) +
          geom_line(data = post_preds_summary, 
                    mapping = aes(x = time, y = ipred, 
                                  group = ID),
                    linetype = 1, linewidth = 1.15) +
          geom_line(data = post_preds_summary, 
                    mapping = aes(x = time, y = pred, 
                                  group = ID),
                    linetype = 2, linewidth = 1.05) +
          geom_point(data = post_preds_summary %>% 
                       filter(!is.na(DV)), 
                     mapping = aes(x = time, y = DV, group = ID),
                     size = 2, color = "red", show.legend = FALSE) +
          scale_y_continuous(name = latex2exp::TeX("$Drug\\;Conc.\\;(\\mu g/mL)$"),
                             limits = c(NA, NA),
                             trans = "identity") +
          scale_x_continuous(name = "Time (h)",
                             breaks = seq(0, max(nonmem_data$time), by = 14),
                             labels = seq(0, max(nonmem_data$time), by = 14),
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
est_ind <- preds_df %>%
  spread_draws(c(CL, VC, KA, 
                 auc_ss, c_max, t_max, t_half,
                 KIN, KOUT, IC50, t_min, r_min)[ID]) %>% 
  mean_qi() %>% 
  select(ID, CL, VC, KA, 
         auc_ss, c_max, t_max, t_half,
         KIN, KOUT, IC50, t_min, r_min) %>% 
  inner_join(post_preds_summary %>% 
               filter(time == 168, cmt == "PK") %>% 
               select(ID, c_trough = "ipred") %>% 
               distinct(),
             by = "ID")

est_ind
