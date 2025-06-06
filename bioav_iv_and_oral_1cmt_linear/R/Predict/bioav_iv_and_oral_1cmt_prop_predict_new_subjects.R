rm(list = ls())
cat("\014")

library(cmdstanr)
library(mrgsolve)
library(tidybayes)
library(posterior)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

fit <- read_rds(
  "bioav_iv_and_oral_1cmt_linear/Stan/Fits/bioav_iv_and_oral_1cmt_prop.rds")

nonmem_data <- read_csv(
  "bioav_iv_and_oral_1cmt_linear/Data/bioav_iv_and_oral_1cmt_prop.csv",
  na = ".") %>% 
  rename_all(tolower) %>% 
  rename(ID = "id",
         DV = "dv") %>% 
  mutate(DV = if_else(is.na(DV), 5555555, DV),    # This value can be anything except NA. It'll be indexed away 
         bloq = if_else(is.na(bloq), -999, bloq)) # This value can be anything except NA. It'll be indexed away 

# For this example, let's simulate 500 mg oral, 250 mg IV, and 1000 mg IV 
# loading dose followed by 500 mg oral
dosing_data_oral <- expand.ev(addl = 6, ii = 24, cmt = 1, amt = 500, ss = 0, 
                              tinf = 0, evid = 1) %>%
  as_tibble() %>% 
  select(ID, time, everything()) %>% 
  mutate(type = "Oral",
         cohort = 1)

dosing_data_iv <- expand.ev(addl = 6, ii = 24, cmt = 2, amt = 250, ss = 0, 
                            tinf = 4, evid = 1) %>%
  as_tibble() %>% 
  select(ID, time, everything()) %>% 
  mutate(cohort = 2,
         ID = 2,
         type = "IV")

dosing_data_iv_and_oral_iv <- expand.ev(addl = 0, ii = 24, cmt = 2, amt = 1000, 
                                        ss = 0, tinf = 4, evid = 1) %>% 
  as.ev()

dosing_data_iv_and_oral_oral <- expand.ev(addl = 5, ii = 24, cmt = 1, amt = 500, 
                                          ss = 0, tinf = 0, evid = 1) %>% 
  as.ev()

dosing_data_iv_and_oral <- ev_seq(dosing_data_iv_and_oral_iv,
                                  dosing_data_iv_and_oral_oral) %>% 
  as_tibble() %>% 
  arrange(ID, time) %>% 
  mutate(cohort = 3,
         ID = 3, 
         type = "IV and Oral")

dosing_data <- dosing_data_oral %>% 
  bind_rows(dosing_data_iv,
            dosing_data_iv_and_oral) %>% 
  mutate(cohort = factor(cohort)) 

t1 <- dosing_data %>% 
  realize_addl() %>% 
  select(time) %>% 
  distinct() %>% 
  deframe()

times_new <- tibble(time = sort(unique(c(t1, seq(0, 168, by = 1)))))

new_data <- bind_rows(replicate(max(dosing_data$ID), times_new, 
                                simplify = FALSE)) %>% 
  mutate(ID = rep(1:max(dosing_data$ID), each = nrow(times_new)),
         amt = 0, 
         evid = 0, 
         rate = 0, 
         addl = 0, 
         ii = 0, 
         cmt = 2, 
         mdv = 1, 
         ss = 0,
         cohort = factor(ID),
         type = case_when(ID == 1 ~ "Oral",
                          ID == 2 ~ "IV",
                          ID == 3 ~ "IV and Oral",
                          TRUE ~ NA_character_)) %>%
  filter(time != 0) %>% 
  select(ID, time, everything()) %>% 
  bind_rows(dosing_data) %>% 
  arrange(ID, time)

# number of individuals in the original dataset
n_subjects <- fit$metadata()$stan_variable_sizes$Z[2]

n_subjects_new <- new_data %>%  # number of new individuals
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
                  n_subjects_new = n_subjects_new,
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
  "bioav_iv_and_oral_1cmt_linear/Stan/Predict/bioav_iv_and_oral_1cmt_prop_predict_new_subjects.stan")

preds <- model$generate_quantities(fit$draws() %>%
                                     thin_draws(1),
                                   data = stan_data,
                                   parallel_chains = 4,
                                   seed = 1234)

preds_df <- preds$draws(format = "draws_df")

post_preds_summary <- preds_df %>%
  spread_draws(epred_stan[i], epred[i]) %>%
  median_qi(epred_stan, epred) %>%
  mutate(ID = new_data$ID[i],
         time = new_data$time[i],
         type = new_data$type[i]) %>%
  select(ID, time, everything(), -i)

tmp <- ggplot(post_preds_summary, aes(x = time, group = ID)) +
  geom_line(aes(y = epred_stan), linetype = 1, size = 1.15) +
  geom_line(aes(y = epred), linetype = 2, size = 1.05) +
  ggforce::facet_wrap_paginate(~ ID, 
                               labeller = label_both,
                               nrow = 2, ncol = 3,
                               page = 1)

for(i in 1:ggforce::n_pages(tmp)){
  print(ggplot(post_preds_summary, aes(x = time, group = ID)) +
          geom_ribbon(aes(ymin = epred.lower, ymax = epred.upper),
                      fill = "blue", alpha = 0.25, show.legend = FALSE) +
          geom_ribbon(aes(ymin = epred.lower, ymax = epred.upper),
                      fill = "blue", alpha = 0.5, show.legend = FALSE) +
          geom_line(aes(y = epred_stan), linetype = 1, size = 1.15) +
          geom_line(aes(y = epred), linetype = 2, size = 1.05) +
          scale_y_continuous(name = latex2exp::TeX("Drug Conc. $(\\mu g/mL)$"),
                             trans = "log10",
                             limits = c(NA, NA)) +
          scale_x_continuous(name = "Time (h)",
                             breaks = seq(0, 168, by = 24),
                             labels = seq(0, 168, by = 24),
                             limits = c(0, 168)) +
          # scale_x_continuous(name = "Time (d)") +
          theme_bw() +
          theme(axis.text = element_text(size = 14, face = "bold"),
                axis.title = element_text(size = 18, face = "bold"),
                legend.position = "bottom") +
          ggforce::facet_wrap_paginate(~ type,
                                       nrow = 1, ncol = 3,
                                       page = i))
  
}

# C_max for the IV groups - the trick in the ODEs doesn't work with IV. Just 
# make sure you've simulated the end of infusion to get C_max
preds_df %>%
  spread_draws(epred_stan[i]) %>% 
  ungroup() %>% 
  mutate(ID = new_data$ID[i],
         time = new_data$time[i],
         type = new_data$type[i]) %>%
  select(ID, time, everything(), -i) %>%
  group_by(ID, .draw) %>% 
  summarize(c_max = max(epred_stan)) %>% 
  filter() %>% 
  group_by(ID) %>% 
  median_qi(c_max) %>% 
  ungroup()

## Individual estimates (posterior mean). Note - for the IV subjects, Cmax and
## Tmax won't be right. You'll need to simulate at the end of infusion and pull
## from that time point
est_ind <- if(stan_data$want_auc_cmax){
  preds_df %>%
    spread_draws(CL[ID], VC[ID], KA[ID], BIOAV[ID],
                 auc_ss[ID], c_max[ID], t_max[ID], t_half[ID]) %>% 
    median_qi() %>% 
    select(ID, CL, VC, KA, BIOAV,
           auc_ss, c_max, t_max, t_half) %>% 
    inner_join(post_preds_summary %>% 
                 filter(time == 168) %>% 
                 select(ID, c_trough = "epred_stan", type) %>% 
                 distinct(),
               by = "ID")
}else{
  preds_df %>%
    spread_draws(CL[ID], VC[ID], KA[ID], BIOAV[ID], t_half[ID]) %>% 
    median_qi() %>% 
    select(ID, CL, VC, KA, BIOAV, t_half) %>% 
    inner_join(post_preds_summary %>% 
                 filter(time == 168) %>% 
                 select(ID, c_trough = "epred_stan", type) %>% 
                 distinct(),
               by = "ID")
}

est_ind

