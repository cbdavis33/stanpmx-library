rm(list = ls())
cat("\014")

library(trelliscopejs)
library(cmdstanr)
library(tidybayes)
library(posterior)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

fit <- read_rds("iv_2cmt_linear/Stan/Fits/iv_2cmt_exp.rds")

nonmem_data <- read_csv("iv_2cmt_linear/Data/iv_2cmt_exp.csv",
                        na = ".") %>% 
  rename_all(tolower) %>% 
  rename(ID = "id",
         DV = "dv") %>% 
  mutate(DV = if_else(is.na(DV), 5555555, DV),    # This value can be anything except NA. It'll be indexed away 
         bloq = if_else(is.na(bloq), -999, bloq), # This value can be anything except NA. It'll be indexed away 
         cmt = 1)

# For this example, let's simulate 50 mg, 100 mg, 200 mg, 400 mg, 800 mg Q2W and
# 400 mg, 800 mg Q4W
dosing_data_q2w <- mrgsolve::expand.ev(addl = 5, ii = 14, 
                                       cmt = 1, amt = c(50, 100, 200, 400, 800), 
                                       ss = 0, tinf = 1/24, evid = 1) %>%
  as_tibble() %>% 
  mutate(ss = 0) %>% 
  select(ID, time, everything())

dosing_data_q4w <- mrgsolve::expand.ev(addl = 2, ii = 28, 
                                       cmt = 1, amt = c(400, 800), ss = 0, 
                                       tinf = 1/24, evid = 1) %>%
  as_tibble() %>% 
  mutate(ss = 0,
         ID = ID + 5) %>% 
  select(ID, time, everything())

dosing_data <- dosing_data_q2w %>% 
  bind_rows(dosing_data_q4w)

end_inf_times <- dosing_data %>% 
  realize_addl() %>% 
  mutate(time = time + tinf) %>% 
  select(-tinf) %>% 
  pull(time)

t1 <- dosing_data %>% 
  realize_addl() %>% 
  select(time) %>% 
  distinct() %>% 
  deframe()

times_new <- tibble(time = sort(unique(c(t1, end_inf_times, 
                                         seq(0, 84, by = 12/24)))))

new_data <- bind_rows(replicate(max(dosing_data$ID), times_new, 
                                simplify = FALSE)) %>% 
  mutate(ID = rep(1:max(dosing_data$ID), each = nrow(times_new)),
         amt = 0, 
         evid = 0, 
         rate = 0, 
         addl = 0, 
         ii = 0, 
         cmt = 1, 
         mdv = 1, 
         ss = 0) %>%
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
                  t_1 = 70,
                  t_2 = 84,
                  want_auc_cmax = 1)

model <- cmdstan_model(
  "iv_2cmt_linear/Stan/Predict/iv_2cmt_exp_predict_new_subjects.stan")

preds <- model$generate_quantities(fit$draws() %>%
                                     thin_draws(1),
                                   data = stan_data,
                                   parallel_chains = 4,
                                   seed = 1234)

preds_df <- preds$draws(format = "draws_df")

regimens <- c(str_c(c(50, 100, 200, 400, 800), " mg Q2W"), 
              str_c(c(400, 800), " mg Q4W"))

post_preds_summary <- preds_df %>%
  spread_draws(epred_stan[i], epred[i]) %>%
  median_qi(epred_stan, epred) %>%
  mutate(ID = new_data$ID[i],
         time = new_data$time[i]) %>%
  select(ID, time, everything(), -i) %>% 
  mutate(regimen = factor(regimens[ID],
                          levels = regimens))

tmp <- ggplot(post_preds_summary, aes(x = time, group = ID)) +
  geom_line(aes(y = epred_stan), linetype = 1, size = 1.15) +
  geom_line(aes(y = epred), linetype = 2, size = 1.05) +
  ggforce::facet_wrap_paginate(~ ID, 
                               labeller = label_both,
                               nrow = 3, ncol = 3,
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
          scale_x_continuous(name = "Time (d)",
                             breaks = seq(0, 84, by = 7),
                             labels = seq(0, 84, by = 7),
                             limits = c(0, 84)) +
          # scale_x_continuous(name = "Time (d)") +
          theme_bw() +
          theme(axis.text = element_text(size = 14, face = "bold"),
                axis.title = element_text(size = 18, face = "bold"),
                legend.position = "bottom") +
          ggforce::facet_wrap_paginate(~ regimen,
                                       nrow = 3, ncol = 3,
                                       page = i))
  
}


data <- read_csv("iv_2cmt_linear/Data/iv_2cmt_exp.csv", na = ".") %>% 
  rename_all(tolower) %>% 
  rename(ID = "id",
         DV = "dv") %>% 
  group_by(ID) %>% 
  mutate(Dose = max(amt, na.rm = TRUE),
         mdv = evid) %>% 
  ungroup() %>% 
  mutate(regimen = str_c(Dose, " mg Q2W"))


post_preds_summary %>% 
  filter(regimen %in% str_c(c(50, 100, 200, 400), " mg Q2W")) %>%
  ggplot(aes(x = time, group = ID)) +
  geom_lineribbon(aes(y = epred, ymin = epred.lower, ymax = epred.upper),
                  fill = "blue", color = "blue", linewidth = 1.05,
                  alpha = 0.25, show.legend = FALSE) +
  geom_lineribbon(aes(y = epred_stan, ymin = epred_stan.lower, ymax = epred_stan.upper),
                  fill = "blue", color = "blue", linewidth = 1.15,
                  alpha = 0.5, show.legend = FALSE) +
  geom_point(data = data %>% 
               mutate(regimen = factor(regimen, levels = regimens)) %>% 
               filter(regimen %in% str_c(c(50, 100, 200, 400), " mg Q2W"), 
                      mdv == 0),
             mapping = aes(x = time, y = DV), color = "red", 
             inherit.aes = FALSE) +
  scale_y_continuous(name = latex2exp::TeX("Drug Conc. $(\\mu g/mL)$"),
                     trans = "log10",
                     limits = c(NA, NA)) +
  scale_x_continuous(name = "Time (d)",
                     breaks = seq(0, max(new_data$time), by = 14),
                     labels = seq(0, max(new_data$time), by = 14),
                     limits = c(0, max(new_data$time))) +
  theme_bw() +
  theme(axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "bottom") +
  coord_cartesian(xlim = c(0, 84)) +
  facet_wrap(~ regimen, scales = "free_y")


## Individual estimates (posterior median)
est_ind <- if(stan_data$want_auc_cmax){
  preds_df %>%
    spread_draws(c(CL, VC, Q, VP, 
                   auc_ss, t_half_alpha, t_half_terminal)[ID]) %>% 
    median_qi() %>% 
    select(ID, CL, VC, Q, VP, 
           auc_ss, t_half_alpha, t_half_terminal) %>% 
    inner_join(post_preds_summary %>% 
                 filter(time %in% c(70 + 1/24, 84)) %>% 
                 mutate(time = round(time, 0)) %>% 
                 select(ID, time, epred_stan) %>%
                 pivot_wider(values_from = epred_stan, names_from = time, 
                             names_prefix = "time") %>%
                 rename(c_trough = time84,
                        c_max = time70) %>% 
                 distinct(),
               by = "ID")
}else{
  preds_df %>%
    spread_draws(c(CL, VC, Q, VP)[ID]) %>% 
    median_qi() %>% 
    select(ID, CL, VC, Q, VP) %>% 
    inner_join(post_preds_summary %>% 
                 filter(time %in% c(70 + 1/24, 84)) %>% 
                 mutate(time = round(time, 0)) %>% 
                 select(ID, time, epred_stan) %>%
                 pivot_wider(values_from = epred_stan, names_from = time, 
                             names_prefix = "time") %>%
                 rename(c_trough = time84,
                        c_max = time70) %>% 
                 distinct(),
               by = "ID")
}

est_ind




