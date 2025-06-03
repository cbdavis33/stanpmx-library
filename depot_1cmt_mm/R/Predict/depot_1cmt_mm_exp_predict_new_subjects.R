rm(list = ls())
cat("\014")

library(cmdstanr)
library(mrgsolve)
library(tidybayes)
library(posterior)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

fit <- read_rds("depot_1cmt_mm/Stan/Fits/depot_1cmt_mm_exp.rds")

nonmem_data <- read_csv("depot_1cmt_mm/Data/depot_1cmt_mm_exp.csv",
                        na = ".") %>% 
  rename_all(tolower) %>% 
  rename(ID = "id",
         DV = "dv") %>% 
  mutate(DV = if_else(is.na(DV), 5555555, DV),    # This value can be anything except NA. It'll be indexed away 
         bloq = if_else(is.na(bloq), -999, bloq)) # This value can be anything except NA. It'll be indexed away 

# For this example, let's simulate 50 mg, 100 mg, 200 mg, 400 mg, 600 mg, 800 mg 
dosing_data <- mrgsolve::expand.ev(addl = 6, ii = 24, cmt = 1, 
                                   amt = c(50, 100, 200, 400, 600, 800), 
                                   tinf = 0, evid = 1, mdv = 1) %>%
  as_tibble() %>% 
  mutate(ss = 0) %>% 
  select(ID, time, everything())

t1 <- dosing_data %>% 
  realize_addl() %>% 
  select(time) %>% 
  distinct() %>% 
  deframe()

times_new <- tibble(time = sort(unique(c(t1, 0.25, seq(0, 168, by = 0.5)))))

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
                  t_1 = 144,
                  t_2 = 168,
                  want_auc_cmax = 0)

model <- cmdstan_model(
  "depot_1cmt_mm/Stan/Predict/depot_1cmt_mm_exp_predict_new_subjects.stan")

preds <- model$generate_quantities(fit$draws() %>%
                                     thin_draws(1),
                                   data = stan_data,
                                   parallel_chains = 4,
                                   seed = 1234)

preds_df <- preds$draws(format = "draws_df")

regimens <- str_c(c(50, 100, 200, 400, 600, 800), " mg")

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
          ggforce::facet_wrap_paginate(~ regimen,
                                       nrow = 2, ncol = 3,
                                       page = i))
  
}

data <- read_csv("depot_1cmt_mm/Data/depot_1cmt_mm_exp.csv", na = ".") %>% 
  rename_all(tolower) %>% 
  rename(ID = "id",
         DV = "dv") %>% 
  group_by(ID) %>% 
  mutate(Dose = max(amt, na.rm = TRUE),
         mdv = evid) %>% 
  ungroup() %>% 
  mutate(regimen = str_c(Dose, " mg"))


post_preds_summary %>% 
  filter(regimen %in% str_c(c(100, 200, 400, 800), " mg")) %>%
  ggplot(aes(x = time, group = ID)) +
  geom_lineribbon(aes(y = epred, ymin = epred.lower, ymax = epred.upper),
                  fill = "blue", color = "blue", linewidth = 1.05,
                  alpha = 0.25, show.legend = FALSE) +
  geom_lineribbon(aes(y = epred_stan, ymin = epred_stan.lower, ymax = epred_stan.upper),
                  fill = "blue", color = "blue", linewidth = 1.15,
                  alpha = 0.5, show.legend = FALSE) +
  geom_point(data = data %>% 
               mutate(regimen = factor(regimen, levels = regimens)) %>% 
               filter(regimen %in% str_c(c(100, 200, 400, 800), " mg"), 
                      mdv == 0),
             mapping = aes(x = time, y = DV), color = "red", 
             inherit.aes = FALSE) +
  scale_y_continuous(name = latex2exp::TeX("Drug Conc. $(\\mu g/mL)$"),
                     trans = "log10",
                     limits = c(NA, NA)) +
  scale_x_continuous(name = "Time (h)",
                     breaks = seq(0, 216, by = 24),
                     labels = seq(0, 216, by = 24),
                     limits = c(0, 168)) +
  theme_bw() +
  theme(axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "bottom") +
  coord_cartesian(xlim = c(0, 168)) +
  facet_wrap(~ regimen, scales = "free_y")


## Individual estimates (posterior mean)
est_ind <- if(stan_data$want_auc_cmax){
  preds_df %>%
    spread_draws(c(VC, VMAX, KM, KA, 
                   auc_ss, c_max, t_max)[ID]) %>% 
    mean_qi() %>% 
    select(ID, VC, VMAX, KM, KA, 
           auc_ss, c_max, t_max) %>% 
    inner_join(post_preds_summary %>% 
                 filter(time == 168) %>% 
                 select(ID, c_trough = "epred_stan") %>% 
                 distinct(),
               by = "ID")
}else{
  preds_df %>%
    spread_draws(c(VC, VMAX, KM, KA)[ID]) %>% 
    mean_qi() %>% 
    select(ID, VC, VMAX, KM, KA) %>% 
    inner_join(post_preds_summary %>% 
                 filter(time == 168) %>% 
                 select(ID, c_trough = "epred_stan") %>% 
                 distinct(),
               by = "ID")
}

est_ind

