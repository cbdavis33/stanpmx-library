rm(list = ls())
cat("\014")

library(trelliscopejs)
library(cmdstanr)
library(mrgsolve)
library(tidybayes)
library(posterior)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

fit <- read_rds(
  "depot_1cmt_linear_covariates/Stan/Fits/depot_1cmt_exp_covariates.rds")

nonmem_data <- read_csv(
  "depot_1cmt_linear_covariates/Data/depot_1cmt_exp_covariates.csv",
  na = ".") %>% 
  rename_all(tolower) %>% 
  rename(ID = "id",
         DV = "dv") %>% 
  mutate(DV = if_else(is.na(DV), 5555555, DV),    # This value can be anything except NA. It'll be indexed away 
         bloq = if_else(is.na(bloq), -999, bloq)) # This value can be anything except NA. It'll be indexed away 

# For this example, let's simulate individuals with various covariates at 400 mg
data_new_covariates <- expand_grid(wt = c(50, 70, 90), 
                                   cmppi = c(0, 1), 
                                   egfr = c(15, 30, 60, 90),
                                   race = 1:4) %>% 
  mutate(ID = 1:n()) %>% 
  relocate(ID, .before = 1)

dosing_data <- mrgsolve::expand.ev(ID = data_new_covariates %>% 
                                     pull(ID),
                                   addl = 6, ii = 24, cmt = 1, 
                                   amt = 400, 
                                   tinf = 0, evid = 1, mdv = 1) %>%
  as_tibble() %>% 
  mutate(ss = 0) %>% 
  select(ID, time, everything())

t1 <- dosing_data %>% 
  realize_addl() %>% 
  select(time) %>% 
  distinct() %>% 
  deframe()

times_new <- tibble(time = sort(unique(c(t1, 0.25, 0.5, seq(0, 168, by = 1)))))

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
  arrange(ID, time) %>% 
  left_join(data_new_covariates, by = "ID")


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

wt <- new_data %>% 
  group_by(ID) %>% 
  distinct(wt) %>% 
  ungroup() %>% 
  pull(wt)

cmppi <- new_data %>% 
  group_by(ID) %>% 
  distinct(cmppi) %>% 
  ungroup() %>% 
  pull(cmppi)

egfr <- new_data %>% 
  group_by(ID) %>% 
  distinct(egfr) %>% 
  ungroup() %>% 
  pull(egfr)

race <- new_data %>% 
  group_by(ID) %>% 
  distinct(race) %>% 
  ungroup() %>% 
  pull(race)

n_races <- nonmem_data %>% 
  distinct(race) %>% 
  drop_na() %>%
  pull(race) %>% 
  length()

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
                  wt = wt,
                  cmppi = cmppi,
                  egfr = egfr,
                  n_races = n_races,
                  race = race,
                  t_1 = 144,
                  t_2 = 168,
                  want_auc_cmax = 1)

model <- cmdstan_model(
  "depot_1cmt_linear_covariates/Stan/Predict/depot_1cmt_exp_predict_new_subjects_covariates.stan")

preds <- model$generate_quantities(fit$draws() %>%
                                     thin_draws(10),
                                   data = stan_data,
                                   parallel_chains = 4,
                                   seed = 1234)

preds_df <- preds$draws(format = "draws_df")

post_preds_summary <- preds_df %>%
  spread_draws(epred_stan[i], epred[i]) %>%
  median_qi(epred_stan, epred) %>%
  mutate(ID = new_data$ID[i],
         time = new_data$time[i],
         wt = factor(new_data$wt[i]),
         cmppi = factor(new_data$cmppi[i]),
         egfr = factor(new_data$egfr[i]),
         race = factor(new_data$race[i])) %>% 
  select(ID, time, everything(), -i) 

(p_wt <- ggplot(post_preds_summary %>% 
                  filter(cmppi == 0, egfr == 90), 
                aes(x = time, group = ID)) +
    geom_lineribbon(aes(y = epred_stan, ymin = epred_stan.lower, 
                        ymax = epred_stan.upper, fill = wt, color = wt),
                    linewidth = 1.15,
                    alpha = 0.25, show.legend = TRUE) +
    scale_y_continuous(name = latex2exp::TeX("Drug Conc. $(\\mu g/mL)$"),
                       trans = "log10",
                       limits = c(NA, NA)) +
    scale_x_continuous(name = "Time (h)",
                       breaks = seq(0, 168, by = 24),
                       labels = seq(0, 168, by = 24),
                       limits = c(0, 168)) +
    theme_bw() +
    theme(axis.text = element_text(size = 14, face = "bold"),
          axis.title = element_text(size = 18, face = "bold"),
          legend.position = "bottom") +
    scale_color_manual(name = "Weight (kg)",
                       values = c("50" = "blue", "70" = "red", 
                                  "90" = "green")) +
    scale_fill_manual(name = "Weight (kg)",
                      values = c("50" = "blue", "70" = "red", 
                                 "90" = "green")) +
    facet_wrap(~race, labeller = label_both))

(p_cmppi <- ggplot(post_preds_summary %>% 
                     filter(wt == 70, egfr == 90), 
                   aes(x = time, group = ID)) +
    geom_lineribbon(aes(y = epred_stan, ymin = epred_stan.lower, 
                        ymax = epred_stan.upper, fill = cmppi, color = cmppi),
                    linewidth = 1.15,
                    alpha = 0.25, show.legend = TRUE) +
    scale_y_continuous(name = latex2exp::TeX("Drug Conc. $(\\mu g/mL)$"),
                       trans = "log10",
                       limits = c(NA, NA)) +
    scale_x_continuous(name = "Time (h)",
                       breaks = seq(0, 168, by = 24),
                       labels = seq(0, 168, by = 24),
                       limits = c(0, 168)) +
    theme_bw() +
    theme(axis.text = element_text(size = 14, face = "bold"),
          axis.title = element_text(size = 18, face = "bold"),
          legend.position = "bottom") +
    scale_color_manual(name = "CMPPI",
                       values = c("0" = "blue", "1" = "red"),
                       labels = c("0" = "No", "1" = "Yes")) +
    scale_fill_manual(name = "CMPPI",
                      values = c("0" = "blue", "1" = "red"),
                      labels = c("0" = "No", "1" = "Yes")) +
    facet_wrap(~race, labeller = label_both))

(p_egfr <- ggplot(post_preds_summary %>% 
                    filter(cmppi == 0, wt == 70), 
                  aes(x = time, group = ID)) +
    geom_lineribbon(aes(y = epred_stan, ymin = epred_stan.lower, 
                        ymax = epred_stan.upper, fill = egfr, color = egfr),
                    linewidth = 1.15,
                    alpha = 0.25, show.legend = TRUE) +
    scale_y_continuous(name = latex2exp::TeX("Drug Conc. $(\\mu g/mL)$"),
                       trans = "log10",
                       limits = c(NA, NA)) +
    scale_x_continuous(name = "Time (h)",
                       breaks = seq(0, 168, by = 24),
                       labels = seq(0, 168, by = 24),
                       limits = c(0, 168)) +
    theme_bw() +
    theme(axis.text = element_text(size = 14, face = "bold"),
          axis.title = element_text(size = 18, face = "bold"),
          legend.position = "bottom") +
    scale_color_manual(name = latex2exp::TeX("$eGFR\\;(mL/min/1.73 m^2)$"),
                       values = c("15" = "blue", "30" = "red",
                                  "60" = "green", "90" = "black")) +
    scale_fill_manual(name = latex2exp::TeX("$eGFR\\;(mL/min/1.73 m^2)$"),
                      values = c("15" = "blue", "30" = "red",
                                 "60" = "green", "90" = "black")) +
    facet_wrap(~race, labeller = label_both))

## Individual estimates (posterior median)
est_ind <- if(stan_data$want_auc_cmax){
  preds_df %>%
    spread_draws(CL[ID], VC[ID], KA[ID], 
                 auc_ss[ID], c_max[ID], t_max[ID], t_half[ID]) %>% 
    median_qi() %>% 
    select(ID, CL, VC, KA, 
           auc_ss, c_max, t_max, t_half) %>% 
    inner_join(post_preds_summary %>% 
                 filter(time == 168) %>% 
                 select(ID, c_trough = "epred_stan") %>% 
                 distinct(),
               by = "ID") %>% 
    left_join(new_data %>% 
                filter(evid == 1) %>%
                distinct(ID, amt, race, wt, cmppi, egfr),
              by = "ID") %>% 
    select(ID, amt, race, wt, cmppi, egfr, everything())
}else{
  preds_df %>%
    spread_draws(CL[ID], VC[ID], KA[ID]) %>% 
    median_qi() %>% 
    select(ID, CL, VC, KA) %>% 
    inner_join(post_preds_summary %>% 
                 filter(time == 168) %>% 
                 select(ID, c_trough = "epred_stan") %>% 
                 distinct(),
               by = "ID") %>% 
    left_join(new_data %>% 
                filter(evid == 1) %>%
                distinct(ID, amt, race, wt, cmppi, egfr),
              by = "ID") %>% 
    select(ID, amt, race, wt, cmppi, egfr, everything())
}

est_ind %>% 
  as_tibble()

