rm(list = ls())
cat("\014")

library(trelliscopejs)
library(cmdstanr)
library(tidybayes)
library(posterior)
library(mrgsolve)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")
locf <- FALSE

file_path <- if_else(isTRUE(locf), 
                     "depot_1cmt_linear_covariates_time_varying/Generic/Data/depot_1cmt_ppa_covariates_time_varying_generic_locf.csv",
                     "depot_1cmt_linear_covariates_time_varying/Generic/Data/depot_1cmt_ppa_covariates_time_varying_generic_nocb.csv")


fit <- if(isTRUE(locf)){
  read_rds("depot_1cmt_linear_covariates_time_varying/Generic/Stan/Fits/ppa_locf.rds")
}else{
  read_rds("depot_1cmt_linear_covariates_time_varying/Generic/Stan/Fits/ppa_nocb.rds")
}


nonmem_data <- read_csv(
  file_path,
  na = ".") %>% 
  rename_all(tolower) %>% 
  rename(ID = "id",
         DV = "dv") %>% 
  mutate(DV = if_else(is.na(DV), 5555555, DV),    # This value can be anything except NA. It'll be indexed away 
         bloq = if_else(is.na(bloq), -999, bloq)) # This value can be anything except NA. It'll be indexed away 


# For this example, let's simulate individuals with various covariates at 400 
# mg. These covariates are constant, as if we can't predict how they change over 
# time. You could give time-varying covariates if you wanted
data_new_covariates <- expand_grid(wt = c(50, 70, 90), 
                                   cmppi = c(0, 1), 
                                   egfr = c(15, 30, 60, 90)) %>% 
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

times_new <- tibble(time = sort(unique(c(t1, 0.25, seq(0, 168, by = 1)))))

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
  pull(wt)

cmppi <- new_data %>% 
  pull(cmppi)

egfr <- new_data %>%  
  pull(egfr)

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
                  t_1 = 144,
                  t_2 = 168,
                  want_auc_cmax = 0)

model <- cmdstan_model(
  "depot_1cmt_linear_covariates_time_varying/Generic/Stan/Predict/depot_1cmt_ppa_predict_new_subjects_covariates_time_varying.stan")

preds <- model$generate_quantities(fit$draws() %>%
                                     thin_draws(20),
                                   data = stan_data,
                                   parallel_chains = 4,
                                   seed = 1234)

preds_df <- preds$draws(format = "draws_df")

post_preds_summary <- preds_df %>%
  spread_draws(ipred[i], pred[i], dv[i]) %>%
  median_qi(ipred, pred, dv) %>%
  mutate(ID = new_data$ID[i],
         time = new_data$time[i],
         wt = factor(new_data$wt[i]),
         cmppi = factor(new_data$cmppi[i]),
         egfr = factor(new_data$egfr[i])) %>% 
  select(ID, time, everything(), -i) 

(p_wt <- ggplot(post_preds_summary %>% 
                  filter(cmppi == 0, egfr == 90), 
                aes(x = time, group = ID)) +
    geom_ribbon(aes(ymin = ipred.lower, ymax = ipred.upper, fill = wt),
                alpha = 0.25, show.legend = FALSE) +
    geom_line(aes(y = ipred, color = wt), linetype = 1, linewidth = 1.15) +
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
                                 "90" = "green")))

(p_cmppi <- ggplot(post_preds_summary %>% 
                     filter(wt == 70, egfr == 90), 
                   aes(x = time, group = ID)) +
    geom_ribbon(aes(ymin = ipred.lower, ymax = ipred.upper, fill = cmppi),
                alpha = 0.25, show.legend = FALSE) +
    geom_line(aes(y = ipred, color = cmppi), linetype = 1, linewidth = 1.15) +
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
                      labels = c("0" = "No", "1" = "Yes")))

(p_egfr <- ggplot(post_preds_summary %>% 
                    filter(cmppi == 0, wt == 70), 
                  aes(x = time, group = ID)) +
    geom_ribbon(aes(ymin = ipred.lower, ymax = ipred.upper, fill = egfr),
                alpha = 0.25, show.legend = FALSE) +
    geom_line(aes(y = ipred, color = egfr), linetype = 1, linewidth = 1.15) +
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
    scale_color_manual(name = latex2exp::TeX("$eGFR \\; (mL/min/1.73 m^2)$"),
                       values = c("15" = "blue", "30" = "red",
                                  "60" = "green", "90" = "black")) +
    scale_fill_manual(name = latex2exp::TeX("$eGFR\\;(mL/min/1.73 m^2)$"),
                      values = c("15" = "blue", "30" = "red",
                                 "60" = "green", "90" = "black")))


