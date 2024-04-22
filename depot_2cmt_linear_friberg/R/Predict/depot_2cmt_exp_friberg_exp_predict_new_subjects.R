rm(list = ls())
cat("\014")

library(trelliscopejs)
library(cmdstanr)
library(tidybayes)
library(posterior)
library(mrgsolve)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

fit <- read_rds("depot_2cmt_linear_friberg/Stan/Fits/depot_2cmt_exp_friberg_exp.rds")

nonmem_data <- read_csv(
  "depot_2cmt_linear_friberg/Data/depot_2cmt_exp_friberg_exp.csv",
  na = ".") %>% 
  rename_all(tolower) %>% 
  rename(ID = "id",
         DV = "dv") %>% 
  mutate(DV = if_else(is.na(DV), 5555555, DV),    # This value can be anything except NA. It'll be indexed away 
         bloq = if_else(is.na(bloq), -999, bloq)) # This value can be anything except NA. It'll be indexed away 

# For this example, let's simulate 10 mg, 20 mg, 40 mg, 80 mg
dosing_data <- mrgsolve::expand.ev(addl = 6, ii = 24, cmt = 1, 
                                   amt = c(10, 20, 40, 80)*1000, 
                                   tinf = 0, evid = 1, mdv = 1) %>%
  as_tibble() %>% 
  mutate(ss = 0) %>% 
  select(ID, time, everything())

t1 <- dosing_data %>% 
  realize_addl() %>% 
  select(time) %>% 
  distinct() %>% 
  deframe()

times_new_pk <- tibble(time = sort(unique(c(t1, 0.25, seq(0, 168, by = 1),
                                            seq(168, 672, by = 24)))))
times_new_pd <- tibble(time = seq(0, 672, by = 12))

new_data_pk <- bind_rows(replicate(max(dosing_data$ID), times_new_pk, 
                                   simplify = FALSE)) %>% 
  mutate(ID = rep(1:max(dosing_data$ID), each = nrow(times_new_pk)),
         amt = 0, 
         evid = 0, 
         rate = 0, 
         addl = 0, 
         ii = 0, 
         cmt = 2, 
         mdv = 1, 
         ss = 0) %>%
  filter(time != 0) %>% 
  select(ID, time, everything())

new_data_pd <- bind_rows(replicate(max(dosing_data$ID), times_new_pd, 
                                   simplify = FALSE)) %>% 
  mutate(ID = rep(1:max(dosing_data$ID), each = nrow(times_new_pd)),
         amt = 0, 
         evid = 0, 
         rate = 0, 
         addl = 0, 
         ii = 0, 
         cmt = 4, 
         mdv = 1, 
         ss = 0) %>%
  select(ID, time, everything())


new_data <- new_data_pk %>% 
  bind_rows(new_data_pd, dosing_data) %>% 
  arrange(ID, time, cmt)


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
                  t_2 = 168)

model <- cmdstan_model(
  "depot_2cmt_linear_friberg/Stan/Predict/depot_2cmt_exp_friberg_exp_predict_new_subjects.stan")

# preds <- model$generate_quantities(fit,
#                                    data = stan_data,
#                                    parallel_chains = 4,
#                                    seed = 1234) 

preds <- model$generate_quantities(fit %>%
                                     as_draws_array() %>%
                                     thin_draws(40),
                                   data = stan_data,
                                   parallel_chains = 4,
                                   seed = 1234)

preds_df <- preds$draws(format = "draws_df")

regimens <- str_c(c(10, 20, 40, 80), " mg")

post_preds_summary <- preds_df %>%
  spread_draws(c(pred, ipred, dv)[i]) %>%
  median_qi(ipred, pred, dv) %>%
  mutate(ID = new_data$ID[i],
         time = new_data$time[i],
         cmt = new_data$cmt[i]) %>%
  select(ID, time, everything(), -i) %>% 
  mutate(regimen = factor(regimens[ID],
                          levels = regimens),
         cmt = case_when(cmt == 2 ~ "PK",
                         cmt == 4 ~ "PD",
                         TRUE ~ NA_character_) %>%
           factor(levels = c("PK", "PD"))) %>% 
  filter(!is.na(cmt))


tmp <- ggplot(post_preds_summary %>% filter(cmt == "PK"), 
              aes(x = time, group = ID)) +
  geom_line(aes(y = ipred), linetype = 1, size = 1.15) +
  geom_line(aes(y = dv), linetype = 2, size = 1.05) +
  ggforce::facet_wrap_paginate(~ ID + cmt, 
                               labeller = label_both,
                               nrow = 1, ncol = length(regimens),
                               page = 1)

p_pk <- p_pd <- vector(mode = "list", length = ggforce::n_pages(tmp))

for(i in 1:ggforce::n_pages(tmp)){
  p_pk[[i]] <- ggplot(post_preds_summary %>% filter(cmt == "PK"), 
                      aes(x = time, group = ID)) +
    geom_ribbon(aes(ymin = dv.lower, ymax = dv.upper),
                fill = "blue", alpha = 0.25, show.legend = FALSE) +
    geom_ribbon(aes(ymin = ipred.lower, ymax = ipred.upper),
                fill = "blue", alpha = 0.5, show.legend = FALSE) +
    geom_line(aes(y = ipred), linetype = 1, size = 1.15) +
    geom_line(aes(y = dv), linetype = 2, size = 1.05) +
    scale_y_continuous(name = latex2exp::TeX("Drug Conc. $(ng/mL)$"),
                       trans = "log10",
                       limits = c(NA, NA)) +
    scale_x_continuous(name = "Time (w)",
                       breaks = seq(0, max(nonmem_data$time), by = 24*7),
                       labels = seq(0, max(nonmem_data$time)/(24*7), by = 1),
                       limits = c(0, max(nonmem_data$time))) +
    # scale_x_continuous(name = "Time (d)") +
    theme_bw() +
    theme(axis.text = element_text(size = 14, face = "bold"),
          axis.title = element_text(size = 18, face = "bold"),
          legend.position = "bottom") +
    ggforce::facet_wrap_paginate(~ regimen,
                                 nrow = 1, ncol = length(regimens),
                                 page = i)
  
}

for(i in 1:ggforce::n_pages(tmp)){
  p_pd[[i]] <- ggplot(post_preds_summary %>% filter(cmt == "PD"), 
                      aes(x = time, group = ID)) +
    geom_ribbon(aes(ymin = dv.lower, ymax = dv.upper),
                fill = "blue", alpha = 0.25, show.legend = FALSE) +
    geom_ribbon(aes(ymin = ipred.lower, ymax = ipred.upper),
                fill = "blue", alpha = 0.5, show.legend = FALSE) +
    geom_line(aes(y = ipred), linetype = 1, size = 1.15) +
    geom_line(aes(y = dv), linetype = 2, size = 1.05) +
    scale_y_continuous(name = latex2exp::TeX("Neutrophils $(\\times 10^9/L)"),
                       trans = "identity") +
    scale_x_continuous(name = "Time (w)",
                       breaks = seq(0, 672, by = 168),
                       labels = seq(0, 672/168, by = 1),
                       limits = c(0, 672)) +
    theme_bw() +
    theme(axis.text = element_text(size = 14, face = "bold"),
          axis.title = element_text(size = 18, face = "bold"),
          legend.position = "bottom") +
    ggforce::facet_wrap_paginate(~ regimen,
                                 nrow = 1, ncol = length(regimens),
                                 page = i)
  
}

for(i in 1:length(p_pk)){
  print(p_pk[[i]] /
          p_pd[[i]])
}

data <- read_csv("depot_2cmt_linear_friberg/Data/depot_2cmt_exp_friberg_exp.csv", 
                 na = ".") %>% 
  rename_all(tolower) %>% 
  rename(ID = "id",
         DV = "dv") %>% 
  group_by(ID) %>% 
  mutate(Dose = max(amt, na.rm = TRUE)/1000,
         mdv = evid) %>% 
  ungroup() %>% 
  mutate(regimen = str_c(Dose, " mg"),
         cmt = case_when(cmt == 2 ~ "PK",
                         cmt == 4 ~ "PD",
                         TRUE ~ NA_character_) %>%
           factor(levels = c("PK", "PD")))


p_pk_check <- post_preds_summary %>% 
  filter(regimen %in% str_c(c(10, 20), " mg"),
         cmt == "PK",
         time <= 168) %>%
  ggplot(aes(x = time, group = ID)) +
  geom_ribbon(aes(ymin = dv.lower, ymax = dv.upper),
              fill = "blue", alpha = 0.25, show.legend = FALSE) +
  geom_ribbon(aes(ymin = ipred.lower, ymax = ipred.upper),
              fill = "blue", alpha = 0.5, show.legend = FALSE) +
  geom_line(aes(y = ipred), linetype = 1, size = 1.15) +
  geom_line(aes(y = dv), linetype = 2, size = 1.05) +
  geom_point(data = data %>% 
               mutate(regimen = factor(regimen, levels = regimens)) %>% 
               filter(regimen %in% str_c(c(10, 20), " mg"), 
                      mdv == 0,
                      cmt == "PK"),
             mapping = aes(x = time, y = DV), color = "red", 
             inherit.aes = FALSE) +
  scale_y_continuous(name = latex2exp::TeX("Drug Conc. $(ng/mL)$"),
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
  facet_wrap(~ regimen + cmt, ncol = 1, scales = "free_y")


p_pd_check <- post_preds_summary %>% 
  filter(regimen %in% str_c(c(10, 20), " mg"),
         cmt == "PD") %>%
  ggplot(aes(x = time, group = ID)) +
  geom_ribbon(aes(ymin = dv.lower, ymax = dv.upper),
              fill = "blue", alpha = 0.25, show.legend = FALSE) +
  geom_ribbon(aes(ymin = ipred.lower, ymax = ipred.upper),
              fill = "blue", alpha = 0.5, show.legend = FALSE) +
  geom_line(aes(y = ipred), linetype = 1, size = 1.15) +
  geom_line(aes(y = dv), linetype = 2, size = 1.05) +
  geom_point(data = data %>% 
               mutate(regimen = factor(regimen, levels = regimens)) %>% 
               filter(regimen %in% str_c(c(10, 20), " mg"), 
                      mdv == 0,
                      cmt == "PD"),
             mapping = aes(x = time, y = DV), color = "red", 
             inherit.aes = FALSE) +
  scale_y_continuous(name = latex2exp::TeX("Neutrophils $(\\times 10^9/L)"),
                     trans = "identity") +
  scale_x_continuous(name = "Time (h)",
                     breaks = seq(0, 672, by = 168),
                     labels = seq(0, 672/168, by = 1),
                     limits = c(0, 672)) +
  theme_bw() +
  theme(axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "bottom") +
  coord_cartesian(xlim = c(0, 672)) +
  facet_wrap(~ regimen + cmt, ncol = 1, scales = "free_y")

ggpubr::ggarrange(p_pk_check, p_pd_check, ncol = 2)

## Summary of distribution of parameters for new subjects 
est_ind <- preds_df %>%
  spread_draws(c(CL, VC, Q, VP, KA, 
                 auc_ss, c_max, t_max, t_half_alpha, t_half_terminal,
                 MTT, CIRC0, GAMMA, ALPHA, r_min)[ID]) %>% 
  median_qi() %>% 
  inner_join(preds_df %>%
               spread_draws(ipred[i]) %>% 
               mutate(ID = new_data$ID[i],
                      time = new_data$time[i],
                      cmt = new_data$cmt[i]) %>% 
               filter(cmt == 2, time == 168) %>% 
               rename(c_trough = ipred) %>% 
               group_by(ID) %>%
               median_qi(c_trough) %>% 
               ungroup() %>% 
               select(ID, c_trough, c_trough.lower = .lower, 
                      c_trough.upper = .upper),
             by = "ID")

est_ind
