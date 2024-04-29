rm(list = ls())
cat("\014")

library(trelliscopejs)
library(cmdstanr)
library(tidybayes)
library(posterior)
library(mrgsolve)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

n_depots <- 3
n_depots_with_delay <- n_depots - 1 # "n_depots - 1" means no delay on the 
                                # fastest depot. "n_depots" means a delay on
                                # fastest depot

with_or_without_delay_string <- if_else(n_depots_with_delay == n_depots, 
                                        "with_initial_delay", 
                                        "no_initial_delay")
prior_string <- "prior_ka1_durn"

fit <- read_rds(
  file.path("parallel_first_order_1cmt_linear",
            "zero_into_depot_for_delays",
            "Stan", 
            "Fits",
            str_c("parallel_", n_depots, "_first_order_1cmt_",
                  with_or_without_delay_string, "_", prior_string,
                  "_prop.rds")))

nonmem_data <- read_csv(file.path("parallel_first_order_1cmt_linear",
                                  "zero_into_depot_for_delays",
                                  "Data", 
                                  str_c("parallel_", n_depots, "_first_order_",
                                        with_or_without_delay_string, 
                                        "_prop.csv")),
                        na = ".") %>% 
  rename_all(tolower) %>% 
  rename(ID = "id",
         DV = "dv") %>% 
  mutate(DV = if_else(is.na(DV), 5555555, DV),    # This value can be anything except NA. It'll be indexed away 
         bloq = if_else(is.na(bloq), -999, bloq)) # This value can be anything except NA. It'll be indexed away 

# For this example, let's simulate 100 mg, 200 mg, 400 mg, 800 mg 
dosing_data <- mrgsolve::expand.ev(addl = 0, ii = 0, cmt = 1, 
                                   amt = c(100, 200, 400, 800), 
                                   tinf = 0, evid = 1, mdv = 1) %>%
  as_tibble() %>% 
  bind_rows(mrgsolve::expand.ev(addl = 0, ii = 0, cmt = 2, 
                                amt = c(100, 200, 400, 800), 
                                tinf = 0, evid = 1, mdv = 1) %>% 
              as_tibble(),
            mrgsolve::expand.ev(addl = 0, ii = 0, cmt = 3, 
                                amt = c(100, 200, 400, 800), 
                                tinf = 0, evid = 1, mdv = 1) %>% 
              as_tibble()) %>% 
  mutate(ss = 0) %>% 
  select(ID, time, everything()) %>% 
  arrange(ID, time, cmt)

t1 <- dosing_data %>% 
  realize_addl() %>% 
  select(time) %>% 
  distinct() %>% 
  deframe()

times_new <- tibble(time = sort(unique(c(t1, seq(0, 24, by = 1), 
                                         seq(24, max(nonmem_data$time), 
                                             by = 24)))))

new_data <- bind_rows(replicate(max(dosing_data$ID), times_new, 
                                simplify = FALSE)) %>% 
  mutate(ID = rep(1:max(dosing_data$ID), each = nrow(times_new)),
         amt = 0, 
         evid = 0, 
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
                  ii = new_data$ii,
                  addl = new_data$addl,
                  ss = new_data$ss,
                  subj_start = subj_start,
                  subj_end = subj_end,
                  n_depots = n_depots,
                  n_depots_with_delay = n_depots_with_delay,
                  t_1 = 0,
                  t_2 = 1500)

model <- cmdstan_model(
  file.path("parallel_first_order_1cmt_linear",
            "zero_into_depot_for_delays",
            "Stan", 
            "Predict",
            "parallel_first_order_1cmt_prop_predict_new_subjects.stan"))

# preds <- model$generate_quantities(fit,
#                                    data = stan_data,
#                                    parallel_chains = 4,
#                                    seed = 1234) 

preds <- model$generate_quantities(fit$draws() %>%
                                     thin_draws(100),
                                   data = stan_data,
                                   parallel_chains = 4,
                                   seed = 1234)

preds_df <- preds$draws(format = "draws_df")

regimens <- str_c(c(100, 200, 400, 800), " mg")

post_preds_summary <- preds_df %>%
  spread_draws(ipred[i], pred[i], dv[i]) %>%
  median_qi(ipred, pred, dv) %>%
  mutate(ID = new_data$ID[i],
         time = new_data$time[i]) %>%
  select(ID, time, everything(), -i) %>% 
  mutate(regimen = factor(regimens[ID],
                          levels = regimens))

tmp <- ggplot(post_preds_summary, aes(x = time, group = ID)) +
  geom_line(aes(y = ipred), linetype = 1, size = 1.15) +
  geom_line(aes(y = dv), linetype = 2, size = 1.05) +
  ggforce::facet_wrap_paginate(~ ID, 
                               labeller = label_both,
                               nrow = 1, ncol = 4,
                               page = 1)

for(i in 1:ggforce::n_pages(tmp)){
  print(ggplot(post_preds_summary, aes(x = time, group = ID)) +
          geom_ribbon(aes(ymin = dv.lower, ymax = dv.upper),
                      fill = "blue", alpha = 0.25, show.legend = FALSE) +
          geom_ribbon(aes(ymin = ipred.lower, ymax = ipred.upper),
                      fill = "blue", alpha = 0.5, show.legend = FALSE) +
          geom_line(aes(y = ipred), linetype = 1, size = 1.15) +
          geom_line(aes(y = dv), linetype = 2, size = 1.05) +
          scale_y_continuous(name = latex2exp::TeX("Drug Conc. $(\\mu g/mL)$"),
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
                                       nrow = 1, ncol = 4,
                                       page = i))
  
}

data <- read_csv(file.path("parallel_first_order_1cmt_linear",
                           "zero_into_depot_for_delays",
                           "Data", 
                           str_c("parallel_", n_depots, "_first_order_",
                                 with_or_without_delay_string, 
                                 "_prop.csv")),
                 na = ".") %>% 
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
  geom_ribbon(aes(ymin = dv.lower, ymax = dv.upper),
              fill = "blue", alpha = 0.25, show.legend = FALSE) +
  geom_ribbon(aes(ymin = ipred.lower, ymax = ipred.upper),
              fill = "blue", alpha = 0.5, show.legend = FALSE) +
  geom_line(aes(y = ipred), linetype = 1, size = 1.15) +
  geom_line(aes(y = dv), linetype = 2, size = 1.05) +
  geom_point(data = data %>% 
               mutate(regimen = factor(regimen, levels = regimens)) %>% 
               filter(regimen %in% str_c(c(100, 200, 400, 800), " mg"), 
                      mdv == 0),
             mapping = aes(x = time, y = DV), color = "red", 
             inherit.aes = FALSE) +
  scale_y_continuous(name = latex2exp::TeX("Drug Conc. $(\\mu g/mL)$"),
                     trans = "log10",
                     limits = c(NA, NA)) +
  scale_x_continuous(name = "Time (w)",
                     breaks = seq(0, max(nonmem_data$time), by = 24*7),
                     labels = seq(0, max(nonmem_data$time)/(24*7), by = 1),
                     limits = c(0, max(nonmem_data$time))) +
  theme_bw() +
  theme(axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        legend.position = "bottom") +
  # coord_cartesian(xlim = c(0, 168)) +
  facet_wrap(~ regimen, scales = "free_y")

est_ind <- preds_df %>%
  spread_draws(c(CL, VC, c_max, t_max, t_half)[ID]) %>% 
  median_qi() %>% 
  inner_join(preds_df %>%
               gather_draws(KA[ID, Depot], DUR[ID, Depot]) %>% 
               mutate(.variable = str_c(.variable, "_", Depot)) %>% 
               ungroup() %>% 
               select(-Depot) %>% 
               pivot_wider(names_from = ".variable", values_from = ".value") %>% 
               group_by(ID) %>% 
               mean_qi(),
             by = "ID") %>% 
  select(!contains("."))

est_ind

