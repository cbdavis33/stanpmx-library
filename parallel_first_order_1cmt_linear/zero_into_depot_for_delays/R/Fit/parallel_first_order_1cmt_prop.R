rm(list = ls())
cat("\014")

library(trelliscopejs)
library(cmdstanr)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

n_depots <- 3
n_depots_with_delay <- n_depots - 1 # "n_depots - 1" means no delay on the 
                                    # fastest depot. "n_depots" means a delay on
                                    # fastest depot

with_or_without_delay_string <- if_else(n_depots_with_delay == n_depots, 
                                        "with_initial_delay", 
                                        "no_initial_delay")

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

## Summary of BLOQ values
nonmem_data %>%
  filter(evid == 0) %>%
  group_by(ID) %>%
  summarize(lloq = unique(lloq),
            n_obs = n(),
            n_bloq = sum(bloq)) %>%
  filter(n_bloq > 0)

(p1 <- ggplot(nonmem_data %>%
                group_by(ID) %>%
                mutate(Dose = factor(max(amt, na.rm = TRUE))) %>%
                ungroup() %>%
                filter(mdv == 0)) +
    geom_line(mapping = aes(x = time, y = DV, group = ID, color = Dose)) +
    geom_point(mapping = aes(x = time, y = DV, group = ID, color = Dose)) +
    scale_color_discrete(name = "Dose (mg)") +
    scale_y_continuous(name = latex2exp::TeX("$Drug\\;Conc.\\;(\\mu g/mL)$"),
                       limits = c(NA, NA),
                       trans = "log10") +
    scale_x_continuous(name = "Time (w)",
                       breaks = seq(0, max(nonmem_data$time), by = 24*7),
                       labels = seq(0, max(nonmem_data$time)/(24*7), by = 1),
                       limits = c(0, max(nonmem_data$time))) +
    theme_bw(18) +
    theme(axis.text = element_text(size = 14, face = "bold"),
          axis.title = element_text(size = 18, face = "bold"),
          axis.line = element_line(linewidth = 2),
          legend.position = "bottom"))

# p1 +
#   facet_wrap(~ID, scales = "free_y", labeller = label_both, ncol = 4)
# 
# p1 +
#   facet_trelliscope(~ID, scales = "free_y", ncol = 2, nrow = 2)


n_subjects <- nonmem_data %>%  # number of individuals
  distinct(ID) %>%
  count() %>%
  deframe()

n_total <- nrow(nonmem_data)   # total number of records

i_obs <- nonmem_data %>%
  mutate(row_num = 1:n()) %>%
  filter(evid == 0) %>%
  select(row_num) %>%
  deframe()

n_obs <- length(i_obs)

subj_start <- nonmem_data %>%
  mutate(row_num = 1:n()) %>%
  group_by(ID) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(row_num) %>%
  deframe()

subj_end <- c(subj_start[-1] - 1, n_total)

# length_tvka_prior = 1 means you're setting a single prior for all depots
# length_tvka_prior = n_depots means you're giving each depot its own prior
# length_tvdur_prior = 1 means you're setting a single prior for all depots
# length_tvdur_prior = n_depots_with_delay means you're giving each depot its own prior
# The user then needs to make sure they input the priors with the correct length
#   - ex. if length_tvka_prior = 1, then scale_tvka = 0.1
#         if length_tvdur_prior = 1, then scale_tvdur = 1000
#         if length_tvka_prior = n_depots (= 3), scale_tvka = c(0.001, 0.01, 0.1) (slowest -> fastest)
#         if length_tvdur_prior = n_depots_with_delay (= 3), scale_tvdur = c(2, 400, 1000) (fastest -> slowest)
#         if length_tvdur_prior = n_depots_with_delay (= 2), scale_tvdur = c(400, 1000) (fastest -> slowest)
stan_data <- list(n_subjects = n_subjects,
                  n_total = n_total,
                  n_obs = n_obs,
                  i_obs = i_obs,
                  ID = nonmem_data$ID,
                  amt = nonmem_data$amt,
                  cmt = nonmem_data$cmt,
                  evid = nonmem_data$evid,
                  ii = nonmem_data$ii,
                  addl = nonmem_data$addl,
                  ss = nonmem_data$ss,
                  time = nonmem_data$time,
                  dv = nonmem_data$DV,
                  subj_start = subj_start,
                  subj_end = subj_end,
                  lloq = nonmem_data$lloq,
                  bloq = nonmem_data$bloq,
                  n_depots = n_depots,
                  n_depots_with_delay = n_depots_with_delay,
                  location_tvcl = 1,
                  location_tvvc = 8,
                  length_tvka_prior = 1,  # = 1 or n_depots
                  length_tvdur_prior = n_depots_with_delay, # = 1 or n_depots_with_delay
                  scale_tvka = 0.1, # must have length 'length_tvka_prior'
                  scale_tvdur = c(400, 1000),               # must have length 'length_tvdur_prior'
                  scale_tvcl = 1,
                  scale_tvvc = 1,
                  scale_omega_ka = rep(0.4, times = n_depots),
                  scale_omega_cl = 0.4,
                  scale_omega_vc = 0.4,
                  scale_omega_dur = rep(0.4, times = n_depots_with_delay),
                  alpha_tvfrac = 1,
                  lkj_df_omega = 2,
                  scale_sigma_p = 0.5,
                  prior_only = 0)

rsimplex <- function(n_depots){
  
  x <- 1 + runif(n_depots, -0.5, 0.5)
  
  return(x/sum(x))
  
}

model <- cmdstan_model(
  file.path("parallel_first_order_1cmt_linear",
            "zero_into_depot_for_delays",
            "Stan", 
            "Fit",
            "parallel_first_order_1cmt_prop.stan"),
  cpp_options = list(stan_threads = TRUE))

fit <- model$sample(
  data = stan_data,
  seed = 235813,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = parallel::detectCores()/4,
  iter_warmup = 500,
  iter_sampling = 1000,
  adapt_delta = 0.8,
  refresh = 100,
  max_treedepth = 10,
  init = function() list(TVKA = sort(rlnorm(n_depots, log(0.01), 1)),
                         TVCL = rlnorm(1, log(1), 0.3),
                         TVVC = rlnorm(1, log(8), 0.3),
                         TVDUR_backward = sort(rlnorm(n_depots_with_delay, 
                                                      log(500), 0.7), 
                                               decreasing = FALSE),
                         TVFRAC = rsimplex(n_depots),
                         omega = rlnorm(n_depots + n_depots + 2, 
                                        log(0.3), 0.3),
                         sigma_p = rlnorm(1, log(0.2), 0.3)))

prior_string_ka <- if_else(stan_data$length_tvka_prior == 1, "ka1", "kan")

prior_string_dur <- if_else(stan_data$length_tvdur_prior == 1, "dur1", "durn")

prior_string <- str_c("prior", prior_string_ka, prior_string_dur, sep = "_")

fit$save_object(
  file.path("parallel_first_order_1cmt_linear",
            "zero_into_depot_for_delays",
            "Stan", 
            "Fits",
            str_c("parallel_", n_depots, "_first_order_1cmt_",
                  with_or_without_delay_string, "_", prior_string,
                  "_prop.rds")))

