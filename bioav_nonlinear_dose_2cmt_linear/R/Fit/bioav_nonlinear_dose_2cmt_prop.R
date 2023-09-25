rm(list = ls())
cat("\014")

library(trelliscopejs)
library(cmdstanr)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

nonmem_data <- read_csv(
  "bioav_nonlinear_dose_2cmt_linear/Data/bioav_nonlinear_dose_2cmt_prop.csv",
  na = ".") %>% 
  rename_all(tolower) %>% 
  rename(ID = "id",
         DV = "dv") %>% 
  mutate(dose = factor(dose),
         DV = if_else(is.na(DV), 5555555, DV),    # This value can be anything except NA. It'll be indexed away 
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
                filter(mdv == 0)) +
    geom_line(mapping = aes(x = time, y = DV, group = ID, color = dose)) +
    geom_point(mapping = aes(x = time, y = DV, group = ID, color = dose)) +
    scale_color_discrete(name = "Dose") +
    scale_y_continuous(name = latex2exp::TeX("$Drug\\;Conc.\\;(\\mu g/mL)$"),
                       limits = c(NA, NA),
                       trans = "log10") +
    scale_x_continuous(name = "Time (d)",
                       breaks = seq(0, 168, by = 24),
                       labels = seq(0, 168/24, by = 24/24),
                       limits = c(0, NA)) +
    theme_bw(18) +
    theme(axis.text = element_text(size = 14, face = "bold"),
          axis.title = element_text(size = 18, face = "bold"),
          axis.line = element_line(linewidth = 2),
          legend.position = "bottom"))

p1 +
  facet_wrap(~dose, scales = "free_y", labeller = label_both, ncol = 3)

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

dose <- nonmem_data %>% 
  group_by(ID) %>%
  distinct(dose) %>% 
  ungroup() %>% 
  pull(dose)

stan_data <- list(n_subjects = n_subjects,
                  n_total = n_total,
                  n_obs = n_obs,
                  i_obs = i_obs,
                  ID = nonmem_data$ID,
                  amt = nonmem_data$amt,
                  cmt = nonmem_data$cmt,
                  evid = nonmem_data$evid,
                  rate = nonmem_data$rate,
                  ii = nonmem_data$ii,
                  addl = nonmem_data$addl,
                  ss = nonmem_data$ss,
                  time = nonmem_data$time,
                  dv = nonmem_data$DV,
                  subj_start = subj_start,
                  subj_end = subj_end,
                  lloq = nonmem_data$lloq,
                  bloq = nonmem_data$bloq,
                  location_tvcl = 0.75,
                  location_tvvc = 18,
                  location_tvq = 3,
                  location_tvvp = 35,
                  location_tvka = 0.8,
                  location_fmax = 0.5,
                  location_f50 = 400,
                  location_hill = 1,
                  scale_tvcl = 1,
                  scale_tvvc = 1,
                  scale_tvq = 1,
                  scale_tvvp = 1,
                  scale_tvka = 1,
                  scale_fmax = 3.5,
                  scale_f50 = 0.6,
                  scale_hill = 0.75,
                  scale_omega_cl = 0.4,
                  scale_omega_vc = 0.4,
                  scale_omega_q = 0.4,
                  scale_omega_vp = 0.4,
                  scale_omega_ka = 0.4,
                  reference_dose = min(as.numeric(as.character(dose))),
                  dose = as.numeric(as.character(dose)),
                  lkj_df_omega = 2,
                  scale_sigma_p = 0.5,
                  prior_only = 0)

model <- cmdstan_model(
  "bioav_nonlinear_dose_2cmt_linear/Stan/Fit/bioav_nonlinear_dose_2cmt_prop.stan",
  cpp_options = list(stan_threads = TRUE))

fit <- 
  model$sample(data = stan_data,
               seed = 11235,
               chains = 4,
               parallel_chains = 4,
               threads_per_chain = parallel::detectCores()/4,
               iter_warmup = 500,
               iter_sampling = 1000,
               adapt_delta = 0.8,
               refresh = 500,
               max_treedepth = 10,
               init = function() list(TVCL = rlnorm(1, log(1), 0.3),
                                      TVVC = rlnorm(1, log(18), 0.3),
                                      TVQ = rlnorm(1, log(4), 0.3),
                                      TVVP = rlnorm(1, log(35), 0.3),
                                      TVKA = rlnorm(1, log(1), 0.3),
                                      FMAX = rbeta(1, 3, 3),
                                      F50 = rlnorm(1, log(500), 0.3),
                                      HILL = runif(1, 0.5, 4),
                                      omega = rlnorm(5, log(0.3), 0.3),
                                      sigma_p = rlnorm(1, log(0.4), 0.3)))

fit$save_object(
  "bioav_nonlinear_dose_2cmt_linear/Stan/Fits/bioav_nonlinear_dose_2cmt_prop.rds")

