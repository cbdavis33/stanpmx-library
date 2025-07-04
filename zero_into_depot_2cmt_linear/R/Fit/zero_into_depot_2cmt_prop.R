rm(list = ls())
cat("\014")

library(cmdstanr)
library(tidyverse)

set_cmdstan_path("~/Torsten/cmdstan")

nonmem_data <- read_csv(
  "zero_into_depot_2cmt_linear/Data/zero_into_depot_2cmt_prop.csv",
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

# (p1 <- ggplot(nonmem_data %>%
#                 group_by(ID) %>%
#                 mutate(Dose = factor(max(amt, na.rm = TRUE))) %>%
#                 ungroup() %>%
#                 filter(mdv == 0)) +
#     geom_line(mapping = aes(x = time, y = DV, group = ID, color = Dose)) +
#     geom_point(mapping = aes(x = time, y = DV, group = ID, color = Dose)) +
#     scale_color_discrete(name = "Dose (mg)") +
#     scale_y_continuous(name = latex2exp::TeX("$Drug\\;Conc.\\;(\\mu g/mL)$"),
#                        limits = c(NA, NA),
#                        trans = "identity") +
#     scale_x_continuous(name = "Time (d)",
#                        breaks = seq(0, 216, by = 24),
#                        labels = seq(0, 216/24, by = 24/24),
#                        limits = c(0, NA)) +
#     theme_bw(18) +
#     theme(axis.text = element_text(size = 14, face = "bold"),
#           axis.title = element_text(size = 18, face = "bold"),
#           axis.line = element_line(linewidth = 2),
#           legend.position = "bottom"))
# 
# p1 +
#   facet_wrap(~ID, scales = "free_y", labeller = label_both, ncol = 4)


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
                  location_tvcl = 1,
                  location_tvvc = 18,
                  location_tvq = 3,
                  location_tvvp = 45,
                  location_tvka = 0.8,
                  location_tvdur = 1,
                  scale_tvcl = 1,
                  scale_tvvc = 1,
                  scale_tvq = 1,
                  scale_tvvp = 1,
                  scale_tvka = 1,
                  scale_tvdur = 1,
                  scale_omega_cl = 0.4,
                  scale_omega_vc = 0.4,
                  scale_omega_q = 0.4,
                  scale_omega_vp = 0.4,
                  scale_omega_ka = 0.4,
                  scale_omega_dur = 0.4,
                  lkj_df_omega = 2,
                  scale_sigma_p = 0.5,
                  prior_only = 0,
                  no_gq_predictions = 0) 

model <- cmdstan_model(
  "zero_into_depot_2cmt_linear/Stan/Fit/zero_into_depot_2cmt_prop.stan",
  cpp_options = list(stan_threads = TRUE))

fit <- model$sample(
  data = stan_data,
  seed = 112358,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = parallel::detectCores()/4,
  iter_warmup = 500,
  iter_sampling = 1000,
  adapt_delta = 0.8,
  refresh = 500,
  max_treedepth = 10,
  output_dir = "zero_into_depot_2cmt_linear/Stan/Fits/Output",
  output_basename = "prop",
  init = function() 
    with(stan_data,
         list(TVCL = rlnorm(1, log(location_tvcl), scale_tvcl/10),
              TVVC = rlnorm(1, log(location_tvvc), scale_tvvc/10),
              TVQ = rlnorm(1, log(location_tvq), scale_tvq/10),
              TVVP = rlnorm(1, log(location_tvvp), scale_tvvp/10),
              TVKA = rlnorm(1, log(location_tvka), scale_tvka/10),
              TVDUR = rlnorm(1, log(location_tvdur), scale_tvdur),
              omega = abs(rnorm(6, 0, c(scale_omega_cl,
                                        scale_omega_vc,
                                        scale_omega_q,
                                        scale_omega_vp,
                                        scale_omega_ka,
                                        scale_omega_dur))),
              sigma_p = abs(rnorm(1, 0, scale_sigma_p)))))

fit$save_object(
  "zero_into_depot_2cmt_linear/Stan/Fits/zero_into_depot_2cmt_prop.rds")

fit$save_data_file(dir = "zero_into_depot_2cmt_linear/Stan/Fits/Stan_Data",
                   basename = "prop", timestamp = FALSE, random = FALSE)
