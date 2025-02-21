rm(list = ls())
cat("\014")

library(gganimate)
library(loo)
library(cmdstanr)
library(tidybayes)
library(bayesplot)
library(patchwork)
library(posterior)
library(tidyverse)

nonmem_data <- read_csv("depot_1cmt_linear/Data/depot_1cmt_prop.csv",
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


## Read in fit
fit <- read_rds("depot_1cmt_linear/Stan/Fits/depot_1cmt_prop.rds")

## Summary of parameter estimates
parameters_to_summarize <- c(str_subset(fit$metadata()$stan_variables, "TV"),
                             str_c("omega_", c("cl", "vc", "ka")),
                             str_subset(fit$metadata()$stan_variables, "cor_"),
                             "sigma_p")

summary <- summarize_draws(fit$draws(parameters_to_summarize), 
                           mean, median, sd, mcse_mean,
                           ~quantile2(.x, probs = c(0.025, 0.975)), rhat,
                           ess_bulk, ess_tail)

summary %>%
  mutate(rse = sd/mean*100) %>% 
  select(variable, mean, sd, q2.5, median, q97.5, rse, rhat, 
         starts_with("ess")) %>% 
  print(n = 50)


## Check the sampler (this is very non-comprehensive)
sampler_diagnostics <- fit$sampler_diagnostics(format = "draws_df")

sampler_diagnostics %>% 
  group_by(.chain) %>% 
  summarize(sum(divergent__)) %>% 
  ungroup()

sampler_diagnostics %>% 
  group_by(.chain) %>% 
  summarize(sum(treedepth__ >= 10)) %>% 
  ungroup()


# Check your largest R-hats. These should be close to 1 (< 1.05 for sure, < 1.02
# ideally)
summary %>%
  arrange(-rhat)

# Density Plots and Traceplots
mcmc_combo(fit$draws(c("TVCL", "TVVC", "TVKA")),
           combo = c("dens_overlay", "trace"))
mcmc_combo(fit$draws(c("omega_cl", "omega_vc", "omega_ka")),
           combo = c("dens_overlay", "trace"))
mcmc_combo(fit$draws(c("sigma_p")),
           combo = c("dens_overlay", "trace"))

mcmc_rank_hist(fit$draws(c("TVCL", "TVVC", "TVKA")))
mcmc_rank_hist(fit$draws(c("omega_cl", "omega_vc", "omega_ka")))
mcmc_rank_hist(fit$draws(c("sigma_p")))

## Check Leave-One-Out Cross-Validation
fit_loo <- fit$loo()
fit_loo
plot(fit_loo, label_points = TRUE)

draws_df <- fit$draws(format = "draws_df")

# Look at predictions
post_preds_summary <- draws_df %>%
  spread_draws(c(pred, epred_stan, epred, ipred, dv_ppc)[i]) %>% 
  mean_qi(pred, epred_stan, epred, ipred, dv_ppc, .width = 0.95) %>%
  mutate(DV = nonmem_data$DV[nonmem_data$evid == 0][i],
         bloq = nonmem_data$bloq[nonmem_data$evid == 0][i],
         lloq = nonmem_data$lloq[nonmem_data$evid == 0][i],
         ID = nonmem_data$ID[nonmem_data$evid == 0][i],
         time = nonmem_data$time[nonmem_data$evid == 0][i]) %>% 
  mutate(DV = if_else(bloq == 1, runif(n(), 0.4*lloq, 0.6*lloq), DV))

post_preds_summary %>% 
  filter(bloq == 0) %>% 
  mutate(in_epred_stan_interval = between(DV, epred_stan.lower, epred_stan.upper),
         in_epred_interval = between(DV, epred.lower, epred.upper),
         in_ipred_interval = between(DV, ipred.lower, ipred.upper),
         in_ppc_interval = between(DV, dv_ppc.lower, dv_ppc.upper)) %>% 
  summarize(across(c(in_epred_stan_interval, in_epred_interval, 
                     in_ipred_interval, in_ppc_interval), mean))

p_dv_vs_pred <- ggplot(post_preds_summary, # %>% filter(bloq == 0), 
                       aes(x = pred, y = DV, color = factor(bloq))) +
  geom_point() +
  scale_color_manual(guide = "none",
                     values = c("0" = "black", "1" = "green")) +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, color = "blue", linewidth = 1.5) +
  geom_smooth(color = "red", se = FALSE, linewidth = 1.5) +
  xlab("Population Predictions (PRED)") +
  ylab("Observed Concentration") +
  theme(axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        axis.line = element_line(linewidth = 2))

p_dv_vs_epred <- ggplot(post_preds_summary, # %>% filter(bloq == 0),  
                        aes(x = epred, y = DV, color = factor(bloq))) +
  geom_point() +
  scale_color_manual(guide = "none",
                     values = c("0" = "black", "1" = "green")) +
  geom_errorbarh(aes(xmin = epred.lower, xmax = epred.upper)) +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, color = "blue", linewidth = 1.5) +
  geom_smooth(color = "red", se = FALSE, linewidth = 1.5) +
  xlab("Population Predictions (EPRED)") +
  ylab("Observed Concentration") +
  theme(axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        axis.line = element_line(linewidth = 2))

p_dv_vs_ipred <- ggplot(post_preds_summary, # %>% filter(bloq == 0),  
                        aes(x = pred, y = DV, color = factor(bloq))) +
  geom_point() +
  scale_color_manual(guide = "none",
                     values = c("0" = "black", "1" = "green")) +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, color = "blue", linewidth = 1.5) +
  geom_smooth(color = "red", se = FALSE, linewidth = 1.5) +
  xlab("Individual Predictions (IPRED)") +
  ylab("Observed Concentration") +
  theme(axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        axis.line = element_line(linewidth = 2))

p_dv_vs_ppc <- ggplot(post_preds_summary, # %>% filter(bloq == 0), 
                      aes(x = dv_ppc, y = DV, color = factor(bloq))) +
  geom_point() +
  scale_color_manual(guide = "none",
                     values = c("0" = "black", "1" = "green")) +
  geom_errorbarh(aes(xmin = dv_ppc.lower, xmax = dv_ppc.upper)) +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, color = "blue", linewidth = 1.5) +
  geom_smooth(color = "red", se = FALSE, linewidth = 1.5) +
  xlab("Individual Predictions (PPC)") +
  ylab("Observed Concentration") +
  theme(axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        axis.line = element_line(linewidth = 2))

p_dv_vs_pred + 
  p_dv_vs_epred +
  p_dv_vs_ipred +
  p_dv_vs_ppc

(p_dv_vs_pred +
    scale_y_log10() +
    scale_x_log10()) +
  (p_dv_vs_epred +
     scale_y_log10() +
     scale_x_log10()) +
  (p_dv_vs_ipred +
     scale_y_log10() +
     scale_x_log10()) +
  (p_dv_vs_ppc +
     scale_y_log10() +
     scale_x_log10())

tmp <- ggplot(post_preds_summary) +
  ggforce::facet_wrap_paginate(~ID, labeller = label_both,
                               nrow = 3, ncol = 4,
                               page = 1, scales = "free_y")

for(i in 1:ggforce::n_pages(tmp)){
  print(ggplot() +
          geom_lineribbon(data = post_preds_summary,
                          mapping = aes(x = time, y = ipred, ymin = ipred.lower,
                                        ymax = ipred.upper, group = ID),
                          fill = "blue", color = "blue", linewidth = 2,
                          alpha = 0.5, show.legend = FALSE) +
          geom_lineribbon(data = post_preds_summary,
                          mapping = aes(x = time, ymin = dv_ppc.lower,
                                        ymax = dv_ppc.upper, group = ID),
                          fill = "blue", linewidth = 2, 
                          alpha = 0.25, show.legend = FALSE) +
          geom_lineribbon(data = post_preds_summary,
                          mapping = aes(x = time, y = epred, ymin = epred.lower,
                                        ymax = epred.upper, group = ID),
                          fill = "red", color = "red", linewidth = 2,
                          alpha = 0.25, show.legend = FALSE) +
          geom_point(data = post_preds_summary, # %>% filter(bloq == 0),
                     mapping = aes(x = time, y = DV, group = ID, color = factor(bloq)),
                     size = 2, show.legend = FALSE) +
          scale_color_manual(values = c("0" = "black", "1" = "green")) +
          scale_y_continuous(name = latex2exp::TeX("$Drug\\;Conc.\\;(\\mu g/mL)$"),
                             limits = c(NA, NA),
                             trans = "identity") +
          scale_x_continuous(name = "Time (h)",
                             breaks = seq(0, max(nonmem_data$time), by = 14),
                             labels = seq(0, max(nonmem_data$time), by = 14),
                             limits = c(0, max(nonmem_data$time))) +
          theme_bw() +
          theme(axis.text = element_text(size = 14, face = "bold"),
                axis.title = element_text(size = 18, face = "bold"),
                legend.position = "bottom") +
          ggforce::facet_wrap_paginate(~ ID, labeller = label_both,
                                       nrow = 3, ncol = 4,
                                       page = i, scales = "free"))
  
}

## Look at residuals and epsilon shrinkage - I think there are two ways to look
## at this, both of which have value:
##   1) by draw 
##   2) using posterior summaries

# Look at it by draw first
residuals <- draws_df %>%
  spread_draws(c(iwres, ipred)[i]) %>% 
  mutate(time = nonmem_data$time[nonmem_data$evid == 0][i],
         ID = nonmem_data$ID[nonmem_data$evid == 0][i],
         bloq = nonmem_data$bloq[nonmem_data$evid == 0][i])

(shrinkage_eps <- residuals %>% 	
    filter(bloq == 0) %>% 
    group_by(.draw) %>% 	
    summarize(sd_iwres = sd(iwres)) %>% 	
    ungroup() %>% 	
    mutate(shrinkage = 1 - sd_iwres) %>% 	
    mean_qi(shrinkage))

residuals %>% 
  filter(bloq == 0) %>%
  sample_draws(4) %>% 
  ggplot(aes(sample = iwres)) + 
  geom_qq() +
  geom_abline() +
  theme_bw() +
  facet_wrap(~.draw, labeller = label_both)

some_residuals_qq <- residuals %>% 
  filter(bloq == 0) %>%
  sample_draws(100) %>% 
  ggplot(aes(sample = iwres)) + 
  geom_qq() +
  geom_abline() +
  theme_bw() +
  transition_manual(.draw)

animate(some_residuals_qq, nframes = 100, width = 384, height = 384, res = 96, 
        dev = "png", type = "cairo")

residuals %>% 
  filter(bloq == 0) %>%
  sample_draws(4) %>%
  ggplot(aes(x = iwres)) + 
  geom_histogram(aes(y = after_stat(density))) +
  geom_density(color = "blue", linewidth = 1.5) +
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1),
                color = "red", size = 1) +
  theme_bw() +
  facet_wrap(~ .draw, labeller = label_both)

residuals %>% 
  filter(bloq == 0) %>%
  sample_draws(9) %>%
  mutate(qn_lower = qnorm(0.025),
         qn_upper = qnorm(0.975)) %>% 
  ggplot(aes(x = time, y = iwres)) + 
  geom_point() +
  geom_hline(yintercept = 0, color = "blue", linewidth = 1.5) +
  geom_hline(aes(yintercept = qn_lower), linetype = "dashed", color = "blue",
             linewidth = 1.25) +
  geom_hline(aes(yintercept = qn_upper), linetype = "dashed", color = "blue",
             linewidth = 1.25) +
  geom_smooth(se = FALSE, color = "red", linewidth = 1.25) +
  theme_bw() +
  xlab("Time (h)") +
  facet_wrap(~ .draw, labeller = label_both)

residuals %>% 
  filter(bloq == 0) %>%
  sample_draws(9) %>%
  mutate(qn_lower = qnorm(0.025),
         qn_upper = qnorm(0.975)) %>% 
  ggplot(aes(x = ipred, y = iwres)) + 
  geom_point() +
  geom_hline(yintercept = 0, color = "blue", linewidth = 1.5) +
  geom_hline(aes(yintercept = qn_lower), linetype = "dashed", color = "blue",
             linewidth = 1.25) +
  geom_hline(aes(yintercept = qn_upper), linetype = "dashed", color = "blue",
             linewidth = 1.25) +
  geom_smooth(se = FALSE, color = "red", linewidth = 1.25) +
  theme_bw() +
  facet_wrap(~ .draw, labeller = label_both)

some_residuals_scatter <- residuals %>% 
  filter(bloq == 0) %>%
  sample_draws(30) %>% 
  mutate(qn_lower = qnorm(0.025),
         qn_upper = qnorm(0.975)) %>%
  ggplot(aes(x = time, y = iwres)) + 
  geom_point() +
  geom_hline(yintercept = 0, color = "blue", linewidth = 1.5) +
  geom_hline(aes(yintercept = qn_lower), linetype = "dashed", color = "blue",
             linewidth = 1.25) +
  geom_hline(aes(yintercept = qn_upper), linetype = "dashed", color = "blue",
             linewidth = 1.25) +
  geom_smooth(se = FALSE, color = "red", linewidth = 1.25) +
  theme_bw() +
  transition_manual(.draw)

animate(some_residuals_scatter, nframes = 30, width = 384, height = 384, 
        res = 96, dev = "png", type = "cairo", fps = 5)

# |iwres| vs. IPRED
residuals %>% 
  filter(bloq == 0) %>%
  sample_draws(9) %>% 
  mutate(abs_iwres = abs(iwres)) %>% 
  ggplot(aes(x = ipred, y = abs_iwres)) +
  geom_point() + 
  geom_smooth(method = "loess", se = FALSE, color = "red", linewidth = 1.25) +
  theme_bw() +
  ylab(("|IWRES|")) +
  facet_wrap(~ .draw, labeller = label_both)

## Some might calculate shrinkage this way based on point estimates, but I don't	
## think this is the right way to go about it, because it uses point estimates	
## rather than the full posterior, but here it is, anyways	
## (1 - sd(posterior means of iwres))	
iwres_mean <- residuals %>%
  filter(bloq == 0) %>%
  group_by(i, time) %>% 
  summarize(iwres_mean = mean(iwres)) %>% 
  ungroup()

iwres_mean %>% 
  summarize(shrinkage = 1 - sd(iwres_mean))

iwres_mean %>% 
  ggplot(aes(x = time, y = iwres_mean)) +
  theme_bw() +
  geom_point() +
  geom_hline(yintercept = 0, color = "blue", linewidth = 1.5) +
  scale_y_continuous(name = "iwres") +
  scale_x_continuous(name = "Time (h)",
                     limits = c(NA, NA))

iwres_mean %>% 
  ggplot(aes(sample = iwres_mean)) + 
  geom_qq() +
  geom_abline() +
  theme_bw()

# # This one is hard to see when each individual has very similar timepoints
# residuals %>%
#   filter(bloq == 0) %>%
#   ggplot(aes(x = time, y = iwres, group = i)) +
#   stat_pointinterval(.width = c(0.95), point_size = 3) +
#   theme_bw() +
#   geom_hline(yintercept = 0, color = "blue", linewidth = 1.5) +
#   scale_x_continuous(name = "Time (h)",
#                      limits = c(NA, NA))

## Look at eta shrinkage 
# We can look at individual posterior densities on top of the density of the 
# population parameter
blah <- draws_df %>%
  gather_draws(CL[ID], VC[ID], KA[ID], TVCL, TVVC, TVKA) %>%
  ungroup() %>%
  mutate(across(c(ID, .variable), as.factor))

ggplot() +
  geom_density(data = blah %>% 
                 filter(is.na(ID)) %>% 
                 mutate(.variable = blah %>% 
                          filter(is.na(ID)) %>% 
                          pull(.variable) %>% 
                          str_remove("TV"),
                        ID = "Population"),
               mapping = aes(x = .value, group = .variable,
                             fill = .variable), color = "black",
               position = "identity") +
  stat_density(data = blah %>% 
                 filter(!is.na(ID)), 
               mapping = aes(x = .value, group = ID, color = ID),
               geom = "line", position = "identity") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~.variable, ncol = 1, scales = "free") +
  theme(legend.position = "bottom") +
  scale_fill_manual(name = "Population Parameter",
                    values = c("red", "blue", "green")) +
  guides(color = "none") +
  ggtitle("Individual Parameter Posterior Densities") 


## Shrinkage
# A function to visualize the shrinkage
plot_shrinkage <- function(.variable, std_dev, ind_params, ...){
  
  dots <- list(...)
  
  x_label <- glue::glue("$\\eta_{[.variable]}$", .open = "[", .close = "]")
  
  p_1 <- ggplot() +
    geom_histogram(aes(x = ind_params$ind_params,
                       y = after_stat(density))) +
    stat_function(fun = dnorm, 
                  args = list(mean = 0, sd = std_dev),
                  color = "red", linewidth = 1) +
    scale_x_continuous(name = latex2exp::TeX(x_label),
                       limits = c(qnorm(c(0.001, 0.995), 
                                        0, std_dev))) +
    theme_bw()
  
  
  if(is.numeric(dots$.draw)){
    
    p_1 <- p_1 +
      labs(subtitle = str_c("draw = ", dots$.draw)) +
      theme(plot.subtitle = element_text(hjust = 0.5))
    
  }
  
  return(p_1)
  
}

# This function stacks multiple shrinkage plots (e.g. by draw)
plot_shrinkage_multiple <- function(draw){
  
  data_shrinkage_by_draw %>% 
    filter(.draw == draw) %>% 
    pmap(plot_shrinkage) %>% 
    wrap_plots(nrow = 1)
  
}

# First look at shrinkage by draw
(shrinkage <- draws_df %>%
    gather_draws(`eta_.*`[ID], regex = TRUE) %>% 
    group_by(.draw, .variable) %>% 
    summarize(sd_eta = sd(.value)) %>% 
    ungroup() %>% 
    mutate(.variable = str_remove(.variable, "eta_")) %>% 
    left_join(draws_df %>%
                gather_draws(`omega_.*`, regex = TRUE) %>% 
                ungroup() %>% 
                arrange(.draw) %>% 
                rename(omega = ".value") %>% 
                mutate(.variable = str_remove(.variable, "omega_")),
              by = c(".draw", ".variable")) %>% 
    mutate(shrinkage = 1 - sd_eta/omega) %>% 
    select(.draw, .variable, shrinkage) %>% 
    group_by(.variable) %>% 
    mean_qi())

draws_for_shrinkage <- draws_df %>%
  sample_draws(4)

data_shrinkage_by_draw <- draws_for_shrinkage %>%
  gather_draws(`omega_.*`, regex = TRUE) %>%
  mutate(.variable = str_remove(.variable, "omega_") %>%
           toupper()) %>%
  rename(std_dev = .value) %>%
  ungroup() %>% 
  inner_join(draws_for_shrinkage %>%
               gather_draws(`eta_.*`[ID], regex = TRUE) %>%
               ungroup()  %>% 
               mutate(.variable = str_remove(.variable, "eta_") %>%
                        toupper()) %>%
               rename(ind_params = `.value`) %>%
               select(-ID) %>%
               group_by(.draw) %>%
               nest(ind_params = ind_params) %>%
               ungroup(),
             by = c(".draw", ".chain", ".iteration", ".variable")) %>%
  arrange(.draw, .variable) %>%
  select(-.chain, -.iteration)

map(data_shrinkage_by_draw %>% 
      distinct(.draw) %>% 
      pull(.draw), 
    .f = plot_shrinkage_multiple) %>% 
  wrap_plots(ncol = 1)

## Some might calculate shrinkage this way based on point estimates, but I don't
## think this is the right way to go about it, because it uses point estimates
## rather than the full posterior, but here it is, anyways
## (1 - sd(individual means)/(posterior mean of omega))
draws_df %>%
  gather_draws(`eta_.*`[ID], regex = TRUE) %>% 
  summarize(estimate = mean(.value)) %>% 
  ungroup() %>% 
  group_by(.variable) %>% 
  summarize(std_dev = sd(estimate)) %>% 
  ungroup() %>% 
  mutate(.variable = str_remove(.variable, "eta_")) %>% 
  left_join(draws_df %>%
              gather_draws(`omega_.*`, regex = TRUE) %>%
              summarize(estimate = mean(.value)) %>% 
              ungroup() %>% 
              mutate(.variable = str_remove(.variable, "omega_")),
            by = ".variable") %>% 
  mutate(shrinkage = 1 - std_dev/estimate)

data_shrinkage_with_point_estimates <- draws_df %>% 
  gather_draws(`omega_.*`, regex = TRUE) %>%
  summarize(std_dev = mean(.value)) %>% 
  mutate(.variable = str_remove(.variable, "omega_") %>% 
           toupper()) %>% 
  inner_join(draws_df %>%
               gather_draws(`eta_.*`[ID], regex = TRUE) %>% 
               summarize(ind_params = mean(.value)) %>% 
               ungroup() %>% 
               mutate(.variable = str_remove(.variable, "eta_") %>%
                        toupper()) %>%
               select(-ID) %>% 
               nest(ind_params = ind_params), 
             by = ".variable")

pmap(data_shrinkage_with_point_estimates, 
     plot_shrinkage) %>% 
  wrap_plots(nrow = 1)


# Individual point estimates (posterior mean)
est_ind <- draws_df %>%
  spread_draws(c(CL, VC, KA, eta_cl, eta_vc, eta_ka)[ID]) %>%
  mean_qi() %>% 
  select(ID, CL, VC, KA,
         eta_cl, eta_vc, eta_ka) %>% 
  mutate(ID = factor(ID))

## Standardized Random Effects (posterior mean)
eta_std <- est_ind %>% 
  select(ID, starts_with("eta")) %>% 
  pivot_longer(cols = starts_with("eta"), 
               names_to = "parameter", 
               values_to = "eta",
               names_prefix = "eta_") %>% 
  group_by(parameter) %>% 
  mutate(eta_std = (eta - mean(eta))/sd(eta)) %>% 
  ungroup() %>% 
  mutate(parameter = toupper(parameter))

(est_ind %>% 
    ggplot(aes(x = eta_cl, y = eta_vc)) + 
    geom_point() +
    theme_bw() + 
    geom_smooth(method = "loess", span = 0.9) +
    ggtitle("Point Estimates")) +
  (draws_df %>%
     sample_draws(10) %>%
     spread_draws(`eta_.*`[ID], regex = TRUE) %>% 
     ungroup()  %>% 
     mutate(ID = factor(ID)) %>% 
     ggplot(aes(x = eta_cl, y = eta_vc)) + 
     geom_point() +
     theme_bw() + 
     geom_smooth(method = "loess", span = 0.9) +
     ggtitle("Samples"))


