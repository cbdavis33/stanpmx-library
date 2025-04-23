# LOO-CV
# LOGO-CV

# 1) Fit PPA data with _prop and _exp and do model comparison
#     a) loo-cv
#     b) logo-cv

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

# Leave-one-out cross-validation
filename_base <- "depot_1cmt_linear/Stan/Fits"

fit_ppa <- read_rds(file.path(filename_base, "depot_1cmt_ppa.rds"))
fit_prop <- read_rds(file.path(filename_base, "ppa_fit_with_prop.rds"))
fit_exp <- read_rds(file.path(filename_base, "ppa_fit_with_exp.rds"))

fit_ppa_loo <- fit_ppa$loo(cores = min(4, parallel::detectCores()/2))
fit_prop_loo <- fit_prop$loo(cores = min(4, parallel::detectCores()/2))
fit_exp_loo <- fit_exp$loo(cores = min(4, parallel::detectCores()/2))

loo_compare(fit_ppa_loo, fit_prop_loo, fit_exp_loo)
