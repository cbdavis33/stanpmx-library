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

fit_no_covariates <- read_rds("depot_1cmt_linear_covariates/Stan/Fits/depot_1cmt_prop_no_covariates.rds")
fit_wt_fixed_allometric <- read_rds("depot_1cmt_linear_covariates/Stan/Fits/depot_1cmt_prop_wt_fixed_allometric.rds")
fit_wt <- read_rds("depot_1cmt_linear_covariates/Stan/Fits/depot_1cmt_prop_wt.rds")
fit_covariates <- read_rds("depot_1cmt_linear_covariates/Stan/Fits/depot_1cmt_prop_covariates.rds")


fit_no_covariates_loo <- fit_no_covariates$loo()
fit_wt_fixed_allometric_loo <- fit_wt_fixed_allometric$loo()
fit_wt_loo <- fit_wt$loo()
fit_covariates_loo <- fit_covariates$loo()

loo_compare(fit_no_covariates_loo, fit_wt_fixed_allometric_loo,
            fit_wt_loo, fit_covariates_loo)
