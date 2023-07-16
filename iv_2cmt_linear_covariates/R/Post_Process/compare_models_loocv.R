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

fit_no_covariates <- read_rds("iv_2cmt_linear_covariates/Stan/Fits/iv_2cmt_prop_no_covariates.rds")
fit_fixed_allometric <- read_rds("iv_2cmt_linear_covariates/Stan/Fits/iv_2cmt_prop_covariates_fixed_allometric.rds")
fit_covariates <- read_rds("iv_2cmt_linear_covariates/Stan/Fits/iv_2cmt_prop_covariates.rds")


fit_no_covariates_loo <- fit_no_covariates$loo()
fit_fixed_allometric_loo <- fit_wt_fixed_allometric$loo()
fit_covariates_loo <- fit_covariates$loo()

loo_compare(fit_no_covariates_loo, fit_fixed_allometric_loo,
            fit_covariates_loo)
