# Bayesian model evaluation of sEBR and flow data

## Packages
library(tidyverse)
library(rstan)
library(brms)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)

## Path to load and save
the_path <- getwd() ## Modify accordingly

## Load models
r1_models <- readr::read_rds(file.path(the_path,"R1_models.rds"))
r2_models <- readr::read_rds(file.path(the_path,"R2_models.rds"))
r3_models <- readr::read_rds(file.path(the_path,"R3_models.rds"))
r4_models <- readr::read_rds(file.path(the_path,"R4_models.rds"))

## Leave-One-Out Cross-Validation
r1_loo <- brms::loo(r1_models$models[[1]], r1_models$models[[2]],
                    r1_models$models[[3]], r1_models$models[[4]],
                    reloo=TRUE)
r2_loo <- brms::loo(r2_models$models[[1]], r2_models$models[[2]],
                    r2_models$models[[3]], r2_models$models[[4]],
                    r2_models$models[[5]], r2_models$models[[6]],
                    r2_models$models[[7]], r2_models$models[[8]],
                    r2_models$models[[9]], r2_models$models[[10]],
                    r2_models$models[[11]], r2_models$models[[12]],
                    reloo=TRUE)
r3_loo <- brms::loo(r3_models$models[[1]], r3_models$models[[2]],
                    reloo=TRUE)
r4_loo <- brms::loo(r4_models$models[[1]], r4_models$models[[2]],
                    reloo=TRUE)

## Save LOO comparisons
readr::write_rds(r1_loo, file.path(the_path,"R1_loo.rds"))
readr::write_rds(r2_loo, file.path(the_path,"R2_loo.rds"))
readr::write_rds(r3_loo, file.path(the_path,"R3_loo.rds"))
readr::write_rds(r4_loo, file.path(the_path,"R4_loo.rds"))
