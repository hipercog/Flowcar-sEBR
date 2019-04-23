# Bayesian data analysis of sEBR and flow data

## Packages
library(tidyverse)
library(rstan)
library(brms)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)

## Path (for large data files to be saved)
the_path <- getwd() ## Modify accordingly

## Data
df_r1 <- readr::read_csv("rq1_between.csv")
blk <- readr::read_rds("BLINK_ANALYSIS/blink_data.RData")
fssl <- readr::read_rds("fss_learning.RData")
str(blk)
str(fssl)

## Wrangling
inc.ssn <- c(1, 5, 6, 7, 8)
brate <- cbind(blk$blinkrate.session1, 
               blk$blinkrate.session5, 
               blk$blinkrate.session6, 
               blk$blinkrate.session7, 
               blk$blinkrate.session8)
long.br <- tidyr::gather(as.data.frame(t(brate)), "Parti", "brate")
ssn.lc <- fssl %>% 
  dplyr::filter(Session %in% inc.ssn) %>% 
  dplyr::group_by(Participant, Session) %>% 
  dplyr::summarise(ssn.lc = as.numeric(lm(ln.duration ~ log(Run))$coefficients[2]), 
                   ssn.in = as.numeric(lm(ln.duration ~ log(Run))$coefficients[1]))
df <- fssl %>% 
  dplyr::filter(Session %in% inc.ssn) %>% 
  dplyr::select(-Run, -cumrun, -slope, -intercept, -flow_z, -demand, -skill, -skilldemand) %>% 
  dplyr::group_by(Participant, Session) %>% 
  dplyr::summarise_all(funs(mean), na.rm = TRUE)
df2 <- df %>%
  dplyr::bind_cols(ssn.lc, long.br, dplyr::filter(fssl, Session %in% inc.ssn & Run == 5) %>% 
              dplyr::select(demand, skill, skilldemand)) %>%
  dplyr::select(1:2, brate, duration, ln.duration, 
                learning_curve, distance, ssn.lc, ssn.in, everything(), 
                -Participant1, -Session1, -Parti) %>%
  dplyr::mutate(skidem = skill - demand)
### Normalisation: centering
df2 <- df2 %>% 
  dplyr::group_by(Participant) %>% 
  dplyr::mutate_at(vars(brate), funs(scale), scale = FALSE)
df_final <- df2 %>% 
  dplyr::mutate_at(vars(-Participant, -Session, -brate), funs(scale), scale = FALSE)
df_n <- df_r1 %>%
  dplyr::mutate_at(vars(-Part, -sEBR), funs(scale), scale = TRUE)

# Model fitting (with Hamiltonian Monte Carlo algorithm via Stan)
## General fitting parameters
chains <- 4 # Four markov chains
itr <- 2000 # Number of iterations per chain
wu <- 1000 # The length of warm-up in each chain
ctrl <- list(adapt_delta = .95)
## Defining weakly informative prior on fixed effects (regulating)
all_priors <- brms::prior(student_t(3, 0, 5), class="Intercept") +
  brms::prior(normal(0, 5), class="b")
null_priors <- brms::prior(student_t(3, 0, 5), class="Intercept")

## Research question 1
### Assuming that the response variable of learning curve is normally distributed
fam <- "lognormal"
### Save list
r1_save <- list(
  models = list(),
  parameters = list(
    chains = chains, iterations = itr,
    warmup = wu, control = ctrl,
    family = fam
  ),
  priors = list(all_priors, null_priors)
)
### Setting up a null model on learning curve with intercept-only
bf_r1_null <- brms::bf(sEBR ~ 1)
fit_r1_null <- brms::brm(bf_r1_null, data=df_n,
                         family=fam, prior = null_priors,
                         chains = chains, iter = itr, warmup = wu,
                         control = ctrl)
r1_save$models[[length(r1_save$models)+1]] <- fit_r1_null
### sEBR as an indenpendent variable
bf_r1_sebr <- brms::bf(sEBR ~ LC)
fit_r1_sebr <- brms::brm(bf_r1_sebr, data=df_n,
                         family=fam, prior = all_priors,
                         chains = chains, iter = itr, warmup = wu,
                         control = ctrl)
r1_save$models[[length(r1_save$models)+1]] <- fit_r1_sebr
### sEBR + flow mean as an indenpendent variable
bf_r1_sf <- brms::bf(sEBR ~ LC + FM)
fit_r1_sf <- brms::brm(bf_r1_sf, data=df_n,
                        family=fam, prior = all_priors,
                        chains = chains, iter = itr, warmup = wu,
                        control = ctrl)
r1_save$models[[length(r1_save$models)+1]] <- fit_r1_sf
### sEBR x flow mean as an indenpendent variable
bf_r1_ia <- brms::bf(sEBR ~ LC * FM)
fit_r1_ia <- brms::brm(bf_r1_ia, data=df_n,
                       family=fam, prior = all_priors,
                       chains = chains, iter = itr, warmup = wu,
                       control = ctrl)
r1_save$models[[length(r1_save$models)+1]] <- fit_r1_ia
#### Saving the model file
readr::write_rds(r1_save, file.path(the_path,"R1_models.rds"))

## Research question 2
### Assuming that the response variable of blink rate is log-normally distributed
fam <- "gaussian"
### Save list
r2_save <- list(
  models = list(),
  parameters = list(
    chains = chains, iterations = itr,
    warmup = wu, control = ctrl,
    family = fam
  ),
  priors = list(all_priors, null_priors)
)
### Setting up a null model on brate with random effect of participant
bf_r2_null <- brms::bf(brate ~ (1 | Participant))
fit_r2_null <- brms::brm(bf_r2_null, data=df_final,
                         family=fam, prior = null_priors,
                         chains = chains, iter = itr, warmup = wu,
                         control = ctrl)
r2_save$models[[length(r2_save$models)+1]] <- fit_r2_null
### Flow as an indenpendent variable
bf_r2_flow <- brms::bf(brate ~ flow + (1 | Participant))
fit_r2_flow <- brms::brm(bf_r2_flow, data=df_final,
                      family=fam, prior = all_priors,
                      chains = chains, iter = itr, warmup = wu,
                      control = ctrl)
r2_save$models[[length(r2_save$models)+1]] <- fit_r2_flow
### Absorption as an indenpendent variable
bf_r2_abs <- brms::bf(brate ~ absorption + (1 | Participant))
fit_r2_abs <- brms::brm(bf_r2_abs, data=df_final,
                         family=fam, prior = all_priors,
                         chains = chains, iter = itr, warmup = wu,
                         control = ctrl)
r2_save$models[[length(r2_save$models)+1]] <- fit_r2_abs
### Fluency as an indenpendent variable
bf_r2_flu <- brms::bf(brate ~ fluency + (1 | Participant))
fit_r2_flu <- brms::brm(bf_r2_flu, data=df_final,
                        family=fam, prior = all_priors,
                        chains = chains, iter = itr, warmup = wu,
                        control = ctrl)
r2_save$models[[length(r2_save$models)+1]] <- fit_r2_flu
### Perceived importance 1 as an indenpendent variable
bf_r2_pi1 <- brms::bf(brate ~ pi1 + (1 | Participant))
fit_r2_pi1 <- brms::brm(bf_r2_pi1, data=df_final,
                        family=fam, prior = all_priors,
                        chains = chains, iter = itr, warmup = wu,
                        control = ctrl)
r2_save$models[[length(r2_save$models)+1]] <- fit_r2_pi1
### Perceived importance 2 as an indenpendent variable
bf_r2_pi2 <- brms::bf(brate ~ pi2 + (1 | Participant))
fit_r2_pi2 <- brms::brm(bf_r2_pi2, data=df_final,
                        family=fam, prior = all_priors,
                        chains = chains, iter = itr, warmup = wu,
                        control = ctrl)
r2_save$models[[length(r2_save$models)+1]] <- fit_r2_pi2
### Perceived importance 3 as an indenpendent variable
bf_r2_pi3 <- brms::bf(brate ~ pi3 + (1 | Participant))
fit_r2_pi3 <- brms::brm(bf_r2_pi3, data=df_final,
                        family=fam, prior = all_priors,
                        chains = chains, iter = itr, warmup = wu,
                        control = ctrl)
r2_save$models[[length(r2_save$models)+1]] <- fit_r2_pi3
### Perceived importance total as an indenpendent variable
bf_r2_pit <- brms::bf(brate ~ pi_total + (1 | Participant))
fit_r2_pit <- brms::brm(bf_r2_pit, data=df_final,
                        family=fam, prior = all_priors,
                        chains = chains, iter = itr, warmup = wu,
                        control = ctrl)
r2_save$models[[length(r2_save$models)+1]] <- fit_r2_pit
### Skill as an indenpendent variable
bf_r2_skill <- brms::bf(brate ~ skill + (1 | Participant))
fit_r2_skill <- brms::brm(bf_r2_skill, data=df_final,
                        family=fam, prior = all_priors,
                        chains = chains, iter = itr, warmup = wu,
                        control = ctrl)
r2_save$models[[length(r2_save$models)+1]] <- fit_r2_skill
### Skill-demand as an indenpendent variable
bf_r2_sd <- brms::bf(brate ~ skidem + (1 | Participant))
fit_r2_sd <- brms::brm(bf_r2_sd, data=df_final,
                          family=fam, prior = all_priors,
                          chains = chains, iter = itr, warmup = wu,
                          control = ctrl)
r2_save$models[[length(r2_save$models)+1]] <- fit_r2_sd
### INTERACTION MDOELS ###
### Skill-demand x Perceived importance 1 as indenpendent variables
bf_r2_ia1 <- brms::bf(brate ~ skidem*pi1 + (1 | Participant))
fit_r2_ia1 <- brms::brm(bf_r2_ia1, data=df_final,
                          family=fam, prior = all_priors,
                          chains = chains, iter = itr, warmup = wu,
                          control = ctrl)
r2_save$models[[length(r2_save$models)+1]] <- fit_r2_ia1
### Skill-demand x Perceived importance 1 as indenpendent variables
bf_r2_ia2 <- brms::bf(brate ~ skidem*pi2 + (1 | Participant))
fit_r2_ia2 <- brms::brm(bf_r2_ia2, data=df_final,
                        family=fam, prior = all_priors,
                        chains = chains, iter = itr, warmup = wu,
                        control = ctrl)
r2_save$models[[length(r2_save$models)+1]] <- fit_r2_ia2
#### Saving the model file
readr::write_rds(r2_save, file.path(the_path,"R2_models.rds"))

## Research Question 4
### Assuming that the response variable of duration is normally distributed
fam <- "gaussian"
### Save list
r3_save <- list(
  models = list(),
  parameters = list(
    chains = chains, iterations = itr,
    warmup = wu, control = ctrl,
    family = fam
  ),
  priors = list(all_priors, null_priors)
)
df_mod <- df_final %>%
  dplyr::filter(!is.na(brate) & !is.na(skill) & !is.na(pi1))
### Setting up a null model on duration with random effect of participant
bf_r3_null <- brms::bf(duration ~ (1 | Participant))
fit_r3_null <- brms::brm(bf_r3_null, data=df_mod,
                         family=fam, prior = null_priors,
                         chains = chains, iter = itr, warmup = wu,
                         control = ctrl)
r3_save$models[[length(r3_save$models)+1]] <- fit_r3_null
### Blink rate x Skill x Perceived Importance 1 as independent variables
bf_r3_ia <- brms::bf(duration ~ brate*skill*pi1 + (1 | Participant))
fit_r3_ia <- brms::brm(bf_r3_ia, data=df_mod,
                         family=fam, prior = null_priors,
                         chains = chains, iter = itr, warmup = wu,
                         control = ctrl)
r3_save$models[[length(r3_save$models)+1]] <- fit_r3_ia
#### Saving the model file
readr::write_rds(r3_save, file.path(the_path,"R3_models.rds"))

## Research Question 4
### Assuming that the response variable of distance is normally distributed
fam <- "gaussian"
### Save list
r4_save <- list(
  models = list(),
  parameters = list(
    chains = chains, iterations = itr,
    warmup = wu, control = ctrl,
    family = fam
  ),
  priors = list(all_priors, null_priors)
)
df_mod2 <- df_final %>%
  dplyr::filter(!is.na(brate) & !is.na(pi2))
### Setting up a null model on distance with random effect of participant
bf_r4_null <- brms::bf(distance ~ (1 | Participant))
fit_r4_null <- brms::brm(bf_r4_null, data=df_mod2,
                         family=fam, prior = null_priors,
                         chains = chains, iter = itr, warmup = wu,
                         control = ctrl)
r4_save$models[[length(r4_save$models)+1]] <- fit_r4_null
### Blink rate x Perceived Importance 2 as independent variables
bf_r4_ia <- brms::bf(distance ~ brate*pi2 + (1 | Participant))
fit_r4_ia <- brms::brm(bf_r4_ia, data=df_mod2,
                       family=fam, prior = null_priors,
                       chains = chains, iter = itr, warmup = wu,
                       control = ctrl)
r4_save$models[[length(r4_save$models)+1]] <- fit_r4_ia
#### Saving the model file
readr::write_rds(r4_save, file.path(the_path,"R4_models.rds"))