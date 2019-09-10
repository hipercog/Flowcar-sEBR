# Flowcar analysis workflow
## Libraries
library(tidyverse)
library(modelr)
library(tidybayes)
library(ggridges)
library(rstan)
library(brms)
library(gridExtra)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()-1)
## Sourcing
source("wf_utils.R")

## Data
df_r1 <- readr::read_csv("rq1_between.csv")
blk <- readr::read_rds("BLINK_ANALYSIS/blink_data.RData")
fssl <- readr::read_rds("fss_learning.RData")

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


############# Workflow #############

# Assign dataset to be investigated
the_data <- df_n
# Name of the response variable
resp_var <- "sEBR"
# Assign the model for investigation
model_formula <- brms::bf(sEBR ~ LC * FM + (1 | Part))
# Assign response variable distribution family
fam <- "lognormal"

# Steo 1: Prior predictive checks
## Check the required parameters
brms::get_prior(model_formula, the_data, family = fam)
## Set priors
prior_info <- list(
  intercept = "t(3, log(10), 1), expectation and variance retrieved from van Slooten 2019",
  b = "Fixed effects (FM, LC, LC:FM) diffuse and weakly informative, N(0, 0.5). Limited variance due to lognormal response family.",
  sd = "Standard deviation of intercept-varying random effects of participants also diffuse and weakly informative but with fat tails, t(3, 0, 0.5)",
  sigma = "Residuals also diffuse and weakly informative but with greater variance and fat tails, t(3, 0, 0.75)"
)
priors <- prior(student_t(3, 2.3, 1), class="Intercept") +
  prior(normal(0, .5), class="b") +
  prior(student_t(3, 0, .5), class="sd") +
  prior(student_t(3, 0, .75), class="sigma")
## Sampling from prior
### Sampling parameters
chains <- 4
iter <- 3000
warmup <- 1500
control <- list(adapt_delta = .99)
config <- list(chains=chains, iter=iter, warmup=warmup, control=control)
### Sampling
prior_fit <- brms::brm(model_formula,
                       family = fam,
                       data = the_data,
                       prior = priors,
                       chains = chains,
                       iter = iter,
                       warmup = warmup,
                       control = control,
                       sample_prior = "only")
### Extracting samples
true_params <- brms::posterior_samples(prior_fit)
samps <- fitted(prior_fit, summary = FALSE)
samps_lin <- fitted(prior_fit, summary = FALSE, scale = "linear")
prior_df <- tibble(Response = as.vector(samps),
                   Linear = as.vector(samps_lin)) %>%
  tidyr::gather(variable, value) %>%
  dplyr::mutate(
    var_label = ifelse(variable=="Response",
                      "Response variable (original scale)",
                      "Linear predictor (link-function scale)")
  )

# Step 2: Computational faithfulness
expdesign <- the_data
## Estimating the dataset for N-simulated data points
the_max <- dim(samps)[1]
### Number of simulations
nsims <- 10
### Picks randomly n-response samples from simulated data
if (nsims<the_max) picks <- sample(1:the_max, nsims, replace = F) else picks <- 1:the_max
expdesign[[resp_var]] <- samps[picks[1],]
brm1 <- brms::brm(model_formula,
                       family = fam,
                       data = expdesign,
                       prior = priors,
                       chains = chains,
                       iter = iter,
                       warmup = warmup,
                       control = control)
mods_gw <- list(brm1)
for (i in 2:length(picks)) {
  expdesign[[resp_var]] <- samps[picks[i],]
  mods_gw[[i]] <- update(brm1, newdata=expdesign, recompile=FALSE)
}
### Extracting diagnostics
logPost <- div_trans <- rhat <- neffRatio <- NA
for (i in 1:length(picks)) {
  logPost[i] <- brms::log_posterior(mods_gw[[i]])
  np <- brms::nuts_params(mods_gw[[i]])
  div_trans[i] <- sum(subset(np, Parameter == "divergent__")$Value) 
  rhat[i] <- brms::rhat(mods_gw[[i]])
  neffRatio[i] <- brms::neff_ratio(mods_gw[[i]])
}
diagnostic_summary <- list(
  log_post = logPost,
  div_trans = div_trans,
  rhat = rhat,
  neff_ratios = neffRatio
)
## Simulation-Based Calibration (SBC)
### Choose thinning step (to remove autocorrelation from HMC samples)
thin_step <- 8
### Name of the parameters
params <- c("b_LC","b_FM","b_LC:FM")
### Generating rankings
sbc_ranks <- sbc_rank_n(mods_gw, params, true_params, thin_step, picks)

# Step 3: Model sensitivity
all_sensitivity <- all_model_sensitivity(mods_gw, params, true_params, picks, priors)

# Step 4: Posterior predictive checks
real_model <- brms::brm(model_formula,
                  family = fam,
                  data = the_data,
                  prior = priors,
                  chains = chains,
                  iter = iter,
                  warmup = warmup,
                  control = control)
post_samples <- brms::posterior_samples(real_model)
ps_tidy <- post_samples %>%
  tidyr::gather(variable, value)

# Saving the workflow
## Compiling
workflow_output <- list(
  model = list(
    formula = model_formula,
    data = the_data,
    family = fam
  ),
  prior_pc_s1 = list(
    info = prior_info,
    priors = priors,
    config = config,
    fit = prior_fit,
    samples = prior_df
  ),
  comp_faithf_s2 = list(
    nsims = nsims,
    diagnostics = diagnostic_summary,
    sbc = list(
      parameters = params,
      thinning = thin_step,
      sbc_ranks = sbc_ranks
    )
  ),
  model_sens_s3 = all_sensitivity,
  post_pc_s4 = list(
    fit = real_model,
    samples = ps_tidy
  )
)
## Location and ID
output_folder <- "/Users/ville/Documents/HY"
wf_id <- paste0(sample(c(letters, LETTERS, 0:9),4), collapse = "")
file_id <- stringr::str_c("workflow_",wf_id,".rds")
## Saving the output
readr::write_rds(workflow_output, file.path(output_folder,file_id))
