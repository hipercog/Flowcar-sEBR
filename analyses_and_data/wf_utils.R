# Utilities for modelling workflow
library(rethinking)
library(LaplacesDemon)


SimFromPrior <- function(priors,class="b",coef="") {
  priors_par <- priors$prior[priors$class==class & priors$coef==coef]
  priors_par <- strsplit(priors_par,"(",fixed=TRUE)[[1]]
  par <- strsplit(priors_par[2],")")[[1]]
  par <- as.numeric(strsplit(par,",")[[1]])
  priors_par <- priors_par[1]
  if (priors_par=="normal") {
    #cat(paste0("rnorm(1,",par[1],",",par[2],")\n"))
    parSamp <- rnorm(1,par[1],par[2])
    if (class%in%c("sigma","sd"))
      while (parSamp<=0) parSamp <- rnorm(1,par[1],par[2])
  }
  if (priors_par=="lkj") {
    parSamp <- rethinking::rlkjcorr(1,2,par[1])[1,2]
  }
  if (priors_par=="student_t") {
    parSamp <- LaplacesDemon::rst(1, par[2], par[3], par[1])
    if (class%in%c("sigma","sd"))
      while (parSamp<=0) parSamp <- LaplacesDemon::rst(1, par[2], par[3], par[1])
  }
  if (priors_par=="exponential") {
    parSamp <- rexp(1, par[1])
  }
  parSamp
}

# F: Plotting Prior Predictive Distributions
## INPUTS
### prior_df: tibble|data frame = Extracted samples
## OUTPUT
### ggplot2 object including response and linear predictor distributions
priorpc_plot <- function(prior_df) {
  prior_summary <- prior_df %>%
    dplyr::group_by(var_label) %>%
    dplyr::summarise(
      q5 = quantile(value, probs = c(.05))[1],
      Mean = mean(value),
      sd = sd(value),
      q95 = quantile(value, probs = c(.95))[1],
      ypoint = hist(value, breaks=30, plot=F)$counts[which.max(hist(value, breaks=30, plot=F)$counts)]
    ) %>% dplyr::ungroup() %>%
    dplyr::mutate(
      mean_lab1 = stringr::str_c("E(y) == ", round(Mean,2)),
      mean_lab2 = stringr::str_c("E(y) == ", signif(Mean,2)),
      mean_lab = ifelse(nchar(mean_lab1)<nchar(mean_lab2), mean_lab1, mean_lab2),
      q5_lab = stringr::str_c("Q[.05] == ", round(q5,2)),
      q95_lab = stringr::str_c("Q[.95] == ", round(q95,2))
    )
  ggplot2::ggplot(prior_df, aes(x=value, fill=var_label)) +
    geom_histogram(color="black") +
    geom_vline(data=prior_summary, aes(xintercept=q5), linetype=2) +
    geom_vline(data=prior_summary, aes(xintercept=q95), linetype=2, color="green") +
    geom_text(data=prior_summary,
              aes(x=q95*1.1, y=ypoint, label=q5_lab), parse=T, hjust=0, color="black") +
    geom_text(data=prior_summary,
              aes(x=q95*1.1, y=ypoint*.875, label=mean_lab), parse=T, hjust=0, color="purple") +
    geom_text(data=prior_summary,
              aes(x=q95*1.1, y=ypoint*.75, label=q95_lab), parse=T, hjust=0, color="green") +
    scale_fill_manual(values=c("red","blue")) +
    scale_x_continuous(breaks = NULL) +
    facet_wrap(~var_label, scales = "free_x") +
    labs(title="Prior predictive distribution of linear term and responses",
         x = "") +
    guides(fill=F) +
    theme_bw()
}

# F: Simulation-based ranking
sbc_rank_param <- function(models, par, true_pars, thin_step, picks) {
  sbc_rank <- NA
  for (i in 1:length(models)) { # Compute SBC rank
    postgw <- brms::posterior_samples(models[[i]], par, exact_match = TRUE)
    idx_so <- seq(1,nrow(postgw),thin_step)
    d <- picks[i]
    sbc_rank[i] <- sum( true_pars[[par]][d] < postgw[[par]][idx_so] )
  }
  sbc_rank
}

# F: Simulation-based ranking on multiple parameters
sbc_rank_n <- function(models, pars, true_pars, thin_step, picks) {
  all_sbc_ranks <- purrr::map(pars, sbc_rank_param, models=models,
                              true_pars=true_pars, thin_step=thin_step,
                              picks=picks)
  sbc_df <- dplyr::bind_cols(all_sbc_ranks)
  names(sbc_df) <- pars
  sbc_df
}

# F: Plot SBC ranks
sbc_plots <- function(sbc_ranks, name) {
  max_rank <- max(purrr::map_int(sbc_ranks, max))
  hbins <- seq(0,max_rank+5,25)-0.5
  binom <- qbinom(c(0.005,0.5,0.995), nrow(sbc_ranks), 1/(length(hbins)-1))
  pdat <- data.frame(min=binom[1],med=binom[2],max=binom[3],x=c(-20,max_rank+20))
  c_light <- "#DCBCBC"; c_light_highlight <- "#C79999"
  c_mid   <- "#B97C7C"; c_mid_highlight   <- "#A25050"
  c_dark  <- "#8F2727"; c_dark_highlight  <- "#7C0000"
  est_pars <- sbc_ranks %>%
    tidyr::gather(variable, value)
  ggplot2::ggplot()+
    geom_ribbon(data=pdat, aes(x=x, ymin=min, ymax=max), fill="grey80") + 
    geom_line( data=pdat, aes(x=x,y=med)) +
    stat_bin(data=est_pars, aes(x=value), breaks=hbins,
             fill=c_dark, colour=c_dark_highlight) + scale_x_continuous(breaks=seq(0,max_rank+5,100))+
    labs(x="Prior Rank",y="", title=name) +
    facet_wrap(~variable, ncol=3) +
    theme_minimal()
}

# F: Model sensitivity scores (z-scores & contraction)
model_sensitivity <- function(models, par, true_pars, picks, prior_sd) {
  z_score <- contraction <- NA
  par2 <- stringr::str_split(par, "_")[[1]][2]
  for (i in 1:length(models)) {
    post_fe <- fixef(models[[i]])
    indx <- which(row.names(post_fe)==par2)
    post_mean <- post_fe[indx, 1]
    post_sd <- post_fe[indx, 2]
    d <- picks[i]
    z_score[i] <- (post_mean - true_pars[[par]][d]) / post_sd
    contraction[i] <- 1 - (post_sd^2 / prior_sd^2)
  }
  df <- tibble(
    par = par,
    z = z_score,
    contraction = contraction
  )
  df
}

# F: Batch model sensitivity scores
all_model_sensitivity <- function(models, pars, true_pars, picks, priors) {
  prior_sds <- rep(NA, length(pars))
  prior_sds <- purrr::map_dbl(
    pars,
    function(x){
      par1 <- stringr::str_split(x, "_")[[1]][1]
      par2 <- stringr::str_split(x, "_")[[1]][2]
      if (par2=="Intercept") par1 <- par2
      ss_df <- priors %>%
        dplyr::filter(class==par1)
      if (nrow(ss_df)==1) {
        the_prior <- ss_df$prior[1]
      } else if (par2%in%ss_df$coef) {
        the_prior <- ss_df$prior[which(ss_df$coef==par2)]
      } else {
        the_prior <- ss_df$prior[which(ss_df$coef=="")]
      }
      priors_par <- strsplit(the_prior,"(",fixed=TRUE)[[1]]
      par <- strsplit(priors_par[2],")")[[1]]
      par <- as.numeric(strsplit(par,",")[[1]])
      if (length(par)==3) {
        final <- par[3]
      } else if (length(par)==2) {
        final <- par[2]
      } else {
        final <- par[1]
      }
      final
    }
  )
  all_sens <- purrr::map2_df(
    pars,
    prior_sds,
    model_sensitivity,
    models=models,
    true_pars=true_pars,
    picks=picks
  )
  all_sens
}

# F: Plotting Model Sensitivity Scores
sensitivity_plot <- function(scores) {
  ggplot(data=scores, aes(x=contraction, y=z))+
    geom_hline(yintercept=0)+
    geom_point()+ xlim(c(0,1))+
    labs(x="Posterior Contraction",y="Posterior z-Score") +
    facet_wrap(~par, ncol=2) +
    theme_minimal()
}

# F: Plotting model adequacy
posterior_adequacy_plot <- function(samples, pars=NULL) {
  if (is.null(pars)) {
    params <- "b_"
    final_samples <- samples %>%
      dplyr::filter(grepl(params,variable))
  } else {
    final_samples <- samples %>%
      dplyr::filter(variable %in% pars)
  }
  p1 <- ggplot2::ggplot(final_samples, aes(y=variable, x=value)) +
    geom_halfeyeh() +
    geom_vline(xintercept = 0, linetype=2) +
    labs(x="Effects (in linear predictor space)",
         y="Parameter") +
    theme_minimal()
  p2 <- brms::pp_check(real_model, nsamples = 600) + 
    labs(x="Response variable") +
    theme_minimal()
  cowplot::plot_grid(p1, p2, ncol=2)
}