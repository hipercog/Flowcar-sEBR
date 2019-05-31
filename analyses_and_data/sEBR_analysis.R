library(plyr)
library(tidyverse)
library(reshape2)
library(haven)
library(nlme)
library(lme4)
library(lmerTest)
library(jtools)
library(huxtable)
library(ggfortify)
library(cowplot)

# setwd("/home/ben/Benslab/FLOW/extra_analyses/Flow_paper2_sEBR/")
setwd("/home/bcowley/Benslab/FLOW/papers/Flowcar-sEBR/analyses_and_data")

getTrialCors <- function(df, X, Y, id = TRUE)
{
  require(data.table)
  tr.cors <- as.matrix(spread(
    setDT(df)[
      ,.(cors = cor(x = get(X), y = get(Y), use="pair", method="k"))
      , by = .(Participant, Session)]
    , Session, `cors`))
  if (!id)
  {
    tr.cors <- tr.cors[,-(1)]
  }
  return(tr.cors)
}

corMat <- function(m1, m2, meth = "k", yous = "pair")
{
  corMat <- seq(1,9)
  for (ix in seq(1, 9))
  {
    idx <- !is.na(m1[ix,]) & !is.na(m2[ix,])
    corMat[ix] <- cor(m1[ix, idx], m2[ix, idx],  use = yous, method = meth)
    print(corMat[ix])
  }
  return(corMat)
}

interact_test <- function(df, dv, iv1, iv2, grp)
{
  fiti <- lmer(get(dv) ~ get(iv1)*get(iv2) + (1 | Participant), data = df)
  summ(fiti)
  ss <- sim_slopes(fiti, pred = iv1, modx = iv2, jnalpha = 0.05, johnson_neyman = TRUE, jnplot = TRUE)
  interact_plot(fiti, pred = iv2, modx = iv1, centered = "none", interval = TRUE, plot.points = TRUE, color.class = "Qual1")
  plot(ss)
  as_huxtable(ss)
}
# interact_test(df, dv, iv1, iv2, grp)

## DATA WRANGLING ----
fssl <- read_rds("fss_learning.RData")
blk <- read_rds("BLINK_ANALYSIS/blink_data.RData")
# runs <- read_rds('run_data_20180123.rds')
# fss <- read_sav("FSS.sav")

inc.ssn <- c(1, 5, 6, 7, 8)

brate <- cbind(blk$blinkrate.session1, blk$blinkrate.session5, blk$blinkrate.session6, blk$blinkrate.session7, blk$blinkrate.session8)
long.br <- gather(as.data.frame(t(brate)), "Parti", "brate")

# brate.ctr <- as.data.frame(scale(t(brate), scale = FALSE))
# long.br.ctr <- gather(brate.ctr, "Parti", "brate")

ssn.lc <- filter(fssl, Session %in% inc.ssn) %>% 
  group_by(Participant, Session) %>% 
  summarise(ssn.lc = as.numeric(lm(ln.duration ~ log(Run))$coefficients[2]), ssn.in = as.numeric(lm(ln.duration ~ log(Run))$coefficients[1]))

df <- filter(fssl, Session %in% inc.ssn) %>% 
  select(-Run, -cumrun, -slope, -intercept, -flow_z, -demand, -skill, -skilldemand) %>% 
  group_by(Participant, Session) %>% 
  summarise_all(funs(mean), na.rm = TRUE) 
df <- df %>%
  bind_cols(ssn.lc, long.br, filter(fssl, Session %in% inc.ssn & Run == 5) %>% select(demand, skill, skilldemand)) %>%
  select(1:2, brate, duration, ln.duration, learning_curve, distance, ssn.lc, ssn.in, everything(), -Participant1, -Session1, -Parti) %>%
  mutate(skidem = skill - demand)
# df <- mutate(df, fabal = fluency - absorption)
df.orig <- df

#### Centering steps: br * Participant, self-rep / perf x all ####
# df <- select(bind_cols(df, long.br.ctr), -Parti)
dfc <- group_by(df, Participant) %>% mutate_at(vars(brate), funs(scale), scale = FALSE)
dfc <- dfc %>% mutate_at(vars(-Participant, -Session, -brate), funs(scale), scale = FALSE)

#### WHICH DATAFRAME TO USE?? ####
df <- dfc
cor(df$skilldemand, df$skidem)

## Check variance of key variables
dfvar <- group_by(df, Participant) %>% 
  summarise_at(vars(brate, skilldemand, pi1, pi2), funs(var), na.rm = TRUE) %>% 
  mutate_at(2:5, funs(scale), center = FALSE) %>%
  mutate_at(2:5, funs(drop))
dfvar$Participant <- factor(dfvar$Participant, levels = dfvar$Participant[order(dfvar$brate)])
# bar plot of key variable variance
dfv.m <- melt(dfvar[1:3], id.vars = "Participant")
ggplot(dfv.m, aes(Participant, value)) + geom_bar(aes(fill = variable), position = "dodge", stat="identity")


## PCA
pca <- prcomp(filter(df[,3:20], !is.na(brate)), scale. = TRUE)
summary(pca)
autoplot(pca, loadings = TRUE, loadings.colour = 'blue', loadings.label = TRUE, loadings.label.size = 3)

## SIMPLE VISUALS ----
# histograms of everything
ggplot(gather(df[-(1:2)]), aes(value)) + 
  geom_histogram() + 
  facet_wrap(~key, scales = 'free_x')


# blk rate by Flow
ggplot(df, aes(x=flow, y=brate, color=Session)) + 
  geom_jitter(alpha=0.4) +
  facet_wrap(~Participant) +
  theme_bw() +
  labs(x = "avg Flow", y = "blink rate (centered)", color = "Session")
# blk rate by fluency
ggplot(df, aes(x=fluency, y=brate, color=Session)) + 
  geom_jitter(alpha=0.4) +
  facet_wrap(~Participant) +
  theme_bw() +
  labs(x = "avg fluency", y = "blink rate (centered)", color = "Session")
# blk rate by absorption
ggplot(df, aes(x=absorption, y=brate, color=Session)) + 
  geom_jitter(alpha=0.4) +
  facet_wrap(~Participant) +
  theme_bw() +
  labs(x = "avg absorption", y = "blink rate (centered)", color = "Session")
# blk rate by perceived importance 1: 'sthg important'
ggplot(df, aes(x=pi1, y=brate, color=Session)) + 
  geom_jitter(alpha=0.4) +
  facet_wrap(~Participant) +
  theme_bw() +
  labs(x = "perceived importance 1: 'sthg important'", y = "blink rate (centered)", color = "Session")
# blk rate by perceived importance 2: no mistakes
ggplot(df, aes(x=pi2, y=brate, color=Session)) + 
  geom_jitter(alpha=0.4) +
  facet_wrap(~Participant) +
  theme_bw() +
  labs(x = "perceived importance 2: no mistakes", y = "blink rate (centered)", color = "Session")
# blk rate by perceived importance 3: fear failure
ggplot(df, aes(x=pi3, y=brate, color=Session)) + 
  geom_jitter(alpha=0.4) +
  facet_wrap(~Participant) +
  theme_bw() +
  labs(x = "perceived importance 3: fear failure", y = "blink rate (centered)", color = "Session")
# blk rate by total perceived importance
ggplot(df, aes(x=pi_total, y=brate, color=Session)) + 
  geom_jitter(alpha=0.4) +
  facet_wrap(~Participant) +
  theme_bw() +
  labs(x = "avg perceived importance", y = "blink rate (centered)", color = "Session")
# blk rate by skill
ggplot(df, aes(x=skill, y=brate, color=Session)) + 
  geom_jitter(alpha=0.4) +
  facet_wrap(~Participant) +
  theme_bw() +
  labs(x = "skill", y = "blink rate (centered)", color = "Session")
# blk rate by demand
ggplot(df, aes(x=demand, y=brate, color=Session)) + 
  geom_jitter(alpha=0.4) +
  facet_wrap(~Participant) +
  theme_bw() +
  labs(x = "demand", y = "blink rate (centered)", color = "Session")
# blk rate by skilldemand
ggplot(df, aes(x=skilldemand, y=brate, color=Session)) + 
  geom_jitter(alpha=0.4) +
  facet_wrap(~Participant) +
  theme_bw() +
  labs(x = "skill v demand", y = "blink rate (centered)", color = "Session")
# blk rate by distance
ggplot(df, aes(x=distance, y=brate, color=Session)) + 
  geom_jitter(alpha=0.4) +
  facet_wrap(~Participant) +
  theme_bw() +
  labs(x = "lc - performance distance", y = "blink rate (centered)", color = "Session")
# blk rate by duration
ggplot(df, aes(x=duration, y=brate, color=Session)) + 
  geom_jitter(alpha=0.4) +
  facet_wrap(~Participant) +
  theme_bw() +
  labs(x = "duration", y = "blink rate (centered)", color = "Session")
# blk rate by session LC
ggplot(df, aes(x=ssn.lc, y=brate, color=Session)) + 
  geom_jitter(alpha=0.4) +
  facet_wrap(~Participant) +
  theme_bw() +
  labs(x = "session LC", y = "blink rate (centered)", color = "Session")
# blk rate by session intercept
ggplot(df, aes(x=ssn.in, y=brate, color=Session)) + 
  geom_jitter(alpha=0.4) +
  facet_wrap(~Participant) +
  theme_bw() +
  labs(x = "session intercept", y = "blink rate (centered)", color = "Session")



# duration by skill
ggplot(df, aes(x=skill, y=duration, color=Session)) + 
  geom_jitter(alpha=0.4) +
  facet_wrap(~Participant) +
  theme_bw() +
  labs(x = "skill", y = "duration", color = "Session")

## ANALYSIS: LINEAR MIXED MODELLING ----
# blink rate predicted by flow components
summary(lmer(brate ~ flow + (1 | Participant), data = df))
summary(lmer(brate ~ fluency + (1 | Participant), data = df))
summary(lmer(brate ~ absorption + (1 | Participant), data = df))
# blink rate predicted by perceived-importance components
summary(lmer(brate ~ pi1 + (1 | Participant), data = df))
summary(lmer(brate ~ pi2 + (1 | Participant), data = df))         # .
summary(lmer(brate ~ pi3 + (1 | Participant), data = df))
summary(lmer(brate ~ pi_total + (1 | Participant), data = df))
# blink rate predicted by skill-demand components
summary(lmer(brate ~ skill + (1 | Participant), data = df))       # *
summary(lmer(brate ~ demand + (1 | Participant), data = df))
summary(lmer(brate ~ skilldemand + (1 | Participant), data = df))
summary(lmer(brate ~ skidem + (1 | Participant), data = df))

# flow components predicted by blink rate
summary(lmer(flow ~ brate + (1 | Participant), data = df))
summary(lmer(fluency ~ brate + (1 | Participant), data = df))
summary(lmer(absorption ~ brate + (1 | Participant), data = df))

summary(lmer(pi1 ~ brate + (1 | Participant), data = df))
summary(lmer(pi2 ~ brate + (1 | Participant), data = df))        # .
summary(lmer(pi3 ~ brate + (1 | Participant), data = df))
summary(lmer(pi_total ~ brate + (1 | Participant), data = df))

summary(lmer(skill ~ brate + (1 | Participant), data = df))      # *
summary(lmer(demand ~ brate + (1 | Participant), data = df))
summary(lmer(skilldemand ~ brate + (1 | Participant), data = df))
summary(lmer(skidem ~ brate + (1 | Participant), data = df))

# interactions of blink rate and flow components
summary(lmer(brate ~ absorption*skill + (1 | Participant), data = df))
summary(lmer(brate ~ absorption*demand + (1 | Participant), data = df))
summary(lmer(brate ~ absorption*skidem + (1 | Participant), data = df))
summary(lmer(brate ~ absorption*skilldemand + (1 | Participant), data = df))
summary(lmer(brate ~ fluency*skill + (1 | Participant), data = df))
summary(lmer(brate ~ fluency*demand + (1 | Participant), data = df))
summary(lmer(brate ~ fluency*skidem + (1 | Participant), data = df))
summary(lmer(brate ~ fluency*skilldemand + (1 | Participant), data = df))
summary(lmer(brate ~ flow*skill + (1 | Participant), data = df))
summary(lmer(brate ~ flow*demand + (1 | Participant), data = df))
summary(lmer(brate ~ flow*skidem + (1 | Participant), data = df))
summary(lmer(brate ~ flow*skilldemand + (1 | Participant), data = df))

summary(lmer(brate ~ absorption*pi1 + (1 | Participant), data = df))       # **
summary(lmer(brate ~ absorption*pi2 + (1 | Participant), data = df))       # **
summary(lmer(brate ~ absorption*pi3 + (1 | Participant), data = df))
summary(lmer(brate ~ absorption*pi_total + (1 | Participant), data = df))  # **
summary(lmer(brate ~ fluency*pi1 + (1 | Participant), data = df))          # *
summary(lmer(brate ~ fluency*pi2 + (1 | Participant), data = df))          # *
summary(lmer(brate ~ fluency*pi3 + (1 | Participant), data = df))
summary(lmer(brate ~ fluency*pi_total + (1 | Participant), data = df))     # *
summary(lmer(brate ~ flow*pi1 + (1 | Participant), data = df))             # *
summary(lmer(brate ~ flow*pi2 + (1 | Participant), data = df))             # **
summary(lmer(brate ~ flow*pi3 + (1 | Participant), data = df))
summary(lmer(brate ~ flow*pi_total + (1 | Participant), data = df))        # **

summary(lmer(brate ~ skill*pi1 + (1 | Participant), data = df))            # *
summary(lmer(brate ~ skill*pi2 + (1 | Participant), data = df))
summary(lmer(brate ~ skill*pi3 + (1 | Participant), data = df))
summary(lmer(brate ~ skill*pi_total + (1 | Participant), data = df))
summary(lmer(brate ~ demand*pi1 + (1 | Participant), data = df))           # .
summary(lmer(brate ~ demand*pi2 + (1 | Participant), data = df))           # *
summary(lmer(brate ~ demand*pi3 + (1 | Participant), data = df))           # .
summary(lmer(brate ~ demand*pi_total + (1 | Participant), data = df))      # *
summary(lmer(brate ~ skidem*pi1 + (1 | Participant), data = df))           # **
summary(lmer(brate ~ skidem*pi2 + (1 | Participant), data = df))           # ***
summary(lmer(brate ~ skidem*pi3 + (1 | Participant), data = df))
summary(lmer(brate ~ skidem*pi_total + (1 | Participant), data = df))      # **
summary(lmer(brate ~ skilldemand*pi1 + (1 | Participant), data = df))      # **
summary(lmer(brate ~ skilldemand*pi2 + (1 | Participant), data = df))      # **
summary(lmer(brate ~ skilldemand*pi3 + (1 | Participant), data = df))
summary(lmer(brate ~ skilldemand*pi_total + (1 | Participant), data = df)) # **

# flipped interactions in blink rate, skilldemand, pi
summary(lmer(skilldemand ~ brate*pi_total + (1 | Participant), data = df)) # **
summary(lmer(pi_total ~ brate*skilldemand + (1 | Participant), data = df)) # *
summary(lmer(skidem ~ brate*pi_total + (1 | Participant), data = df))      # ***
summary(lmer(pi_total ~ brate*skidem + (1 | Participant), data = df))      # **

# Messing with random factors
# summary(lmer(brate ~ skilldemand + (skilldemand | Participant), data = df))
# summary(lmer(brate ~ skilldemand*pi_total + (1 | Participant) + (1 | skilldemand:Participant) + (1 | pi_total:Participant), data = df))

## INTERACTION TESTING ----

# BLINK RATE ~ SKILLDEMAND * PI-1
cor(df$skidem, df$pi1)
fiti <- lmer(brate ~ skidem*pi1 + (1 + Session | Participant), data = df)
summ(fiti)
ss <- sim_slopes(fiti, pred = "pi1", modx = "skidem", jnalpha = 0.05, johnson_neyman = TRUE, jnplot = TRUE, cond.int = TRUE)
interact_plot(fiti, pred = "pi1", modx = "skidem", centered = "none", interval = TRUE, plot.points = TRUE, color.class = "Qual2")
plot(ss)
as_huxtable(ss)
ss

# BLINK RATE ~ skidem * PI-2
cor(df$skidem, df$pi2)
fiti <- lmer(brate ~ skidem*pi2 + (1 | Participant), data = dfc)
summ(fiti)
ss <- sim_slopes(fiti, pred = "pi2", modx = "skidem", jnalpha = 0.05, johnson_neyman = TRUE, jnplot = TRUE, cond.int = TRUE)
interact_plot(fiti, pred = "pi2", modx = "skidem", centered = "none", interval = TRUE, plot.points = TRUE, color.class = "Qual2")
plot(ss)
as_huxtable(ss)
ss


# BLINK RATE ~ ABSORPTION * PI-1
cor(df$absorption, df$pi1)
lmm.br.ap1 <- lmer(brate ~ absorption*pi1 + (1 | Participant), data = df) # **
summary(lmm.br.ap1)
ss <- sim_slopes(lmm.br.ap1, pred = "pi1", modx = "absorption", jnalpha = 0.05, johnson_neyman = TRUE, jnplot = TRUE, cond.int = TRUE)
interact_plot(lmm.br.ap1, pred = "pi1", modx = "absorption", centered = "none", interval = TRUE, plot.points = TRUE, color.class = "Qual2")
plot(ss)
as_huxtable(ss)
ss

# BLINK RATE ~ ABSORPTION * PI-2
cor(df$absorption, df$pi2)
lmm.br.ap2 <- lmer(brate ~ absorption*pi2 + (1 | Participant), data = df) # ***
summary(lmm.br.ap2)
ss <- sim_slopes(lmm.br.ap2, pred = "pi2", modx = "absorption", jnalpha = 0.05, johnson_neyman = TRUE, jnplot = TRUE, cond.int = TRUE)
interact_plot(lmm.br.ap2, pred = "pi2", modx = "absorption", centered = "none", interval = TRUE, plot.points = TRUE, color.class = "Qual2")
plot(ss)
as_huxtable(ss)
ss


## INTERACTIONS - PERFORMANCE ----

# avg run duration predicted by blink rate, flow, skilldemand, pi
summary(lmer(duration ~ brate + (1 | Participant), data = df))
summary(lmer(duration ~ flow + (1 | Participant), data = df))
summary(lmer(duration ~ pi1 + (1 | Participant), data = df))
summary(lmer(duration ~ pi2 + (1 | Participant), data = df))
summary(lmer(duration ~ pi3 + (1 | Participant), data = df))
summary(lmer(duration ~ pi_total + (1 | Participant), data = df))
summary(lmer(duration ~ skill + (1 | Participant), data = df))       # ***
summary(lmer(duration ~ demand + (1 | Participant), data = df))
summary(lmer(duration ~ skidem + (1 | Participant), data = df))      # *
summary(lmer(duration ~ skilldemand + (1 | Participant), data = df))

# avg run duration predicted by interaction of brate * skill-demand 
summary(lmer(duration ~ brate*skill + (1 | Participant), data = df))
summary(lmer(duration ~ brate*demand + (1 | Participant), data = df))
summary(lmer(duration ~ brate*skilldemand + (1 | Participant), data = df))
summary(lmer(duration ~ brate*skidem + (1 | Participant), data = df))
# avg run duration predicted by interaction of brate * pi
summary(lmer(duration ~ brate*pi1 + (1 | Participant), data = df))
summary(lmer(duration ~ brate*pi2 + (1 | Participant), data = df))      # .
summary(lmer(duration ~ brate*pi3 + (1 | Participant), data = df))
summary(lmer(duration ~ brate*pi_total + (1 | Participant), data = df))
# check inverts
summary(lmer(skidem ~ duration*brate + (1 | Participant), data = df))
summary(lmer(pi_total ~ duration*brate + (1 | Participant), data = df))

# Three-way interactions
summary(lmer(duration ~ brate*skill*pi1 + (1 | Participant), data = df))       # **
summary(lmer(duration ~ brate*skill*pi2 + (1 | Participant), data = df))
summary(lmer(duration ~ brate*skill*pi3 + (1 | Participant), data = df))
summary(lmer(duration ~ brate*skill*pi_total + (1 | Participant), data = df))  # .
summary(lmer(duration ~ brate*skidem*pi1 + (1 | Participant), data = df))      # .
summary(lmer(duration ~ brate*skidem*pi2 + (1 | Participant), data = df))      # .
summary(lmer(duration ~ brate*skidem*pi3 + (1 | Participant), data = df))      # *
summary(lmer(duration ~ brate*skidem*pi_total + (1 | Participant), data = df))
summary(lmer(pi_total ~ duration*brate*skidem + (1 | Participant), data = df))

## INTERACTION TESTING - PERFORMANCE ----

# DURATION ~ BRATE * SKILL/DEMAND * PI_TOTAL
fiti <- lmer(duration ~ brate*skill*pi_total + (1 | Participant), data = df)
summ(fiti)
ss <- sim_slopes(fiti, pred = "brate", modx = "skill", mod2 = "pi_total", jnalpha = 0.05, johnson_neyman = TRUE, jnplot = TRUE, cond.int = TRUE)
interact_plot(fiti, pred = "brate", modx = "skill", mod2 = "pi_total", centered = "none", interval = TRUE, plot.points = TRUE, color.class = "Qual2")
plot(ss)
as_huxtable(ss)
ss

fiti <- lmer(duration ~ brate*demand*pi_total + (1 | Participant), data = df)
summ(fiti)
ss <- sim_slopes(fiti, pred = "brate", modx = "demand", mod2 = "pi_total", jnalpha = 0.05, johnson_neyman = TRUE, jnplot = TRUE, cond.int = TRUE)
interact_plot(fiti, pred = "brate", modx = "demand", mod2 = "pi_total", centered = "none", interval = TRUE, plot.points = TRUE, color.class = "Qual2")
plot(ss)
as_huxtable(ss)
ss

fiti <- lmer(duration ~ brate*skidem*pi_total + (1 | Participant), data = df)
summ(fiti)
ss <- sim_slopes(fiti, pred = "brate", modx = "pi_total", mod2 = "skidem", jnalpha = 0.05, johnson_neyman = TRUE, jnplot = TRUE, cond.int = TRUE)
interact_plot(fiti, pred = "brate", modx = "pi_total", mod2 = "skidem", centered = "none", interval = TRUE, plot.points = TRUE, color.class = "Qual2")
plot(ss)
as_huxtable(ss)
ss



# avg distance per session predicted by blink rate, flow, skilldemand, pi
# NB! avg distance per session tends to zero, not a good metric
summary(lmer(distance ~ brate + (1 | Participant), data = df))

summary(lmer(distance ~ flow + (1 | Participant), data = df))         # **
summary(lmer(distance ~ fluency + (1 | Participant), data = df))      # ***
summary(lmer(distance ~ absorption + (1 | Participant), data = df))   # .

summary(lmer(distance ~ pi1 + (1 | Participant), data = df))
summary(lmer(distance ~ pi2 + (1 | Participant), data = df))          # *
summary(lmer(distance ~ pi3 + (1 | Participant), data = df))
summary(lmer(distance ~ pi_total + (1 | Participant), data = df))

summary(lmer(distance ~ skill + (1 | Participant), data = df))
summary(lmer(distance ~ demand + (1 | Participant), data = df))
summary(lmer(distance ~ skidem + (1 | Participant), data = df))
summary(lmer(distance ~ skilldemand + (1 | Participant), data = df))  # **

summary(lmer(distance ~ brate*flow + (1 | Participant), data = df))   # .
summary(lmer(distance ~ brate*fluency + (1 | Participant), data = df))# *
summary(lmer(distance ~ brate*absorption + (1 | Participant), data = df))

summary(lmer(distance ~ brate*pi1 + (1 | Participant), data = df))
summary(lmer(distance ~ brate*pi2 + (1 | Participant), data = df))    # **
summary(lmer(distance ~ brate*pi3 + (1 | Participant), data = df))
summary(lmer(distance ~ brate*pi_total + (1 | Participant), data = df))

summary(lmer(distance ~ brate*skill + (1 | Participant), data = df))
summary(lmer(distance ~ brate*demand + (1 | Participant), data = df))
summary(lmer(distance ~ brate*skilldemand + (1 | Participant), data = df))
summary(lmer(distance ~ brate*skidem + (1 | Participant), data = df))

summary(lmer(distance ~ brate*fluency*skidem + (1 | Participant), data = df)) # *
summary(lmer(distance ~ brate*fluency*skilldemand + (1 | Participant), data = df)) # *

## INTERACTION TESTING - PERCEIVED PERFORMANCE ----

# DISTANCE ~ BRATE * SKILL/DEMAND * PI
fiti <- lmer(distance ~ brate*pi2 + (1 | Participant), data = df)    # **
summ(fiti)
ss <- sim_slopes(fiti, pred = "brate", modx = "pi2", jnalpha = 0.05, johnson_neyman = TRUE, jnplot = TRUE, cond.int = TRUE)
interact_plot(fiti, pred = "brate", modx = "pi2", centered = "none", interval = TRUE, plot.points = TRUE, color.class = "Qual2")
plot(ss)
as_huxtable(ss)
ss

fiti <- lmer(distance ~ brate*fluency*skidem + (1 | Participant), data = df)
summ(fiti)
ss <- sim_slopes(fiti, pred = "brate", modx = "fluency", mod2 = "skidem", jnalpha = 0.05, johnson_neyman = TRUE, jnplot = TRUE, cond.int = TRUE)
interact_plot(fiti, pred = "brate", modx = "fluency", mod2 = "skidem", centered = "none", interval = TRUE, plot.points = TRUE, color.class = "Qual2")
plot(ss)
as_huxtable(ss)
ss

fiti <- lmer(distance ~ brate*fluency*skilldemand + (1 | Participant), data = df)
summ(fiti)
ss <- sim_slopes(fiti, pred = "brate", modx = "fluency", mod2 = "skilldemand", jnalpha = 0.05, johnson_neyman = TRUE, jnplot = TRUE, cond.int = TRUE)
interact_plot(fiti, pred = "brate", modx = "fluency", mod2 = "skilldemand", centered = "none", interval = TRUE, plot.points = TRUE, color.class = "Qual2")
plot(ss)
as_huxtable(ss)
ss

# per-session log-log model fit SLOPE predicted by interaction of blink rate, FSS non-Flow components
summary(lmer(ssn.lc ~ brate + (1 | Participant), data = df))
summary(lmer(ssn.lc ~ brate*skilldemand + (1 | Participant), data = df)) # skilldemand = dot
summary(lmer(ssn.lc ~ brate*pi1 + (1 | Participant), data = df)) # pi1 = one star
summary(lmer(ssn.lc ~ brate*pi2 + (1 | Participant), data = df)) # pi2 = one star
summary(lmer(ssn.lc ~ pi1 + (1 | Participant), data = df)) # pi1 = one star
summary(lmer(ssn.lc ~ pi2 + (1 | Participant), data = df)) # pi2 = one star
summary(lmer(ssn.lc ~ skilldemand*pi1 + (1 | Participant), data = df))
summary(lmer(ssn.lc ~ skilldemand*pi2 + (1 | Participant), data = df))
summary(lmer(ssn.lc ~ brate*skilldemand*pi2 + (1 | Participant), data = df))

# per-session log-log model fit intercept predicted by interaction of blink rate, FSS non-Flow components
summary(lmer(ssn.in ~ brate + (1 | Participant), data = df))
summary(lmer(ssn.in ~ brate*skilldemand + (1 | Participant), data = df)) # skilldemand = dot
summary(lmer(ssn.in ~ brate*pi1 + (1 | Participant), data = df)) # pi1 = one star
summary(lmer(ssn.in ~ brate*pi2 + (1 | Participant), data = df)) # pi2 = one star
summary(lmer(ssn.in ~ pi1 + (1 | Participant), data = df)) # pi1 = one star
summary(lmer(ssn.in ~ pi2 + (1 | Participant), data = df)) # pi2 = one star
summary(lmer(ssn.in ~ skilldemand*pi1 + (1 | Participant), data = df))
summary(lmer(ssn.in ~ skilldemand*pi2 + (1 | Participant), data = df))
summary(lmer(ssn.in ~ brate*skilldemand*pi2 + (1 | Participant), data = df))


## ----
# ebrXflocp <- corMat(brate.ctr, as.matrix(spread(mn.Flo.cmp[, c(1, 2, 5)], "Session", "pi")[, 2:6]))


## ANALYSIS: CORRELATIONS
# BLINK RATE PER SESSION VS CORRELATION OF TRIAL-TIME AND ORDER PER SESSION
tr.cors <- getTrialCors(fssl, 'duration', 'Run', id=FALSE)
ebrXtrCor <- corMat(brate.ctr, tr.cors[, inc.ssn])
# BLINK RATE PER SESSION VS CORRELATION OF FLOW AND 'DISTANCE' (DIFF OF PERFORMANCE AND POWER LAW L.C.) PER SESSION
dF.cors <- getTrialCors(fssl, 'flow_z', 'distance', id=FALSE)
ebrXdfCor <- corMat(brate.ctr, dF.cors[, inc.ssn])

