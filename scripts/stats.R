## Statistics 1: between group-session analysis
# library(BayesFactor)
library(lme4)
library(brms)
Sys.setenv(PATH = paste("C:/Rtools/bin", Sys.getenv("PATH"), sep=";")) # Needed or there will be a pop-up everytime compiling C models.
Sys.setenv(BINPREF = "C:/Rtools/mingw_$(WIN)/bin/")

# Define paths
wrkdir <- "Z://PD_motor//rest_ec//groupanalysis//"
# wrkdir <- "C://Users//Mikkel//Documents//PD-proj_betaEvent//data"
setwd(wrkdir)
load(file = 'workspace.Rdata')
# load(file = '.Rdata')

################################################################################
## N event analysis
# Maximal Likelihood analysis
mod3p <- glmer(nevent ~ group*session+(1|subs), data = neve.data, family='poisson')
mod2p <- update(mod3p, ~. -group:session)
mod1p <- update(mod2p, ~. -group)
mod0p <- update(mod1p, ~. -session)

anova(mod0p,mod1p,mod2p,mod3p)

# b <- anovaBF(nevent ~ session*group+subs, whichRandom = "subs", data = neve.data)
# summary(b)
# plot(b)

# Make BRMS models
br.nev3 <- brm(bf(nevent ~ group*session+(1|subs)), data = neve.data, family = poisson, 
               save_all_pars = TRUE, iter = 5000)
br.nev2 <- brm(bf(nevent ~ group+session+(1|subs)), data = neve.data, family = poisson, 
               save_all_pars = TRUE, iter = 5000)
br.nev1 <- brm(bf(nevent ~ session+(1|subs)), data = neve.data, family = poisson, 
               save_all_pars = TRUE, iter = 5000)
br.nev0 <- brm(bf(nevent ~ 1+(1|subs)), data = neve.data, family = poisson, 
               save_all_pars = TRUE, iter = 5000)

# Save
setwd(wrkdir)
save(br.nev3,br.nev2,br.nev1,br.nev0, file='neve_analysis.RData')

# Model comparison
n.bf10 <- bayes_factor(br.nev1,br.nev0)
n.bf21 <- bayes_factor(br.nev2,br.nev1)
n.bf32 <- bayes_factor(br.nev3,br.nev2)
# n.bf31 <- bayes_factor(br.nev3,br.nev1)

n.bf10
n.bf21
n.bf32
# n.bf31

# Hypothesis testing
# ... 

summary(br.nev3)

################################################################################
## Time to event
load(file = 'itieve.RData')

# Maximal Likelihood analysis
iti.mod3 <- lmer(log.eve.iti ~ group*session+(1|subs), data = itieve.data, REML = F)
iti.mod2 <- lmer(log.eve.iti ~ group+session+(1|subs), data = itieve.data, REML = F)
iti.mod1 <- lmer(log.eve.iti ~ session+(1|subs), data = itieve.data, REML = F)
iti.mod0 <- lmer(log.eve.iti ~ 1+(1|subs), data = itieve.data, REML = F)

anova(iti.mod0,iti.mod1,iti.mod2,iti.mod3)
summary(iti.mod3)

# Make BRMS models
br.iti3 <- brm(bf(eve.iti.ms ~ group*session+(1|subs)), 
               data = itieve.data, family = lognormal,
               save_all_pars = TRUE, iter = 5000)
br.iti2 <- brm(bf(eve.iti.ms ~ group+session+(1|subs)), 
               data = itieve.data, family = lognormal,
               save_all_pars = TRUE, iter = 5000)
br.iti1 <- brm(bf(eve.iti.ms ~ session+(1|subs)),
               data = itieve.data, family = lognormal,
               save_all_pars = TRUE, iter = 5000)
br.iti0 <- brm(bf(eve.iti.ms ~ 1+(1|subs)),
               data = itieve.data, family = lognormal,
               save_all_pars = TRUE, iter = 5000)

# Save
setwd(wrkdir)
save(br.iti3,br.iti2,br.iti1,br.iti0, file='itieve_analysis.RData')

# Model comparison
iti.bf10 <- bayes_factor(br.iti1,br.iti0)
iti.bf21 <- bayes_factor(br.iti2,br.iti1)
iti.bf32 <- bayes_factor(br.iti3,br.iti2)

iti.bf10
iti.bf21
iti.bf32

# Hypothesis testing
# ... 


################################################################################
## Event length
load(file = 'leneve.RData')

# Maximal Likelihood analysis
len.mod3 <- lmer(log.eve.len ~ group*session+(1|subs), data = leneve.data, REML = F)
len.mod2 <- lmer(log.eve.len ~ group+session+(1|subs), data = leneve.data, REML = F)
len.mod1 <- lmer(log.eve.len ~ session+(1|subs), data = leneve.data, REML = F)
len.mod0 <- lmer(log.eve.len ~ 1+(1|subs), data = leneve.data, REML = F)

anova(len.mod0,len.mod1,len.mod2,len.mod3)
summary(len.mod3)

# Make BRMS models
br.len3 <- brm(bf(eve.len.ms ~ group*session+(1|subs)),
               data = leneve.data, family = shifted_lognormal,
               save_all_pars = TRUE, iter = 5000)
br.len2 <- brm(bf(eve.len.ms ~ group+session+(1|subs)),
               data = leneve.data, family = shifted_lognormal,
               save_all_pars = TRUE, iter = 5000)
br.len1 <- brm(bf(eve.len.ms ~ session+(1|subs)),
               data = leneve.data, family = shifted_lognormal,
               save_all_pars = TRUE, iter = 5000)
br.len0 <- brm(bf(eve.len.ms ~ 1+(1|subs)),
               data = leneve.data, family = shifted_lognormal, 
               save_all_pars = TRUE, iter = 5000)

# Save
setwd(wrkdir)
save(br.len3,br.len2,br.len1,br.len0, file='leneve_analysis.RData')

len.bf10 <- bayes_factor(br.len1,br.len0)
len.bf21 <- bayes_factor(br.len2,br.len1)
len.bf32 <- bayes_factor(br.len3,br.len2)

len.bf10
len.bf21
len.bf32

################################################################################
## Max power
load(file = 'maxeve.RData')

# Maximal Likelihood analysis
max.mod3 <- lmer(log.eve.max ~ group*session+(1|subs), data = maxeve.data, REML = F)
max.mod2 <- lmer(log.eve.max ~ group+session+(1|subs), data = maxeve.data, REML = F)
max.mod1 <- lmer(log.eve.max ~ session+(1|subs), data = maxeve.data, REML = F)
max.mod0 <- lmer(log.eve.max ~ 1+(1|subs), data = maxeve.data, REML = F)

anova(max.mod0,max.mod1,max.mod2,max.mod3)
summary(max.mod3)

# BRMS analysis
br.max3 <- brm(bf(eve.max ~ group*session+(1|subs)), 
               data = maxeve.data, family = shifted_lognormal, 
               save_all_pars = TRUE, iter = 5000)
br.max2 <- brm(bf(eve.max ~ group+session+(1|subs)), 
               data = maxeve.data, family = shifted_lognormal, 
               save_all_pars = TRUE, iter = 5000)
br.max1 <- brm(bf(eve.max ~ session+(1|subs)), 
               data = maxeve.data, family = shifted_lognormal,
               save_all_pars = TRUE, iter = 5000)
br.max0 <- brm(bf(eve.max ~ 1+(1|subs)),
               data = maxeve.data, family = shifted_lognormal, 
               save_all_pars = TRUE, iter = 5000)

# Save
setwd(wrkdir)
save(br.max3,br.max2,br.max1,br.max0, file='maxeve_analysis.RData')

max.bf10 <- bayes_factor(br.max1,br.max0)
max.bf21 <- bayes_factor(br.max2,br.max1)
max.bf32 <- bayes_factor(br.max3,br.max2)

# Not good model
# br.max3gau <- brm(bf(eve.max ~ group*session+(1|subs)), 
#                data = maxeve.data, family = gaussian, save_all_pars=T)
# br.max2gau <- brm(bf(eve.max ~ group+session+(1|subs)), 
#                data = maxeve.data, family = gaussian, save_all_pars=T)
# br.max1gau <- brm(bf(eve.max ~ session+(1|subs)), 
#                data = maxeve.data, family = gaussian, save_all_pars=T)
# br.max0gau <- brm(bf(eve.max ~ 1+(1|subs)), 
#                data = maxeve.data, family = gaussian, save_all_pars=T)
# 
# max.bf10.gau <- bayes_factor(br.max1gau,br.max0gau)
# max.bf21.gau <- bayes_factor(br.max2gau,br.max1gau)
# max.bf32.gau <- bayes_factor(br.max3gau,br.max2gau)

# br.max3ln <- brm(bf(log.eve.max ~ group*session+(1|subs)), 
#                data = maxeve.data, family = gaussian, save_all_pars=T)
# br.max2ln <- brm(bf(log.eve.max ~ group+session+(1|subs)), 
#                  data = maxeve.data, family = gaussian, save_all_pars=T)
# br.max1ln <- brm(bf(log.eve.max ~ session+(1|subs)), 
#                  data = maxeve.data, family = gaussian, save_all_pars=T)
# br.max0ln <- brm(bf(log.eve.max ~ 1+(1|subs)), 
#                  data = maxeve.data, family = gaussian, save_all_pars=T)
# 
# max.bf10.ln <- bayes_factor(br.max1ln,br.max0ln)
# max.bf21.ln <- bayes_factor(br.max2ln,br.max1ln)
# max.bf32.ln <- bayes_factor(br.max3ln,br.max2ln)

save.image(file='workspace.Rdata')
# END