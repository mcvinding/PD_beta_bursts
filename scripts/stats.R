## Statistics 1: between group-session analysis
# library(BayesFactor)
library(lme4)
library(brms)
Sys.setenv(PATH = paste("C:/Rtools/bin", Sys.getenv("PATH"), sep=";")) # Needed or there will be a pop-up everytime compiling C models.
Sys.setenv(BINPREF = "C:/Rtools/mingw_$(WIN)/bin/")
options(mc.cores=parallel::detectCores)                   # Try run with multicores !!!

# Define paths
wrkdir <- "Z://PD_motor//rest_ec//groupanalysis//"
# wrkdir <- "C://Users//Mikkel//Documents//PD-proj_betaEvent//data"
setwd(wrkdir)
# load(file = '.Rdata')

################################################################################
## N event analysis
load(file = 'neve.RData')

# Maximal Likelihood analysis
mod3p <- glmer(nevent ~ group*session+(1|subs), data = neve.data, family='poisson')
mod2p <- update(mod3p, ~. -group:session)
mod1p <- update(mod2p, ~. -group)
mod0p <- update(mod1p, ~. -session)

anova(mod0p,mod1p,mod2p,mod3p)

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

# re-load
load(file='neve_analysis.RData')

# Model comparison
n.bf10 <- bayes_factor(br.nev1,br.nev0)
n.bf21 <- bayes_factor(br.nev2,br.nev1)
n.bf32 <- bayes_factor(br.nev3,br.nev2)
n.bf31 <- bayes_factor(br.nev3,br.nev3)

n.bf10
n.bf21
n.bf32
n.bf31

# Hypothesis testing
h <- hypothesis(br.nev3, "groupptns<0")
h <- hypothesis(br.nev3, "session2+groupptns:session2>0")
h <- hypothesis(br.nev3, "session2<0")

# Summaries
summary(br.nev3)

sam.nev3 <- posterior_samples(br.nev3, "^b")

exp(fixef(br.nev3)[2])-1                                          # Ptns session 1 (rel change)
quantile(exp(sam.nev3[,2]), probs=c(0.025,0.975))-1               # 95%CI
mean(exp(sam.nev3[,3]+sam.nev3[,4]))-1                            # Ptns session 2 (rel change)
quantile(exp(sam.nev3[,3]+sam.nev3[,4]), probs=c(0.025,0.975))-1  # 95%CI
mean(exp(sam.nev3[,3]))-1                                         # Ctrls session 2 (rel change)
quantile(exp(sam.nev3[,3]), probs=c(0.025,0.975))-1               # 95%CI

# Estimated number per min
exp(fixef(br.nev3)[1])/3                                          # Ctrl session 1
quantile(exp(sam.nev3[,1]), probs=c(0.025,0.975))/3               # 95%CI
exp(fixef(br.nev3)[1]+fixef(br.nev3)[2])/3                        # Ptns session 1
quantile(exp(sam.nev3[,1]+sam.nev3[,2]), probs=c(0.025,0.975))/3  # 95%CI
exp(fixef(br.nev3)[1]+fixef(br.nev3)[3])/3                        # Ctrl session 2
quantile(exp(sam.nev3[,1]+sam.nev3[,3]), probs=c(0.025,0.975))/3  # 95%CI
exp(fixef(br.nev3)[1]+fixef(br.nev3)[2]+fixef(br.nev3)[3]+fixef(br.nev3)[4])/3                        # Ctrl session 2
quantile(exp(sam.nev3[,1]+sam.nev3[,2]+sam.nev3[,3]+sam.nev3[,4]), probs=c(0.025,0.975))/3  # 95%CI

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

# Re-load
load('itieve_analysis.RData')

# Model comparison
iti.bf10 <- bayes_factor(br.iti1,br.iti0)
iti.bf21 <- bayes_factor(br.iti2,br.iti1)
iti.bf32 <- bayes_factor(br.iti3,br.iti2)

iti.bf10
iti.bf21
iti.bf32

# Hypothesis testing
h1 <- hypothesis(br.iti3, "groupptns>0")
h2 <- hypothesis(br.iti3, "session2+groupptns:session2<0")
h3 <- hypothesis(br.iti3, "session2>0")
hypothesis(br.iti3, "groupptns:session2>0")

#Summaries
summary(br.iti3)

sam.iti3 <- posterior_samples(br.iti3, "^b")

exp(fixef(br.iti3))                                               # All
mean(exp(sam.iti3[,1]+sam.iti3[,2]))                                # Ptns session 1 (ms)
quantile(sam.iti3[,1]+sam.iti3[,2], probs=c(0.025,0.975))           # 95%CI
mean(exp(sam.iti3[,1]+sam.iti3[,2]+sam.iti3[,3]+sam.iti3[,4]))      # Ptns session 2 (ms)
quantile(exp(sam.iti3[,1]+sam.iti3[,2]+sam.iti3[,3]+sam.iti3[,4]),probs=c(0.025,0.975))  # 95%CI
mean(exp(sam.iti3[,3]+sam.iti3[,4]))-1                               # Ptns session 2 (rel change)
quantile(exp(sam.iti3[,3]+sam.iti3[,4]),probs=c(0.025,0.975))-1  # 95%CI

mean(exp(sam.iti3[,1]+sam.iti3[,3]), probs=c(0.025,0.975))-1                   # Ctrls session 2 (rel change)
quantile(exp(sam.iti3[,1]+sam.iti3[,3]), probs=c(0.025,0.975))-1               # 95%CI
mean(exp(sam.iti3[,3]), probs=c(0.025,0.975))-1                   # Ctrls session 2 (ms)
quantile(exp(sam.iti3[,3]), probs=c(0.025,0.975))-1               # 95%CI


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

# Re-load
load(file='leneve_analysis.RData')

len.bf10 <- bayes_factor(br.len1,br.len0)
len.bf21 <- bayes_factor(br.len2,br.len1)
len.bf32 <- bayes_factor(br.len3,br.len2)

len.bf10
len.bf21
len.bf32

# Hypothesis testing
h1 <- hypothesis(br.len3, "groupptns<0")
h2 <- hypothesis(br.len3, "groupptns+session2+groupptns:session2>groupptns")
h3 <- hypothesis(br.len3, "session2>0")

#Summaries
summary(br.len3)

sam.len3 <- posterior_samples(br.len3, "^b")

exp(fixef(br.len3))                                               # All
mean(exp(sam.len3[,1]+sam.len3[,2]))                              # Ptns session 1 (ms)
quantile(exp(sam.len3[,1]+sam.len3[,2]), probs=c(0.025,0.975))    # 95%CI
mean(exp(sam.len3[,1]+sam.len3[,2]+sam.len3[,3]+sam.len3[,4]))                               # Ptns session 2 (abs change)
quantile(exp(sam.len3[,1]+sam.len3[,2]+sam.len3[,3]+sam.len3[,4]),probs=c(0.025,0.975))  # 95%CI

mean(exp(sam.len3[,1]+sam.len3[,3]))                               # Ptns session 2 (abs change)
quantile(exp(sam.len3[,1]+sam.len3[,3]),probs=c(0.025,0.975))  # 95%CI


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

# Re-load
load(file='maxeve_analysis.RData')
     
max.bf10 <- bayes_factor(br.max1,br.max0)
max.bf21 <- bayes_factor(br.max2,br.max1)
max.bf32 <- bayes_factor(br.max3,br.max2)

max.bf10
max.bf21
max.bf32

# Hypothesis testing
h1 <- hypothesis(br.max3, "groupptns<0")
h2 <- hypothesis(br.max3, "session2+groupptns:session2<0")
h3 <- hypothesis(br.max3, "session2<0")
h4 <- hypothesis(br.max3, "groupptns:session2<0")

#Summaries
summary(br.max3)

sam.max3 <- posterior_samples(br.max3, "^b")

mean(exp(sam.max3[,3]))-1                                       # Ctrl session 2 (rel)
quantile(exp(sam.max3[,3]), probs=c(0.025,0.975))-1             # 95%CI
mean(exp(sam.max3[,3]+sam.max3[,4]))-1                          # Ptns session 1 (rel)
quantile(exp(sam.max3[,3]+sam.max3[,4]), probs=c(0.025,0.975))-1           # 95%CIquantile(exp(sam.len3[,1]+sam.len3[,2]+sam.len3[,3]+sam.len3[,4]),probs=c(0.025,0.975))  # 95%CI

mean(exp(sam.max3[,1]))                                         # Ctrl session 1 (F)
quantile(exp(sam.max3[,1]), probs=c(0.025,0.975))               # 95%CI
mean(exp(sam.max3[,1]+sam.max3[,2]))                            # Ptns session 1 (F)
quantile(exp(sam.max3[,1]+sam.max3[,2]), probs=c(0.025,0.975))  # 95% CI
mean(exp(sam.max3[,1]+sam.max3[,3]))                            # Ctrl session 2 (F)
quantile(exp(sam.max3[,1]+sam.max3[,3]), probs=c(0.025,0.975))  # 95% CI
mean(exp(sam.max3[,1]+sam.max3[,2]+sam.max3[,3]+sam.max3[,4]))  # Ptns session 2 (F)
quantile(exp(sam.max3[,1]+sam.max3[,2]+sam.max3[,3]+sam.max3[,4]), probs=c(0.025,0.975))  # 95% CI

# Clean-up
save.image(file='workspace.Rdata')
# END