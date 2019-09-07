## Statistics 1: between group-session analysis
library(lme4)
library(brms)
Sys.setenv(PATH = paste("C:/Rtools/bin", Sys.getenv("PATH"), sep=";")) # Needed or there will be a pop-up everytime compiling C models.
Sys.setenv(BINPREF = "C:/Rtools/mingw_$(WIN)/bin/")
options(mc.cores=parallel::detectCores)                   # Try run with multicores !!!

# Define paths
wrkdir <- "Z://PD_motor//rest_ec//groupanalysis//"
setwd(wrkdir)

################################################################################
## Init.
newdata <- data.frame(group = factor(c("ctrl","ctrl","ptns","ptns")),
                      session = factor(c("1","2","1","2")))

################################################################################
## N event analysis
################################################################################
load(file = 'neve.RData')

# Maximal Likelihood analysis
mod3p <- glmer(nevent.min ~ group*session+(1|subs), data = neve.data, family='poisson')
mod2p <- update(mod3p, ~. -group:session)
mod1p <- update(mod2p, ~. -group)
mod0p <- update(mod1p, ~. -session)

anova(mod0p,mod1p,mod2p,mod3p)

# Make BRMS models
br.nev3 <- brm(nevent.min ~ group*session+(1|subs), data = neve.data, family = poisson, 
               save_all_pars = TRUE, iter = 5000)
br.nev2 <- brm(nevent.min ~ group+session+(1|subs), data = neve.data, family = poisson, 
               save_all_pars = TRUE, iter = 5000)
br.nev1 <- brm(nevent.min ~ session+(1|subs), data = neve.data, family = poisson, 
               save_all_pars = TRUE, iter = 5000)
br.nev0 <- brm(nevent.min ~ 1+(1|subs), data = neve.data, family = poisson, 
               save_all_pars = TRUE, iter = 5000)

# Save
setwd(wrkdir)
save(br.nev3,br.nev2,br.nev1,br.nev0, file='neve_analysis.RData')

# re-load
load(file='neve_analysis.RData')

# Model comparison
n.bf10 <- bayes_factor(br.nev1,br.nev0) #n.bf10
n.bf21 <- bayes_factor(br.nev2,br.nev1) #n.bf21
n.bf32 <- bayes_factor(br.nev3,br.nev2) #n.bf32

# Hypothesis testing
h1 <- hypothesis(br.nev3, "groupptns>0")
h2 <- hypothesis(br.nev3, "session2+groupptns:session2<0")
h3 <- hypothesis(br.nev3, "session2>0")

h1$hypothesis$Post.Prob*2
h2$hypothesis$Post.Prob*2
h3$hypothesis$Post.Prob*2

# Summaries
summary(br.nev3)

# Predict
pp.neve <- predict(br.nev3, newdata=newdata, summary=T, re_formula=NA, robust=T)

sam.nev3 <- posterior_samples(br.nev3, "^b")

exp(fixef(br.nev3)[2])-1                                          # Ptns session 1 (rel change)
quantile(exp(sam.nev3[,2]), probs=c(0.025,0.975))-1               # 95%CI
mean(exp(sam.nev3[,3]+sam.nev3[,4]))-1                            # Ptns session 2 (rel change)
quantile(exp(sam.nev3[,3]+sam.nev3[,4]), probs=c(0.025,0.975))-1  # 95%CI
mean(exp(sam.nev3[,3]))-1                                         # Ctrls session 2 (rel change)
quantile(exp(sam.nev3[,3]), probs=c(0.025,0.975))-1               # 95%CI

################################################################################
## Inter-event interval
################################################################################
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

# Hypothesis testing
h1 <- hypothesis(br.iti3, "groupptns<0")
h2 <- hypothesis(br.iti3, "session2+groupptns:session2>0")
h3 <- hypothesis(br.iti3, "session2<0")
h4 <- hypothesis(br.iti3, "groupptns:session2>0")

h1$hypothesis$Post.Prob*2
h2$hypothesis$Post.Prob*2
h3$hypothesis$Post.Prob*2
h4$hypothesis$Post.Prob*2

#Summaries
summary(br.iti3)

set.seed(100)
pp.iti.mean <- predict(br.iti3, newdata = newdata, re_formula=NA, probs=c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975), robust=F, summary=T)
pp.iti.median <- predict(br.iti3, newdata = newdata, re_formula=NA, probs=c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975), robust=T, summary=T)

# Make data frame for plotting
pp2.iti <- predict(br.iti3, newdata = newdata, re_formula=NA, summary=F)

pp.iti.long <- data.frame(
  pp = c(pp2.iti[,1],pp2.iti[,2],pp2.iti[,3],pp2.iti[,4]),
  label = rep(c('Controls 1', 'Controls 2', 'PD 1/OFF', 'PD 2/ON'), each=dim(pp2.iti)[1]),
  group = rep(c('ctrl','ctrl','ptns','ptns'), each=dim(pp2.iti)[1]),
  session = rep(c("1","2","1","2"), each=dim(pp2.iti)[1])
)
save(pp.iti.long, file='itieve_pp.RData')

# Posterior samples for summary
sam.iti3 <- posterior_samples(br.iti3, "^b")

mean(exp(sam.iti3[,3]+sam.iti3[,4]))-1                             # Ptns session 2 (rel change)
quantile(exp(sam.iti3[,3]+sam.iti3[,4]),probs=c(0.025,0.975))-1  # 95%CI

mean(exp(sam.iti3[,3]))-1                             # Ctrl session 2 (rel change)
quantile(exp(sam.iti3[,3]),probs=c(0.025,0.975))-1  # 95%CI

################################################################################
## Event duration
################################################################################
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

# Hypothesis testing
h1 <- hypothesis(br.len3, "groupptns<0")
h2 <- hypothesis(br.len3, "session2>0")
h3 <- hypothesis(br.len3, "session2+groupptns:session2<0")

#Summaries
summary(br.len3)

set.seed(100)
pp.len.mean <- predict(br.len3, newdata = newdata, re_formula= ~ndt, probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975), robust=F, summary=T)
pp.len.median <- predict(br.len3, newdata = newdata, re_formula= ~ndt, probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975), robust=T, summary=T)

# predict
pp2.len <- predict(br.len3) #, newdata = newdata, re_formula=NA, summary=F)

# Make data frame for plotting
pp.len.long <- data.frame(
  pp = c(pp2.len[,1],pp2.len[,2],pp2.len[,3],pp2.len[,4]),
  label = rep(c('Controls 1', 'Controls 2', 'PD 1/OFF', 'PD 2/ON'), each=dim(pp2.len)[1]),
  group = rep(c('ctrl','ctrl','ptns','ptns'), each=dim(pp2.len)[1]),
  session = rep(c("1","2","1","2"), each=dim(pp2.len)[1])
)
save(pp.len.long, file='leneve_pp.RData')

# Posterior samples
sam.len3 <- posterior_samples(br.len3, "^b")

mean(exp(sam.len3[,3]+sam.len3[,4]))-1                               # Ptns session 2 (%-change)
quantile(exp(sam.len3[,3]+sam.len3[,4]),probs=c(0.025,0.975))-1      # 95%CI

mean(exp(sam.len3[,3]))-1                               # Ptns session 2 (abs change)
quantile(exp(sam.len3[,3]),probs=c(0.025,0.975))-1      # 95%CI


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
h4 <- hypothesis(br.max3, "groupptns:session2>0")

#Summaries
summary(br.max3)
  
# predict
set.seed(100)
pp.max.mean <- predict(br.max3, newdata = newdata, re_formula= ~ndt, probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975), robust=F, summary=T)
pp.max.median <- predict(br.max3, newdata = newdata, re_formula= ~ndt, probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975), robust=T, summary=T)

pp2.max <- predict(br.max3, newdata = newdata, re_formula=NA, summary=F)

# Make data frame for plotting
pp.max.long <- data.frame(
  pp = c(pp2.max[,1],pp2.max[,2],pp2.max[,3],pp2.max[,4]),
  label = rep(c('Controls 1', 'Controls 2', 'PD 1/OFF', 'PD 2/ON'), each=dim(pp2.max)[1]),
  group = rep(c('ctrl','ctrl','ptns','ptns'), each=dim(pp2.max)[1]),
  session = rep(c("1","2","1","2"), each=dim(pp2.max)[1])
)
save(pp.max.long, file='maxeve_pp.RData')

# Posterior samples for summary
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