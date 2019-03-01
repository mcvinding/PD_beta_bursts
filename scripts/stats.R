## Statistics 1: between group-session analysis
library(BayesFactor)
library(lme4)
library(nlme)
# library(multcomp)
library(brms)

# Define paths
wrkdir <- "Z://PD_motor//rest_ec//groupanalysis//"
# wrkdir <- "C://Users//Mikkel//Documents//PD-proj_betaEvent//data"
setwd(wrkdir)
load(file = 'workspace.Rdata')
load(file = '.Rdata')

# Stats N event data: ANOVA
lme.neve <- lme(nevent ~ group*session, data=neve.data, random = ~1|subs)
anova(lme.neve)

summary(glht(lme.neve, linfct = mcp(list(session="Tukey",group="Tukey")), test = adjusted(type = "bonferroni")))
        
aov.neve <- aov(nevent ~ group*session+Error(subs),data=neve.data)
summary(aov.neve)
post.hoc <- TukeyHSD(aov.neve)

mod3p <- glmer(nevent ~ group*session+(1|subs), data = neve.data, family='poisson')
mod2p <- update(mod3p, ~. -group:session)
mod1p <- update(mod2p, ~. -group)
mod0p <- update(mod1p, ~. -session)

anova(mod0p,mod1p,mod2p,mod3p)

mod3g <- lmer(nevent ~ group*session+(1|subs), data = neve.data)
mod2g <- update(mod3g, ~. -group:session)
mod1g <- update(mod2g, ~. -group)
mod0g <- update(mod1g, ~. -session)

anova(mod0g,mod1g,mod2g,mod3g)

b <- anovaBF(nevent ~ session*group+subs, whichRandom = "subs", data = neve.data)
summary(b)
plot(b)

br.nev3 <- brm(bf(nevent ~ group*session+(1|subs)), data = neve.data, family = poisson, save_all_pars = TRUE)
br.nev2 <- brm(bf(nevent ~ group+session+(1|subs)), data = neve.data, family = poisson, save_all_pars = TRUE)
br.nev1 <- brm(bf(nevent ~ session+(1|subs)), data = neve.data, family = poisson, save_all_pars = TRUE)
br.nev0 <- brm(bf(nevent ~ 1+(1|subs)), data = neve.data, family = poisson, save_all_pars = TRUE)

bf10 <- bayes_factor(br.nev1,br.nev0)
bf21 <- bayes_factor(br.nev2,br.nev1)
bf32 <- bayes_factor(br.nev3,br.nev2)

bf10
bf21
bf32

br.nev3ga <- brm(bf(nevent ~ group*session+(1|subs)), 
                 data = neve.data, family = gaussian, save_all_pars = TRUE)
br.nev2ga <- brm(bf(nevent ~ group+session+(1|subs)), 
                 data = neve.data, family = gaussian, save_all_pars = TRUE)
br.nev1ga <- brm(bf(nevent ~ session+(1|subs)), 
                 data = neve.data, family = gaussian, save_all_pars = TRUE)
br.nev0ga <- brm(bf(nevent ~ 1+(1|subs)), 
                 data = neve.data, family = gaussian, save_all_pars = TRUE)

bf10ga <- bayes_factor(br.nev1ga,br.nev0ga)
bf21ga <- bayes_factor(br.nev2ga,br.nev1ga)
bf32ga <- bayes_factor(br.nev3ga,br.nev2ga)

bf10ga
bf21ga
bf32ga

summary(br.nev3)
summary(br.nev3ga)

## Time to event [NB: somehow there are zero values in this data]
itieve.data$log.eve.iti <- log(itieve.data$eve.iti.ms)
iti.mod3 <- lmer(eve.iti.ms ~ group*session+(1|subs), data = itieve.data, REML = F)
iti.mod2 <- lmer(eve.iti.ms ~ group+session+(1|subs), data = itieve.data, REML = F)
iti.mod1 <- lmer(eve.iti.ms ~ session+(1|subs), data = itieve.data, REML = F)
iti.mod0 <- lmer(eve.iti.ms ~ 1+(1|subs), data = itieve.data, REML = F)

anova(iti.mod0,iti.mod1,iti.mod2,iti.mod3)
summary(iti.mod3)

br.iti <- brm(bf(eve.iti.ms ~ group*session+(1|subs)), data = itieve.data, family = lognormal())

## Event length
leneve.data$log.eve.len <- log(leneve.data$eve.len)
len.mod3 <- lmer(eve.len.ms ~ group*session+(1|subs), data = leneve.data, REML = F)
len.mod2 <- lmer(eve.len.ms ~ group+session+(1|subs), data = leneve.data, REML = F)
len.mod1 <- lmer(eve.len.ms ~ session+(1|subs), data = leneve.data, REML = F)
len.mod0 <- lmer(eve.len.ms ~ 1+(1|subs), data = leneve.data, REML = F)

anova(len.mod0,len.mod1,len.mod2,len.mod3)
summary(len.mod3)

br.len3 <- brm(bf(eve.len.ms ~ group*session+(1|subs)), 
               data = leneve.data, family = lognormal, save_all_pars=T)
br.len2 <- brm(bf(eve.len.ms ~ group+session+(1|subs)), 
               data = leneve.data, family = lognormal, save_all_pars=T)
br.len1 <- brm(bf(eve.len.ms ~ session+(1|subs)), 
               data = leneve.data, family = lognormal, save_all_pars=T)
br.len0 <- brm(bf(eve.len.ms ~ 1+(1|subs)), 
               data = leneve.data, family = lognormal, save_all_pars=T)

len.bf10 <- bayes_factor(br.len1,br.len0)
len.bf21 <- bayes_factor(br.len2,br.len1)
len.bf32 <- bayes_factor(br.len3,br.len2)
len.bf10
len.bf21
len.bf32

##
br.max3 <- brm(bf(eve.max ~ group*session+(1|subs)), 
               data = maxeve.data, family = lognormal, save_all_pars=T)
br.max2 <- brm(bf(eve.max ~ group+session+(1|subs)), 
               data = maxeve.data, family = lognormal, save_all_pars=T)
br.max1 <- brm(bf(eve.max ~ session+(1|subs)), 
               data = maxeve.data, family = lognormal, save_all_pars=T)
br.max0 <- brm(bf(eve.max ~ 1+(1|subs)), 
               data = maxeve.data, family = lognormal, save_all_pars=T)

max.bf10 <- bayes_factor(br.max1,br.max0)
max.bf21 <- bayes_factor(br.max2,br.max1)
max.bf32 <- bayes_factor(br.max3,br.max2)
max.bf10
max.bf21
max.bf32

save.image(file='workspace.Rdata')


# Misc
qqnorm(resid(br.max3))
qqline(resid(br.max3))
       

save.image(file='workspace.Rdata')
