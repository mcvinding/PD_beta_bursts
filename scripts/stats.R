## Statistics 1: between group-session analysis
library(BayesFactor)
library(lme4)
library(nlme)
library(multcomp)
library(brms)

# Define paths
wrkdir <- "Z://PD_motor//rest_ec//groupanalysis//"
wrkdir <- "C://Users//Mikkel//Documents//PD-proj_betaEvent//data"
setwd(wrkdir)
load(file = 'workspace.Rdata')

# Stats N event data: ANOVA
lme.neve <- lme(nevent ~ group*session, data=neve.data, random = ~1|subs)
anova(lme.neve)

summary(glht(lme.neve, linfct = mcp(list(session="Tukey",group="Tukey")), test = adjusted(type = "bonferroni")))
        
aov.neve <- aov(nevent ~ group*session+Error(subs),data=neve.data)
summary(aov.neve)
post.hoc <- TukeyHSD(aov.neve)

mod1 <- glmer(nevent ~ group*session+(1|subs), data = neve.data, family='poisson')
mod2 <- update(mod1, ~. -group:session)
mod3 <- update(mod2, ~. -group)
mod4 <- update(mod3, ~. -session)

anova(mod1,mod2,mod3,mod4)

mod1 <- lmer(nevent ~ group*session+(1|subs), data = neve.data)
mod2 <- update(mod1, ~. -group:session)
mod3 <- update(mod2, ~. -group)
mod4 <- update(mod3, ~. -session)

anova(mod1,mod2,mod3,mod4)

b <- anovaBF(nevent ~ group*session+subs, whichRandom = "subs", data = neve.data)
summary(b)
plot(b)

br.nev <- brm(bf(nevent ~ group*session+(1|subs)), data = neve.data, family = poisson)

##
iti.mod1 <- lmer(eve.iti.ms ~ group*session+(1|subs), data = itieve.data)

br.iti <- brm(bf(eve.iti.ms ~ group*session+(1|subs)), data = itieve.data, family = lognormal())
##
br.len <- brm(bf(~ group*session+(1|subs)), data = itieve.data, family = poisson)
##
#br.max <- brm(bf(~ group*session+(1|subs)), data = itieve.data, family = poisson)





