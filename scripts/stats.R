## Statistics 1: between group-session analysis
library(BayesFactor)
library(lme4)
library(nlme)
library(multcomp)

# Define paths
wrkdir <- "Z://PD_motor//rest_ec//groupanalysis//"
setwd(wrkdir)
load(file = 'workspace.Rdata')

# Stats N event data: ANOVA
lme.neve <- lme(nevent ~ group*session, data=neve.data, random = ~1|subs)
anova(lme.neve)

summary(glht(lme.neve, linfct = mcp(list(session="Tukey",group="Tukey")), test = adjusted(type = "bonferroni"))
        
aov.neve <- aov(nevent ~ group*session+Error(subs),data=neve.data))
summary(aov.neve)
post.hoc <- TukeyHSD(aov.neve)

mod <- glmer(nevent ~ group*session+(1|subs), data = neve.data, family='poisson')

b <- anovaBF(nevent ~ group*session+subs, whichRandom = "subs", data = neve.data)

