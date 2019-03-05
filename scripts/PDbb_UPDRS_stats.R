## Statistics 2: within group correlation with clinical scores (MSD-UPDRS-III)
# Define paths
wrkdir <- "Z://PD_motor//rest_ec//groupanalysis//"
setwd(wrkdir)
load(file = 'workspace.Rdata')
# load(file = '.Rdata')

#################### Stats (mean) ##############################
load(file='uData.Rdata')

# # Not used
# library(lme4)
# 
# all.mod1 <- lmer(nevent~Total + (1|subs), data=u.neve.data, REML=F, subset = u.neve.data$group.x=="ptns")
# summary(all.mod1)
# all.mod0 <- lmer(nevent~1 + (1|subs), data=u.neve.data, REML=F, subset = u.neve.data$group.x=="ptns")
# anova(all.mod0,all.mod1)
# all.mod2 <- lmer(nevent~Total+session + (1|subs), data=u.neve.data, REML=F, subset = u.neve.data$group.x=="ptns")
# anova(all.mod0,all.mod1,all.mod2)
# all.mod3 <- lmer(nevent~Total*session + (1|subs), data=u.neve.data, REML=F, subset = u.neve.data$group.x=="ptns")
# anova(all.mod0,all.mod1,all.mod2,all.mod3)

# # Summary
# # pam1patients <- subset(pam1data,Sub_type=='patient')
# t.test(PD.data$Total~PD.data$session, paired=T)
# ttestBF(x=PD.data$Total[PD.data$session=="1"],y=PD.data$Total[PD.data$session=="2"], data=u.neve.data, paired=T, rscale='wide')
# 
# t.test(x=u.neve.data$UPDRS_off[u.neve.data$Sub_type=='control'] , y=u.neve.data$UPDRS_on, paired=F)
# ttestBF(x=pam1data$UPDRS_off[pam1data$Sub_type=='control'], y=PD.data$Total[PD.data$session=="2"], paired=F, rscale='wide')
# ttestBF(x=pam1data$UPDRS_off[pam1data$Sub_type=='control'], y=PD.data$Total[PD.data$session=="1"], paired=F, rscale='wide')
# 
# aggregate(u.neve.data$Total, by=list(u.neve.data$group.x,u.neve.data$session), mean)
# aggregate(u.neve.data$Total, by=list(u.neve.data$group.x,u.neve.data$session), sd)


# Bayes factor (TEST)
load(file='uData.Rdata')
library(BayesFactor)

bf0 <- lmBF(nevent~subs, data=PD.data, whichRandom="subs")
bf1 <- lmBF(nevent~session+subs, data=PD.data, whichRandom="subs")
bf1/bf0

# PD.data.x <- PD.data[complete.cases(PD.data),]
bf0 <- lmBF(nevent~subs, data=PD.data, whichRandom="subs")
bf1 <- lmBF(nevent~session*subs, data=PD.data, whichRandom=c("subs","session"))
bfT <- lmBF(nevent~Total+session+subs, data=PD.data, whichRandom=c("subs","session"))

bfT/bf1

bf.f1 <- lmBF(nevent~F1+session*subs, whichRandom=c("subs","session"),data=PD.data)
bf.f1/bf1
bf.f2 <- lmBF(nevent~F2+session*subs, whichRandom=c("subs","session"),data=PD.data)
bf.f2/bf1
bf.f3 <- lmBF(nevent~F3+session*subs, whichRandom=c("subs","session"),data=PD.data)
bf.f3/bf1
bf.f4 <- lmBF(nevent~F4+session*subs, whichRandom=c("subs","session"),data=PD.data)
bf.f4/bf1
bf.f5 <- lmBF(nevent~F5+session*subs, whichRandom=c("subs","session"),data=PD.data)
bf.f5/bf1
bf.f6 <- lmBF(nevent~F6+session*subs, whichRandom=c("subs","session"),data=PD.data)
bf.f6/bf1
bf.f7 <- lmBF(nevent~F7+session*subs, whichRandom=c("subs","session"),data=PD.data)
bf.f7/bf1

bf1 <- lmBF(nevent~session+subs, data=PD.data, whichRandom=c("subs"))
bfT <- lmBF(nevent~Total+session+subs, data=PD.data, whichRandom=c("subs"))
bfT/bf1
bf.f1 <- lmBF(nevent~F1+session+subs, whichRandom=c("subs"),data=PD.data)
bf.f1/bf1
bf.f2 <- lmBF(nevent~F2+session+subs, whichRandom=c("subs"),data=PD.data)
bf.f2/bf1
bf.f3 <- lmBF(nevent~F3+session+subs, whichRandom=c("subs"),data=PD.data)
bf.f3/bf1
bf.f4 <- lmBF(nevent~F4+session+subs, whichRandom=c("subs"),data=PD.data)
bf.f4/bf1
bf.f5 <- lmBF(nevent~F5+session+subs, whichRandom=c("subs"),data=PD.data)
bf.f5/bf1
bf.f6 <- lmBF(nevent~F6+session+subs, whichRandom=c("subs"),data=PD.data)
bf.f6/bf1
bf.f7 <- lmBF(nevent~F7+session+subs, whichRandom=c("subs"),data=PD.data)
bf.f7/bf1

mF1.post <- posterior(bf.f1, iterations = 1000)
summary(mF1.post)
mF2.post <- posterior(bf.f2, iterations = 1000)
summary(mF2.post)
mF3.post <- posterior(bf.f3, iterations = 1000)
summary(mF3.post)
mF4.post <- posterior(bf.f4, iterations = 1000)
summary(mF4.post)
mF5.post <- posterior(bf.f5, iterations = 1000)
summary(mF5.post)
mF6.post <- posterior(bf.f6, iterations = 1000)
summary(mF6.post)
mF7.post <- posterior(bf.f7, iterations = 1000)
summary(mF7.post)

mT.post <- posterior(bfT, iterations = 2000)
summary(mT.post)

# # Conventional correlation (just for show)
# cor.test(PD.data$F1,PD.data$nevent)
# cor.test(PD.data$F2,PD.data$nevent)
# cor.test(PD.data$F3,PD.data$nevent)
# cor.test(PD.data$F4,PD.data$nevent)
# cor.test(PD.data$F5,PD.data$nevent)
# cor.test(PD.data$F6,PD.data$nevent)
# cor.test(PD.data$F7,PD.data$nevent)
# cor.test(PD.data$Total,PD.data$nevent)

##
save.image(".RData")

# ### New division ###
# bf0 <- lmBF(nevent~subs, data=PD.data, whichRandom="subs")
# bf1 <- lmBF(nevent~session*subs, data=PD.data, whichRandom=c("subs"))
# 
# bf.tremor <- lmBF(nevent~tremor+session*subs, whichRandom=c("subs"),data=PD.data)
# bf.tremor/bf1
# bf.tremor.post <- posterior(bf.tremor, iterations = 2000)
# 
# bf.rigid <- lmBF(nevent~rigid+session*subs, whichRandom=c("subs","session"),data=PD.data)
# bf.rigid/bf1
# bf.rigid.post <- posterior(bf.rigid, iterations = 2000)
# 
# bf.axial <- lmBF(nevent~axial+session*subs, whichRandom=c("subs"),data=PD.data)
# bf.axial/bf1
# bf.axial.post <- posterior(bf.axial, iterations = 2000)
# 
# bf.brady <- lmBF(nevent~brady+session*subs, whichRandom=c("subs"),data=PD.data)
# bf.brady/bf1
# bf.brady.post <- posterior(bf.brady, iterations = 2000)
# 
# # Score
# bf0 <- lmBF(score~subs, data=u.long2, whichRandom="subs")
# bf1 <- lmBF(score ~ factor+session*subs, data=u.long2, whichRandom=c("subs"))
# bf.F <- lmBF(score ~ factor+nevent+session*subs, data=u.long2, whichRandom=c("subs"))
# bf.F/bf1
# bf.Fx <- lmBF(score ~ factor*rebound+session*id, data=u.long2, whichRandom=c("id","session"))
# bf.Fx/bf1
# 
# mT.post <- posterior(bf.Fx, iterations = 2000)
# summary(mT.post)

# # Conventional correlation (just for show)
# cor.test(PD.data$tremor,PD.data$nevent)
# cor.test(PD.data$rigid,PD.data$nevent)
# cor.test(PD.data$axial,PD.data$nevent)
# cor.test(PD.data$brady,PD.data$nevent)
# 
# plot(PD.data$tremor~PD.data$nevent)
# plot(PD.data$rigid~PD.data$nevent)
# plot(PD.data$axial,PD.data$nevent)
# plot(PD.data$brady,PD.data$nevent)

### BRMS (TEST)
library(brms)
Sys.setenv(PATH = paste("C:/Rtools/bin", Sys.getenv("PATH"), sep=";")) # Needed or there will be a pop-up everytime compiling C models.
Sys.setenv(BINPREF = "C:/Rtools/mingw_$(WIN)/bin/")

# N events
nmod0 <- brm(bf(nevent ~ 1+(session|subs)), data=PD.data, save_all_pars=TRUE)
nmodF1 <- brm(bf(nevent ~ F1+(session|subs)), data=PD.data,  save_all_pars=TRUE)
nmodF2 <- brm(bf(nevent ~ F2+(session|subs)), data=PD.data,  save_all_pars=TRUE)
nmodF3 <- brm(bf(nevent ~ F3+(session|subs)), data=PD.data,  save_all_pars=TRUE)
nmodF4 <- brm(bf(nevent ~ F4+(session|subs)), data=PD.data,  save_all_pars=TRUE)
nmodF5 <- brm(bf(nevent ~ F5+(session|subs)), data=PD.data,  save_all_pars=TRUE)
nmodF6 <- brm(bf(nevent ~ F6+(session|subs)), data=PD.data,  save_all_pars=TRUE)
nmodF7 <- brm(bf(nevent ~ F7+(session|subs)), data=PD.data,  save_all_pars=TRUE)
nmodT <- brm(bf(nevent ~ Total+(session|subs)), data=PD.data,  save_all_pars=TRUE)  

bayes_factor(nmodF1,nmod0)
bayes_factor(nmodF2,nmod0)
bayes_factor(nmodF3,nmod0)
bayes_factor(nmodF4,nmod0)
bayes_factor(nmodF5,nmod0)
bayes_factor(nmodF6,nmod0)
bayes_factor(nmodF7,nmod0)
bayes_factor(nmodT,nmod0)


############################## MAX POWER #####################################
# TEST TEST TEST
pmod0 <- brm(bf(nevent ~ 1+(session|subs)), data=PD.data, save_all_pars=TRUE)
pmodF1 <- brm(bf(nevent ~ F1+(session|subs)), data=PD.data, family=poisson,  save_all_pars=TRUE)
pmodF2 <- brm(bf(nevent ~ F2+(session|subs)), data=PD.data, family=poisson, save_all_pars=TRUE)
pmodF3 <- brm(bf(nevent ~ F3+(session|subs)), data=PD.data, family=poisson,  save_all_pars=TRUE)
pmodF4 <- brm(bf(nevent ~ F4+(session|subs)), data=PD.data, family=poisson,  save_all_pars=TRUE)
pmodF5 <- brm(bf(nevent ~ F5+(session|subs)), data=PD.data, family=poisson,  save_all_pars=TRUE)
pmodF6 <- brm(bf(nevent ~ F6+(session|subs)), data=PD.data, family=poisson,  save_all_pars=TRUE)
pmodF7 <- brm(bf(nevent ~ F7+(session|subs)), data=PD.data, family=poisson,  save_all_pars=TRUE)
pmodT <- brm(bf(nevent ~ Total+(session|subs)), data=PD.data, family=poisson,  save_all_pars=TRUE)  

bayes_factor(nmodF1p,nmod0p)
bayes_factor(nmodF2p,nmod0p)
bayes_factor(nmodF3p,nmod0p)
bayes_factor(nmodF4p,nmod0p)
bayes_factor(nmodF5p,nmod0p)
bayes_factor(nmodF6p,nmod0p)
bayes_factor(nmodF7p,nmod0p)
bayes_factor(nmodTp,nmod0p)

# Max power
nmod0 <- brm(bf(eve.max ~ 1+(session|subs)), data=PD.data, save_all_pars=TRUE)
nmodF1 <- brm(bf(nevent ~ F1+(session|subs)), data=PD.data,  save_all_pars=TRUE)
nmodF2 <- brm(bf(nevent ~ F2+(session|subs)), data=PD.data,  save_all_pars=TRUE)
nmodF3 <- brm(bf(nevent ~ F3+(session|subs)), data=PD.data,  save_all_pars=TRUE)
nmodF4 <- brm(bf(nevent ~ F4+(session|subs)), data=PD.data,  save_all_pars=TRUE)
nmodF5 <- brm(bf(nevent ~ F5+(session|subs)), data=PD.data,  save_all_pars=TRUE)
nmodF6 <- brm(bf(nevent ~ F6+(session|subs)), data=PD.data,  save_all_pars=TRUE)
nmodF7 <- brm(bf(nevent ~ F7+(session|subs)), data=PD.data,  save_all_pars=TRUE)
nmodT <- brm(bf(nevent ~ Total+(session|subs)), data=PD.data,  save_all_pars=TRUE)  


br.max0 <- brm(bf(eve.max ~ 1+(1|subs)), 
               data = maxeve.data, family = lognormal, save_all_pars=T)
