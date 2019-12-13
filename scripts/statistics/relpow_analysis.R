## Statistics 1b: between group-session analysis of relative beta power
library(R.matlab)
library(BayesFactor)
library(brms)
Sys.setenv(PATH = paste("C:/Rtools/bin", Sys.getenv("PATH"), sep=";")) # Needed or there will be a pop-up everytime compiling C models.

# Define paths
# wrkdir <- "Z://PD_motor//rest_ec//groupanalysis//"
wrkdir <- "C://Users//Mikkel//Documents//betabursts//groupanalysis"
setwd(wrkdir)

# Load data
temp <- readMat("B_relpow.mat")
ptns.relpow1 <- temp$ptns.relpow1
ptns.relpow2 <- temp$ptns.relpow2
ctrl.relpow1 <- temp$ctrl.relpow1
ctrl.relpow2 <- temp$ctrl.relpow2

# Arrange data
relpow <- c(ptns.relpow1, ptns.relpow2, ctrl.relpow1, ctrl.relpow2)
subs <- as.factor(c(rep(1:length(ptns.relpow1),2), rep(1:length(ctrl.relpow1)+20,2)))
session <- as.factor(c(rep(c(1,2), each=length(ptns.relpow1)), rep(c(1,2), each=length(ctrl.relpow1)*2)))
group <- as.factor(c(rep(1, length(ptns.relpow1)*2), rep(2, length(ctrl.relpow1)*2)))

rel.dat <- data.frame(relpow = c(ptns.relpow1, ptns.relpow2, ctrl.relpow1, ctrl.relpow2),
                      subj = as.factor(c(rep(1:length(ptns.relpow1),2), rep(1:length(ctrl.relpow1)+20,2))),
                      session = as.factor(c(rep(c(1,2), each=length(ptns.relpow1)), rep(c(1,2), each=length(ctrl.relpow1)))),
                      group = as.factor(c(rep(1, length(ptns.relpow1)*2), rep(2, length(ctrl.relpow1)*2)))
                      )

rel.dat$group <- revalue(rel.dat$group, c("1"="ptns", "2"="ctrl"))
save(rel.dat, file = 'reldat.RData')

# # Make BRMS models
# br.rel3 <- brm(bf(relpow ~ group*session+(1|subj)), data = rel.dat, family = gaussian, 
#                save_all_pars = TRUE, iter = 5000)
# br.rel2 <- brm(bf(relpow ~ group+session+(1|subj)), data = rel.dat, family = gaussian, 
#                save_all_pars = TRUE, iter = 5000)
# br.rel1 <- brm(bf(relpow ~ session+(1|subj)), data = rel.dat, family = gaussian, 
#                save_all_pars = TRUE, iter = 5000)
# br.rel0 <- brm(bf(relpow ~ 1+(1|subj)), data = rel.dat, family = gaussian, 
#                save_all_pars = TRUE, iter = 5000)
# 
# # Save
# setwd(wrkdir)
# save(br.rel3,br.rel2,br.rel1,br.rel0, file='relpow_analysis.RData')
# 
# # Model comparison
# r.bf10 <- bayes_factor(br.rel1,br.rel0)
# r.bf21 <- bayes_factor(br.rel2,br.rel1)
# r.bf32 <- bayes_factor(br.rel3,br.rel2)
# r.bf31 <- bayes_factor(br.rel3,br.rel1)

# Summary
# summary(br.rel3)


# Bayesian "ANOVA"
# av <- anovaBF(relpow~group*session+subj, data=rel.dat, whichRandom="subj")

# Bayesian "T-tests"
ttestBF(ptns.relpow1, ptns.relpow2, paird = TRUE)
ttestBF(ctrl.relpow1, ctrl.relpow2, paird = TRUE)
ttestBF(ptns.relpow1, ctrl.relpow1, paird = FALSE)
ttestBF(ptns.relpow2, ctrl.relpow2, paird = FALSE)

ptns.diff <- ptns.relpow1-ptns.relpow2
ctrl.diff <- ctrl.relpow1-ctrl.relpow2

ttestBF(ptns.diff, ctrl.diff, paird = FALSE)

# Student's T-test
t.test(ptns.relpow1, ptns.relpow2, paird = TRUE)
t.test(ctrl.relpow1, ctrl.relpow2, paird = TRUE)
t.test(ptns.relpow1, ctrl.relpow1, paird = FALSE)
t.test(ptns.relpow2, ctrl.relpow2, paird = FALSE)

t.test(ptns.diff, ctrl.diff, paird = FALSE)
