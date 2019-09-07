## Statistics 1: between group-session analysis
# library(BayesFactor)
library(lme4)
library(brms)
library(BayesFactor)
Sys.setenv(PATH = paste("C:/Rtools/bin", Sys.getenv("PATH"), sep=";")) # Needed or there will be a pop-up everytime compiling C models.
Sys.setenv(BINPREF = "C:/Rtools/mingw_$(WIN)/bin/")

# Define paths
wrkdir <- "Z://PD_motor//rest_ec//groupanalysis//"
setwd(wrkdir)
load(file = 'neve_ext.RData')
neve.data$nevent.min <- neve.data$nevent/3

################################################################################
steps <- sort(unique(neve.data$steps))
BF10 <- rep(0, length(steps))
BF21 <- rep(0, length(steps))
BF32 <- rep(0, length(steps))

# for (i in 1:length(steps)){
#   tempdat <- subset(neve.data, steps==steps[i])
#   Bmod3 <- lmBF(nevent.min ~ group*session+subs, data=tempdat, whichRandom='subs')
#   Bmod2 <- lmBF(nevent.min ~ group+session+subs, data=tempdat, whichRandom='subs')
#   Bmod1 <- lmBF(nevent.min ~ session+subs, data=tempdat, whichRandom='subs')
#   Bmod0 <- lmBF(nevent.min ~ subs, data=tempdat, whichRandom='subs')
#   
#   bf10 <- Bmod1/Bmod0
#   bf21 <- Bmod2/Bmod1
#   bf32 <- Bmod3/Bmod2
#   
#   BF10[i] <- bf10@bayesFactor$bf
#   BF21[i] <- bf21@bayesFactor$bf
#   BF32[i] <- bf32@bayesFactor$bf
# }
# 
# BFdata <- data.frame(steps=steps,
#                      BF10=BF10,
#                      BF21=BF21,
#                      BF32=BF32)

# Make BRMS models
steps <- sort(unique(neve.data$steps))
steps <- steps[steps>=0.5 & steps <= 4.0]
BF10 <- rep(0, length(steps))
BF21 <- rep(0, length(steps))
BF32 <- rep(0, length(steps))

for (i in 1:length(steps)){
  tempdat <- subset(neve.data, steps==steps[i])
  br.nev3 <- brm(bf(nevent ~ group*session+(1|subs)), data = tempdat, family = poisson, 
                 save_all_pars = TRUE, iter = 1000, cores = 1)
  br.nev2 <- brm(bf(nevent ~ group+session+(1|subs)), data = tempdat, family = poisson, 
                 save_all_pars = TRUE, iter = 1000, cores = 1)
  br.nev1 <- brm(bf(nevent ~ session+(1|subs)), data = tempdat, family = poisson, 
                 save_all_pars = TRUE, iter = 1000, cores = 1)
  br.nev0 <- brm(bf(nevent ~ 1+(1|subs)), data = tempdat, family = poisson, 
                 save_all_pars = TRUE, iter = 1000, cores = 1)

  n.bf10 <- bayes_factor(br.nev1,br.nev0)
  n.bf21 <- bayes_factor(br.nev2,br.nev1)
  n.bf32 <- bayes_factor(br.nev3,br.nev2)
  
  BF10[i] <- n.bf10$bf
  BF21[i] <- n.bf21$bf
  BF32[i] <- n.bf32$bf
}

# Save
setwd(wrkdir)
save(BF10,BF21,BF32, file='range_bf.RData')
