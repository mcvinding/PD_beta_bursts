## Statistics 1: between group-session analysis
# library(BayesFactor)
library(brms)
Sys.setenv(PATH = paste("C:/Rtools/bin", Sys.getenv("PATH"), sep=";")) # Needed or there will be a pop-up everytime compiling C models.
Sys.setenv(BINPREF = "C:/Rtools/mingw_$(WIN)/bin/")

# Define paths and 
wrkdir <- "Z://PD_motor//rest_ec//groupanalysis//"
setwd(wrkdir)

################################################################################
load(file = 'neve_ext.RData')
neve.data$nevent.min <- round(neve.data$nevent/3)

steps <- sort(unique(neve.data$steps))
steps <- steps[steps>=0.5 & steps <= 4.0]

BF10 <- rep(0, length(steps))
BF21 <- rep(0, length(steps))
BF32 <- rep(0, length(steps))

for (i in 1:length(steps)){
  print(paste('################ step:', steps[i], ' ##################'))
  tempdat <- neve.data[neve.data$steps==steps[i],]
  br.nev3 <- brm(nevent.min ~ group*session+(1|subs), data = tempdat, family = poisson, 
                 save_all_pars = TRUE, iter = 1000, cores = 1)
  br.nev2 <- brm(nevent.min ~ group+session+(1|subs), data = tempdat, family = poisson, 
                 save_all_pars = TRUE, iter = 1000, cores = 1)
  br.nev1 <- brm(nevent.min ~ session+(1|subs), data = tempdat, family = poisson, 
                 save_all_pars = TRUE, iter = 1000, cores = 1)
  br.nev0 <- brm(nevent.min ~ 1+(1|subs), data = tempdat, family = poisson, 
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
