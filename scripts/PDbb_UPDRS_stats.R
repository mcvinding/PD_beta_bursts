## Statistics 2: within group correlation with clinical scores (MSD-UPDRS-III)
library(brms)
# library(BayesFactor)
Sys.setenv(PATH = paste("C:/Rtools/bin", Sys.getenv("PATH"), sep=";")) # Needed or there will be a pop-up everytime compiling C models.
Sys.setenv(BINPREF = "C:/Rtools/mingw_$(WIN)/bin/")
options(mc.cores=parallel::detectCores)                   # Try run with multicores !!!

# Define paths
wrkdir <- "C://Users//Mikkel//Documents//betabursts//groupanalysis"

setwd(wrkdir)
load(file='uData.Rdata')

################################################################################
## N events
u.neve.data.PD <- subset(u.neve.data, group=="ptns")

br.nev.uF1 <- brm(bf(nevent ~ F1+(1|subs:session)), data = u.neve.data.PD, family = poisson,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.nev.uF2 <- brm(bf(nevent ~ F2+(1|subs:session)), data = u.neve.data.PD, family = poisson,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.nev.uF3 <- brm(bf(nevent ~ F3+(1|subs:session)), data = u.neve.data.PD, family = poisson,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.nev.uF4 <- brm(bf(nevent ~ F4+(1|subs:session)), data = u.neve.data.PD, family = poisson,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.nev.uF5 <- brm(bf(nevent ~ F5+(1|subs:session)), data = u.neve.data.PD, family = poisson,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.nev.uF6 <- brm(bf(nevent ~ F6+(1|subs:session)), data = u.neve.data.PD, family = poisson,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.nev.uF7 <- brm(bf(nevent ~ F7+(1|subs:session)), data = u.neve.data.PD, family = poisson,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.nev.uFT <- brm(bf(nevent ~ Total+(1|subs:session)), data = u.neve.data.PD, family = poisson,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.nev.0 <- brm(bf(nevent ~ 1+(1|subs:session)), data = u.neve.data.PD, family = poisson,
                  save_all_pars = TRUE, iter = 5000, cores = 1)

setwd(wrkdir)
save(br.nev.uF1,br.nev.uF2,br.nev.uF3,br.nev.uF4,br.nev.uF5,br.nev.uF6,br.nev.uF7,br.nev.uFT,br.nev.0,
     file = 'updrs_neve_mods.R')

# hypothesis testing
load(file = 'updrs_neve_mods.R')

h1 <- hypothesis(br.nev.uF1, "F1>0")
h2 <- hypothesis(br.nev.uF2, "F2>0")
h3 <- hypothesis(br.nev.uF3, "F3>0")
h4 <- hypothesis(br.nev.uF4, "F4>0")
h5 <- hypothesis(br.nev.uF5, "F5>0")
h6 <- hypothesis(br.nev.uF6, "F6>0")
h7 <- hypothesis(br.nev.uF7, "F7<0")
hT <- hypothesis(br.nev.uFT, "Total>0")

# Bayes Factor
n.bf1 <- bayes_factor(br.nev.uF1,br.nev.0)
n.bf2 <- bayes_factor(br.nev.uF2,br.nev.0)
n.bf3 <- bayes_factor(br.nev.uF3,br.nev.0)
n.bf4 <- bayes_factor(br.nev.uF4,br.nev.0)
n.bf5 <- bayes_factor(br.nev.uF5,br.nev.0)
n.bf6 <- bayes_factor(br.nev.uF6,br.nev.0)
n.bf7 <- bayes_factor(br.nev.uF7,br.nev.0)
n.bfT <- bayes_factor(br.nev.uFT,br.nev.0)


################################################################################
## ITI 
u.iti.data.PD <- subset(u.iti.data, group=="ptns")

br.iti.uF1 <- brm(bf(eve.iti.ms ~ F1+(1|subs:session)), data = u.iti.data.PD, family = lognormal,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.iti.uF2 <- brm(bf(eve.iti.ms ~ F2+(1|subs:session)), data = u.iti.data.PD, family = lognormal,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.iti.uF3 <- brm(bf(eve.iti.ms ~ F3+(1|subs:session)), data = u.iti.data.PD, family = lognormal,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.iti.uF4 <- brm(bf(eve.iti.ms ~ F4+(1|subs:session)), data = u.iti.data.PD, family = lognormal,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.iti.uF5 <- brm(bf(eve.iti.ms ~ F5+(1|subs:session)), data = u.iti.data.PD, family = lognormal,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.iti.uF6 <- brm(bf(eve.iti.ms ~ F6+(1|subs:session)), data = u.iti.data.PD, family = lognormal,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.iti.uF7 <- brm(bf(eve.iti.ms ~ F7+(1|subs:session)), data = u.iti.data.PD, family = lognormal,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.iti.uFT <- brm(bf(eve.iti.ms ~ Total+(1|subs:session)), data = u.iti.data.PD, family = lognormal,
                  save_all_pars = TRUE, iter = 5000, cores = 1)

save(br.iti.uF1,br.iti.uF2,br.iti.uF3,br.iti.uF4,br.iti.uF5,br.iti.uF6,br.iti.uF7,br.iti.uFT,
     file = 'updrs_iti_mods.R')

# hypothesis testing
load(file = 'updrs_iti_mods.R')
h1 <- hypothesis(br.iti.uF1, "F1<0")
h2 <- hypothesis(br.iti.uF2, "F2<0")
h3 <- hypothesis(br.iti.uF3, "F3<0")
h4 <- hypothesis(br.iti.uF4, "F4<0")
h5 <- hypothesis(br.iti.uF5, "F5<0")
h6 <- hypothesis(br.iti.uF6, "F6<0")
h7 <- hypothesis(br.iti.uF7, "F7>0")
hT <- hypothesis(br.iti.uFT, "Total<0")

################################################################################
## Length of events
u.len.data.PD <- subset(u.len.data, group=="ptns")

br.len.uF1 <- brm(bf(eve.len.ms ~ F1+(1|subs:session)), data = u.len.data.PD, family = shifted_lognormal,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.len.uF2 <- brm(bf(eve.len.ms ~ F2+(1|subs:session)), data = u.len.data.PD, family = shifted_lognormal,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.len.uF3 <- brm(bf(eve.len.ms ~ F3+(1|subs:session)), data = u.len.data.PD, family = shifted_lognormal,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.len.uF4 <- brm(bf(eve.len.ms ~ F4+(1|subs:session)), data = u.len.data.PD, family = shifted_lognormal,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.len.uF5 <- brm(bf(eve.len.ms ~ F5+(1|subs:session)), data = u.len.data.PD, family = shifted_lognormal,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.len.uF6 <- brm(bf(eve.len.ms ~ F6+(1|subs:session)), data = u.len.data.PD, family = shifted_lognormal,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.len.uF7 <- brm(bf(eve.len.ms ~ F7+(1|subs:session)), data = u.len.data.PD, family = shifted_lognormal,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.len.uFT <- brm(bf(eve.len.ms ~ Total+(1|subs:session)), data = u.len.data.PD, family = shifted_lognormal,
                  save_all_pars = TRUE, iter = 5000, cores = 1)

save(br.len.uF1,br.len.uF2,br.len.uF3,br.len.uF4,br.len.uF5,br.len.uF6,br.len.uF7,br.len.uFT,
     file = 'updrs_len_mods.R')

# hypothesis testing
load(file = 'updrs_len_mods.R')
h1 <- hypothesis(br.len.uF1, "F1<0")
h2 <- hypothesis(br.len.uF2, "F2>0")
h3 <- hypothesis(br.len.uF3, "F3<0")
h4 <- hypothesis(br.len.uF4, "F4>0")
h5 <- hypothesis(br.len.uF5, "F5<0")
h6 <- hypothesis(br.len.uF6, "F6<0")
h7 <- hypothesis(br.len.uF7, "F7>0")
hT <- hypothesis(br.len.uFT, "Total>0")

################################################################################
## Max power
u.max.data.PD <- subset(u.max.data, group=="ptns")

br.max.uF1 <- brm(bf(eve.max ~ F1+(1|subs:session)), data = u.max.data.PD, family = shifted_lognormal,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.max.uF2 <- brm(bf(eve.max ~ F2+(1|subs:session)), data = u.max.data.PD, family = shifted_lognormal,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.max.uF3 <- brm(bf(eve.max ~ F3+(1|subs:session)), data = u.max.data.PD, family = shifted_lognormal,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.max.uF4 <- brm(bf(eve.max ~ F4+(1|subs:session)), data = u.max.data.PD, family = shifted_lognormal,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.max.uF5 <- brm(bf(eve.max ~ F5+(1|subs:session)), data = u.max.data.PD, family = shifted_lognormal,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.max.uF6 <- brm(bf(eve.max ~ F6+(1|subs:session)), data = u.max.data.PD, family = shifted_lognormal,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.max.uF7 <- brm(bf(eve.max ~ F7+(1|subs:session)), data = u.max.data.PD, family = shifted_lognormal,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.max.uFT <- brm(bf(eve.max ~ Total+(1|subs:session)), data = u.max.data.PD, family = shifted_lognormal,
                  save_all_pars = TRUE, iter = 5000, cores = 1)

save(br.max.uF1,br.max.uF2,br.max.uF3,br.max.uF4,br.max.uF5,br.max.uF6,br.max.uF7,br.max.uFT,
     file = 'updrs_max_mods.R')

# hypothesis testing
load(file = 'updrs_max_mods.R')
h1 <- hypothesis(br.max.uF1, "F1>0")
h2 <- hypothesis(br.max.uF2, "F2>0")
h3 <- hypothesis(br.max.uF3, "F3<0")
h4 <- hypothesis(br.max.uF4, "F4>0")
h5 <- hypothesis(br.max.uF5, "F5>0")
h6 <- hypothesis(br.max.uF6, "F6>0")
h7 <- hypothesis(br.max.uF7, "F7>0")
hT <- hypothesis(br.max.uFT, "Total>0")

# END