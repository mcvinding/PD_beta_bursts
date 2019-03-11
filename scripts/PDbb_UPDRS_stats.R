## Statistics 2: within group correlation with clinical scores (MSD-UPDRS-III)
library(brms)
# library(BayesFactor)
Sys.setenv(PATH = paste("C:/Rtools/bin", Sys.getenv("PATH"), sep=";")) # Needed or there will be a pop-up everytime compiling C models.
Sys.setenv(BINPREF = "C:/Rtools/mingw_$(WIN)/bin/")
options(mc.cores=parallel::detectCores)                   # Try run with multicores !!!

# Define paths
wrkdir <- "Z://PD_motor//rest_ec//groupanalysis//"
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

seted(wrkdir)
save(br.nev.uF1,br.nev.uF2,br.nev.uF3,br.nev.uF4,br.nev.uF5,br.nev.uF6,br.nev.uF7,br.nev.uFT,
     file = 'updrs_neve_mods.R')

# hypothesis testing
# ...

################################################################################
## ITI 
u.iti.data.PD <- subset(u.iti.data, group=="ptns")

br.iti.uF1 <- brm(bf(eve.iti ~ F1+(1|subs:session)), data = u.iti.data.PD, family = lognormal,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.iti.uF2 <- brm(bf(eve.iti ~ F2+(1|subs:session)), data = u.iti.data.PD, family = lognormal,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.iti.uF3 <- brm(bf(eve.iti ~ F3+(1|subs:session)), data = u.iti.data.PD, family = lognormal,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.iti.uF4 <- brm(bf(eve.iti ~ F4+(1|subs:session)), data = u.iti.data.PD, family = lognormal,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.iti.uF5 <- brm(bf(eve.iti ~ F5+(1|subs:session)), data = u.iti.data.PD, family = lognormal,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.iti.uF6 <- brm(bf(eve.iti ~ F6+(1|subs:session)), data = u.iti.data.PD, family = lognormal,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.iti.uF7 <- brm(bf(eve.iti ~ F7+(1|subs:session)), data = u.iti.data.PD, family = lognormal,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.iti.uFT <- brm(bf(eve.iti ~ Total+(1|subs:session)), data = u.iti.data.PD, family = lognormal,
                  save_all_pars = TRUE, iter = 5000, cores = 1)

save(br.iti.uF1,br.iti.uF2,br.iti.uF3,br.iti.uF4,br.iti.uF5,br.iti.uF6,br.iti.uF7,br.iti.uFT,
     file = 'updrs_iti_mods.R')

# hypothesis testing
# ...


################################################################################
## Length of events
u.len.data.PD <- subset(u.len.data, group=="ptns")

br.len.uF1 <- brm(bf(eve.len ~ F1+(1|subs:session)), data = u.len.data.PD, family = shifted_lognormal,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.len.uF2 <- brm(bf(eve.len ~ F2+(1|subs:session)), data = u.len.data.PD, family = shifted_lognormal,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.len.uF3 <- brm(bf(eve.len ~ F3+(1|subs:session)), data = u.len.data.PD, family = shifted_lognormal,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.len.uF4 <- brm(bf(eve.len ~ F4+(1|subs:session)), data = u.len.data.PD, family = shifted_lognormal,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.len.uF5 <- brm(bf(eve.len ~ F5+(1|subs:session)), data = u.len.data.PD, family = shifted_lognormal,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.len.uF6 <- brm(bf(eve.len ~ F6+(1|subs:session)), data = u.len.data.PD, family = shifted_lognormal,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.len.uF7 <- brm(bf(eve.len ~ F7+(1|subs:session)), data = u.len.data.PD, family = shifted_lognormal,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.len.uFT <- brm(bf(eve.len ~ Total+(1|subs:session)), data = u.len.data.PD, family = shifted_lognormal,
                  save_all_pars = TRUE, iter = 5000, cores = 1)

save(br.len.uF1,br.len.uF2,br.len.uF3,br.len.uF4,br.len.uF5,br.len.uF6,br.len.uF7,br.len.uFT,
     file = 'updrs_len_mods.R')

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


# END