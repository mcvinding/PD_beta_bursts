## Statistics 2: within group mixed-models with clinical scores (MSD-UPDRS-III)
library(brms)
Sys.setenv(PATH = paste("C:/Rtools/bin", Sys.getenv("PATH"), sep=";")) # Needed or there will be a pop-up everytime compiling C models.
Sys.setenv(BINPREF = "C:/Rtools/mingw_$(WIN)/bin/")
# options(mc.cores=parallel::detectCores)  # Try run with multicores. For some reason, doing this means bridge sampling does not work later!

# Define paths
wrkdir <- "C://Users//Mikkel//Documents//betabursts//groupanalysis"
setwd(wrkdir)
load(file='uData.Rdata')

################################################################################
## N events
u.neve.data.PD <- subset(u.neve.data, group=="ptns")

br.nev.uF1 <- brm(F1 ~ nevent.min+(1|subs:session), data = u.neve.data.PD, family = poisson,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.nev.uF2 <- brm(F2 ~ nevent.min+(1|subs:session), data = u.neve.data.PD, family = poisson,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.nev.uF3 <- brm(F3 ~ nevent.min+(1|subs:session), data = u.neve.data.PD, family = poisson,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.nev.uF45 <- brm(F45 ~ nevent.min+(1|subs:session), data = u.neve.data.PD, family = poisson,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.nev.uF6 <- brm(F6 ~ nevent.min+(1|subs:session), data = u.neve.data.PD, family = poisson,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.nev.uF7 <- brm(F7 ~ nevent.min+(1|subs:session), data = u.neve.data.PD, family = poisson,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.nev.uFT <- brm(Total ~ nevent.min+(1|subs:session), data = u.neve.data.PD, family = poisson,
                  save_all_pars = TRUE, iter = 5000, cores = 1)

# Save
setwd(wrkdir)
save(br.nev.uF1,br.nev.uF2,br.nev.uF3,br.nev.uF45,br.nev.uF6,br.nev.uF7, br.nev.uFT,
     file = 'updrs_neve_mods.R')

# hypothesis testing
load(file = 'updrs_neve_mods.R')

h1 <- hypothesis(br.nev.uF1, "nevent.min>0")
h2 <- hypothesis(br.nev.uF2, "nevent.min<0")
h3 <- hypothesis(br.nev.uF3, "nevent.min>0")
h45 <- hypothesis(br.nev.uF45, "nevent.min>0")
h6 <- hypothesis(br.nev.uF6, "nevent.min>0")
h7 <- hypothesis(br.nev.uF7, "nevent.min<0")
hT <- hypothesis(br.nev.uFT, "nevent.min>0")

P1 <- h1$hypothesis$Post.Prob*2
P2 <- h2$hypothesis$Post.Prob*2
P3 <- h3$hypothesis$Post.Prob*2
P45 <- h45$hypothesis$Post.Prob*2
P6 <- h6$hypothesis$Post.Prob*2
P7 <- h7$hypothesis$Post.Prob*2
PT <- hT$hypothesis$Post.Prob*2

# Pediction (how much % increase/decrease per N extra bursts)
N <- 10

sam.nev1 <- posterior_samples(br.nev.uF1, "^b")
quantile(exp(sam.nev1$b_nevent*N), probs=c(0.025,0.5,0.975))-1

sam.nev2 <- posterior_samples(br.nev.uF2, "^b")
quantile(exp(sam.nev2$b_nevent*N), probs=c(0.025,0.5,0.975))-1

sam.nev3 <- posterior_samples(br.nev.uF3, "^b")
quantile(exp(sam.nev3$b_nevent*N), probs=c(0.025,0.5,0.975))-1

sam.nev45 <- posterior_samples(br.nev.uF45, "^b")
quantile(exp(sam.nev45$b_nevent*N), probs=c(0.025,0.5,0.975))-1

sam.nev6 <- posterior_samples(br.nev.uF6, "^b")
quantile(exp(sam.nev6$b_nevent*N), probs=c(0.025,0.5,0.975))-1

sam.nev7 <- posterior_samples(br.nev.uF7, "^b")
quantile(exp(sam.nev7$b_nevent*N), probs=c(0.025,0.5,0.975))-1

#END