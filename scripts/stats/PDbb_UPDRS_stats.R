## Statistics 2: within group correlation with clinical scores (MSD-UPDRS-III)
library(brms)
# library(BayesFactor)
Sys.setenv(PATH = paste("C:/Rtools/bin", Sys.getenv("PATH"), sep=";")) # Needed or there will be a pop-up everytime compiling C models.
Sys.setenv(BINPREF = "C:/Rtools/mingw_$(WIN)/bin/")
# options(mc.cores=parallel::detectCores)  # Try run with multicores. For some reason, doing this means bridge sampling does not work later!

# Define paths
wrkdir <- "Z://PD_motor//rest_ec//groupanalysis//"
setwd(wrkdir)
load(file='uData.Rdata')

################################################################################
## N events (TEST1)
u.neve.data.PD <- subset(u.neve.data, group=="ptns")

br.nev.uF1 <- brm(bf(F1 ~ nevent+(nevent|subs)), data = u.neve.data.PD, family = poisson,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.nev.uF2 <- brm(bf(F2 ~ nevent+(nevent|subs)), data = u.neve.data.PD, family = poisson,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.nev.uF3 <- brm(bf(F3 ~ nevent+(nevent|subs)), data = u.neve.data.PD, family = poisson,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.nev.uF4 <- brm(bf(F4 ~ nevent+(nevent|subs)), data = u.neve.data.PD, family = poisson,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.nev.uF5 <- brm(bf(F5 ~ nevent+(nevent|subs)), data = u.neve.data.PD, family = poisson,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.nev.uF6 <- brm(bf(F6 ~ nevent+(nevent|subs)), data = u.neve.data.PD, family = poisson,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.nev.uF7 <- brm(bf(F7 ~ nevent+(nevent|subs)), data = u.neve.data.PD, family = poisson,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
# br.nev.uFT <- brm(bf(Total ~ nevent+(1|subs:session)), data = u.neve.data.PD, family = poisson,
                  # save_all_pars = TRUE, iter = 5000, cores = 1)

# Null models
br.nev0.uF1 <- brm(bf(F1 ~ 1+(nevent|subs)), data = u.neve.data.PD, family = poisson,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.nev0.uF2 <- brm(bf(F2 ~ 1+(nevent|subs)), data = u.neve.data.PD, family = poisson,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.nev0.uF3 <- brm(bf(F3 ~ 1+(nevent|subs)), data = u.neve.data.PD, family = poisson,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.nev0.uF4 <- brm(bf(F4 ~ 1+(nevent|subs)), data = u.neve.data.PD, family = poisson,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.nev0.uF5 <- brm(bf(F5 ~ 1+(nevent|subs)), data = u.neve.data.PD, family = poisson,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.nev0.uF6 <- brm(bf(F6 ~ 1+(nevent|subs)), data = u.neve.data.PD, family = poisson,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
br.nev0.uF7 <- brm(bf(F7 ~ 1+(nevent|subs)), data = u.neve.data.PD, family = poisson,
                  save_all_pars = TRUE, iter = 5000, cores = 1)
# br.nev0.uFT <- brm(bf(Total ~ 1+(1|subs:session)), data = u.neve.data.PD, family = poisson,
#                   save_all_pars = TRUE, iter = 5000, cores = 1)

# Model comparison
n.bf.F1 <- bayes_factor(br.nev.uF1, br.nev0.uF1)
n.bf.F2 <- bayes_factor(br.nev.uF2, br.nev0.uF2)
n.bf.F3 <- bayes_factor(br.nev.uF3, br.nev0.uF3)
n.bf.F4 <- bayes_factor(br.nev.uF4, br.nev0.uF4)
n.bf.F5 <- bayes_factor(br.nev.uF5, br.nev0.uF5)
n.bf.F6 <- bayes_factor(br.nev.uF6, br.nev0.uF6)
n.bf.F7 <- bayes_factor(br.nev.uF7, br.nev0.uF7)
# n.bf.FT <- bayes_factor(br.nev.uFT, br.nev0.uFT)

# Save
setwd(wrkdir)
save(br.nev.uF1,br.nev.uF2,br.nev.uF3,br.nev.uF4,br.nev.uF5,br.nev.uF6,br.nev.uF7, #br.nev.uFT,
     br.nev0.uF1,br.nev0.uF2,br.nev0.uF3,br.nev0.uF4,br.nev0.uF5,br.nev0.uF6,br.nev0.uF7, #br.nev0.uFT,
     file = 'updrs_neve_mods2.R')

# hypothesis testing
load(file = 'updrs_neve_mods.R')

h1 <- hypothesis(br.nev.uF1, "nevent>0")
h2 <- hypothesis(br.nev.uF2, "nevent>0")
h3 <- hypothesis(br.nev.uF3, "nevent>0")
h4 <- hypothesis(br.nev.uF4, "nevent>0")
h5 <- hypothesis(br.nev.uF5, "nevent>0")
h6 <- hypothesis(br.nev.uF6, "nevent>0")
h7 <- hypothesis(br.nev.uF7, "nevent<0")
hT <- hypothesis(br.nev.uFT, "nevent>0")

P1 <- h1$hypothesis$Post.Prob*2
P2 <- h2$hypothesis$Post.Prob*2
P3 <- h3$hypothesis$Post.Prob*2
P4 <- h4$hypothesis$Post.Prob*2
P5 <- h5$hypothesis$Post.Prob*2
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

sam.nev4 <- posterior_samples(br.nev.uF4, "^b")
quantile(exp(sam.nev4$b_nevent*N), probs=c(0.025,0.5,0.975))-1

sam.nev5 <- posterior_samples(br.nev.uF5, "^b")
quantile(exp(sam.nev5$b_nevent*N), probs=c(0.025,0.5,0.975))-1

sam.nev6 <- posterior_samples(br.nev.uF6, "^b")
quantile(exp(sam.nev6$b_nevent*N), probs=c(0.025,0.5,0.975))-1

sam.nev7 <- posterior_samples(br.nev.uF7, "^b")
quantile(exp(sam.nev7$b_nevent*N), probs=c(0.025,0.5,0.975))-1

#END