# Ideal observer analysis.
# * Beta burst rate
# * Inter-burst interval (mean, median, mode)
# * Burst duration (mean, median, mode)
# * Max burst amplitude (mean, median, mode)
# * Relative beta PSD
library(R.matlab)
# Import function
source("C:/Users/Mikkel/Documents/betabursts/scripts/functions/ideal_obs_fun.R")
source("C:/Users/Mikkel/Documents/betabursts/scripts/functions/misc_funs.R")

# Set work dir
wrkdir <- "C://Users//Mikkel//Documents//betabursts//groupanalysis"
setwd(wrkdir)

###################################################################
# Burst rate
load(file = 'neve.RData')

# Ptns vs. ctrl - session 1
id.obs(neve.data$nevent[neve.data$group=="ptns" & neve.data$session=="1"],
       neve.data$nevent[neve.data$group=="ctrl" & neve.data$session=="1"])

# Ptns vs. ctrl - session 2
id.obs(neve.data$nevent[neve.data$group=="ptns" & neve.data$session=="2"],
       neve.data$nevent[neve.data$group=="ctrl" & neve.data$session=="2"])

###################################################################
# Inter-burst interval
load(file = 'itieve.RData')

iti_mean <- aggregate(itieve.data$eve.iti, by=list(subs=itieve.data$subs, group=itieve.data$group, session=itieve.data$session), mean)
iti_medi <- aggregate(itieve.data$eve.iti, by=list(subs=itieve.data$subs, group=itieve.data$group, session=itieve.data$session), median)
iti_mode <- aggregate(itieve.data$eve.iti.ms, by=list(subs=itieve.data$subs, group=itieve.data$group, session=itieve.data$session), Mode)

# MEAN
# Ptns vs. ctrl - session 1
id.obs(iti_mean$x[iti_mean$group=="ptns" & iti_mean$session=="1"],
       iti_mean$x[iti_mean$group=="ctrl" & iti_mean$session=="1"])

# Ptns vs. ctrl - session 2
id.obs(iti_mean$x[iti_mean$group=="ptns" & iti_mean$session=="2"],
       iti_mean$x[iti_mean$group=="ctrl" & iti_mean$session=="2"])

#MEDIAN
# Ptns vs. ctrl - session 1
id.obs(iti_medi$x[iti_medi$group=="ptns" & iti_medi$session=="1"],
       iti_medi$x[iti_medi$group=="ctrl" & iti_medi$session=="1"])

# Ptns vs. ctrl - session 2
id.obs(iti_medi$x[iti_medi$group=="ptns" & iti_medi$session=="2"],
       iti_medi$x[iti_medi$group=="ctrl" & iti_medi$session=="2"])

# MODE
# Ptns vs. ctrl - session 1
id.obs(iti_mode$x[iti_mode$group=="ptns" & iti_mode$session=="1"],
       iti_mode$x[iti_mode$group=="ctrl" & iti_mode$session=="1"])

# Ptns vs. ctrl - session 2
id.obs(iti_mode$x[iti_mode$group=="ptns" & iti_mode$session=="2"],
       iti_mode$x[iti_mode$group=="ctrl" & iti_mode$session=="2"])


###################################################################
# Event duration
load(file = 'leneve.RData')

len_mean <- aggregate(leneve.data$eve.len, by=list(subs=leneve.data$subs, group=leneve.data$group, session=leneve.data$session), mean)
len_medi <- aggregate(leneve.data$eve.len, by=list(subs=leneve.data$subs, group=leneve.data$group, session=leneve.data$session), median)
len_mode <- aggregate(leneve.data$eve.len.ms, by=list(subs=leneve.data$subs, group=leneve.data$group, session=leneve.data$session), Mode)

# MEAN
# Ptns vs. ctrl - session 1
id.obs(len_mean$x[len_mean$group=="ptns" & len_mean$session=="1"],
       len_mean$x[len_mean$group=="ctrl" & len_mean$session=="1"])

# Ptns vs. ctrl - session 2
id.obs(len_mean$x[len_mean$group=="ptns" & len_mean$session=="2"],
       len_mean$x[len_mean$group=="ctrl" & len_mean$session=="2"])

# MEDIAN
# Ptns vs. ctrl - session 1
id.obs(len_medi$x[len_medi$group=="ptns" & len_medi$session=="1"],
       len_medi$x[len_medi$group=="ctrl" & len_medi$session=="1"])

# Ptns vs. ctrl - session 2
id.obs(len_medi$x[len_medi$group=="ptns" & len_medi$session=="2"],
       len_medi$x[len_medi$group=="ctrl" & len_medi$session=="2"])

# MODE
# Ptns vs. ctrl - session 1
id.obs(len_mode$x[len_mode$group=="ptns" & len_mode$session=="1"],
       len_mode$x[len_mode$group=="ctrl" & len_mode$session=="1"])

# Ptns vs. ctrl - session 2
id.obs(len_mode$x[len_mode$group=="ptns" & len_mode$session=="2"],
       len_mode$x[len_mode$group=="ctrl" & len_mode$session=="2"])


###################################################################
# Max amplitude
load(file = 'maxeve.RData')

amp_mean <- aggregate(maxeve.data$eve.max, by=list(subs=maxeve.data$subs, group=maxeve.data$group, session=maxeve.data$session), mean)
amp_medi <- aggregate(maxeve.data$eve.max, by=list(subs=maxeve.data$subs, group=maxeve.data$group, session=maxeve.data$session), median)
amp_mode <- aggregate(maxeve.data$eve.max, by=list(subs=maxeve.data$subs, group=maxeve.data$group, session=maxeve.data$session), Mode, digits=2)

# MEAN
# Ptns vs. ctrl - session 1
id.obs(amp_mean$x[amp_mean$group=="ptns" & amp_mean$session=="1"],
       amp_mean$x[amp_mean$group=="ctrl" & amp_mean$session=="1"])

# Ptns vs. ctrl - session 2
id.obs(amp_mean$x[amp_mean$group=="ptns" & amp_mean$session=="2"],
       amp_mean$x[amp_mean$group=="ctrl" & amp_mean$session=="2"])

# MEDIAN
# Ptns vs. ctrl - session 1
id.obs(amp_medi$x[amp_medi$group=="ptns" & amp_medi$session=="1"],
       amp_medi$x[amp_medi$group=="ctrl" & amp_medi$session=="1"])

# Ptns vs. ctrl - session 2
id.obs(amp_medi$x[amp_medi$group=="ptns" & amp_medi$session=="2"],
       amp_medi$x[amp_medi$group=="ctrl" & amp_medi$session=="2"])

# MODE
# Ptns vs. ctrl - session 1
id.obs(amp_mode$x[amp_mode$group=="ptns" & amp_mode$session=="1"],
       amp_mode$x[amp_mode$group=="ctrl" & amp_mode$session=="1"])

# Ptns vs. ctrl - session 2
id.obs(amp_mode$x[amp_mode$group=="ptns" & amp_mode$session=="2"],
       amp_mode$x[amp_mode$group=="ctrl" & amp_mode$session=="2"])


###################################################################
# Relative beta PSD
load(file = 'reldat.RData')

# Ptns vs. ctrl - session 1
id.obs(rel.dat$relpow[rel.dat$group=="1" & rel.dat$session=="1"],
       rel.dat$relpow[rel.dat$group=="2" & rel.dat$session=="1"])

# Ptns vs. ctrl - session 2
id.obs(rel.dat$relpow[rel.dat$group=="1" & rel.dat$session=="2"],
       rel.dat$relpow[rel.dat$group=="2" & rel.dat$session=="2"])

###################################################################
# Beta power substracted 1/f spectrum





###################################################################

# # Plot distribution
# library(ggplot2)
# d1 <- data.frame(val=r0)
# d2 <- data.frame(val=r1)
# d1$var <- 'r1'
# d2$var <- 'r2'
# dat <- rbind(d1, d2)
# ggplot(dat, aes(val, fill = var)) + geom_histogram(bins=20, alpha=.8, position = 'identity')


#END