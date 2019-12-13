## Statistics 1c: between group-session analysis of "fooof" analysis 
# * beta power - 1/f characteristics
# * 1/f Slope / Intercept
library(plyr)
library(BayesFactor)

# Define paths
#wrkdir <- "Z://PD_motor//rest_ec//groupanalysis//"
wrkdir <- "C://Users//Mikkel//Documents//betabursts//groupanalysis"
setwd(wrkdir)

########################################################################################
# 1/f regression parameters
f.dat <- read.csv("df_aper.csv", header=TRUE, sep=";")
f.dat$group <- as.factor(f.dat$group)
f.dat$group <- revalue(f.dat$group, c("1"="ptns", "2"="ctrl"))
f.dat$session <- as.factor(f.dat$session)
f.dat$subj <- as.factor(rep(1:(length(f.dat$intercept)/2),2))

# Summaries: intercept
aggregate(f.dat$intercept, by=list(group=f.dat$group, session=f.dat$session), mean)
aggregate(f.dat$intercept, by=list(group=f.dat$group, session=f.dat$session), median)
aggregate(f.dat$intercept, by=list(group=f.dat$group, session=f.dat$session), sd)

# BF stats: intercept
ttestBF(f.dat$intercept[f.dat$group=="ptns" & f.dat$session=="1"], f.dat$intercept[f.dat$group=="ctrl" & f.dat$session=="1"])
ttestBF(f.dat$intercept[f.dat$group=="ptns" & f.dat$session=="2"], f.dat$intercept[f.dat$group=="ctrl" & f.dat$session=="2"])
ttestBF(f.dat$intercept[f.dat$group=="ptns" & f.dat$session=="1"], f.dat$intercept[f.dat$group=="ptns" & f.dat$session=="2"], paired=TRUE)
ttestBF(f.dat$intercept[f.dat$group=="ctrl" & f.dat$session=="1"], f.dat$intercept[f.dat$group=="ctrl" & f.dat$session=="2"], paired=TRUE)

# NHST stats: intercept
t.test(f.dat$intercept[f.dat$group=="ptns" & f.dat$session=="1"], f.dat$intercept[f.dat$group=="ctrl" & f.dat$session=="1"])
t.test(f.dat$intercept[f.dat$group=="ptns" & f.dat$session=="2"], f.dat$intercept[f.dat$group=="ctrl" & f.dat$session=="2"])
t.test(f.dat$intercept[f.dat$group=="ptns" & f.dat$session=="1"], f.dat$intercept[f.dat$group=="ptns" & f.dat$session=="2"], paired=TRUE)
t.test(f.dat$intercept[f.dat$group=="ctrl" & f.dat$session=="1"], f.dat$intercept[f.dat$group=="ctrl" & f.dat$session=="2"], paired=TRUE)

# Summaries: slope
aggregate(f.dat$slope, by=list(group=f.dat$group, session=f.dat$session), mean)
aggregate(f.dat$slope, by=list(group=f.dat$group, session=f.dat$session), median)
aggregate(f.dat$slope, by=list(group=f.dat$group, session=f.dat$session), sd)

# BF stats: slope
ttestBF(f.dat$slope[f.dat$group=="ptns" & f.dat$session=="1"], f.dat$slope[f.dat$group=="ctrl" & f.dat$session=="1"])
ttestBF(f.dat$slope[f.dat$group=="ptns" & f.dat$session=="2"], f.dat$slope[f.dat$group=="ctrl" & f.dat$session=="2"])
ttestBF(f.dat$slope[f.dat$group=="ptns" & f.dat$session=="1"], f.dat$slope[f.dat$group=="ptns" & f.dat$session=="2"], paired=TRUE)
ttestBF(f.dat$slope[f.dat$group=="ctrl" & f.dat$session=="1"], f.dat$slope[f.dat$group=="ctrl" & f.dat$session=="2"], paired=TRUE)

# NHST stats: slpow
t.test(f.dat$slope[f.dat$group=="ptns" & f.dat$session=="1"], f.dat$slope[f.dat$group=="ctrl" & f.dat$session=="1"])
t.test(f.dat$slope[f.dat$group=="ptns" & f.dat$session=="2"], f.dat$slope[f.dat$group=="ctrl" & f.dat$session=="2"])
t.test(f.dat$slope[f.dat$group=="ptns" & f.dat$session=="1"], f.dat$slope[f.dat$group=="ptns" & f.dat$session=="2"], paired=TRUE)
t.test(f.dat$slope[f.dat$group=="ctrl" & f.dat$session=="1"], f.dat$slope[f.dat$group=="ctrl" & f.dat$session=="2"], paired=TRUE)

########################################################################################
# BETA
# Load data
b.dat <- read.csv("df_beta.csv", header=TRUE, sep=";")
b.dat$group <- as.factor(b.dat$group)
b.dat$group <- revalue(b.dat$group, c("1"="ptns", "2"="ctrl"))
b.dat$session <- as.factor(b.dat$session)
b.dat$subj <- as.factor(rep(1:(length(b.dat$peak_freq)/2),2))

# Summaries: peak freq.
aggregate(b.dat$peak_freq, by=list(group=b.dat$group, session=b.dat$session), mean)
aggregate(b.dat$peak_freq, by=list(group=b.dat$group, session=b.dat$session), median)
aggregate(b.dat$peak_freq, by=list(group=b.dat$group, session=b.dat$session), sd)

# BF stats: peak freq.
ttestBF(b.dat$peak_freq[b.dat$group=="ptns" & b.dat$session=="1"], b.dat$peak_freq[b.dat$group=="ctrl" & b.dat$session=="1"])
ttestBF(b.dat$peak_freq[b.dat$group=="ptns" & b.dat$session=="2"], b.dat$peak_freq[b.dat$group=="ctrl" & b.dat$session=="2"])
ttestBF(b.dat$peak_freq[b.dat$group=="ptns" & b.dat$session=="1"], b.dat$peak_freq[b.dat$group=="ptns" & b.dat$session=="2"], paired=TRUE)
ttestBF(b.dat$peak_freq[b.dat$group=="ctrl" & b.dat$session=="1"], b.dat$peak_freq[b.dat$group=="ctrl" & b.dat$session=="2"], paired=TRUE)

# NHST stats: peak freq.
t.test(b.dat$peak_freq[b.dat$group=="ptns" & b.dat$session=="1"], b.dat$peak_freq[b.dat$group=="ctrl" & b.dat$session=="1"])
t.test(b.dat$peak_freq[b.dat$group=="ptns" & b.dat$session=="2"], b.dat$peak_freq[b.dat$group=="ctrl" & b.dat$session=="2"])
t.test(b.dat$peak_freq[b.dat$group=="ptns" & b.dat$session=="1"], b.dat$peak_freq[b.dat$group=="ptns" & b.dat$session=="2"], paired=TRUE)
t.test(b.dat$peak_freq[b.dat$group=="ctrl" & b.dat$session=="1"], b.dat$peak_freq[b.dat$group=="ctrl" & b.dat$session=="2"], paired=TRUE)

# Summaries: peak pow
aggregate(b.dat$peak_pow, by=list(group=b.dat$group, session=b.dat$session), mean)
aggregate(b.dat$peak_pow, by=list(group=b.dat$group, session=b.dat$session), median)
aggregate(b.dat$peak_pow, by=list(group=b.dat$group, session=b.dat$session), sd)

# BF stats: peak pow
ttestBF((b.dat$peak_pow[b.dat$group=="ptns" & b.dat$session=="1"]), (b.dat$peak_pow[b.dat$group=="ctrl" & b.dat$session=="1"]))
ttestBF((b.dat$peak_pow[b.dat$group=="ptns" & b.dat$session=="2"]), (b.dat$peak_pow[b.dat$group=="ctrl" & b.dat$session=="2"]))
ttestBF((b.dat$peak_pow[b.dat$group=="ptns" & b.dat$session=="1"]), (b.dat$peak_pow[b.dat$group=="ptns" & b.dat$session=="2"]), paired=TRUE)
ttestBF((b.dat$peak_pow[b.dat$group=="ctrl" & b.dat$session=="1"]), (b.dat$peak_pow[b.dat$group=="ctrl" & b.dat$session=="2"]), paired=TRUE)

# NHST stats: peak pow
t.test((b.dat$peak_pow[b.dat$group=="ptns" & b.dat$session=="1"]), (b.dat$peak_pow[b.dat$group=="ctrl" & b.dat$session=="1"]))
t.test((b.dat$peak_pow[b.dat$group=="ptns" & b.dat$session=="2"]), (b.dat$peak_pow[b.dat$group=="ctrl" & b.dat$session=="2"]))
t.test((b.dat$peak_pow[b.dat$group=="ptns" & b.dat$session=="1"]), (b.dat$peak_pow[b.dat$group=="ptns" & b.dat$session=="2"]), paired=TRUE)
t.test((b.dat$peak_pow[b.dat$group=="ctrl" & b.dat$session=="1"]), (b.dat$peak_pow[b.dat$group=="ctrl" & b.dat$session=="2"]), paired=TRUE)

hist((b.dat$peak_pow[b.dat$group=="ptns" & b.dat$session=="1"]), 10)
hist((b.dat$peak_pow[b.dat$group=="ptns" & b.dat$session=="2"]), 10)
hist((b.dat$peak_pow[b.dat$group=="ctrl" & b.dat$session=="1"]), 10)
hist((b.dat$peak_pow[b.dat$group=="ctrl" & b.dat$session=="2"]), 10)

########################################################################################
# Alpha
# Load data
a.dat <- read.csv("df_alpha.csv", header=TRUE, sep=";")
a.dat$group <- as.factor(a.dat$group)
a.dat$group <- revalue(a.dat$group, c("1"="ptns", "2"="ctrl"))
a.dat$session <- as.factor(a.dat$session)
a.dat$subj <- as.factor(rep(1:(length(a.dat$peak_freq)/2),2))
a.dat$is.pk <- !is.na(a.dat$peak_freq)

# How many subjects had alpha peak
aggregate(a.dat$is.pk, by=list(group=a.dat$group, session=a.dat$session), sum)

# Remove missing cases
a.dat0 <- subset(a.dat, a.dat$is.pk==1)
for (i in 1:length(unique(a.dat$subj))){
  tmp <- subset(a.dat, a.dat$subj==unique(a.dat$subj)[i])
  if (any(!tmp$is.pk)){
    a.dat$rm[a.dat$subj==unique(a.dat$subj)[i]] <- 1
  } else {
    a.dat$rm[a.dat$subj==unique(a.dat$subj)[i]] <- 0
  }
}
a.dat1 <- subset(a.dat, a.dat$rm==0)

# Summaries: peak freq
aggregate(a.dat$peak_freq, by=list(group=a.dat$group, session=a.dat$session), mean, na.rm=TRUE)
aggregate(a.dat$peak_freq, by=list(group=a.dat$group, session=a.dat$session), median, na.rm=TRUE)
aggregate(a.dat$peak_freq, by=list(group=a.dat$group, session=a.dat$session), sd, na.rm=TRUE)

# BF stats: peak freq.
ttestBF(a.dat0$peak_freq[a.dat0$group=="ptns" & a.dat0$session=="1"], a.dat0$peak_freq[a.dat0$group=="ctrl" & a.dat0$session=="1"])
ttestBF(a.dat0$peak_freq[a.dat0$group=="ptns" & a.dat0$session=="2"], a.dat0$peak_freq[a.dat0$group=="ctrl" & a.dat0$session=="2"])
ttestBF(a.dat1$peak_freq[a.dat1$group=="ptns" & a.dat1$session=="1"], a.dat1$peak_freq[a.dat1$group=="ptns" & a.dat1$session=="2"], paired=TRUE)
ttestBF(a.dat1$peak_freq[a.dat1$group=="ctrl" & a.dat1$session=="1"], a.dat1$peak_freq[a.dat1$group=="ctrl" & a.dat1$session=="2"], paired=TRUE)

# NHST stats: peak freq.
t.test(a.dat0$peak_freq[a.dat0$group=="ptns" & a.dat0$session=="1"], a.dat0$peak_freq[a.dat0$group=="ctrl" & a.dat0$session=="1"])
t.test(a.dat0$peak_freq[a.dat0$group=="ptns" & a.dat0$session=="2"], a.dat0$peak_freq[a.dat0$group=="ctrl" & a.dat0$session=="2"])
t.test(a.dat1$peak_freq[a.dat1$group=="ptns" & a.dat1$session=="1"], a.dat1$peak_freq[a.dat1$group=="ptns" & a.dat1$session=="2"], paired=TRUE)
t.test(a.dat1$peak_freq[a.dat1$group=="ctrl" & a.dat1$session=="1"], a.dat1$peak_freq[a.dat1$group=="ctrl" & a.dat1$session=="2"], paired=TRUE)

# Summaries: peak pow
aggregate(a.dat$peak_pow, by=list(group=a.dat$group, session=a.dat$session), mean, na.rm=TRUE)
aggregate(a.dat$peak_pow, by=list(group=a.dat$group, session=a.dat$session), median, na.rm=TRUE)
aggregate(a.dat$peak_pow, by=list(group=a.dat$group, session=a.dat$session), sd, na.rm=TRUE)

# BF stats: peak pow
ttestBF(a.dat0$peak_pow[a.dat0$group=="ptns" & a.dat0$session=="1"], a.dat0$peak_pow[a.dat0$group=="ctrl" & a.dat0$session=="1"])
ttestBF(a.dat0$peak_pow[a.dat0$group=="ptns" & a.dat0$session=="2"], a.dat0$peak_pow[a.dat0$group=="ctrl" & a.dat0$session=="2"])
ttestBF(a.dat1$peak_pow[a.dat1$group=="ptns" & a.dat1$session=="1"], a.dat1$peak_pow[a.dat1$group=="ptns" & a.dat1$session=="2"], paired=TRUE)
ttestBF(a.dat1$peak_pow[a.dat1$group=="ctrl" & a.dat1$session=="1"], a.dat1$peak_pow[a.dat1$group=="ctrl" & a.dat1$session=="2"], paired=TRUE)

# NHST stats: peak pow
t.test(a.dat0$peak_pow[a.dat0$group=="ptns" & a.dat0$session=="1"], a.dat0$peak_pow[a.dat0$group=="ctrl" & a.dat0$session=="1"])
t.test(a.dat0$peak_pow[a.dat0$group=="ptns" & a.dat0$session=="2"], a.dat0$peak_pow[a.dat0$group=="ctrl" & a.dat0$session=="2"])
t.test(a.dat1$peak_pow[a.dat1$group=="ptns" & a.dat1$session=="1"], a.dat1$peak_pow[a.dat1$group=="ptns" & a.dat1$session=="2"], paired=TRUE)
t.test(a.dat1$peak_pow[a.dat1$group=="ctrl" & a.dat1$session=="1"], a.dat1$peak_pow[a.dat1$group=="ctrl" & a.dat1$session=="2"], paired=TRUE)

########################################################################################
# Theta
# Load data
t.dat <- read.csv("df_theta.csv", header=TRUE, sep=";")
t.dat$group <- as.factor(t.dat$group)
t.dat$group <- revalue(t.dat$group, c("1"="ptns", "2"="ctrl"))
t.dat$session <- as.factor(t.dat$session)
t.dat$subj <- as.factor(rep(1:(length(t.dat$peak_freq)/2),2))
t.dat$is.pk <- !is.na(t.dat$peak_freq)

# How many subjects had alpha peak
aggregate(t.dat$is.pk, by=list(group=t.dat$group, session=t.dat$session), sum)

# Remove missing cases
t.dat0 <- subset(t.dat, t.dat$is.pk==1)
for (i in 1:length(unique(t.dat$subj))){
  tmp <- subset(t.dat, t.dat$subj==unique(t.dat$subj)[i])
  if (any(!tmp$is.pk)){
    t.dat$rm[t.dat$subj==unique(t.dat$subj)[i]] <- 1
  } else {
    t.dat$rm[t.dat$subj==unique(t.dat$subj)[i]] <- 0
  }
}
t.dat1 <- subset(t.dat, t.dat$rm==0)

# Summaries: peak freq
aggregate(t.dat$peak_freq, by=list(group=t.dat$group, session=t.dat$session), mean, na.rm=TRUE)
aggregate(t.dat$peak_freq, by=list(group=t.dat$group, session=t.dat$session), median, na.rm=TRUE)
aggregate(t.dat$peak_freq, by=list(group=t.dat$group, session=t.dat$session), sd, na.rm=TRUE)

# BF stats: peak freq.
ttestBF(t.dat0$peak_freq[t.dat0$group=="ptns" & t.dat0$session=="1"], t.dat0$peak_freq[t.dat0$group=="ctrl" & t.dat0$session=="1"])
ttestBF(t.dat0$peak_freq[t.dat0$group=="ptns" & t.dat0$session=="2"], t.dat0$peak_freq[t.dat0$group=="ctrl" & t.dat0$session=="2"])
ttestBF(t.dat1$peak_freq[t.dat1$group=="ptns" & t.dat1$session=="1"], t.dat1$peak_freq[t.dat1$group=="ptns" & t.dat1$session=="2"], paired=TRUE)
ttestBF(t.dat1$peak_freq[t.dat1$group=="ctrl" & t.dat1$session=="1"], t.dat1$peak_freq[t.dat1$group=="ctrl" & t.dat1$session=="2"], paired=TRUE)

# NHST stats: peak freq.
t.test(t.dat0$peak_freq[t.dat0$group=="ptns" & t.dat0$session=="1"], t.dat0$peak_freq[t.dat0$group=="ctrl" & t.dat0$session=="1"])
t.test(t.dat0$peak_freq[t.dat0$group=="ptns" & t.dat0$session=="2"], t.dat0$peak_freq[t.dat0$group=="ctrl" & t.dat0$session=="2"])
t.test(t.dat1$peak_freq[t.dat1$group=="ptns" & t.dat1$session=="1"], t.dat1$peak_freq[t.dat1$group=="ptns" & t.dat1$session=="2"], paired=TRUE)
t.test(t.dat1$peak_freq[t.dat1$group=="ctrl" & t.dat1$session=="1"], t.dat1$peak_freq[t.dat1$group=="ctrl" & t.dat1$session=="2"], paired=TRUE)

# Summaries: peak pow
aggregate(t.dat$peak_pow, by=list(group=t.dat$group, session=t.dat$session), mean, na.rm=TRUE)
aggregate(t.dat$peak_pow, by=list(group=t.dat$group, session=t.dat$session), median, na.rm=TRUE)
aggregate(t.dat$peak_pow, by=list(group=t.dat$group, session=t.dat$session), sd, na.rm=TRUE)

# BF stats: peak pow
ttestBF(t.dat0$peak_pow[t.dat0$group=="ptns" & t.dat0$session=="1"], t.dat0$peak_pow[t.dat0$group=="ctrl" & t.dat0$session=="1"])
ttestBF(t.dat0$peak_pow[t.dat0$group=="ptns" & t.dat0$session=="2"], t.dat0$peak_pow[t.dat0$group=="ctrl" & t.dat0$session=="2"])
ttestBF(t.dat1$peak_pow[t.dat1$group=="ptns" & t.dat1$session=="1"], t.dat1$peak_pow[t.dat1$group=="ptns" & t.dat1$session=="2"], paired=TRUE)
ttestBF(t.dat1$peak_pow[t.dat1$group=="ctrl" & t.dat1$session=="1"], t.dat1$peak_pow[t.dat1$group=="ctrl" & t.dat1$session=="2"], paired=TRUE)

# NHST stats: peak pow
t.test(t.dat0$peak_pow[t.dat0$group=="ptns" & t.dat0$session=="1"], t.dat0$peak_pow[t.dat0$group=="ctrl" & t.dat0$session=="1"])
t.test(t.dat0$peak_pow[t.dat0$group=="ptns" & t.dat0$session=="2"], t.dat0$peak_pow[t.dat0$group=="ctrl" & t.dat0$session=="2"])
t.test(t.dat1$peak_pow[t.dat1$group=="ptns" & t.dat1$session=="1"], t.dat1$peak_pow[t.dat1$group=="ptns" & t.dat1$session=="2"], paired=TRUE)
t.test(t.dat1$peak_pow[t.dat1$group=="ctrl" & t.dat1$session=="1"], t.dat1$peak_pow[t.dat1$group=="ctrl" & t.dat1$session=="2"], paired=TRUE)

###################################################################################
# Save data for later
save("f.dat", "t.dat", "b.dat", "a.dat",
     file="fooof_dat.R")

#END