# Prepare MSD-UPDRS-III data.
library(xlsx)
library(reshape2)

# Define paths
wrkdir <- "Z://PD_motor//rest_ec//groupanalysis//"
# wrkdir <- "C://Users//Mikkel//Documents//PD-proj_betaEvent//data"
setwd(wrkdir)
load(file = 'workspace.Rdata')
load(file = 'maxeve.RData')
load(file = 'itieve.RData')
load(file = 'leneve.RData')
load(file = 'neve.RData')

#Load data
raw.UPDRS <- read.xlsx('Z://PD_motor//subj_data//UPDRS_raw.xlsx',1,header=T)
raw.UPDRS$n.id = as.factor(raw.UPDRS$id)
raw.UPDRS$session <- as.factor(raw.UPDRS$session)

load(file='Z://PD_motor//subj_data//alldata.RData')

# Prepare data
updrs.long <- data.frame(id = rep(alldata$MEG_ID,2),
                         UPDRS.old = c(alldata$UPDRS_off,alldata$UPDRS_on),
                         group = rep(alldata$Sub_type,2),
                         session = as.factor(rep(c(1,2),each=length(alldata$MEG_ID))),
                         n.id = rep(alldata$ID,2),
                         hand = rep(alldata$PAM_hand,2),
                         HY_stage = rep(alldata$HY_stage,2))
levels(updrs.long$group) <- c("ctrl","ptns")
updrs.data <- merge(updrs.long,raw.UPDRS,by.x=c('n.id','session'), by.y=c('n.id','session'))
updrs.data$id.x <- paste('0',updrs.data$id.x, sep="")

u.neve.data <- merge(neve.data, updrs.data, by.x=c('subs','session','group'), by.y=c('id.x','session','group'))
u.len.data <- merge(leneve.data, updrs.data, by.x=c('subs','session','group'), by.y=c('id.x','session','group')) 
u.iti.data <- merge(itieve.data, updrs.data, by.x=c('subs','session','group'), by.y=c('id.x','session','group')) 
u.max.data <- merge(maxeve.data, updrs.data, by.x=c('subs','session','group'), by.y=c('id.x','session','group')) 


setwd(wrkdir)
save(u.neve.data,u.len.data,u.iti.data,u.max.data,
     file='uData.Rdata')


## Make difference data
updrs.data.PD <- subset(updrs.data, group=="ptns")

subs <- unique(updrs.data.PD$id.x)
diff.dat <- data.frame(subs=subs, F1=0, F2=0, F3=0, F4=0, F5=0, F6=0, F7=0, Ft = 0, dN = 0)

for (i in 1:length(updrs.data.PD$id.x)) {
  s <- updrs.data.PD$id.x[i]
  tmp <- subset(updrs.data, id.x==s)
  tmp.N <- subset(neve.data, subs==s)
  diff.dat$F1[diff.dat$subs==s] <- tmp$F1[1]-tmp$F1[2]
  diff.dat$F2[diff.dat$subs==s] <- tmp$F2[1]-tmp$F2[2]
  diff.dat$F3[diff.dat$subs==s] <- tmp$F3[1]-tmp$F3[2]
  diff.dat$F4[diff.dat$subs==s] <- tmp$F4[1]-tmp$F4[2]
  diff.dat$F5[diff.dat$subs==s] <- tmp$F5[1]-tmp$F5[2]
  diff.dat$F6[diff.dat$subs==s] <- tmp$F6[1]-tmp$F6[2]
  diff.dat$F7[diff.dat$subs==s] <- tmp$F7[1]-tmp$F7[2]
  diff.dat$Ft[diff.dat$subs==s] <- tmp$Total[1]-tmp$Total[2]
  diff.dat$dN[diff.dat$subs==s] <- tmp.N$nevent[2]-tmp.N$nevent[1]
}

#END