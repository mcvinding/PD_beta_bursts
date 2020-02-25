# Prepare MSD-UPDRS-III data.
library(xlsx)
library(reshape2)

# Define paths
wrkdir <- "C://Users//Mikkel//Documents//betabursts//groupanalysis"
setwd(wrkdir)
load(file = 'maxeve.RData')
load(file = 'itieve.RData')
load(file = 'leneve.RData')
load(file = 'neve.RData')

#Load data
raw.UPDRS <- read.xlsx('UPDRS_raw.xlsx',1,header=T)
raw.UPDRS$n.id = as.factor(raw.UPDRS$id)
raw.UPDRS$session <- as.factor(raw.UPDRS$session)

load(file='alldata.RData')

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
updrs.data$F45 <- updrs.data$F4 + updrs.data$F5     # Combine L+R bradykinesia
updrs.data$id.x <- paste('0',updrs.data$id.x, sep="")

u.neve.data <- merge(neve.data, updrs.data, by.x=c('subs','session','group'), by.y=c('id.x','session','group'))
u.len.data <- merge(leneve.data, updrs.data, by.x=c('subs','session','group'), by.y=c('id.x','session','group')) 
u.iti.data <- merge(itieve.data, updrs.data, by.x=c('subs','session','group'), by.y=c('id.x','session','group')) 
u.max.data <- merge(maxeve.data, updrs.data, by.x=c('subs','session','group'), by.y=c('id.x','session','group')) 

#save
setwd(wrkdir)
save(u.neve.data,u.len.data,u.iti.data,u.max.data,
     file='uData.Rdata')

#END