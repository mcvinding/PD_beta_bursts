### Read data from Matlab. Save for analysis and plotting ###
library(R.matlab)

# Define paths
wrkdir <- "Z://PD_motor//rest_ec//groupanalysis//"
setwd(wrkdir)

# Read N event data
temp <- readMat("nevent.mat")
PDn1 <- temp$PDn1
PDn2 <- temp$PDn2
ctrln1 <- temp$ctrln1
ctrln2 <- temp$ctrln2
nevent <- c(temp$PDn1,temp$PDn2,temp$ctrln1,temp$ctrln2)
subs <- unlist(c(temp$PD.subs,temp$PD.subs,temp$ctrl.subs,temp$ctrl.subs))
group <- c(rep("ptns",2*length(temp$PD.subs)),rep("ctrl",2*length(temp$ctrl.subs)))
session <- c(rep(c(rep("1",length(temp$PD.subs)), rep("2",length(temp$PD.subs))),2))

neve.data <- data.frame(nevent=nevent,
                        subs=subs,
                        group=group,
                        session=session)

save(neve.data,PDn1,PDn2,ctrln1,ctrln2, file = 'neve.RData')

ptns.subs <- unlist(temp$PD.subs)
ctrl.subs <- unlist(temp$ctrl.subs)

# Read event duration data
temp <- readMat("lenevent.mat")
len <- c(temp$len1,temp$len2)
subs <- c(temp$sub1,temp$sub2)
group <- ifelse(subs %in% ptns.subs,'ptns','ctrl')
session <- c(rep("1",length(temp$len1)),rep("2",length(temp$len2)))

leneve.data <- data.frame(eve.len=len,
                        subs=subs,
                        group=group,
                        session=session)
leneve.data$eve.len.ms <- leneve.data$eve.len*1000
leneve.data$log.eve.len <- log(leneve.data$eve.len)

save(leneve.data, file = 'leneve.RData')

# Read time to event data
temp <- readMat("toevent.mat")
iti <- c(temp$alltoe1,temp$alltoe2)
subs <- c(temp$sub1,temp$sub2)
group <- ifelse(subs %in% ptns.subs,'ptns','ctrl')
session <- c(rep("1",length(temp$alltoe1)),rep("2",length(temp$alltoe2)))

itieve.data <- data.frame(eve.iti=iti,
                          subs=subs,
                          group=group,
                          session=session)
itieve.data <- subset(itieve.data, eve.iti != 0) # This should be fixed in the Matlab scripts! Update: it is, but the old data has not meed overwritten!
itieve.data$eve.iti.ms <- itieve.data$eve.iti*1000
itieve.data$log.eve.iti <- log(itieve.data$eve.iti)
save(itieve.data, file = 'itieve.RData')

# NB THIS IS WORK IN PROGRESS
data <- data.frame()
for (s in unique(itieve.data$subs)){
  temp <- subset(itieve.data, subs==s)
  temp$eve.iti.n1 <- c(temp$eve.iti[-1],NA)
  
  data <- rbind(data,temp)
  
}

# Read max peak in event data
temp <- readMat("pkmaxevent.mat")
maxpk <- c(temp$max1,temp$max2)
subs <- c(temp$sub1,temp$sub2)
group <- ifelse(subs %in% ptns.subs,'ptns','ctrl')
session <- c(rep("1",length(temp$max1)),rep("2",length(temp$max2)))

maxeve.data <- data.frame(eve.max=maxpk,
                          subs=subs,
                          group=group,
                          session=session)
maxeve.data$log.eve.max <- log(maxeve.data$eve.max)

save(maxeve.data, file = 'maxeve.RData')

# Clean-up and save entire workspace
rm(temp,iti,group,len,maxpk,nevent,session,subs)
save.image(file='workspace.Rdata')

#END