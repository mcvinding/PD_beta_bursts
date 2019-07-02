### Read data from Matlab. Save for analysis and plotting ###
library(R.matlab)

# Define paths
wrkdir <- "Z://PD_motor//rest_ec//groupanalysis//"
setwd(wrkdir)

# Read N event data
temp <- readMat("nevent_extended.mat")
subs <- unlist(temp$subs.long)
group <- unlist(temp$grup.long)
session <- as.factor(drop(unlist(temp$sesi.long)))
nevent <- drop(unlist(temp$neve.long))
steps <- drop(unlist(temp$step.long))

neve.data <- data.frame(nevent=nevent,
                        steps=steps,
                        subs=subs,
                        group=group,
                        session=session)

save(neve.data, file = 'neve_ext.RData')

#END