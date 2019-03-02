# Prepare MSD-UPDRS-III data.
library(xlsx)
library(reshape2)

# Define paths
wrkdir <- "Z://PD_motor//rest_ec//groupanalysis//"
wrkdir <- "C://Users//Mikkel//Documents//PD-proj_betaEvent//data"
setwd(wrkdir)
load(file = 'workspace.Rdata')
load(file = '.Rdata')

#Load data
raw.UPDRS <- read.xlsx('Z://PD_motor//subj_data//UPDRS_raw.xlsx',1,header=T)
raw.UPDRS$n.id = as.factor(raw.UPDRS$id)
raw.UPDRS$session <- as.factor(raw.UPDRS$session)
# levels(raw.UPDRS$session) <- c('off','on')

# raw.UPDRS$F1 <- as.numeric(as.character(raw.UPDRS$F1))
# raw.UPDRS$F2 <- as.numeric(as.character(raw.UPDRS$F2))
# raw.UPDRS$F3 <- as.numeric(as.character(raw.UPDRS$F3))
# raw.UPDRS$F4 <- as.numeric(as.character(raw.UPDRS$F4))
# raw.UPDRS$F5 <- as.numeric(as.character(raw.UPDRS$F5))
# raw.UPDRS$F6 <- as.numeric(as.character(raw.UPDRS$F6))
# raw.UPDRS$F7 <- as.numeric(as.character(raw.UPDRS$F7))

load(file='Z://PD_motor//subj_data//alldata.RData')

# Prepare data
updrs.long <- data.frame(id = rep(alldata$MEG_ID,2),
                         UPDRS.old = c(alldata$UPDRS_off,alldata$UPDRS_on),
                         group = rep(alldata$Sub_type,2),
                         session = rep(c(1,2),each=length(alldata$MEG_ID)),
                         n.id = rep(alldata$ID,2),
                         hand = rep(alldata$PAM_hand,2),
                         HY_stage = rep(alldata$HY_stage,2))

updrs.data <- merge(updrs.long,raw.UPDRS,by.x=c('n.id','session'), by.y=c('n.id','session'))
updrs.data$id.x <- paste('0',updrs.data$id.x, sep="")

u.neve.data <- merge(neve.data,updrs.data, by.x=c('subs','session','group'), by.y=c('id.x','session','group'))

# Flip lateralized factors for left hand subjects
F4.flip <- ifelse(u.neve.data$hand=="left",u.neve.data$F5,u.neve.data$F4)
F5.flip <- ifelse(u.neve.data$hand=="left",u.neve.data$F4,u.neve.data$F5)

# u.neve.data$F4 <- F4.flip
# u.neve.data$F5 <- F5.flip

save(u.neve.data, file='uNeveData.Rdata')

### New division ###
# Another way to divide the MDS-UPDRS-III is: 
# .	tremor (sum of items 15-18)
u.neve.data$tremor = u.neve.data$X3.15_postTremorHand_R+
  u.neve.data$X3.15_postTremorHand_L+
  u.neve.data$X3.16_kinTremorHands_L+
  u.neve.data$X3.16_kinTremorHands_R+
  u.neve.data$X3.17_restTremorAmp_L.J+
  u.neve.data$X3.17_restTremorAmp_RUE+
  u.neve.data$X3.17_restTremorAmp_LUE+
  u.neve.data$X3.17_restTremorAmp_RLE+
  u.neve.data$X3.17_restTremorAmp_LLE+
  u.neve.data$X3.18_consRestTremor
# .	rigidity (item 3)
u.neve.data$rigid = u.neve.data$X3.3_rigidity_neck+
  u.neve.data$X3.3_rigidity_RUE+
  u.neve.data$X3.3_rigidity_LUE+
  u.neve.data$X3.3_rigidity_RLE+
  u.neve.data$X3.3_rigidity_LLE
# .	bradykinesia (sum of items 2, 4-9 and 14) 
u.neve.data$brady = u.neve.data$X3.2_facial_exp+
  u.neve.data$X3.4_fingerTap_L+
  u.neve.data$X3.4_fingerTap_R+
  u.neve.data$X3.5_moveHands_L+
  u.neve.data$X3.5_movemHands_R+
  u.neve.data$X3.6_pro.supMoveHands_L+
  u.neve.data$X3.6_pro.supMoveHands_R+
  u.neve.data$X3.7_toeTap_R+
  u.neve.data$X3.7_toeTap_L+
  u.neve.data$X3.8_legAgil_L+
  u.neve.data$X3.8_legAgil_R+
  u.neve.data$X3.9_chair+
  u.neve.data$X3.14_globalSpont
# .	axial (sum of items 1 and 9-13). 
u.neve.data$axial = u.neve.data$X3.11_freeze+
  u.neve.data$X3.9_chair+
  u.neve.data$X3.10_gait+
  u.neve.data$X3.11_freeze+
  u.neve.data$X3.12_postSablil+
  u.neve.data$X3.13_posture

PD.data <- subset(u.neve.data, u.neve.data$group.x=="ptns")
PD.data$group.x <- factor(PD.data$group.x)
PD.data$id <- factor(PD.data$id)

# Super-long data
u.long <- melt(PD.data,
               # ID variables - all the variables to keep but not split apart on
               id.vars=c("subs","nevent","session"),
               # The source columns
               measure.vars=c("F1", "F2", "F3", "F4", "F5", "F6", "F7" ),
               # Name of the destination column that will identify the original
               # column that the measurement came from
               variable.name="factor",
               value.name="score"
)

u.long2 <- melt(PD.data,
                # ID variables - all the variables to keep but not split apart on
                id.vars=c("subs","nevent","session"),
                # The source columns
                measure.vars=c("tremor", "rigid", "brady", "axial" ),
                # Name of the destination column that will identify the original
                # column that the measurement came from
                variable.name="factor",
                value.name="score"
)

# Get summary
aggregate(u.long2$score, by=list(session=u.long2$session, factor=u.long2$factor),mean)
aggregate(u.long2$score, by=list(session=u.long2$session, factor=u.long2$factor),median)
aggregate(u.long2$score, by=list(session=u.long2$session, factor=u.long2$factor),range)