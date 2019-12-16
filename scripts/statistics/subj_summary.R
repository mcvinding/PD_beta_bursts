#### Data on participant for PAM1 (rebound) part ####
library(BayesFactor)

load(file='Z://PD_motor//subj_data///alldata.RData')

# Descriptive summary: all subjects recruited
aggregate(alldata$age,list(alldata$Sub_type),mean)
aggregate(alldata$age,list(alldata$Sub_type),range)
gendertab <- xtabs(~sex+Sub_type,alldata)
gendertab

# Reject subjects
bbdata <- subset(alldata, alldata$MEG_ID != '345')
bbdata <- subset(bbdata, bbdata$MEG_ID != '364')
bbdata$MEG_ID <- factor(bbdata$MEG_ID)

# Summary and comparison: subjects in analysis
## Age
aggregate(bbdata$age,list(bbdata$Sub_type),mean)
aggregate(bbdata$age,list(bbdata$Sub_type),range)
t.test(age~Sub_type, bbdata)
ttestBF(formula=age~Sub_type, data=bbdata, paired=F)

## Sex
gendertab <- xtabs(~sex+Sub_type,bbdata)
fisher.test(gendertab)
contingencyTableBF(gendertab, sampleType = 'indepMulti',fixedMargin='cols')   # Testbetween colums (i.e group)

## MoCA
aggregate(bbdata$MoCA, by=list(bbdata$Sub_type), mean)
aggregate(bbdata$MoCA, by=list(bbdata$Sub_type), sd)
t.test(MoCA~Sub_type, bbdata)
ttestBF(formula=MoCA~Sub_type, data=bbdata, paired=F)

## HADS
aggregate(bbdata$HADS_angst, by=list(bbdata$Sub_type), mean)
aggregate(bbdata$HADS_angst, by=list(bbdata$Sub_type), sd)
t.test(HADS_angst~Sub_type, bbdata)
ttestBF(formula=HADS_angst~Sub_type, data=bbdata, paired=F)

aggregate(bbdata$HADS_depression, by=list(bbdata$Sub_type), mean)
aggregate(bbdata$HADS_depression, by=list(bbdata$Sub_type), sd)
t.test(HADS_depression~Sub_type, bbdata)
ttestBF(formula=HADS_depression~Sub_type, data=bbdata, paired=F)

## UPDRS (Ptns vs Ctrl)
aggregate(bbdata$UPDRS_off, by=list(bbdata$Sub_type), mean)
aggregate(bbdata$UPDRS_off, by=list(bbdata$Sub_type), sd)
aggregate(bbdata$UPDRS_off, by=list(bbdata$Sub_type), range)
aggregate(bbdata$UPDRS_off, by=list(bbdata$Sub_type), median)
t.test(UPDRS_off~Sub_type, bbdata)
ttestBF(formula=UPDRS_off~Sub_type, data=bbdata, paired=F)

bbpatients <- subset(bbdata, bbdata$Sub_type == 'patient')
range(bbpatients$UPDRS_on, na.rm=T)
median(bbpatients$UPDRS_on, na.rm=T)

bbpatients$UPDRS_diff <- bbpatients$UPDRS_off - bbpatients$UPDRS_on
bbpatients$UPDRS_pctdiff <- bbpatients$UPDRS_diff*100 / bbpatients$UPDRS_off

mean(bbpatients$UPDRS_diff)
median(bbpatients$UPDRS_diff)
range(bbpatients$UPDRS_diff)

mean(bbpatients$UPDRS_pctdiff)
median(bbpatients$UPDRS_pctdiff)
range(bbpatients$UPDRS_pctdiff)

ttestBF(bbpatients$UPDRS_off,bbpatients$UPDRS_on, paired = T)

## Disesae duration and medication
median(bbpatients$disease_dur, na.rm=T)
range(bbpatients$disease_dur, na.rm=T)
median(bbpatients$LEDD, na.rm=T)
mean(bbpatients$LEDD, na.rm=T)
sd(bbpatients$LEDD, na.rm=T)
range(bbpatients$LEDD, na.rm=T)
