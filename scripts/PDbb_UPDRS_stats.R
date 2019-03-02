## Statistics 2: within group correlation with clinical scores (MSD-UPDRS-III)
# Define paths
wrkdir <- "Z://PD_motor//rest_ec//groupanalysis//"
setwd(wrkdir)
load(file = 'workspace.Rdata')
# load(file = '.Rdata')

#################### Stats (mean) ##############################
load(file='uData.Rdata')

# Not used
library(lme4)

all.mod1 <- lmer(nevent~Total + (1|subs), data=u.neve.data, REML=F, subset = u.neve.data$group.x=="ptns")
summary(all.mod1)
all.mod0 <- lmer(nevent~1 + (1|subs), data=u.neve.data, REML=F, subset = u.neve.data$group.x=="ptns")
anova(all.mod0,all.mod1)
all.mod2 <- lmer(nevent~Total+session + (1|subs), data=u.neve.data, REML=F, subset = u.neve.data$group.x=="ptns")
anova(all.mod0,all.mod1,all.mod2)
all.mod3 <- lmer(nevent~Total*session + (1|subs), data=u.neve.data, REML=F, subset = u.neve.data$group.x=="ptns")
anova(all.mod0,all.mod1,all.mod2,all.mod3)

# # Summary
# # pam1patients <- subset(pam1data,Sub_type=='patient')
# t.test(PD.data$Total~PD.data$session, paired=T)
# ttestBF(x=PD.data$Total[PD.data$session=="1"],y=PD.data$Total[PD.data$session=="2"], data=u.neve.data, paired=T, rscale='wide')
# 
# t.test(x=u.neve.data$UPDRS_off[u.neve.data$Sub_type=='control'] , y=u.neve.data$UPDRS_on, paired=F)
# ttestBF(x=pam1data$UPDRS_off[pam1data$Sub_type=='control'], y=PD.data$Total[PD.data$session=="2"], paired=F, rscale='wide')
# ttestBF(x=pam1data$UPDRS_off[pam1data$Sub_type=='control'], y=PD.data$Total[PD.data$session=="1"], paired=F, rscale='wide')
# 
# aggregate(u.neve.data$Total, by=list(u.neve.data$group.x,u.neve.data$session), mean)
# aggregate(u.neve.data$Total, by=list(u.neve.data$group.x,u.neve.data$session), sd)

# Bayes factor (used)
load(file='uData.Rdata')
library(BayesFactor)

bf0 <- lmBF(nevent~subs, data=PD.data, whichRandom="subs")
bf1 <- lmBF(nevent~session+subs, data=PD.data, whichRandom="subs")
bf1/bf0

PD.data.x <- PD.data[complete.cases(PD.data),]
bf0 <- lmBF(nevent~subs, data=PD.data, whichRandom="subs")
bf1 <- lmBF(nevent~session*subs, data=PD.data, whichRandom=c("subs","session"))
bfT <- lmBF(nevent~Total+session+subs, data=PD.data, whichRandom=c("subs","session"))

bfT/bf1

bf.f1 <- lmBF(nevent~F1+session*subs, whichRandom=c("subs","session"),data=PD.data)
bf.f1/bf1
bf.f2 <- lmBF(nevent~F2+session*subs, whichRandom=c("subs","session"),data=PD.data)
bf.f2/bf1
bf.f3 <- lmBF(nevent~F3+session*subs, whichRandom=c("subs","session"),data=PD.data)
bf.f3/bf1
bf.f4 <- lmBF(nevent~F4+session*subs, whichRandom=c("subs","session"),data=PD.data)
bf.f4/bf1
bf.f5 <- lmBF(nevent~F5+session*subs, whichRandom=c("subs","session"),data=PD.data)
bf.f5/bf1
bf.f6 <- lmBF(nevent~F6+session*subs, whichRandom=c("subs","session"),data=PD.data)
bf.f6/bf1
bf.f7 <- lmBF(nevent~F7+session*subs, whichRandom=c("subs","session"),data=PD.data)
bf.f7/bf1

mF1.post <- posterior(bf.f1, iterations = 2000)
summary(mF1.post)
mF2.post <- posterior(bf.f2, iterations = 2000)
summary(mF2.post)
mF3.post <- posterior(bf.f3, iterations = 2000)
summary(mF3.post)
mF4.post <- posterior(bf.f4, iterations = 2000)
summary(mF4.post)
mF5.post <- posterior(bf.f5, iterations = 2000)
summary(mF5.post)
mF6.post <- posterior(bf.f6, iterations = 2000)
summary(mF6.post)
mF7.post <- posterior(bf.f7, iterations = 2000)
summary(mF7.post)

mT.post <- posterior(bfT, iterations = 2000)
summary(mT.post)

# Conventional correlation (just for show)
cor.test(PD.data$F1,PD.data$nevent)
cor.test(PD.data$F2,PD.data$nevent)
cor.test(PD.data$F3,PD.data$nevent)
cor.test(PD.data$F4,PD.data$nevent)
cor.test(PD.data$F5,PD.data$nevent)
cor.test(PD.data$F6,PD.data$nevent)
cor.test(PD.data$F7,PD.data$nevent)
cor.test(PD.data$Total,PD.data$nevent)

##
save.image(".RData")

### New division ###
bf0 <- lmBF(nevent~subs, data=PD.data, whichRandom="subs")
bf1 <- lmBF(nevent~session*subs, data=PD.data, whichRandom=c("subs"))

bf.tremor <- lmBF(nevent~tremor+session*subs, whichRandom=c("subs"),data=PD.data)
bf.tremor/bf1
bf.tremor.post <- posterior(bf.tremor, iterations = 2000)

bf.rigid <- lmBF(nevent~rigid+session*subs, whichRandom=c("subs","session"),data=PD.data)
bf.rigid/bf1
bf.rigid.post <- posterior(bf.rigid, iterations = 2000)

bf.axial <- lmBF(nevent~axial+session*subs, whichRandom=c("subs"),data=PD.data)
bf.axial/bf1
bf.axial.post <- posterior(bf.axial, iterations = 2000)

bf.brady <- lmBF(nevent~brady+session*subs, whichRandom=c("subs"),data=PD.data)
bf.brady/bf1
bf.brady.post <- posterior(bf.brady, iterations = 2000)

# Score
bf0 <- lmBF(score~subs, data=u.long2, whichRandom="subs")
bf1 <- lmBF(score ~ factor+session*subs, data=u.long2, whichRandom=c("subs"))
bf.F <- lmBF(score ~ factor+nevent+session*subs, data=u.long2, whichRandom=c("subs"))
bf.F/bf1
bf.Fx <- lmBF(score ~ factor*rebound+session*id, data=u.long2, whichRandom=c("id","session"))
bf.Fx/bf1

mT.post <- posterior(bf.Fx, iterations = 2000)
summary(mT.post)

# Conventional correlation (just for show)
cor.test(PD.data$tremor,PD.data$nevent)
cor.test(PD.data$rigid,PD.data$nevent)
cor.test(PD.data$axial,PD.data$nevent)
cor.test(PD.data$brady,PD.data$nevent)


plot(PD.data$tremor,PD.data$nevent)
plot(PD.data$rigid,PD.data$nevent)
plot(PD.data$axial,PD.data$nevent)
plot(PD.data$brady,PD.data$nevent)


# ## Each item
# bf0 <- lmBF(rebound~session*id, data=PD.data, whichRandom=c("id","session"))
# 
# bf.x1 <- lmBF(rebound~X3.1+session*id, whichRandom=c("id","session"),data=PD.data)
# bf.x1/bf0
# bf.x2 <- lmBF(rebound~X3.2+session*id, whichRandom=c("id","session"),data=PD.data)
# bf.x2/bf0
# bf.x3 <- lmBF(rebound~X3.3+session*id, whichRandom=c("id","session"),data=PD.data)
# bf.x3/bf0
# bf.x4 <- lmBF(rebound~X3.4+session*id, whichRandom=c("id","session"),data=PD.data)
# bf.x4/bf0
# bf.x5 <- lmBF(rebound~X3.5+session*id, whichRandom=c("id","session"),data=PD.data)
# bf.x5/bf0
# bf.x6 <- lmBF(rebound~X3.6+session*id, whichRandom=c("id","session"),data=PD.data)
# bf.x6/bf0
# bf.x7 <- lmBF(rebound~X3.7+session*id, whichRandom=c("id","session"),data=PD.data)
# bf.x7/bf0
# bf.x8 <- lmBF(rebound~X3.8+session*id, whichRandom=c("id","session"),data=PD.data)
# bf.x8/bf0
# bf.x9 <- lmBF(rebound~X3.9+session*id, whichRandom=c("id","session"),data=PD.data)
# bf.x9/bf0
# bf.x10 <- lmBF(rebound~X3.10+session*id, whichRandom=c("id","session"),data=PD.data)
# bf.x10/bf0
# bf.x11 <- lmBF(rebound~X3.11+session*id, whichRandom=c("id","session"),data=PD.data)
# bf.x11/bf0
# bf.x12 <- lmBF(rebound~X3.12+session*id, whichRandom=c("id","session"),data=PD.data)
# bf.x12/bf0
# bf.x13 <- lmBF(rebound~X3.13+session*id, whichRandom=c("id","session"),data=PD.data)
# bf.x13/bf0
# bf.x14 <- lmBF(rebound~X3.14+session*id, whichRandom=c("id","session"),data=PD.data)
# bf.x14/bf0
# bf.x15 <- lmBF(rebound~X3.15+session*id, whichRandom=c("id","session"),data=PD.data)
# bf.x15/bf0
# bf.x16 <- lmBF(rebound~X3.16+session*id, whichRandom=c("id","session"),data=PD.data)
# bf.x16/bf0
# bf.x17 <- lmBF(rebound~X3.17+session*id, whichRandom=c("id","session"),data=PD.data)
# bf.x17/bf0
# bf.x18 <- lmBF(rebound~X3.18+session*id, whichRandom=c("id","session"),data=PD.data)
# bf.x18/bf0

############################ PLOTS #######################################
library(ggplot2)

### Plot F1
N=dim(mF1.post)[1]
x = 0:max(PD.data$F1)
y = vector("list",length=N)
for(i in 1:N) {y[[i]] = mF1.post[i,1]+mF1.post[i,2]*x}
df = data.frame(F1=rep(x,N),nevent=unlist(y))
df$f = rep(1:N,each=max(PD.data$F1)+1)

g1 <- ggplot(df, aes(F1, nevent)) + 
  geom_line(alpha=1/25,col="grey",aes(group=f)) + 
  stat_summary(geom="line", fun.y=mean, color="black", lwd=0.5)+
  theme_bw()+
  geom_point(data=PD.data, aes(x=F1, y=nevent, color=session), size=1)+
  scale_color_manual(values=c('red','blue'))+
  xlim(0,max(PD.data$F1))+ylim(250,400)+
  xlab('Factor score') + ylab('N events')+
  ggtitle('F1: Midline function')
g1

### Plot F2
N=dim(mF2.post)[1]
x = 0:max(PD.data$F2)
y = vector("list",length=N)
for(i in 1:N) {y[[i]] = mF2.post[i,1]+mF2.post[i,2]*x}
df = data.frame(F2=rep(x,N),nevent=unlist(y))
df$f = rep(1:N,each=max(PD.data$F2)+1)

g2 <- ggplot(df, aes(F2, nevent)) + 
  geom_line(alpha=1/25,col="grey",aes(group=f)) + 
  stat_summary(geom="line", fun.y=mean, color="black", lwd=0.5)+
  theme_bw()+
  geom_point(data=PD.data, aes(x=F2, y=nevent, color=session), size=1)+
  scale_color_manual(values=c('red','blue'))+
  xlim(0,max(PD.data$F2))+ylim(250,400)+
  xlab('Factor score') + ylab('N events')+
  ggtitle('F2: Rest tremor')
g2

### Plot F3
N=dim(mF3.post)[1]
x = 0:max(PD.data$F3)
y = vector("list",length=N)
for(i in 1:N) {y[[i]] = mF3.post[i,1]+mF3.post[i,2]*x}
df = data.frame(F3=rep(x,N),nevent=unlist(y))
df$f = rep(1:N,each=max(PD.data$F3)+1)

g3 <- ggplot(df, aes(F3, nevent)) + 
  geom_line(alpha=1/25,col="grey",aes(group=f)) + 
  stat_summary(geom="line", fun.y=mean, color="black", lwd=0.5)+
  theme_bw()+
  geom_point(data=PD.data, aes(x=F3, y=nevent, color=session), size=1)+
  scale_color_manual(values=c('red','blue'))+
  xlim(0,max(PD.data$F3))+ylim(250,400)+
  xlab('Factor score') + ylab('N events')+
  ggtitle('F3: Rigidity')
g3

### Plot F4
N=dim(mF4.post)[1]
x = 0:max(PD.data$F4)
y = vector("list",length=N)
for(i in 1:N) {y[[i]] = mF4.post[i,1]+mF4.post[i,2]*x}
df = data.frame(F4=rep(x,N),nevent=unlist(y))
df$f = rep(1:N,each=max(PD.data$F4)+1)

g4 <- ggplot(df, aes(F4, nevent)) + 
  geom_line(alpha=1/25,col="grey",aes(group=f)) + 
  stat_summary(geom="line", fun.y=mean, color="black", lwd=0.5)+
  theme_bw()+
  geom_point(data=PD.data, aes(x=F4, y=nevent, color=session), size=1)+
  scale_color_manual(values=c('red','blue'))+
  xlim(0,max(PD.data$F4))+ylim(250,400)+
  xlab('Factor score') + ylab('N events')+
  ggtitle('F4: Bradykinesia right')
g4

### Plot F5
N=dim(mF5.post)[1]
x = 0:max(PD.data$F5)
y = vector("list",length=N)
for(i in 1:N) {y[[i]] = mF5.post[i,1]+mF5.post[i,2]*x}
df = data.frame(F5=rep(x,N),nevent=unlist(y))
df$f = rep(1:N,each=max(PD.data$F5)+1)

g5 <- ggplot(df, aes(F5, nevent)) + 
  geom_line(alpha=1/25,col="grey",aes(group=f)) + 
  stat_summary(geom="line", fun.y=mean, color="black", lwd=0.5)+
  theme_bw()+
  geom_point(data=PD.data, aes(x=F5, y=nevent, color=session), size=1)+
  scale_color_manual(values=c('red','blue'))+
  xlim(0,max(PD.data$F5))+ylim(250,400)+
  xlab('Factor score') + ylab('N events')+
  ggtitle('F5: Bradykinesia left')
g5

### Plot F6
N=dim(mF6.post)[1]
x = 0:max(PD.data$F6)
y = vector("list",length=N)
for(i in 1:N) {y[[i]] = mF6.post[i,1]+mF6.post[i,2]*x}
df = data.frame(F6=rep(x,N),nevent=unlist(y))
df$f = rep(1:N,each=max(PD.data$F6)+1)

g6 <- ggplot(df, aes(F6, nevent)) + 
  geom_line(alpha=1/25,col="grey",aes(group=f)) + 
  stat_summary(geom="line", fun.y=mean, color="black", lwd=0.5)+
  theme_bw()+
  geom_point(data=PD.data, aes(x=F6, y=nevent, color=session), size=1)+
  scale_color_manual(values=c('red','blue'))+
  xlim(0,max(PD.data$F6))+ylim(250,400)+
  xlab('Factor score') + ylab('N events')+
  ggtitle('F6: Postural and kinetic tremor')
g6

### Plot F7
N=dim(mF7.post)[1]
x = 0:max(PD.data$F7)
y = vector("list",length=N)
for(i in 1:N) {y[[i]] = mF7.post[i,1]+mF7.post[i,2]*x}
df = data.frame(F7=rep(x,N),nevent=unlist(y))
df$f = rep(1:N,each=max(PD.data$F7)+1)

g7 <- ggplot(df, aes(F7, nevent)) + 
  geom_line(alpha=1/25,col="grey",aes(group=f)) + 
  stat_summary(geom="line", fun.y=mean, color="black", lwd=0.5)+
  theme_bw()+
  geom_point(data=PD.data, aes(x=F7, y=nevent, color=session), size=1)+
  scale_color_manual(values=c('red','blue'))+
  xlim(0,max(PD.data$F7))+ylim(250,400)+
  xlab('Factor score') + ylab('N events')+
  ggtitle('F7: Lower limb bradykinesia')
g7

### Plot Total
N=dim(mT.post)[1]
x = 0:max(PD.data$Total)
y = vector("list",length=N)
for(i in 1:N) {y[[i]] = mT.post[i,1]+mT.post[i,2]*x}
df = data.frame(Total=rep(x,N),nevent=unlist(y))
df$f = rep(1:N,each=max(PD.data$Total)+1)

gT <- ggplot(df, aes(Total, nevent)) + 
  geom_line(alpha=1/25,col="grey",aes(group=f)) + 
  stat_summary(geom="line", fun.y=mean, color="black", lwd=0.5)+
  theme_bw()+
  geom_point(data=PD.data, aes(x=Total, y=nevent, color=session), size=1)+
  scale_color_manual(values=c('red','blue'))+
  xlim(0,max(PD.data$Total))+ylim(250,400)+
  ggtitle('Total UPDRS-III')+
  xlab('Total score') + ylab('N events')+
  ggtitle('Total UPDRS-III')
gT


# Add layout
lay <- theme_bw() + theme(legend.position = "none",
                          text = element_text(size = 11, family="Helvetica"),
                          title = element_text(size = 12, vjust = 1.5, face="bold",lineheight = NULL),
                          # axis.text = element_text(size=10),
                          axis.title = element_text(size = 11, vjust = .5, face="plain"),
                          # axis.title.x = element_text(face="plain", size=4),
                          # axis.title.y = element_text(face="plain", size=4),
                          axis.text = element_text(face="bold", size=11),
                          panel.grid = element_blank())
  
  
g1 <- g1 + lay
g2 <- g2 + lay
g3 <- g3 + lay
g4 <- g4 + lay
g5 <- g5 + lay
g6 <- g6 + lay
g7 <- g7 + lay
gT <- gT + lay
  
ggsave("F1.jpeg", plot=g1, device="jpeg", units="mm", width=35, height=35, dpi=500, scale=2.5)
ggsave("F2.jpeg", plot=g2, device="jpeg", units="mm", width=35, height=35, dpi=500, scale=2.5)
ggsave("F3.jpeg", plot=g3, device="jpeg", units="mm", width=35, height=35, dpi=500, scale=2.5)
ggsave("F4.jpeg", plot=g4, device="jpeg", units="mm", width=35, height=35, dpi=500, scale=2.5)
ggsave("F5.jpeg", plot=g5, device="jpeg", units="mm", width=35, height=35, dpi=500, scale=2.5)
ggsave("F6.jpeg", plot=g6, device="jpeg", units="mm", width=35, height=35, dpi=500, scale=2.5)
ggsave("F7.jpeg", plot=g7, device="jpeg", units="mm", width=35, height=35, dpi=500, scale=2.5)
ggsave("Total.jpeg", plot=gT, device="jpeg", units="mm", width=35, height=35, dpi=500, scale=2.5)

save.image(".RData")


## Plot new division
### Plot F1
N=dim(bf.tremor.post)[1]
x = 0:max(PD.data$tremor)
y = vector("list",length=N)
for(i in 1:N) {y[[i]] = bf.tremor.post[i,1]+bf.tremor.post[i,2]*x}
df = data.frame(tremor=rep(x,N),nevent=unlist(y))
df$f = rep(1:N,each=max(PD.data$tremor)+1)

ggplot(df, aes(tremor, nevent)) + 
  geom_line(alpha=1/25,col="grey",aes(group=f)) + 
  stat_summary(geom="line", fun.y=mean, color="black", lwd=0.5)+
  theme_bw()+
  geom_point(data=PD.data, aes(x=tremor, y=nevent, color=session), size=1)+
  scale_color_manual(values=c('red','blue'))+
  xlim(0,max(PD.data$tremor))+ylim(-0.2,0.4)+
  xlab('Factor score') + ylab('Beta nevent %-change')+
  ggtitle('Tremor')

N=dim(bf.rigid.post)[1]
x = 0:max(PD.data$rigid)
y = vector("list",length=N)
for(i in 1:N) {y[[i]] = bf.rigid.post[i,1]+bf.rigid.post[i,2]*x}
df = data.frame(rigid=rep(x,N),nevent=unlist(y))
df$f = rep(1:N,each=max(PD.data$rigid)+1)

ggplot(df, aes(rigid, nevent)) + 
  geom_line(alpha=1/25,col="grey",aes(group=f)) + 
  stat_summary(geom="line", fun.y=mean, color="black", lwd=0.5)+
  theme_bw()+
  geom_point(data=PD.data, aes(x=rigid, y=nevent, color=session), size=1)+
  scale_color_manual(values=c('red','blue'))+
  xlim(0,max(PD.data$rigid))+ylim(-0.2,0.4)+
  xlab('Factor score') + ylab('Beta nevent %-change')+
  ggtitle('Rigid')




