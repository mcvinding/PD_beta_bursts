# Make plots of UPDRS ~ MEG summary statistics based on brms models
library(ggplot2)
library(brms)

## Load data
wrkdir <- "Z://PD_motor//rest_ec//groupanalysis//"
wrkdir <- "C:\\Users\\Mikkel\\Documents\\PD-proj_betaEvent\\data"
setwd(wrkdir)
load(file='uData.Rdata')

## Setting
N <- 100      # Number of "error" lines


######### 
# REDO EVERYTHIN WITH buildin function marginal_effects and predict for individual values

# Make prediction like this (buy not for N event plot )
newdat <- PD.data.uF1
newdat$F4 <- 0
newdat$pred <- predict(br.nev.uF4, newdata = newdat)

mef <- plot(marginal_effects(br.nev.uF4), line_args=list(color="black"))
plt <- mef$F4 + 
  geom_point(data=PD.data.uF4, aes(x=F4, y=nevent, color=session), size=1,inherit.aes = FALSE) + 
  theme_bw() +
  scale_color_manual(values=c('red','blue'))+
  ylim(250,420)+
  xlab('Factor score') + ylab('N events')+
  ggtitle('F4: Stuff')
plt

############################################################################
## Plot N event
load(file = 'updrs_neve_mods.RData')
pd.nev.uF1 <- posterior_samples(br.nev.uF1, "^b")
PD.data.uF1 <- br.nev.uF1$data
pd.nev.uF2 <- posterior_samples(br.nev.uF2, "^b")
PD.data.uF2 <- br.nev.uF2$data
pd.nev.uF3 <- posterior_samples(br.nev.uF3, "^b")
PD.data.uF3 <- br.nev.uF3$data
pd.nev.uF4 <- posterior_samples(br.nev.uF4, "^b")
PD.data.uF4 <- br.nev.uF4$data
pd.nev.uF5 <- posterior_samples(br.nev.uF5, "^b")
PD.data.uF5 <- br.nev.uF5$data
pd.nev.uF6 <- posterior_samples(br.nev.uF6, "^b")
PD.data.uF6 <- br.nev.uF6$data
pd.nev.uF7 <- posterior_samples(br.nev.uF7, "^b")
PD.data.uF7 <- br.nev.uF7$data
pd.nev.uFT <- posterior_samples(br.nev.uFT, "^b")
PD.data.uFT <- br.nev.uFT$data



### Plot F1
x <- seq(0,max(PD.data.uF1$F1)+1,0.5)
y <- vector("list",length=N)
for(i in 1:N) {y[[i]] <- exp(pd.nev.uF1[i,1]+pd.nev.uF1[i,2]*x)}
df1 = data.frame(F1=rep(x,N),nevent=unlist(y))
df1$f = rep(1:N,each=length(x))

g1 <- ggplot(df1, aes(F1, nevent)) + 
  geom_line(alpha=1/2,col="grey",aes(group=f)) + 
  stat_summary(geom="line", fun.y=mean, color="black", lwd=0.5)+      # WRONG!!!!
  theme_bw()+
  geom_point(data=PD.data.uF1, aes(x=F1, y=nevent, color=session), size=1)+
  scale_color_manual(values=c('red','blue'))+
  xlim(0,max(x))+ylim(250,420)+
  xlab('Factor score') + ylab('N events')+
  ggtitle('F1: Midline function')
g1

### Plot F2
x <- seq(0,max(PD.data.uF2$F2)+1,0.5)
y <- vector("list",length=N)
for(i in 1:N) {y[[i]] <- exp(pd.nev.uF2[i,1]+pd.nev.uF2[i,2]*x)}
df2 = data.frame(F2=rep(x,N),nevent=unlist(y))
df2$f = rep(1:N,each=length(x))

g2 <- ggplot(df2, aes(F2, nevent)) + 
  geom_line(alpha=1/2,col="grey",aes(group=f)) + 
  stat_summary(geom="line", fun.y=mean, color="black", lwd=0.5)+
  theme_bw()+
  geom_point(data=PD.data.uF2, aes(x=F2, y=nevent, color=session), size=1)+
  scale_color_manual(values=c('red','blue'))+
  xlim(0,max(x))+ylim(250,420)+
  xlab('Factor score') + ylab('N events')+
  ggtitle('F2: Rest tremor')
g2

### Plot F3
x <- seq(0,max(PD.data.uF3$F3)+1,0.5)
y <- vector("list",length=N)
for(i in 1:N) {y[[i]] <- exp(pd.nev.uF3[i,1]+pd.nev.uF3[i,2]*x)}
df3 = data.frame(F3=rep(x,N),nevent=unlist(y))
df3$f = rep(1:N,each=length(x))

g3 <- ggplot(df3, aes(F3, nevent)) + 
  geom_line(alpha=1/2,col="grey",aes(group=f)) + 
  stat_summary(geom="line", fun.y=mean, color="black", lwd=0.5)+
  theme_bw()+
  geom_point(data=PD.data.uF3, aes(x=F3, y=nevent, color=session), size=1)+
  scale_color_manual(values=c('red','blue'))+
  xlim(0,max(x))+ylim(250,420)+
  xlab('Factor score') + ylab('N events')+
  ggtitle('F3: Rigidity')
g3

### Plot F4
x <- seq(0,max(PD.data.uF4$F4)+1,0.5)
y <- vector("list",length=N)
for(i in 1:N) {y[[i]] <- exp(pd.nev.uF4[i,1]+pd.nev.uF4[i,2]*x)}
df4 = data.frame(F4=rep(x,N),nevent=unlist(y))
df4$f = rep(1:N,each=length(x))

g4 <- ggplot(df4, aes(F4, nevent)) + 
  geom_line(alpha=1/2,col="grey",aes(group=f)) + 
  stat_summary(geom="line", fun.y=mean, color="black", lwd=0.5)+
  theme_bw()+
  geom_point(data=PD.data.uF4, aes(x=F4, y=nevent, color=session), size=1)+
  scale_color_manual(values=c('red','blue'))+
  xlim(0,max(x))+ylim(250,420)+
  xlab('Factor score') + ylab('N events')+
  ggtitle('F4: Bradykinesia right')
g4

### Plot F5
x <- seq(0,max(PD.data.uF5$F5)+1,0.5)
y <- vector("list",length=N)
for(i in 1:N) {y[[i]] <- exp(pd.nev.uF5[i,1]+pd.nev.uF5[i,2]*x)}
df5 = data.frame(F5=rep(x,N),nevent=unlist(y))
df5$f = rep(1:N,each=length(x))

g5 <- ggplot(df5, aes(F5, nevent)) + 
  geom_line(alpha=1/2,col="grey",aes(group=f)) + 
  stat_summary(geom="line", fun.y=mean, color="black", lwd=0.5)+
  theme_bw()+
  geom_point(data=PD.data.uF5, aes(x=F5, y=nevent, color=session), size=1)+
  scale_color_manual(values=c('red','blue'))+
  xlim(0,max(x))+ylim(250,420)+
  xlab('Factor score') + ylab('N events')+
  ggtitle('F5: Bradykinesia left')
g5

### Plot F6
x <- seq(0,max(PD.data.uF6$F6)+1,0.5)
y <- vector("list",length=N)
for(i in 1:N) {y[[i]] <- exp(pd.nev.uF6[i,1]+pd.nev.uF6[i,2]*x)}
df6 = data.frame(F6=rep(x,N),nevent=unlist(y))
df6$f = rep(1:N,each=length(x))

g6 <- ggplot(df6, aes(F6, nevent)) + 
  geom_line(alpha=1/2,col="grey",aes(group=f)) + 
  stat_summary(geom="line", fun.y=mean, color="black", lwd=0.5)+
  theme_bw()+
  geom_point(data=PD.data.uF6, aes(x=F6, y=nevent, color=session), size=1)+
  scale_color_manual(values=c('red','blue'))+
  xlim(0,max(x))+ylim(250,420)+
  xlab('Factor score') + ylab('N events')+
  ggtitle('F6: Postural and kinetic tremor')
g6

### Plot F7
x <- seq(0,max(PD.data.uF7$F7)+1,0.5)
y <- vector("list",length=N)
for(i in 1:N) {y[[i]] <- exp(pd.nev.uF7[i,1]+pd.nev.uF7[i,2]*x)}
df7 = data.frame(F7=rep(x,N),nevent=unlist(y))
df7$f = rep(1:N,each=length(x))

g7 <- ggplot(df7, aes(F7, nevent)) + 
  geom_line(alpha=1/2,col="grey",aes(group=f)) + 
  stat_summary(geom="line", fun.y=mean, color="black", lwd=0.5)+
  theme_bw()+
  geom_point(data=PD.data.uF7, aes(x=F7, y=nevent, color=session), size=1)+
  scale_color_manual(values=c('red','blue'))+
  xlim(0,max(x))+ylim(250,420)+
  xlab('Factor score') + ylab('N events')+
  ggtitle('F7: Lower limb bradykinesia')
g7

### Plot Total
x <- seq(0,max(PD.data.uFT$Total)+1,0.5)
y <- vector("list",length=N)
for(i in 1:N) {y[[i]] <- exp(pd.nev.uFT[i,1]+pd.nev.uFT[i,2]*x)}
dfT = data.frame(FT=rep(x,N),nevent=unlist(y))
dfT$f = rep(1:N,each=length(x))

gT <- ggplot(dfT, aes(FT, nevent)) + 
  geom_line(alpha=1/2,col="grey",aes(group=f)) + 
  stat_summary(geom="line", fun.y=mean, color="black", lwd=0.5)+
  theme_bw()+
  geom_point(data=PD.data.uFT, aes(x=Total, y=nevent, color=session), size=1)+
  scale_color_manual(values=c('red','blue'))+
  xlim(0,max(x))+ylim(250,420)+
  xlab('Factor score') + ylab('N events')+
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

ggsave("neve_F1.jpeg", plot=g1, device="jpeg", units="mm", width=40, height=35, dpi=500, scale=2.5)
ggsave("neve_F2.jpeg", plot=g2, device="jpeg", units="mm", width=40, height=35, dpi=500, scale=2.5)
ggsave("neve_F3.jpeg", plot=g3, device="jpeg", units="mm", width=40, height=35, dpi=500, scale=2.5)
ggsave("neve_F4.jpeg", plot=g4, device="jpeg", units="mm", width=40, height=35, dpi=500, scale=2.5)
ggsave("neve_F5.jpeg", plot=g5, device="jpeg", units="mm", width=40, height=35, dpi=500, scale=2.5)
ggsave("neve_F6.jpeg", plot=g6, device="jpeg", units="mm", width=40, height=35, dpi=500, scale=2.5)
ggsave("neve_F7.jpeg", plot=g7, device="jpeg", units="mm", width=40, height=35, dpi=500, scale=2.5)
ggsave("neve_Total.jpeg", plot=gT, device="jpeg", units="mm", width=40, height=35, dpi=500, scale=2.5)

save.image(".RData")

############################################################################
## Plot time between events
getItiDat <- function(mod,Fx,Scl=1){
  b <- posterior_samples(mod, "^b")
  temp1 <- posterior_samples(mod, "^r")
  temp1 <- exp(temp1+b$b_Intercept)*Scl
  temp2 <- rapply(temp1,mean)
  session <- as.factor(ifelse(grepl("_1",names(temp2)),"1","2"))
  subs <- vector()
  for (i in 1:length(names(temp2))){subs[i] <- substr(names(temp2)[i],16,21)}
  dd <- data.frame(subs=subs,iti=temp2,session=session)
  tempF <- aggregate(get(Fx,mod$data), by=list(mod$data$`subs:session`,mod$data$session), median)
  outdat <- merge(dd,tempF, by.x=c("subs","session"),by.y=c("Group.1","Group.2"))
  names(outdat) <- c("ss","session","iti",Fx)
  
  x <- seq(0,max(get(Fx,outdat))+1,0.5)
  y <- vector("list",length=N)
  for(i in 1:N) {y[[i]] <- exp(b[i,1]+b[i,2]*x)*Scl}
  df2 = data.frame(x=rep(x,N),y=unlist(y))
  names(df2) <- c(Fx,"iti")
  df2$f = rep(1:N,each=length(x))
  yp <- exp(fixef(mod)[1]+fixef(mod)[2]*x)*Scl
  df2p <- data.frame(y=yp,x=x)
  names(df2p) <- c("iti",Fx)
  
  output <- list(df=outdat, allpred=df2, fixpred=df2p)
  return(output)
}

load(file = 'updrs_iti_mods.RData')

### Plot F1 iti
d <- getItiDat(br.iti.uF1,"F1")

giti1 <- ggplot(d$allpred, aes(F1, iti)) + 
  geom_line(alpha=1/2,col="grey",aes(group=f))+ 
  geom_line(data = d$fixpred, aes(x=F1, y=iti))+
  geom_point(data=d$df, aes(x=F1, y=iti, color=session), size=1)+
  scale_color_manual(values=c('red','blue'))+
  theme_bw()+
  xlim(0,max(d$fixpred$F1))+
  xlab('Factor score') + ylab('Time between events (ms)')+
  ggtitle('F1: Midline function')
giti1

### Plot F2 iti
d <- getItiDat(br.iti.uF2,"F2",1000)

giti2 <- ggplot(d$allpred, aes(F2, iti)) + 
  geom_line(alpha=1/2,col="grey",aes(group=f)) + 
  geom_line(data = d$fixpred, aes(x=F2, y=iti))+
  geom_point(data=d$df, aes(x=F2, y=iti, color=session), size=1)+
  scale_color_manual(values=c('red','blue'))+
  theme_bw()+
  xlim(0,max(d$fixpred$F2))+
  xlab('Factor score') + ylab('Time between events (ms)')+
  ggtitle('F2: Rest tremor')
giti2

### Plot F3 iti
d <- getItiDat(br.iti.uF3,"F3",1000)
giti3 <- ggplot(d$allpred, aes(F3, iti)) + 
  geom_line(alpha=1/2,col="grey",aes(group=f))+
  geom_line(data=d$fixpred, aes(x=F3, y=iti))+
  geom_point(data=d$df, aes(x=F3, y=iti, color=session), size=1)+
  scale_color_manual(values=c('red','blue'))+
  xlim(0,max(d$fixpred$F3))+
  theme_bw()+
  xlab('Factor score') + ylab('Time between events (ms)')+
  ggtitle('F3: Rigidity')
giti3

### Plot F4 iti
d <- getItiDat(br.iti.uF4,"F4",1000)
giti4 <- ggplot(d$allpred, aes(F4, iti)) + 
  geom_line(alpha=1/2,col="grey",aes(group=f))+
  geom_line(data=d$fixpred, aes(x=F4, y=iti))+
  geom_point(data=d$df, aes(x=F4, y=iti, color=session), size=1)+
  scale_color_manual(values=c('red','blue'))+
  xlim(0,max(d$fixpred$F4))+
  theme_bw()+
  xlab('Factor score') + ylab('Time between events (ms)')+
  ggtitle('F4: Bradykinesia right')
giti4

### Plot F5 iti
d <- getItiDat(br.iti.uF5,"F5")
giti5 <- ggplot(d$allpred, aes(F5, iti)) + 
  geom_line(alpha=1/2,col="grey",aes(group=f))+
  geom_line(data=d$fixpred, aes(x=F5, y=iti))+
  geom_point(data=d$df, aes(x=F5, y=iti, color=session), size=1)+
  scale_color_manual(values=c('red','blue'))+
  xlim(0,max(d$fixpred$F5))+
  theme_bw()+
  xlab('Factor score') + ylab('Time between events (ms)')+
  ggtitle('F5: Bradykinesia left')
giti5

### Plot F6 iti
d <- getItiDat(br.iti.uF6,"F6")
giti6 <- ggplot(d$allpred, aes(F6, iti)) + 
  geom_line(alpha=1/2,col="grey",aes(group=f))+
  geom_line(data=d$fixpred, aes(x=F6, y=iti))+
  geom_point(data=d$df, aes(x=F6, y=iti, color=session), size=1)+
  scale_color_manual(values=c('red','blue'))+
  xlim(0,max(d$fixpred$F6))+
  theme_bw()+
  xlab('Factor score') + ylab('Time between events (ms)')+
  ggtitle('F6: Postural and kinetic tremor')
giti6

### Plot F7 iti
d <- getItiDat(br.iti.uF7,"F7")
giti7 <- ggplot(d$allpred, aes(F7, iti)) + 
  geom_line(alpha=1/2,col="grey",aes(group=f))+
  geom_line(data=d$fixpred, aes(x=F7, y=iti))+
  geom_point(data=d$df, aes(x=F7, y=iti, color=session), size=1)+
  scale_color_manual(values=c('red','blue'))+
  xlim(0,max(d$fixpred$F7))+
  theme_bw()+
  xlab('Factor score') + ylab('Time between events (ms)')+
  ggtitle('F7: Lower limb bradykinesia')
giti7

### Plot Total iti
d <- getItiDat(br.iti.uFT,"Total")
gitiT <- ggplot(d$allpred, aes(Total, iti)) + 
  geom_line(alpha=1/2,col="grey",aes(group=f))+
  geom_line(data=d$fixpred, aes(x=Total, y=iti))+
  geom_point(data=d$df, aes(x=Total, y=iti, color=session), size=1)+
  scale_color_manual(values=c('red','blue'))+
  xlim(0,max(d$fixpred$Total))+
  theme_bw()+
  xlab('Factor score') + ylab('Time between events (ms)')+
  ggtitle('Total UPDRS-III')
gitiT

############################################################################
## Plot event duration
getLenDat <- function(mod,Fx,Scl=1){
  b <- posterior_samples(mod, "^b")
  temp1 <- posterior_samples(mod, "^r")
  temp1 <- exp(temp1+b$b_Intercept)*Scl
  temp2 <- rapply(temp1,mean)
  session <- as.factor(ifelse(grepl("_1",names(temp2)),"1","2"))
  subs <- vector()
  for (i in 1:length(names(temp2))){subs[i] <- substr(names(temp2)[i],16,21)}
  dd <- data.frame(subs=subs,iti=temp2,session=session)
  tempF <- aggregate(get(Fx,mod$data), by=list(mod$data$`subs:session`,mod$data$session), median)
  outdat <- merge(dd,tempF, by.x=c("subs","session"),by.y=c("Group.1","Group.2"))
  names(outdat) <- c("ss","session","iti",Fx)
  
  x <- seq(0,max(get(Fx,outdat))+1,0.5)
  y <- vector("list",length=N)
  for(i in 1:N) {y[[i]] <- exp(b[i,1]+b[i,2]*x)*Scl}
  df2 = data.frame(x=rep(x,N),y=unlist(y))
  names(df2) <- c(Fx,"iti")
  df2$f = rep(1:N,each=length(x))
  yp <- exp(fixef(mod)[1]+fixef(mod)[2]*x)*Scl
  df2p <- data.frame(y=yp,x=x)
  names(df2p) <- c("iti",Fx)
  
  output <- list(df=outdat, allpred=df2, fixpred=df2p)
  return(output)
}

load(file = 'updrs_len_mods.RData')



#R code for the above (fewer simulations than I did, but enough)
#
z <- rnorm(1000000,0,1) # 1 million normals
delta <- 2*exp(1/2)
y <- exp(z)             # 2-parameter lognormal
x <- y + delta          # shifted (i.e. 3-parameter) lognormal
par(mfrow=c(1,2)) 
hist(z,n=200,col="skyblue",bord="skyblue",freq=FALSE)
hist(log(x),n=200,col="lightgreen",bord="lightgreen",xlim=c(1,5),freq=FALSE)



