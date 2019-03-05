# Make plots of UPDRS ~ MEG summary statistics based on brms models
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


############################################################################
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