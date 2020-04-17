# Plot extended analysis
library(ggplot2)

wrkdir <- "C://Users//Mikkel//Documents//betabursts//groupanalysis"
outdir <- "C://Users//Mikkel//Documents//betabursts//Figures//"
setwd(wrkdir)

#######################################################################################
# Plot BF across steps
load(file='range_bf.RData')

## Arragne data
bf.dat <- data.frame("steps"=seq(0.5,4,by=0.1),
                    "test"=rep(c('bf10','bf21','bf32'), each=36),
                    "BF"=c(BF10,BF21,BF32))
bf.dat$logBF <- log(bf.dat$BF)

## Plot BF across steps
plt <- ggplot(bf.dat, aes(x=steps, y=logBF, fill=test))+
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=log(1/3), ymax=log(3), alpha=0.1, fill="red") +
  geom_hline(yintercept=log(10),lty=2, color='blue') +
  geom_hline(yintercept=log(1/10),lty=2, color='blue') +
  geom_hline(yintercept=0,lty=2, alpha=0.5) +
  geom_vline(xintercept=1.3, lty=3) +
  scale_x_continuous(limits=c(0.45, 4.05), breaks=seq(1,4,1))+
  geom_line()+
  geom_point(shape=21, size=3)+theme_bw()+
  scale_fill_discrete(labels = c("Session", "Group","Session x Group"),
                      name = "Main effect") +
  labs(title="Main effects across thresholds",
       x='Threshold (median + median * x)')+
  ylim(c(-5,10))+
  theme(legend.position=c(0,1),
        legend.justification=c(0,1),
        legend.background=element_blank(),
        legend.title = element_text(face="bold", size=13),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size=16, face="bold"),
        axis.title = element_text(face="bold", size=13),
        axis.text = element_text(color="black", size=12),
        plot.margin = unit(c(0, 0, 0, 0), "cm"))

plt

## Save
ggsave(paste(outdir,"extended_bf.png", sep=""), plt, 
       dpi=600, width=10, height=5, units="cm", scale=3)


#######################################################################################
# Plot N across steps
load(file = 'neve_ext.RData')

neve.mean <- data.frame(aggregate(neve.data$nevent.min, by=list(neve.data$steps,neve.data$group,neve.data$session),mean))
neve.sd <- data.frame(aggregate(neve.data$nevent.min, by=list(neve.data$steps,neve.data$group,neve.data$session),sd))
neve.summary <- merge(neve.mean, neve.sd, by=c("Group.1","Group.2","Group.3"))
colnames(neve.summary) <- c("steps","group","session","mean","sd")
  
n.plt <- ggplot(neve.summary, aes(x=steps, y=mean, color=group, shape=session))+
  geom_vline(xintercept=1.3, lty=3)+
  geom_line(position=position_dodge(0.05))+
  geom_errorbar(aes(x=steps, ymin=mean-sd, ymax=mean+sd), width=0.05, size=0.5,position=position_dodge(0.05)) + 
  geom_point(size=2, position=position_dodge(0.05))+
  scale_shape_manual(values=c(16,17), label=c("1/OFF","2/ON"))+
  scale_color_manual(values=c("red","blue"), label=c("Control","PD"))+
  scale_x_continuous(limits=c(0.45, 4.05), breaks=seq(1,4,1))+
  scale_y_continuous(limits=c(0.5, 300), breaks=seq(0,300,50))+
  theme_bw()+
  labs(title="Burst rate across thresholds",
       x='Threshold (median + median * x)',
       y='Burst/min',
       color="Group", shape="Session")+
  theme(legend.position=c(.99,.99),
        legend.justification=c(1,1),
        legend.background=element_blank(), #element_rect(fill='white',color=NA),
        legend.title = element_text(face="bold", size=13),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size=16, face="bold"),
        axis.title = element_text(face="bold", size=13),
        axis.text = element_text(color="black", size=12),
        plot.margin = unit(c(0, 0, 0, 0), "cm"))
n.plt  

# Save
ggsave(paste(outdir,"extended_neve.png", sep=""), n.plt, 
       dpi=600, width=10, height=5, units="cm", scale=3)

#######################################################################################
# Plot ideal observer analysis across steps
load(file='range_roc.RData')

## Arragne data
roc.dat <- data.frame("steps"=seq(0.5,4,by=0.1),
                      "test"=rep(c('Session 1/OFF','Session 2/ON'), each=36),
                      "ROC"=c(roc1,roc2))

ios.plt <- ggplot(roc.dat, aes(x=steps, y=ROC, color=test))+
  geom_vline(xintercept=1.3, lty=3)+
  geom_hline(yintercept=0.5, lty=2)+
  geom_line(position=position_dodge(0.05))+
  geom_point(size=2, position=position_dodge(0.05))+
  scale_color_manual(values=c("red","blue"), label=c("1/OFF","2/ON"))+
  scale_x_continuous(limits=c(0.45, 4.05), breaks=seq(1,4,1))+
  scale_y_continuous(limits=c(0.49, 1), breaks=seq(0.5,1,0.1))+
  theme_bw()+
  labs(title="Ideal observer analysis",
       x='Threshold (median + median * x)',
       y='AUROC',
       color="Session")+
  theme(legend.position=c(0.01,1),
        legend.justification=c(0,1),
        legend.background=element_blank(), #element_rect(fill='white',color=NA),
        legend.title = element_text(face="bold", size=13),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size=16, face="bold"),
        axis.title = element_text(face="bold", size=13),
        axis.text = element_text(color="black", size=12),
        plot.margin = unit(c(0, 0, 0, 0), "cm"))
ios.plt  

ggsave(paste(outdir,"extended_ios.png", sep=""), ios.plt, 
       dpi=600, width=6, height=3, units="cm", scale=3)

# END  