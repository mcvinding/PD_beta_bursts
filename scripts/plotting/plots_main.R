## Pretty plots
library(ggplot2)

# Load data
wrkdir <- "Z://PD_motor//rest_ec//groupanalysis//"
outdir <- "Z://PD_motor//rest_ec//figures//"
setwd(wrkdir)

#####################################################################################
# Plot N events
load(file = 'neve.RData')

## Make summary
neve.data$task <- paste(neve.data$group, neve.data$session)
n.summary <- aggregate(neve.data$nevent.min, list(neve.data$task), mean)
names(n.summary) <- c("task","mean")
n.summary.sd <- aggregate(sub.summary$x, list(sub.summary$Group.1), sd)
n.summary$sd <- n.summary.sd$x
n.summary$se <- n.summary$sd/sqrt(19)*2

## Plot
set.seed(9000)
nplt <- ggplot(neve.data, aes(x=task, y=nevent.min))+
  geom_crossbar(data=n.summary, aes(x=task,y=mean, ymin=mean, ymax=mean), width = 0.5) +
  geom_errorbar(data=n.summary, aes(x=task,y=mean, ymin=mean-sd, ymax=mean+sd), width=0.2) +
  geom_point(aes(fill=subs),position=position_jitter(width = 0.15),color='black',shape=21,size=2) +
  scale_y_continuous(limits=c(75,160), breaks=seq(80,140,20))+
  labs(title="Beta burst rate",
       x='Group/Session',
       y = "Burst/min")+
  scale_x_discrete(labels=c("Control/1","Control/2","PD/1 (OFF)","PD/2 (ON)"))+
  theme_bw() +
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5, size=rel(2), face="bold"),
        axis.title = element_text(face="bold", size=rel(1.5)),
        axis.text = element_text(color="black", size=rel(1.2)),
        panel.grid = element_blank())
nplt

## Save
ggsave(paste(outdir,"neve_group.jpeg",sep=""), plot=nplt, device="png", units="mm", width=60, height=40, dpi=600, scale=3)


#####################################################################################
# Plot iti (confetti plots; not used in publication)
load(file = 'itieve.RData')

# Make summary
itieve.data$task <- paste(itieve.data$group,itieve.data$session)
i.pltdata <- aggregate(itieve.data$log.eve.iti, list(itieve.data$task, itieve.data$subs), mean)
names(i.pltdata) <- c("task","subs","log.mean")
i.pltdata$median <- exp(i.pltdata$log.mean)*1000
i.summary <- aggregate(i.pltdata$log.mean, list(i.pltdata$task), mean)
names(i.summary) <- c("task","log.mean")
i.summary.sd <- aggregate(i.pltdata$log.mean, list(i.pltdata$task), sd)
i.summary$log.sd.l <- i.summary$log.mean-i.summary.sd$x
i.summary$log.sd.u <- i.summary$log.mean+i.summary.sd$x
i.summary$log.ci.l <- i.summary$log.mean-i.summary.sd$x/sqrt(19)*2
i.summary$log.ci.u <- i.summary$log.mean+i.summary.sd$x/sqrt(19)*2
i.summary$median <- exp(i.summary$log.mean)*1000
i.summary$sd.l <- exp(i.summary$log.sd.l)*1000
i.summary$sd.u <- exp(i.summary$log.sd.u)*1000

i.plt <- ggplot(i.pltdata, aes(x=task, y=median))+
  geom_point(aes(fill=subs),position=position_jitter(width = 0.2),color='black',shape=21,size=2) +
  geom_crossbar(data=i.summary, aes(x=task,y=median, ymin=median,ymax=median), width = 0.5) +
  geom_errorbar(data=i.summary, aes(x=task,y=median, ymin=sd.l, ymax=sd.u), width=0.2) +
  guides(fill=F)+
  # scale_y_continuous(limits=c(225,475), breaks=seq(200,500,50))+
  labs(title="Time between events",
       x='Group/Session',
       y = "Time (ms)")+
  scale_x_discrete(labels=c("Ctrl/1","Ctrl/2","PD/1 (OFF)","PD/2 (ON)"))+
  theme_bw() +
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5, size=rel(2), face="bold"),
        axis.title = element_text(face="bold", size=rel(1.5)),
        axis.text = element_text(color="black", size=rel(1.2)),
        panel.grid = element_blank())
i.plt

ggsave(paste(outdir,"iti_group.jpeg",sep=""),
       plot=i.plt, device="jpeg", units="mm", width=60, height=40, dpi=600, scale=3)

#####################################################################################
# Plot max (confetti plot; not used in publication)
load(file = 'maxeve.RData')
maxeve.data$task <- paste(maxeve.data$group,maxeve.data$session)
m.pltdata <- aggregate(maxeve.data$log.eve.max, list(maxeve.data$task, maxeve.data$subs), mean)
names(m.pltdata) <- c("task","subs","log.mean")
m.pltdata$median <- exp(m.pltdata$log.mean)
m.summary <- aggregate(m.pltdata$log.mean, list(m.pltdata$task), mean)
names(m.summary) <- c("task","log.mean")
m.summary.sd <- aggregate(i.pltdata$log.mean, list(i.pltdata$task), sd)
m.summary$log.sd.l <- m.summary$log.mean-m.summary.sd$x
m.summary$log.sd.u <- m.summary$log.mean+m.summary.sd$x
m.summary$log.ci.l <- m.summary$log.mean-m.summary.sd$x/sqrt(19)*2
m.summary$log.ci.u <- m.summary$log.mean+m.summary.sd$x/sqrt(19)*2
m.summary$median <- exp(m.summary$log.mean)
m.summary$sd.l <- exp(m.summary$log.sd.l)
m.summary$sd.u <- exp(m.summary$log.sd.u)

m.plt <- ggplot(m.pltdata, aes(x=task, y=median))+
  geom_point(aes(fill=subs),position=position_jitter(width = 0.2),color='black',shape=21,size=2) +
  geom_crossbar(data=m.summary, aes(x=task,y=median, ymin=median,ymax=median), width = 0.5) +
  geom_errorbar(data=m.summary, aes(x=task,y=median, ymin=sd.l, ymax=sd.u), width=0.2) +
  guides(fill=F)+
  xlab("Group/session") + ylab("Amplitude (F-score)")
m.plt

ggsave(paste(outdir,"max_group.jpeg",sep=""), plot=m.plt, device="jpeg", units="mm", width=40, height=35, dpi=500, scale=4)

#####################################################################################
# Inspection of distributions (pooled data)
#####################################################################################
#####################################################################################
## Time between events 
load(file = 'itieve.RData')
load(file='itieve_pp.RData')

itieve.data$label <-  as.factor(paste(itieve.data$group, itieve.data$session, sep = " ", collapse = NULL))
levels(itieve.data$label) <- c('Controls 1', 'Controls 2', 'PD 1/OFF', 'PD 2/ON')
itieve.data$label <- factor(itieve.data$label, 
                            levels =c('PD 1/OFF','Controls 1','PD 2/ON','Controls 2'))

i.den <- ggplot(itieve.data, aes(x=eve.iti.ms, fill=group))+
  facet_wrap(.~label) +
  geom_density(bw=20, alpha=.5) +
  geom_line(stat="density",data=pp.iti.long, aes(x=pp, fill=NA), bw=20, linetype = "dashed", size=0.5, color='black') +
  xlim(0,1500)+
  labs(x='Time between events (ms)', y = "Density")+
  theme_bw() +
  scale_fill_manual(values=c('red','blue'), guide=FALSE)+
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(face="bold", lineheight=12),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, lineheight=12, face="bold"),
        axis.title = element_text(face="bold", lineheight=11),
        axis.text = element_text(color="black", lineheight=9),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())+
  ggtitle('Inter-event interval')
i.den
ggsave(paste(outdir,"iti_denst.jpeg",sep=""),
       plot=i.den, device="png", units="cm", width=8, height=5, dpi=500, scale=2)

#####################################################################################
## Event duration
load(file = 'leneve.RData')
load(file='leneve_pp.RData')

leneve.data$label <-  as.factor(paste(leneve.data$group, leneve.data$session, sep = " ", collapse = NULL))
levels(leneve.data$label) <- c('Controls 1', 'Controls 2', 'PD 1/OFF', 'PD 2/ON')
leneve.data$label <- factor(leneve.data$label, 
                            levels =c('PD 1/OFF','Controls 1','PD 2/ON','Controls 2'))

l.den <- ggplot(leneve.data, aes(x=eve.len.ms, fill=group))+
  # geom_histogram(binwidth=10)+
  geom_density(bw=2.5,alpha=.5) +
  geom_line(stat="density",data=pp.len.long, aes(x=pp, fill=NA), bw=2.5, linetype = "dashed", size=0.5, color='black') +
  facet_wrap(~label) + 
  xlim(0,300)+
  labs(x='Duration (ms)', y = "Density")+
  theme_bw()+
  scale_fill_manual(values=c('red','blue'), guide=FALSE)+
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(face="bold", lineheight=12),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, lineheight=12, face="bold"),
        axis.title = element_text(face="bold", lineheight=11),
        axis.text = element_text(color="black", lineheight=9),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) +
  ggtitle('Beta event duration')
l.den

ggsave(paste(outdir,"len_dens.jpeg",sep=""), 
       plot=l.den, device="png", units="cm", width=8, height=4, dpi=500, scale=2)

#####################################################################################
load(file = 'maxeve.RData')

maxeve.data$label <-  as.factor(paste(maxeve.data$group, maxeve.data$session, sep = " ", collapse = NULL))
levels(maxeve.data$label) <- c('Controls 1', 'Controls 2', 'PD 1/OFF', 'PD 2/ON')
maxeve.data$label <- factor(maxeve.data$label, 
                            levels =c('PD 1/OFF','Controls 1','PD 2/ON','Controls 2'))


## Max power
m.den <- ggplot(maxeve.data, aes(x=eve.max, fill=group))+
  geom_density(alpha=.5) +
  geom_line(stat="density", data=pp.max.long, aes(x=pp, fill=NA), linetype = "dashed", size=0.5, color='black') +
  facet_wrap(~label) + 
  xlim(0,5)+
  xlab("Max peak") +
  labs(x='Amplitude (F-score)', y = "Density")+
  theme_bw()+
  scale_fill_manual(values=c('red','blue'), guide=FALSE)+
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(face="bold", lineheight=12),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, lineheight=12, face="bold"),
        axis.title = element_text(face="bold", lineheight=11),
        axis.text = element_text(color="black", lineheight=9),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())+
  ggtitle('Maximum amplitude')
m.den
ggsave(paste(outdir,"max_dens.jpeg",sep=""),
       plot=m.den, device="jpeg", units="cm", width=8, height=4, dpi=500, scale=2)
