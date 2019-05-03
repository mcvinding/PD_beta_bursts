## Pretty plots
library(ggplot2)

# Load data
wrkdir <- "Z://PD_motor//rest_ec//groupanalysis//"
outdir <- "Z://PD_motor//rest_ec//figures//"
setwd(wrkdir)

# Plot N events
load(file = 'neve.RData')
neve.data$task <- paste(neve.data$group,neve.data$session)
n.summary <- aggregate(neve.data$nevent, list(neve.data$task), mean)
names(n.summary) <- c("task","mean")
n.summary.sd <- aggregate(neve.data$nevent, list(neve.data$task), sd)
n.summary$sd <- n.summary.sd$x
n.summary$se <- n.summary$sd/sqrt(19)*2

## Plot
nplt <- ggplot(neve.data, aes(x=task, y=nevent))+
  geom_crossbar(data=n.summary, aes(x=task,y=mean, ymin=mean, ymax=mean), width = 0.5) +
  geom_errorbar(data=n.summary, aes(x=task,y=mean, ymin=mean-sd, ymax=mean+sd), width=0.2) +
  geom_point(aes(fill=subs),position=position_jitter(width = 0.15),color='black',shape=21,size=2) +
  scale_y_continuous(limits=c(225,475), breaks=seq(200,500,50))+
  labs(title="Number of events",
       x='Group/Session',
       y = "N events")+
  scale_x_discrete(labels=c("Ctrl/1","Ctrl/2","PD/1 (OFF)","PD/2 (ON)"))+
  theme_bw() +
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5, size=rel(2), face="bold"),
        axis.title = element_text(face="bold", size=rel(1.5)),
        axis.text = element_text(color="black", size=rel(1.2)),
        panel.grid = element_blank())
  
nplt

## Save
ggsave(paste(outdir,"neve_group.jpeg",sep=""), plot=nplt, device="png", units="mm", width=60, height=40, dpi=600, scale=3)

# Plot iti
load(file = 'itieve.RData')
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

# Plot max
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


# Inspection of distributions (pooled data)
## Time between events 
i.den <- ggplot(itieve.data, aes(x=eve.iti.ms, fill=group))+
  # geom_histogram(binwidth=10)+
  geom_density(alpha=.3) +
  facet_wrap(.~session, labeller=as_labeller(c("1"="1/OFF","2"="2/ON"))) +
  xlim(0,3000)+
  theme_bw()+
  scale_fill_manual(values=c('red','blue'), label=c("Ctrl","PD"))+
  labs(x='Time between events (ms)',
       y = "Density",
       fill = "Group")+
  theme(legend.position=c(.99,.99),
        legend.justification=c(1,1),
        legend.background=element_rect(fill='white',color=NA),
        legend.title = element_text(face="bold"),
        strip.background = element_rect(fill="white"),
        strip.text = element_text(face="bold",size=rel(1.5)),
        plot.title = element_text(hjust = 0.5, size=rel(2), face="bold"),
        axis.title = element_text(face="bold", size=rel(1.5)),
        axis.text = element_text(color="black", size=rel(1.2)),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())
i.den
ggsave(paste(outdir,"iti_denst.jpeg",sep=""),
       plot=i.den, device="jpeg", units="mm", width=60, height=40, dpi=600, scale=4)

## Event duration
l.den <- ggplot(leneve.data, aes(x=eve.len.ms, fill=group))+
  # geom_histogram(binwidth=10)+
  geom_density(alpha=.3) +
  facet_wrap(~session, labeller=as_labeller(c("1"="1/OFF","2"="2/ON"))) + 
  xlim(0,350)+xlab("Event duration") +
  theme_bw()+
  scale_fill_manual(values=c('red','blue'), label=c("Ctrl","PD"))+
  labs(x='Event duration (ms)',
       y = "Density",
       fill = "Group")+
  theme(legend.position=c(.99,.99),
        legend.justification=c(1,1),
        legend.background=element_rect(fill='white',color=NA),
        legend.title = element_text(face="bold"),
        strip.background = element_rect(fill="white"),
        strip.text = element_text(face="bold",size=rel(1.5)),
        plot.title = element_text(hjust = 0.5, size=rel(2), face="bold"),
        axis.title = element_text(face="bold", size=rel(1.5)),
        axis.text = element_text(color="black", size=rel(1.2)),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())
l.den
ggsave(paste(outdir,"len_dens.jpeg",sep=""), plot=l.den, device="png", units="mm", width=60, height=40, dpi=600, scale=4)

## Max power
m.den <- ggplot(maxeve.data, aes(x=eve.max, fill=group))+
  # geom_histogram(binwidth=10)+
  geom_density(alpha=.3) +
  facet_wrap(~session) + 
  xlab("Max peak") +
  theme_bw()
m.den
ggsave("max_dens.jpeg", plot=i.den, device="jpeg", units="mm", width=40, height=35, dpi=500, scale=4)
