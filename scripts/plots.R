## Pretty plots
library(ggplot2)

# Plot N events
neve.data$task <- paste(neve.data$group,neve.data$session)
n.summary <- aggregate(neve.data$nevent, list(neve.data$task), mean)
names(n.summary) <- c("task","mean")
n.summary.sd <- aggregate(neve.data$nevent, list(neve.data$task), sd)
n.summary$sd <- n.summary.sd$x
n.summary$se <- n.summary$sd/sqrt(19)*2

nplt <- ggplot(neve.data, aes(x=task, y=nevent))+
  geom_point(aes(fill=subs),position=position_jitter(width = 0.2),color='black',shape=21,size=2) +
  geom_crossbar(data=n.summary, aes(x=task,y=mean, ymin=mean,ymax=mean), width = 0.5) +
  geom_errorbar(data=n.summary, aes(x=task,y=mean, ymin=mean-sd, ymax=mean+sd), width=0.2) +
  guides(fill=F) +
  ylim(225, 475) +
  xlab("Group/session") + ylab("N events")
nplt

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
nplt <- nplt + lay

ggsave("neve_group.jpeg", plot=nplt, device="jpeg", units="mm", width=40, height=35, dpi=500, scale=4)

# Plot iti
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
  xlab("Group/session") + ylab("Time between events (ms)") +
  lay
i.plt

ggsave("iti_group.jpeg", plot=i.plt, device="jpeg", units="mm", width=40, height=35, dpi=500, scale=4)

# Plot iti
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
  xlab("Group/session") + ylab("Amplitude (Z-score[?])") +
  lay
m.plt

ggsave("max_group.jpeg", plot=m.plt, device="jpeg", units="mm", width=40, height=35, dpi=500, scale=4)


## Inspection of distributions (pooled data)
i.den <- ggplot(itieve.data, aes(x=eve.iti.ms, fill=group))+
  # geom_histogram(binwidth=10)+
  geom_density(alpha=.3) +
  facet_wrap(~session) + 
  xlim(0,4000)+xlab("Time between events") +
  theme_bw()
ggsave("iti_denst.jpeg", plot=i.den, device="jpeg", units="mm", width=40, height=35, dpi=500, scale=4)


l.den <- ggplot(leneve.data, aes(x=eve.len.ms, fill=group))+
  # geom_histogram(binwidth=10)+
  geom_density(alpha=.3) +
  facet_wrap(~session) + 
  xlim(0,500)+xlab("Event duration") +
  theme_bw()
ggsave("len_dens.jpeg", plot=l.den, device="jpeg", units="mm", width=40, height=35, dpi=500, scale=4)

m.den <- ggplot(maxeve.data, aes(x=eve.max, fill=group))+
  # geom_histogram(binwidth=10)+
  geom_density(alpha=.3) +
  facet_wrap(~session) + 
  xlab("Max peak") +
  theme_bw()
ggsave("max_dens.jpeg", plot=i.den, device="jpeg", units="mm", width=40, height=35, dpi=500, scale=4)


## Cumulative plots of burst iti
ggplot(itieve.data, aes(eve.iti.ms, colour = group:session)) + 
  stat_ecdf(pad = FALSE)+
  xlim(0,5000) + theme_bw() +
  facet_wrap(~subs)
