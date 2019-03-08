## Pretty plots
library(ggplot2)

# Plot N events
neve.data$task <- paste(neve.data$group,neve.data$session)
n.summary <- aggregate(neve.data$nevent, list(neve.data$task), mean)
names(n.summary) <- c("task","mean")
n.summary.sd <- aggregate(neve.data$nevent, list(neve.data$task), sd)
n.summary$sd <- n.summary.sd$x
n.summary$se <- n.summary$sd/sqrt(20)*2

nplt <- ggplot(neve.data, aes(x=task, y=nevent))+
  geom_point(aes(fill=subs),position=position_jitter(width = 0.2),color='black',shape=21,size=2) +
  geom_crossbar(data=n.summary, aes(x=task,y=mean, ymin=mean,ymax=mean), width = 0.5) +
  geom_errorbar(data=n.summary,aes(x=task,y=mean, ymin=mean-sd, ymax=mean+sd), width=0.2) +
  guides(fill=F) + theme_bw()
  # ylim(225, 475)
nplt

## Inspection of distributions (pooled data)
ggplot(itieve.data, aes(x=eve.iti.ms, fill=group))+
  # geom_histogram(binwidth=10)+
  geom_density(alpha=.3) +
  facet_wrap(~session) + 
  xlim(0,4000)+xlab("Time between events") +
  theme_bw()

ggplot(leneve.data, aes(x=eve.len.ms, fill=group))+
  # geom_histogram(binwidth=10)+
  geom_density(alpha=.3) +
  facet_wrap(~session) + 
  xlim(0,500)+xlab("Event duration") +
  theme_bw()

ggplot(maxeve.data, aes(x=eve.max, fill=group))+
  # geom_histogram(binwidth=10)+
  geom_density(alpha=.3) +
  facet_wrap(~session) + 
  xlab("Max peak") +
  theme_bw()

## Cumulative plots of burst iti
ggplot(itieve.data, aes(eve.iti.ms, colour = group:session)) + 
  stat_ecdf(pad = FALSE)+
  xlim(0,5000) + theme_bw() +
  facet_wrap(~subs)
