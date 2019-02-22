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
