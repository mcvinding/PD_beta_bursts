# Plot extended analysis
library(ggplot2)

wrkdir <- "Z://PD_motor//rest_ec//groupanalysis//"
outdir <- "Z://PD_motor//rest_ec//figures//"
setwd(wrkdir)
load(file='range_bf.RData')

# Arragne data
bf.dat <- data.frame("steps"=seq(0.1,5,by=0.1),
                    "test"=rep(c('bf10','bf21','bf32'), each=50),
                    "BF"=c(BF10,BF21,BF32))
bf.dat$logBF <- log(bf.dat$BF)

# Plot BF across steps
plt <- ggplot(bf.dat, aes(x=steps, y=logBF, fill=test))+
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=log(1/3), ymax=log(3), alpha=0.2, fill="red") +
  # geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=log(1/3), ymax=log(3)), fill='red', alpha=0.1)+
  geom_hline(yintercept=log(10),lty=2, color='red')+
  geom_hline(yintercept=log(1/10),lty=2, color='red')+
  geom_hline(yintercept=0,lty=2, alpha=0.5)+
  geom_vline(xintercept=1.3)+
  geom_line()+
  geom_point(shape=21, size=2)+theme_bw()+
  scale_fill_discrete(labels = c("Session > Intercept", "Session+Group > Session","Session*Group > Session+Group"),
                      name = "Model comparison") +
  labs(title="Comparison across thresholds",
       x='Threshold (Med+x*Med)')+
  ylim(c(-10,19))+
  theme(legend.position=c(0.01,0.99),
        legend.justification=c(0,1),
        legend.background=element_rect(fill='white',color=NA),
        legend.title = element_text(face="bold"),
        plot.title = element_text(hjust = 0.5, size=rel(2), face="bold"),
        axis.title = element_text(face="bold", size=rel(1.5)),
        axis.text = element_text(color="black", size=rel(1.2)))
plt
# Save
ggsave(paste(outdir,"extended_bf.png", sep=""), plt, dpi=600, width=10 ,height=5, units="cm", scale=3)

# Plot N across steps
load(file = 'neve_ext.RData')

neve.mean <- data.frame(aggregate(neve.data$nevent, by=list(neve.data$steps,neve.data$group,neve.data$session),mean))
neve.sd <- data.frame(aggregate(neve.data$nevent, by=list(neve.data$steps,neve.data$group,neve.data$session),sd))
neve.summary <- merge(neve.mean,neve.sd,by=c("Group.1","Group.2","Group.3"))
colnames(neve.summary) <- c("steps","group","session","mean","sd")
  
n.plt <- ggplot(neve.summary, aes(x=steps, y=mean,color=group, shape=session))+
  geom_vline(xintercept=1.3)+
  geom_line(position=position_dodge(0.05))+
  geom_errorbar(aes(x=steps, ymin=mean-sd, ymax=mean+sd), width=0.05, size=0.5,position=position_dodge(0.05)) + 
  geom_point(size=2, position=position_dodge(0.05))+
  # scale_shape_manual(values=c(21,23), label=c("1/OFF","2/ON"))+
  scale_shape_manual(values=c(16,17), label=c("1/OFF","2/ON"))+
  scale_color_manual(values=c("red","blue"))+
  theme_bw()+
  labs(title="N events across thresholds",
       x='Threshold (Med+x*Med)',
       y='N events (Mean±sd)',
       fill="Group", shape="Session")+
  theme(legend.position=c(.99,.99),
        legend.justification=c(1,1),
        legend.background=element_rect(fill='white',color=NA),
        legend.title = element_text(face="bold"),
        plot.title = element_text(hjust = 0.5, size=rel(2), face="bold"),
        axis.title = element_text(face="bold", size=rel(1.5)),
        axis.text = element_text(color="black", size=rel(1.2)))
n.plt  

# Save
ggsave(paste(outdir,"extended_neve.png", sep=""), n.plt, dpi=600, width=10 ,height=5, units="cm", scale=3)

# END  