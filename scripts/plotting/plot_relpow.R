# Plot relative beta power in confetti plot
library(ggplot2)

# Define paths
wrkdir <- "Z://PD_motor//rest_ec//groupanalysis//"
outdir <- "Z://PD_motor//rest_ec//figures//"

setwd(wrkdir)

# Load data
load(file = 'reldat.RData')
rel.dat$task <- paste(rel.dat$group,rel.dat$session)
b.summary <- aggregate(rel.dat$relpow, list(rel.dat$task), mean)
names(b.summary) <- c("task","mean")
b.summary.sd <- aggregate(rel.dat$relpow, list(rel.dat$task), sd)
b.summary$sd <- b.summary.sd$x
b.summary$se <- b.summary$sd/sqrt(19)*2

## Plot
set.seed(999)
rel.plt <- ggplot(rel.dat, aes(x=task, y=relpow))+
  geom_crossbar(data=b.summary, aes(x=task, y=mean, ymin=mean, ymax=mean), width = 0.5) +
  geom_errorbar(data=b.summary, aes(x=task, y=mean, ymin=mean-sd, ymax=mean+sd), width=0.2) +
  geom_point(aes(fill=subs), position=position_jitter(width = 0.15), color='black', shape=21, size=2) +
  # scale_y_continuous(limits=c(225,475), breaks=seq(200,500,50))+
  labs(title="Relative power",
       x="",
       y = "")+
  scale_x_discrete(labels=c("PD/1 (OFF)","PD/2 (ON)","Ctrl/1","Ctrl/2"))+
  theme_bw() +
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5, size=rel(2), face="bold"),
        axis.title = element_blank(),
        axis.text = element_text(color="black", face="bold", size=rel(1)),
        # panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line())
rel.plt

# Save
ggsave(paste(outdir,"relBpow.jpeg",sep=""), plot=rel.plt,
       device="png", units="mm", width=30, height=20, dpi=600, scale=4)

# END