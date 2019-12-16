# Plot 1/f params and beta power in confetti plot
library(ggplot2)

# Define paths
#wrkdir <- "Z://PD_motor//rest_ec//groupanalysis//"
wrkdir <- "C://Users//Mikkel//Documents//betabursts//groupanalysis"
outdir <- "C://Users//Mikkel//Documents//betabursts//Figures//"
setwd(wrkdir)

# Load data
f.dat <- read.csv("df_aper.csv", header=TRUE, sep=";")
f.dat$group <- as.factor(f.dat$group)
f.dat$group <- revalue(f.dat$group, c("1"="ptns", "2"="ctrl"))
f.dat$session <- as.factor(f.dat$session)
f.dat$subj <- as.factor(rep(1:(length(f.dat$intercept)/2),2))

# Prepare data
f.dat$task <- paste(f.dat$group, f.dat$session)

fi.summary <- aggregate(f.dat$intercept, list(f.dat$task), mean)
names(fi.summary) <- c("task","mean")
fi.summary.sd <- aggregate(f.dat$intercept, list(f.dat$task), sd)
fi.summary$sd <- fi.summary.sd$x

fs.summary <- aggregate(f.dat$slope, list(f.dat$task), mean)
names(fs.summary) <- c("task","mean")
fs.summary.sd <- aggregate(f.dat$slope, list(f.dat$task), sd)
fs.summary$sd <- fs.summary.sd$x

##############################################################################
## Plot intercept
set.seed(666)
icp.plt <- ggplot(f.dat, aes(x=task, y=intercept))+
  geom_crossbar(data=fi.summary, aes(x=task, y=mean, ymin=mean, ymax=mean), width = 0.5) +
  geom_errorbar(data=fi.summary, aes(x=task, y=mean, ymin=mean-sd, ymax=mean+sd), width=0.2) +
  geom_point(aes(fill=subj), position=position_jitter(width = 0.15), color='black', shape=21, size=2) +
  scale_y_continuous(limits=c(-3.2,-0.8), breaks=seq(-3.0,-1.0,0.5))+
  labs(title="1/f intercept",
       x="",
       y = "")+
  scale_x_discrete(labels=c("Ctrl/1","Ctrl/2","PD/1 (OFF)","PD/2 (ON)"))+
  theme_bw() +
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5, size=rel(2), face="bold"),
        axis.title = element_blank(),
        axis.text = element_text(color="black", face="bold", size=rel(1)),
        # panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line())
icp.plt

# Save
ggsave(paste(outdir,"oof_intercept.jpeg",sep=""), plot=icp.plt,
       device="png", units="mm", width=30, height=20, dpi=600, scale=4)

##############################################################################
## Plot slope
set.seed(666)
slp.plt <- ggplot(f.dat, aes(x=task, y=slope))+
  geom_crossbar(data=fs.summary, aes(x=task, y=mean, ymin=mean, ymax=mean), width = 0.5) +
  geom_errorbar(data=fs.summary, aes(x=task, y=mean, ymin=mean-sd, ymax=mean+sd), width=0.2) +
  geom_point(aes(fill=subj), position=position_jitter(width = 0.15), color='black', shape=21, size=2) +
  scale_y_continuous(limits=c(0.3,1.3), breaks=seq(0.4,1.3,0.2))+
  labs(title="1/f slope",
       x="",
       y = "")+
  scale_x_discrete(labels=c("Ctrl/1","Ctrl/2","PD/1 (OFF)","PD/2 (ON)"))+
  theme_bw() +
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5, size=rel(2), face="bold"),
        axis.title = element_blank(),
        axis.text = element_text(color="black", face="bold", size=rel(1)),
        # panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line())
slp.plt

# Save
ggsave(paste(outdir,"oof_slope.jpeg",sep=""), plot=slp.plt,
       device="png", units="mm", width=30, height=20, dpi=600, scale=4)


##############################################################################
## Plot beta
# Load data
b.dat <- read.csv("df_beta.csv", header=TRUE, sep=";")
b.dat$group <- as.factor(b.dat$group)
b.dat$group <- revalue(b.dat$group, c("1"="ptns", "2"="ctrl"))
b.dat$session <- as.factor(b.dat$session)
b.dat$subj <- as.factor(rep(1:(length(b.dat$peak_pow)/2),2))

# Prepare data
b.dat$task <- paste(b.dat$group, b.dat$session)

b.summary <- aggregate(b.dat$peak_pow, list(b.dat$task), mean)
names(b.summary) <- c("task","mean")
b.summary.sd <- aggregate(b.dat$peak_pow, list(b.dat$task), sd)
b.summary$sd <- b.summary.sd$x

set.seed(666)
bpk.plt <- ggplot(b.dat, aes(x=task, y=peak_pow))+
  geom_crossbar(data=b.summary, aes(x=task, y=mean, ymin=mean, ymax=mean), width = 0.5) +
  geom_errorbar(data=b.summary, aes(x=task, y=mean, ymin=mean-sd, ymax=mean+sd), width=0.2) +
  geom_point(aes(fill=subj), position=position_jitter(width = 0.15), color='black', shape=21, size=2) +
  # scale_y_continuous(limits=c(0.1,0.85), breaks=seq(0.2,0.8,0.1))+
  labs(title="Beta peak power (- 1/f)",
       x="",
       y = "")+
  scale_x_discrete(labels=c("Ctrl/1","Ctrl/2","PD/1 (OFF)","PD/2 (ON)"))+
  theme_bw() +
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5, size=rel(2), face="bold"),
        axis.title = element_blank(),
        axis.text = element_text(color="black", face="bold", size=rel(1)),
        # panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line())
bpk.plt

# Save
ggsave(paste(outdir,"oof_beta.jpeg",sep=""), plot=bpk.plt,
       device="png", units="mm", width=30, height=20, dpi=600, scale=4)

#END