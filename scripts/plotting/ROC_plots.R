# ROC plots
library(ROCR)

# Load data
load(file="C:/Users/Mikkel/Documents/betabursts/groupanalysis/roc_results.Rdata")
setwd("C:/Users/Mikkel/Documents/betabursts/Figures")

#################################################################################
# ROC plot function
ROCplot <- function(roc1, roc2, heading, col.1="black", col.2="black", leg.cex=1){
  
  # Plot 1
  plot(roc1$roc, colorize=F,  main = heading, lwd=1, col=col.1)
  legend(.4, .3, 
         paste("\u2014 AUC: ", round(roc1$auc,2)), 
         text.font=2, text.col=col.1,
         bty="n", cex = leg.cex, xjust=0, yjust=.5)

  par(new=TRUE)
  
  # Plot 2
  plot(roc2$roc, colorize=F, lwd=1, col=col.2, lty=2)
  legend(.4, .2, 
         paste("- - AUC: ", round(roc2$auc,2)), 
         text.font=2, text.col=col.2,
         bty="n", cex = leg.cex, xjust=0, yjust=.5)
  
  # Add diagonal line
  abline(a=0, b=1, lty=3)
}

#################################################################################
# Multiple plots: all plots 3x4
jpeg("roc_4x4_plot.jpg", units="in", width=6, height=6.5, res=600)
par(mfrow=c(4,4),
    cex.lab=.5, cex.main=.7, cex.axis=.6,
    xaxt='s', yaxt='s',
    oma=c(2,2,0,0), mar=c(1,1,1,0.5), mgp=c(3, 0.2, 0))

ROCplot(n1.roc, n2.roc, "Burst rate", leg.cex=.7)
ROCplot(i1mean.roc, i2mean.roc, "Inter-burst interval (mean)", leg.cex=.7)
ROCplot(l1mean.roc, l2mean.roc, "Burst duration (mean)", leg.cex=.7)
ROCplot(a1mean.roc, a2mean.roc, "Burst peak amplitude (mean)", leg.cex=.7)

ROCplot(r1.roc, r2.roc, "Relative beta power", leg.cex=.7)
ROCplot(i1medi.roc, i2medi.roc, "Inter-burst interval (median)", leg.cex=.7)
ROCplot(l1medi.roc, l2medi.roc, "Burst duration (median)", leg.cex=.7)
ROCplot(a1medi.roc, a2medi.roc, "Burst peak amplitude (median)", leg.cex=.7)

ROCplot(b1.roc, b2.roc, "Beta peak power", leg.cex=.7)
ROCplot(i1mode.roc, i2mode.roc, "Inter-burst interval (mode)", leg.cex=.7)
ROCplot(l1mode.roc, l2mode.roc, "Burst duration (mode)", leg.cex=.7)
ROCplot(a1mode.roc, a2mode.roc, "Burst peak amplitude (mode)", leg.cex=.7)

ROCplot(fi1.roc, fi2.roc, "1/f intercept", leg.cex=.7)
ROCplot(fs1.roc, fs2.roc, "1/f slope", leg.cex=.7)

mtext(text="Sensitivity (TPR)", side=2, line=.5, outer=TRUE, font=2, cex=.8)
mtext(text="1-specificity (FPR)", side=1, line=1, outer=TRUE, font=2, cex=.8)
dev.off()

#################################################################################
##### Plot without 3xIBI, Dur and Amp.
jpeg("roc_4x2_plot.jpg", units="in", width=4, height=6.5, res=600)

par(oma=c(3,3,0,0),mar=c(2,2,2,2),mfrow=c(4,2),
    cex.lab=.6, cex.main=1.2, cex.axis=.8,
    oma=c(3,3,0,0),mar=c(2,1,2,2))

ROCplot(n1.roc, n2.roc, "Burst rate", leg.cex=.85)
ROCplot(r1.roc, r2.roc, "Relative beta power", leg.cex=.85)
ROCplot(i1medi.roc, i2medi.roc, "Inter-burst interval", leg.cex=.85)
ROCplot(b1.roc, b2.roc, "Beta peak power", leg.cex=.85)
ROCplot(l1medi.roc, l2medi.roc, "Burst duration", leg.cex=.85)
ROCplot(fi1.roc, fi2.roc, "1/f intercept", leg.cex=.85)
ROCplot(a1medi.roc, a2medi.roc, "Burst peak amplitude", leg.cex=.85)
ROCplot(fs1.roc, fs2.roc, "1/f slope", leg.cex=.85)

mtext(text="Sensitivity (TPR)", side=2, line=1.5, outer=TRUE, font=2)
mtext(text="1-specificity (FPR)", side=1, line=1.5, outer=TRUE, font=2)

dev.off()

#END