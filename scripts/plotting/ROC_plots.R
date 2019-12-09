# ROC plots
library(ROCR)



# ROC plot function
ROCplot <- function(roc, auc, heading){
  
  jpeg("test_plot.jpg", units="in", width=3, height=3, res=600)
  par(cex.lab=.6, cex.main=.8, cex.axis=.5,
      oma=c(3,3,0,0),mar=c(1,1,1,1)# ps = 12,
  )
  plot(roc, colorize=F, 
       main = heading,
       ylab = "Sensitivity (TPR)",
       xlab = "1-specificity (FPR)",
       lwd=2)
  abline(a=0, b=1, lty=3)
  legend(.65, .3, 
         paste("AUC: ", round(auc,2)), 
         text.font=2,
         bty="n", 
         cex = 1,
         xjust=.5, yjust=.5)
  dev.off()
  
  # Double plot
  plot(n1.roc)
  par(new=TRUE)
  plot(n2.roc, lty=2)
  
  # M;ultiple plots
  par(oma=c(3,3,0,0),mar=c(2,2,2,2),mfrow=c(3,4))
  
  plot(1,1,ylab="",xlab="",type="n"); title('lololol')
  plot(1,1,ylab="",xlab="",type="n")
  plot(1,1,ylab="",xlab="",type="n")
  plot(1,1,ylab="",xlab="",type="n")
  plot(1,1,ylab="",xlab="",type="n"); title('lololol')
  plot(1,1,ylab="",xlab="",type="n")
  plot(1,1,ylab="",xlab="",type="n")
  plot(1,1,ylab="",xlab="",type="n")  
  plot(1,1,ylab="",xlab="",type="n"); title('lololol')
  plot(1,1,ylab="",xlab="",type="n")
  plot(1,1,ylab="",xlab="",type="n")
  plot(1,1,ylab="",xlab="",type="n")  
  
  
  mtext(text="Sensitivity (TPR)", side=1, line=1, outer=TRUE)
  mtext(text="1-specificity (FPR)", side=2, line=1, outer=TRUE)
  
  
}

saveROcplt <- function(ROCplot) {
  jpeg("test_plot.jpg", units="in", width=5, height=5, res=300)
  dev.off()
}