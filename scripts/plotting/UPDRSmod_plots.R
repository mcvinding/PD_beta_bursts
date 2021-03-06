# Make plots of UPDRS ~ MEG summary statistics based on brms models
library(ggplot2)
library(brms)

## Load data
wrkdir <- "C://Users//Mikkel//Documents//betabursts//groupanalysis"
outdir <- "C://Users//Mikkel//Documents//betabursts//Figures//"
setwd(wrkdir)
load(file = 'updrs_neve_mods.R')

############################################################################
# Make a generic plot function (errors based on k-fold)
pred.plot.kfold <- function(bms.mod, K=10){
  # Read names
  data <- bms.mod$data
  hdr <- names(data)
  F.name <- names(data)[1]
  
  # kfold
  kf <- kfold(bms.mod, save_fits=TRUE, K=K)
  
  # Arrange data
  x <- seq(min(data$nevent.min)-1, max(data$nevent.min)+1,0.5)
  yk <- vector("list",length=K)
  for(i in 1:K) {yk[[i]] <- exp(fixef(kf$fits[i,]$fit)[1]+fixef(kf$fits[i,]$fit)[2]*x)}
  df.yk = data.frame(x=rep(x,K), y=unlist(yk))
  df.yk$f = rep(1:K, each=length(x))
  y <- exp(fixef(bms.mod)[1]+fixef(bms.mod)[2]*x)
  df.y <- data.frame(x=x, y=y)
  
  #Plot
  plt <- ggplot(df.y, aes(y=y, x=x)) + 
    geom_line(data=df.yk, aes(y=y, x=x, group=f), alpha=1, col="grey") +
    geom_line(data=df.y, aes(y=y, x=x), alpha=1, col="black", size=1) +
    
    geom_point(data=data, aes_string(y=F.name, x="nevent.min", color="session"), size=1, inherit.aes = FALSE) +
    theme_bw() +
    scale_color_manual(values=c('red','blue')) +
    xlab('Burst/min') + ylab('Factor score') +
    theme_bw() + 
    theme(legend.position = "none",
          text = element_text(size = 12),
          title = element_text(size = 16, vjust = 1.5, face="bold"),
          axis.text = element_text(size=12),
          axis.title = element_text(size = 13, vjust = .5, face="bold"),
          plot.margin = unit(c(0, 15, 0, 0), "pt"))
  
  return(plt)
}  

# Make a generic plot function (errors based on marginal effects)
pred.plot <- function(bms.mod){
  data <- bms.mod$data
  hdr <- names(data)
  F.name <- names(data)[1]
  plt <- plot(marginal_effects(bms.mod), line_args=list(color="black"))
  plt <- plt$nevent +
    geom_point(data=u.neve.data.PD, aes_string(y=F.name, x="nevent.min", color="session"), size=1, inherit.aes = FALSE) +
    theme_bw() +
    scale_color_manual(values=c('red','blue'))+
    xlab('Burst/min') + ylab('Factor score')+
    theme_bw() + 
    theme(legend.position = "none",
          text = element_text(size = 11),
          title = element_text(size = 16, vjust = 1.5, face="bold",lineheight = NULL),
          axis.text = element_text(size=12),
          axis.title = element_text(size = 13, vjust = .5, face="bold"),
          axis.text = element_text(face="bold", size=11),
          panel.grid = element_blank())
  return(plt)
}

############################################################################
# Plot N event
plt.f1 <- pred.plot.kfold(br.nev.uF1) + ggtitle('Midline function')
plt.f2 <- pred.plot.kfold(br.nev.uF2) + ggtitle('Rest tremor')
plt.f3 <- pred.plot.kfold(br.nev.uF3) + ggtitle('Rigidity')
plt.f45 <- pred.plot.kfold(br.nev.uF45) + ggtitle('Bradykinesia')
plt.f6 <- pred.plot.kfold(br.nev.uF6) + ggtitle('Postural and kinetic tremor')
plt.f7 <- pred.plot.kfold(br.nev.uF7) + ggtitle('Lower limb bradykinesia')
plt.fT <- pred.plot.kfold(br.nev.uFT) + ggtitle('Total UPDRS-III') + ylab('Score')

plt.f1 <- plt.f1 + theme(plot.margin= unit(c(0, 15, 0, 0), "pt"))
plt.f3 <- plt.f3 + theme(plot.margin= unit(c(0, 15, 0, 0), "pt"))
plt.f45 <- plt.f45 + theme(plot.margin= unit(c(0, 15, 0, 0), "pt"))
plt.f6 <- plt.f6 + theme(plot.margin= unit(c(0, 15, 0, 0), "pt"))
plt.f7 <- plt.f7 + theme(plot.margin= unit(c(0, 15, 0, 0), "pt"))
plt.fT <- plt.fT + theme(plot.margin= unit(c(0, 15, 0, 0), "pt"))

# Save
setwd(outdir)
ggsave("new_neve_F1.jpeg", plot=plt.f1, device="jpeg", units="mm", width=50, height=45, dpi=600, scale=2)
ggsave("new_neve_F2.jpeg", plot=plt.f2, device="jpeg", units="mm", width=50, height=45, dpi=600, scale=2)
ggsave("new_neve_F3.jpeg", plot=plt.f3, device="jpeg", units="mm", width=50, height=45, dpi=600, scale=2)
ggsave("new_neve_F45.jpeg", plot=plt.f45, device="jpeg", units="mm", width=50, height=45, dpi=600, scale=2)
ggsave("new_neve_F6.jpeg", plot=plt.f6, device="jpeg", units="mm", width=50, height=45, dpi=600, scale=2)
ggsave("new_neve_F7.jpeg", plot=plt.f7, device="jpeg", units="mm", width=50, height=45, dpi=600, scale=2)
ggsave("new_neve_Total.jpeg", plot=plt.fT, device="jpeg", units="mm", width=50, height=45, dpi=600, scale=2)

#END