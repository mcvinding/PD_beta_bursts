# Make plots of UPDRS ~ MEG summary statistics based on brms models
library(ggplot2)
library(brms)

## Load data
wrkdir <- "Z://PD_motor//rest_ec//groupanalysis//"
outdir <- "Z://PD_motor//rest_ec//figures//"
setwd(wrkdir)
load(file = 'updrs_neve_mods.R')

############################################################################
# Make a generic plot function
pred.plot <- function(bms.mod){
  data <- bms.mod$data
  hdr <- names(data)
  F.name <- names(data)[1]
  plt <- plot(marginal_effects(bms.mod), line_args=list(color="black"))
  plt <- plt$nevent +
    geom_point(data=u.neve.data.PD, aes_string(y=F.name, x="nevent", color="session"), size=1, inherit.aes = FALSE) +
    theme_bw() +
    scale_color_manual(values=c('red','blue'))+
    xlab('N events') + ylab('Factor score')+
    theme_bw() + 
    theme(legend.position = "none",
          text = element_text(size = 11),
          title = element_text(size = 12, vjust = 1.5, face="bold",lineheight = NULL),
          # axis.text = element_text(size=10),
          axis.title = element_text(size = 11, vjust = .5, face="plain"),
          # axis.title.x = element_text(face="plain", size=4),
          # axis.title.y = element_text(face="plain", size=4),
          axis.text = element_text(face="bold", size=11),
          panel.grid = element_blank())
  return(plt)
}

############################################################################
# Plot N event
plt.f1 <- pred.plot(br.nev.uF1) + ggtitle('Midline function')
plt.f2 <- pred.plot(br.nev.uF2) + ggtitle('Rest tremor')
plt.f3 <- pred.plot(br.nev.uF3) + ggtitle('Rigidity')
plt.f4 <- pred.plot(br.nev.uF4) + ggtitle('Bradykinesia (right side)')
plt.f5 <- pred.plot(br.nev.uF5) + ggtitle('Bradykinesia (left side)')
plt.f6 <- pred.plot(br.nev.uF6) + ggtitle('Postural and kinetic tremor')
plt.f7 <- pred.plot(br.nev.uF7) + ggtitle('Lower limb bradykinesia')
# plt.fT <- pred.plot(br.nev.uFT) + ggtitle('Total UPDRS-III')

# Save
setwd(outdir)
ggsave("new_neve_F1.jpeg", plot=plt.f1, device="jpeg", units="mm", width=40, height=35, dpi=500, scale=2.5)
ggsave("new_neve_F2.jpeg", plot=plt.f2, device="jpeg", units="mm", width=40, height=35, dpi=500, scale=2.5)
ggsave("new_neve_F3.jpeg", plot=plt.f3, device="jpeg", units="mm", width=40, height=35, dpi=500, scale=2.5)
ggsave("new_neve_F4.jpeg", plot=plt.f4, device="jpeg", units="mm", width=40, height=35, dpi=500, scale=2.5)
ggsave("new_neve_F5.jpeg", plot=plt.f5, device="jpeg", units="mm", width=40, height=35, dpi=500, scale=2.5)
ggsave("new_neve_F6.jpeg", plot=plt.f6, device="jpeg", units="mm", width=40, height=35, dpi=500, scale=2.5)
ggsave("new_neve_F7.jpeg", plot=plt.f7, device="jpeg", units="mm", width=40, height=35, dpi=500, scale=2.5)
# ggsave("new_neve_Total.jpeg", plot=plt.fT, device="jpeg", units="mm", width=40, height=35, dpi=500, scale=2.5)

#END