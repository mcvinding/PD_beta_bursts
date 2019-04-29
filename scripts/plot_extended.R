# Plot extended analysis
library(ggplot2)

wrkdir <- "Z://PD_motor//rest_ec//groupanalysis//"
setwd(wrkdir)
load(file='range_bf.RData')

bf.dat <- data.frame("steps"=seq(0.1,5,by=0.1),
                    "test"=rep(c('bf10','bf21','bf32'), each=50),
                    "BF"=c(BF10,BF21,BF32))
bf.dat$logBF <- log(bf.dat$BF)

ggplot(bf.dat, aes(x=steps, y=logBF, fill=test))+
  geom_line()+
  geom_point(shape=21, size=2)+theme_bw()

# + horisontal line or shaded area to indicate n.s.
# + vertical line to indicate max value