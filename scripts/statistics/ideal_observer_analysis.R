# Ideal observer analysis (cf. Purushothaman & Casagrande, 2013, Front.Psych.).

# Make function 
id.obs <- function(r0, r1){
  sig1 <- rep(0,length(r1))
  for (i in 1:length(r1)) {
    sig1[i] <- sum(r0 <= r1[i])
  }
  sig2 <- sum(sig1)
  
  Phat <- sig2/(length(r1)*length(r0))
  if (Phat < 0.5){
    Phat <- 1-Phat
  }
  return(Phat)
}

# Define distributions
r0 <- rnorm(50, 20, 5)
r1 <- rnorm(50, 10, 5)

r0 <- neve.data$nevent[neve.data$group=="ptns" & neve.data$session=="1"]
r1 <- neve.data$nevent[neve.data$group=="ctrl" & neve.data$session=="1"]



id.obs(r0, r1)


# Plot distribution
library(ggplot2)
d1 <- data.frame(val=r0)
d2 <- data.frame(val=r1)
d1$var <- 'r1'
d2$var <- 'r2'
dat <- rbind(d1, d2)
ggplot(dat, aes(val, fill = var)) + geom_histogram(bins=20, alpha=.8, position = 'identity')
