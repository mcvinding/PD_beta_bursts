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

#END