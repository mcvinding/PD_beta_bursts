# Misc functions for ROC/ideal observer analysis
library(ROCR)
library(lme4)

############################################################################
# Get area under ROC curve
cut.val <- function(x, mod){
  val <- (logit(x)-coef(mod)[1])/coef(mod)[2]
  # val <- (log(x/(1-x))-coef(mod)[1])/coef(mod)[2]
  return(val)
}

# Funtion that do the analysis
auroc.fun <- function(r0, r1){
  
  df <- data.frame(x <- c(r0, r1),
                   group <- as.factor(c(rep(0,length(r0)), rep(1,length(r1))))
  )
  
  lmod <- glm(group~x, data=df, family = 'binomial')
  
  pred <- predict(lmod, df, type="response")
  pred <- prediction(pred, df$group)
  
  auc <- performance(pred, "auc")
  auc <- unlist(slot(auc, "y.values"))
  # print(c(AUC=auc))
  
  return(auc)
}

############################################################################
# Ideal observer analysis (cf. Purushothaman & Casagrande, 2013, Front.Psych.). Not used.
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