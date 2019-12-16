# Misc. functions for PD-BB proj.

## MODE
Mode <- function(x, digits=0) {
  if  (!is.na(digits)){
    x <- round(x, digits=digits)
  }
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
