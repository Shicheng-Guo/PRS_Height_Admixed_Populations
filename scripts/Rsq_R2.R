
###########
rsq.R2 <- function(formula1, formula2, data, indices) {
  d <- data[indices,] # allows boot to select sample
  fit <- lm(formula1, data=d)
  fit2<- lm(formula2, data=d)
  res<-partial.R2(fit, fit2)
  return(res)
}
##############
