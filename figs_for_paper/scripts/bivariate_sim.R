library(mvtnorm)
library(ggplot2)
library(wesanderson)
expected.or <- function(x, rho){
val <- qnorm(1-x)
prob <- pmvnorm(lower=c(val, val), upper=Inf, sigma=matrix(c(1,rho,rho,1),2,2))[1]
OR <- prob/x/x
return(OR)
}

x.vals <- seq(0.01,0.5, l=50)
rho15 <- sapply(x.vals, expected.or, sqrt(0.15))
rho034 <- sapply(x.vals, expected.or, sqrt(0.034))
rho036 <- sapply(x.vals, expected.or, sqrt(0.036))
rho041 <- sapply(x.vals, expected.or, sqrt(0.041))


my_col<-wes_palette("Darjeeling1",4)

pdf('OR_expected.pdf')
plot(1-x.vals, rho15, type="l", xlim=c(0.5,1), col=my_col[1], bty="n", xlab="Quantile", ylab="Bivariate normal OR", lwd=2, cex.lab=1.5)
lines(1-x.vals, rho041, type="l",col=my_col[2], lwd=2)
lines(1-x.vals, rho036, type="l",col=my_col[3], lwd=2)
lines(1-x.vals, rho034, type="l",col=my_col[4], lwd=2)
legend("topleft", c(expression(paste(R^2, "=0.15",sep="")), expression(paste(R^2,"=0.041", sep="")), expression(paste(R^2,"=0.036", sep="")), expression(paste(R^2,"=0.034", sep=""))), col=my_col, bty="n", lty=1, lwd=2)
dev.off()
