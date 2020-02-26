library(gamlss)
library(data.table)
library(ggplot2)

data <- fread('~/height_prediction/output/WHI_JHS_UKB_HRS_summ.txt')

data <- data[!is.na(data$HEIGHTX),]

model<-lm(HEIGHTX~sex+age+age2+Dataset, data=data)

data$RES <- model$residuals

gamlss.mod1 <- gamlss(RES~EUR_ANC, family=NO2, data=data)
gamlss.mod2 <- gamlss(RES~EUR_ANC, sigma.formula=~EUR_ANC, family=NO2(sigma.link=identity), data=data)

pdf("~/height_prediction/figs_for_paper/figs/Variance_model1.pdf")
plot(data$EUR_ANC, data$RES, pch=16, cex=0.3, xlab="European Ancestry", ylab="Residual Height", bty="n" )

c1 <-  gamlss.mod2$mu.coefficients
c2 <- gamlss.mod2$sigma.coefficients
pd <- c(c2[1], -0.24*c2[1])

## abline(c1[1], c1[2], col="red", lwd=3)
x <- seq(0, 1, l=101)
sd.diff <- sqrt(c2[1]+x*c2[2])
pd.diff <- sqrt(pd[1]+x*pd[2])

mn <- c1[1] + x*c1[2]
lines(x, mn, col="red", lwd=2)


legend("topright", c("Mean"), col=c("red" ), lwd=2, bty="n")
dev.off()

pdf("~/height_prediction/figs_for_paper/figs/Variance_model2.pdf")
plot(data$EUR_ANC, data$RES, pch=16, cex=0.3, xlab="European Ancestry", ylab="Residual Height", bty="n" )

c1 <-  gamlss.mod2$mu.coefficients
c2 <- gamlss.mod2$sigma.coefficients
pd <- c(c2[1], -0.24*c2[1])

## abline(c1[1], c1[2], col="red", lwd=3)
x <- seq(0, 1, l=101)
sd.diff <- sqrt(c2[1]+x*c2[2])
pd.diff <- sqrt(pd[1]+x*pd[2])

mn <- c1[1] + x*c1[2]
lines(x, mn, col="red", lwd=2)
lines(x, mn+pd.diff, col="blue", lwd=2)
lines(x, mn-pd.diff, col="blue", lwd=2)


legend("topright", c("Mean",  "Expected +/- 1 SD"), col=c("red", "blue"), lwd=2, bty="n")
dev.off()

pdf("~/height_prediction/figs_for_paper/figs/Variance_model3.pdf")
plot(data$EUR_ANC, data$RES, pch=16, cex=0.3, xlab="European Ancestry", ylab="Residual Height", bty="n" )

c1 <-  gamlss.mod2$mu.coefficients
c2 <- gamlss.mod2$sigma.coefficients
pd <- c(c2[1], -0.24*c2[1])

## abline(c1[1], c1[2], col="red", lwd=3)
x <- seq(0, 1, l=101)
sd.diff <- sqrt(c2[1]+x*c2[2])
pd.diff <- sqrt(pd[1]+x*pd[2])

mn <- c1[1] + x*c1[2]
lines(x, mn, col="red", lwd=2)
lines(x, mn+sd.diff, col="orange", lwd=2)
lines(x, mn-sd.diff, col="orange", lwd=2)
lines(x, mn+pd.diff, col="blue", lwd=2)
lines(x, mn-pd.diff, col="blue", lwd=2)

legend("topright", c("Mean height", "Expected +/- 1 SD",  "Fitted +/- 1 SD"), col=c("red", "blue", "orange"), lwd=2, bty="n")
dev.off()

