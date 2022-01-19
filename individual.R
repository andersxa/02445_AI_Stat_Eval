library(lawstat)
data = read.table("WR_2021.txt",header=TRUE)
data$sex = factor(data$sex)

pdf("DATA_VIS.pdf",12,4)
par(mfrow=c(1,3))
interaction.plot(data$sex, data$distance, data$time, fun=identity, type="p", xlab="Gender", ylab="Record Time", trace.label="Distance")
plot(data$distance,data$time,col=data$sex, xlab="Distance", ylab="Record Time")
legend("topleft",legend=c("men","women"), col=1:2, pch=1)
plot(log(data$distance),log(data$time),col=data$sex, xlab="Log Distance", ylab="Log Record Time")
legend("topleft",legend=c("men","women"), col=1:2, pch=1)
dev.off()

pdf("INTERACTION_PLOT.pdf",12,6)
par(mfrow=c(1,2))
interaction.plot(data$distance, data$sex, data$time, fun=identity, type="l", xlab="Distance", ylab="Record Time", trace.label="Gender")
interaction.plot(log(data$distance), data$sex, log(data$time), fun=identity, type="l", xlab="Log Distance", ylab="Log Record Time", trace.label="Gender")
dev.off()

levene.test(log(data$time), data$sex, location="mean")
lapply(split(data, data$sex), function(x) shapiro.test(log(x$time)))

levene.test(data$time, data$sex, location="mean")
lapply(split(data, data$sex), function(x) shapiro.test(x$time))

#OVERALL NORMALITY
shapiro.test(log(data$time))
#TEST FOR HOMOGENUITY - I.E. EQUAL VARIANCE IN GROUPS

L = lm(log(time) ~ sex*log(distance), data=data)
L2 = lm(log(time) ~ sex + log(distance), data=data)
ks.test(residuals(L),"pnorm",sd=sd(residuals(L)))
ks.test(residuals(L2),"pnorm",sd=sd(residuals(L2)))

anova(L)

pdf("FITTED_MODELS.pdf",12,12)
par(mfrow=c(1,2))

plot(log(data$distance),log(data$time),col=data$sex,xlab="Log Distance", ylab="Log Record Time")
lines(log(data$distance)[data$sex == "men"],fitted(L)[data$sex == "men"],col=3,lty=1)
lines(log(data$distance)[data$sex == "women"],fitted(L)[data$sex == "women"],col=3,lty=2)
legend("topleft",legend=c("men","women"),lty=1:2,col=3,title="model")
legend("bottomright",legend=c("men","women"),pch=1,col=1:2,title="data")

plot(log(data$distance),log(data$time),col=data$sex,xlab="Log Distance", ylab="Log Record Time")
lines(log(data$distance)[data$sex == "men"],fitted(L2)[data$sex == "men"],col=4,lty=1)
lines(log(data$distance)[data$sex == "women"],fitted(L2)[data$sex == "women"],col=4,lty=2)
legend("topleft",legend=c("men","women"),lty=1:2,col=4,title="model")
legend("bottomright",legend=c("men","women"),pch=1,col=1:2,title="data")

plot(data$distance,data$time,col=data$sex,xlab="Distance", ylab="Record Time")
lines(data$distance[data$sex == "men"],exp(fitted(L))[data$sex == "men"],col=3,lty=1)
lines(data$distance[data$sex == "women"],exp(fitted(L))[data$sex == "women"],col=3,lty=2)
legend("topleft",legend=c("men","women"),lty=1:2,col=3,title="model")
legend("bottomright",legend=c("men","women"),pch=1,col=1:2,title="data")

plot(data$distance,data$time,col=data$sex,xlab="Distance", ylab="Record Time")
lines(data$distance[data$sex == "men"],exp(fitted(L2))[data$sex == "men"],col=4,lty=1)
lines(data$distance[data$sex == "women"],exp(fitted(L2))[data$sex == "women"],col=4,lty=2)
legend("topleft",legend=c("men","women"),lty=1:2,col=4,title="model")
legend("bottomright",legend=c("men","women"),pch=1,col=1:2,title="data")
dev.off()









