branches <- read.table("child_parent.txt", as.is=T)$V1[-1]

pdf("Dist2ConservedHMR_OE.pdf", width=3.5, height=3.5, pointsize=8)
par(mfrow=c(3,4), mar=c(4,3,3,1))
for (child in branches) {
x <- read.table(paste("branch_", child, "_gain_dist_conserved.txt", sep=""))$V1+1
xs <- read.table(paste("branch_", child,"_gain_dist_conserved_shuffled.txt", sep=""))$V1+1
xh <- hist(log(x)/log(10), breaks=seq(0, 7, 1), plot=F)
xsh <- hist(log(xs)/log(10), breaks=seq(0, 7, 1), plot=F)
f <- (xh$counts)/(xsh$counts) * sum(xsh$counts)/sum(xh$counts)
xh$counts <- f
plot(xh, ylim=range(0,xh$counts, na.rm=T), ylab="Observed/expected",
    xlab="Distance to conserved HMR (log 10)", main=child, xlim=c(0,7), col="black")
abline(h=1, col="red", lty="dashed")
}

par(mfrow=c(3,4), mar=c(4,3,3,1))
for (child in branches) {
x <- read.table(paste("branch_", child, "_loss_dist_conserved.txt", sep=""))$V1+1
xs <- read.table(paste("branch_", child,"_loss_dist_conserved_shuffled.txt", sep=""))$V1+1
xh <- hist(log(x)/log(10), breaks=seq(0, 7, 1), plot=F)
xsh <- hist(log(xs)/log(10), breaks=seq(0, 7, 1), plot=F)
f <- (xh$counts)/(xsh$counts) * (sum(xsh$counts)/sum(xh$counts))
xh$counts <- f
plot(xh, ylim=range(0, xh$counts, na.rm=T), ylab="Observed/expected",
    xlab="Distance to conserved HMR (log 10)", main=child, xlim=c(0,7), col="black")
abline(h=1, col="red", lty="dashed")
}
dev.off()

