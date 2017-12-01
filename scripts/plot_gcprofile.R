#### Human-Mouse

hg1 <- read.table("Human_Mouse/hg19_profile/I-II_ref_HSBD.bed.max500bp.aligned.GCprofile.1bp", header=T)
hg2 <- read.table("Human_Mouse/hg19_profile/II-III_ref_MSBD.bed.max500bp.aligned.GCprofile.1bp", header=T)
hg3 <- read.table("Human_Mouse/hg19_profile/III-IV_ref_HEBD.bed.max500bp.aligned.GCprofile.1bp", header=T)
hg1$Group <- 1 ; hg1$Group[hg1$Dist>0] <- 2
hg2$Group <- 3 ; hg2$Group[hg2$Dist>0] <- 4
hg3$Group <- 5 ; hg3$Group[hg3$Dist>0] <- 6
hg1$Dist <- hg1$Dist - 1000
hg3$Dist <- hg3$Dist + 1000

mm1 <- read.table("Human_Mouse/mm10_profile/I-II_ref_HSBD.bed.max500bp.aligned.GCprofile.1bp", header=T)
mm2 <- read.table("Human_Mouse/mm10_profile/II-III_ref_MSBD.bed.max500bp.aligned.GCprofile.1bp", header=T)
mm3 <- read.table("Human_Mouse/mm10_profile/III-IV_ref_HEBD.bed.max500bp.aligned.GCprofile.1bp", header=T)
mm1$Group <- 1 ; mm1$Group[mm1$Dist>0] <- 2
mm2$Group <- 3 ; mm2$Group[mm2$Dist>0] <- 4
mm3$Group <- 5 ; mm3$Group[mm3$Dist>0] <- 6
mm1$Dist <- mm1$Dist - 1000
mm3$Dist <- mm3$Dist + 1000

hg <- rbind(hg1, hg2, hg3); hg$Dist <- hg$Dist/100
mm <- rbind(mm1, mm2, mm3); mm$Dist <- mm$Dist/100
hg <- hg[hg$Total>200, ]
mm <- mm[mm$Total>200, ]


pdf("Human_Mouse_gcprofile.pdf", width=2.3, height=3.7, pointsize=8)
par(mfrow=c(2,1), mar=c(3,3,1,1),las=1)
plot(hg$Dist, hg$CpGoe, ylim=c(0,1), pch=".", cex=0.1, col="#127d09")
for (i in 1:6) {
  points(loess.smooth(x = hg$Dist[hg$Group==i], y = hg$CpGoe[hg$Group==i], span = 0.1), type="l", lwd=1)
}
abline(v=c(-10,0,10), lty="dashed")
plot(mm$Dist, mm$CpGoe, ylim=c(0,1), pch=".", cex=0.1, col="#127d09")
for (i in 1:6) {
  points(loess.smooth(x=mm$Dist[mm$Group==i], y=mm$CpGoe[mm$Group==i], span = 0.1), type="l", lwd=1)
}
abline(v=c(-10,0,10), lty="dashed")
dev.off()


#### Dog-Mouse

hg1 <- read.table("Dog_Mouse/canFam3_profile/I-II_ref_HSBD.bed.max500bp.aligned.GCprofile.1bp", header=T)
hg2 <- read.table("Dog_Mouse/canFam3_profile/II-III_ref_MSBD.bed.max500bp.aligned.GCprofile.1bp", header=T)
hg3 <- read.table("Dog_Mouse/canFam3_profile/III-IV_ref_HEBD.bed.max500bp.aligned.GCprofile.1bp", header=T)
hg1$Group <- 1 ; hg1$Group[hg1$Dist>0] <- 2
hg2$Group <- 3 ; hg2$Group[hg2$Dist>0] <- 4
hg3$Group <- 5 ; hg3$Group[hg3$Dist>0] <- 6
hg1$Dist <- hg1$Dist - 1000
hg3$Dist <- hg3$Dist + 1000

mm1 <- read.table("Dog_Mouse/mm10_profile/I-II_ref_HSBD.bed.max500bp.aligned.GCprofile.1bp", header=T)
mm2 <- read.table("Dog_Mouse/mm10_profile/II-III_ref_MSBD.bed.max500bp.aligned.GCprofile.1bp", header=T)
mm3 <- read.table("Dog_Mouse/mm10_profile/III-IV_ref_HEBD.bed.max500bp.aligned.GCprofile.1bp", header=T)
mm1$Group <- 1 ; mm1$Group[mm1$Dist>0] <- 2
mm2$Group <- 3 ; mm2$Group[mm2$Dist>0] <- 4
mm3$Group <- 5 ; mm3$Group[mm3$Dist>0] <- 6
mm1$Dist <- mm1$Dist - 1000
mm3$Dist <- mm3$Dist + 1000

hg <- rbind(hg1, hg2, hg3); hg$Dist <- hg$Dist/100
mm <- rbind(mm1, mm2, mm3); mm$Dist <- mm$Dist/100
hg <- hg[hg$Total>100, ]
mm <- mm[mm$Total>100, ]
pdf("Dog_Mouse_gcprofile.pdf", width=2.3, height=3.7, pointsize=8)
par(mfrow=c(2,1), mar=c(3,3,1,1),las=1)
plot(hg$Dist, hg$CpGoe, ylim=c(0,1.2), pch=".", cex=0.1, col="#CB4D42")
for (i in 1:6) {
  points(loess.smooth(x = hg$Dist[hg$Group==i], y = hg$CpGoe[hg$Group==i], span = 0.2), type="l", lwd=1)
}
abline(v=c(-10,0,10), lty="dashed")
plot(mm$Dist, mm$CpGoe, ylim=c(0,1.2), pch=".", cex=0.1, col="#00989D")
for (i in 1:6) {
  points(loess.smooth(x=mm$Dist[mm$Group==i], y=mm$CpGoe[mm$Group==i], span = 0.2), type="l", lwd=1)
}
abline(v=c(-10,0,10), lty="dashed")
dev.off()

