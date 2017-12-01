
ProjDir="/home/cmb-panasas2/jqu/sperm_methylome_evolution/"

p2sig <- function(p){
  sig <- "  "
  if (p< 0.05) sig <- "*";
  if (p< 0.01) sig <- paste(rep("*", floor(min(-log(p)/log(10), 3))), collapse="")
  return(sig)
}

types  <- c("Boreoeutheria", "Euarchontoglires", "Catarrhini", "Homininae", "Hominini", "Human")
x <- list()
m <- list()
for (i in 1:6){
  y <- read.table(paste("Human_HYPO_since_",types[i],".masked.cpgobex.filtered", sep=""))
  z <- read.table(paste("Human_HYPO_since_",types[i],".human_sperm_roi.filtered", sep=""))
  x[[i]] <- y$V5 
  m[[i]] <- z$V5 
} # x & m [1-5]:all to human

pdf("HMRage_CpGdensity_bymethgroup.pdf", width=6.6, height=3, pointsize=8)
LC <- c(0, 0.05)
HC <- c(0.05, 0.1)
par(mfrow=c(length(LC),2), mar=c(2,5,1,1))
for (part in 1:length(LC)){
  lc <- LC[part]
  hc <- HC[part]
  lowmethx <- list()
  lowmethm <- list()
  for(i in 1:5) {
    ind <- which(m[[i]] >= lc & m[[i]] < hc )
    lowmethx[[6-i]] <- x[[i]][ind] # lowmethx[1-5]: human to all
    lowmethm[[6-i]] <- m[[i]][ind]
  }
  names(lowmethx) <- rev(types)
  names(lowmethm) <- rev(types)
  
  boxplot(lowmethm, las=1, main="Methylation level (human)", outline=F, horizontal = T)
  for (j in 1:4) {
    s <- wilcox.test(lowmethm[[j]], lowmethm[[j+1]], alt="less")$p.val
    sg <- wilcox.test(lowmethm[[j]], lowmethm[[j+1]], alt="greater")$p.val
    text(y=j+0.3, x=hc-0.1*(hc-lc), pos = 4, labels=p2sig(s), col="blue")
    text(y=j+0.3, x=hc-0.1*(hc-lc), pos = 4, labels=p2sig(sg), col="red")
  } 
  
  boxplot(lowmethx, las=1, main="CpG enrichment", outline=F, horizontal = T, ylim=c(0,1.1))
  for (j in 1:4) {
    p <- wilcox.test(lowmethx[[j]], lowmethx[[j+1]], alt="less")$p.val
    pg <- wilcox.test(lowmethx[[j]], lowmethx[[j+1]], alt="greater")$p.val
    text(y=j+0.3, x=0.9, pos = 4, labels=p2sig(p), col="red")
    text(y=j+0.3, x=0.9, pos = 4, labels=p2sig(pg), col="blue")
  }
}
dev.off()

# violinplot
library(ggplot2)
types  <- c("Boreoeutheria", "Euarchontoglires", "Catarrhini", "Homininae", "Hominini", "Human")
for (i in 1:6){
  y <- read.table(paste("Human_HYPO_since_",types[i],".masked.cpgobex.filtered", sep=""))
  z <- read.table(paste("Human_HYPO_since_",types[i],".human_sperm_roi.filtered", sep=""))
  indy <- which(as.numeric(matrix(unlist(strsplit(as.vector(y$V4), ":")), byrow=T, ncol=4)[,2])>0)
  indz <- which(as.numeric(matrix(unlist(strsplit(as.vector(z$V4), ":")), byrow=T, ncol=5)[,2])>0)
  if (i==1) {
    x <- data.frame(cpg_density=y$V5[indy],group=rep(types[i], length(indy)), stringsAsFactors=F)
    w <- data.frame(meth=z$V5[indz],group=rep(types[i], length(indz)), stringsAsFactors=F)
    
  } else {
    x <- rbind(x, data.frame(cpg_density=y$V5[indy], group=rep(types[i], length(indy)), stringsAsFactors=F))
    w <- rbind(w, data.frame(meth=z$V5[indz], group=rep(types[i], length(indz)), stringsAsFactors=F))
  }
}
x$group <- factor(x$group, levels = types)
w$group <- factor(w$group, levels = types)


bp_density <- ggplot(x, aes(x=group, y=cpg_density, fill=group)) +
geom_violin(trim=TRUE, scale="width", size=0.37)+
geom_boxplot(width=0.15, fill="white", outlier.shape=NA, size=0.37)+
coord_cartesian(ylim = c(0, 1.3)) +
labs(title="",x="", y = "CpG o/e ratio") +
scale_fill_brewer(palette="Blues", name="HMR age\n(Myr)",
breaks=levels(x$group),
labels=c(">90", "30-90", "9-30", "7-9", "<7")) +
theme_classic()+ theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
theme(panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
panel.border = element_rect(color="black", size=0.37),
axis.line = element_line(colour = "black", size=0.37),
axis.ticks = element_line(size=0.37),
line=element_line(size=0.37),
panel.background = element_blank(),
legend.key.size = unit(0.12, "in"),
legend.key = element_blank(),
legend.title=element_text(size=8),
legend.text=element_text(size=8),
axis.text=element_text(size=8),
axis.title=element_text(size=10))

wilcox.test(x$cpg_density[x$group=="Root"], x$cpg_density[x$group=="Primate"], alt="greater" )$p.val #[1] 9.12629e-82
wilcox.test(x$cpg_density[x$group=="Primate"], x$cpg_density[x$group=="HCG"], alt="greater" )$p.val #[1] 1.357827e-47
wilcox.test(x$cpg_density[x$group=="HCG"], x$cpg_density[x$group=="HC"], alt="greater" )$p.val  #[1] 0.7208867
wilcox.test(x$cpg_density[x$group=="HC"], x$cpg_density[x$group=="Human"], alt="greater" )$p.val #[1] 2.177688e-15
wilcox.test(x$cpg_density[x$group=="HCG"], x$cpg_density[x$group=="Human"], alt="greater" )$p.val #[1] 5.924737e-22

bp_meth <- ggplot(w, aes(x=group, y=meth, fill=group)) +
geom_violin(trim=TRUE, scale="width", size=0.37)+
geom_boxplot(width=0.15, fill="white", outlier.shape=NA, size=0.37)+
coord_cartesian(ylim = c(0, 0.5)) +
labs(title="", x="", y = "Methylation level") +
scale_fill_brewer(palette="Blues", name="HMR age\n(Myr)",
breaks=levels(x$group),
labels=c(">90", "30-90", "9-30", "7-9", "<7")) +
theme_classic()+ theme_bw() +
theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
theme(panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
panel.border = element_rect(color="black", size=0.37),
axis.line = element_line(colour = "black", size=0.37),
axis.ticks = element_line(size=0.37),
line=element_line(size=0.37),
panel.background = element_blank(),
legend.key.size = unit(0.12, "in"),
legend.key = element_blank(),
legend.title=element_text(size=8),
legend.text=element_text(size=8),
axis.text=element_text(size=8),
axis.title=element_text(size=10))


pdf("hypo_history_CpG_density.pdf", width=6, height=2.5, pointsize=8)
source(paste(ProjDir, "scripts/multiplot.R", sep="/")
multiplot(bp_density,bp_meth, cols=2)
dev.off()
 

