files <- c("Human_Sperm_HMR.since_Boreoeutheria",
"Human_Sperm_HMR.since_Euarchontoglire",
"Human_Sperm_HMR.since_Catarrhini",
"Human_Sperm_HMR.since_Homininae",
"Human_Sperm_HMR.since_Hominini",
"Human_Sperm_HMR.since_Human")

mat <- data.frame(region=character(0), context = character(0), median=numeric(0), l95ci=numeric(0), h95ci=numeric(0))
for (f in files){
  for (i in c("promoter", "nonpromoter")) {
    x <- read.table(paste(f,i, sep="."))
    x <- x[x$V5>=3, ]
    m <- median(x$V3-x$V2)
    n <- nrow(x)
    l <- as.numeric(quantile(x$V3-x$V2, 1/2 -1.96/2/sqrt(n), na.rm=T))
    h <- as.numeric(quantile(x$V3-x$V2, 1/n + 1/2 + 1.96/2/sqrt(n), na.rm=T))
    mat <- rbind(mat, data.frame(region=f, context = i, median=m, l95ci=l, h95ci=h))
  }
}
write.table(mat, file="Human_median_stat.txt", col.names=T, row.names=F, sep="\t")

mat <- data.frame(region=character(0), context = character(0), mean=numeric(0), l95ci=numeric(0), h95ci=numeric(0))
for (f in files){
  for (i in c("promoter", "nonpromoter")) {
    x <- read.table(paste(f,i, sep="."))
    x <- x[x$V5>=3, ]
    m <- mean(x$V3-x$V2)
    n <- nrow(x)
    l <- m-sd(x$V3-x$V2)/sqrt(n)*1.96
    h <- m+sd(x$V3-x$V2)/sqrt(n)*1.96
    mat <- rbind(mat, data.frame(region=f, context = i, mean=m, l95ci=l, h95ci=h))
  }
}
write.table(mat, file="Human_mean_stat.txt", col.names=T, row.names=F, sep="\t")



files <- c("Mouse_Sperm_HMR.since_Boreoeutheria",
"Mouse_Sperm_HMR.since_Euarchontoglire",
"Mouse_Sperm_HMR.since_Murinae",
"Mouse_Sperm_HMR.since_Mouse")
mat <- data.frame(region=character(0), context = character(0), median=numeric(0), l95ci=numeric(0), h95ci=numeric(0))
for (f in files){
  for (i in c("promoter", "nonpromoter")) {
    x <- read.table(paste(f,i, sep="."))
    x <- x[x$V5>=3, ]
    m <- median(x$V3-x$V2)
    n <- nrow(x)
    l <- as.numeric(quantile(x$V3-x$V2, 1/2 -1.96/2/sqrt(n), na.rm=T))
    h <- as.numeric(quantile(x$V3-x$V2, 1/n + 1/2 + 1.96/2/sqrt(n), na.rm=T))
    mat <- rbind(mat, data.frame(region=f, context = i, median=m, l95ci=l, h95ci=h))
  }
}
write.table(mat, file="Mouse_median_stat.txt", col.names=T, row.names=F, sep="\t")

mat <- data.frame(region=character(0), context = character(0), mean=numeric(0), l95ci=numeric(0), h95ci=numeric(0))
for (f in files){
  for (i in c("promoter", "nonpromoter")) {
    x <- read.table(paste(f,i, sep="."))
    x <- x[x$V5>=3, ]
    m <- mean(x$V3-x$V2)
    n <- nrow(x)
    l <- m-sd(x$V3-x$V2)/sqrt(n)*1.96
    h <- m+sd(x$V3-x$V2)/sqrt(n)*1.96
    mat <- rbind(mat, data.frame(region=f, context = i, mean=m, l95ci=l, h95ci=h))
  }
}
write.table(mat, file="Mouse_mean_stat.txt", col.names=T, row.names=F, sep="\t")


library(ggplot2)
x <- read.table("Human_mean_stat.txt", sep="\t", header=T, as.is=T)
y <- read.table("Mouse_mean_stat.txt", sep="\t", header=T, as.is=T)
x <- rbind(x,y)
x$region <- factor(x$region, levels= c(paste("Human_Sperm_HMR.since_",
c("Boreoeutheria", "Euarchontoglire", "Catarrhini", "Homininae", "Hominini", "Human"), sep=""),
paste("Mouse_Sperm_HMR.since_", c("Boreoeutheria","Euarchontoglire","Murinae","Mouse"), sep="")) )
x$context <- factor(x$context, levels = c("promoter", "nonpromoter"))
p <- ggplot(x, aes(x=region, y=mean, fill=context)) +
      geom_bar(position=position_dodge(), stat="identity",
               colour="black", # Use black outlines,
               size=.3) +      # Thinner lines
      geom_errorbar(aes(ymin=l95ci, ymax=h95ci),
                    size=.3,    # Thinner lines
                    width=.2,
                    position=position_dodge(.9)) +
      xlab("Age") +
      ylab("Mean HMR size") +
      scale_fill_manual(values=c("#0071bc", "#39b54a"),
                        name="Context", # Legend label, use darker colors
                        breaks=c("promoter", "nonpromoter"),
                        labels=c("Promoter", "Nonpromoter")) +
      ggtitle("HMR size by hypomethylation age group") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 30, hjust=1, vjust=1, size=8 ),
            axis.title.x=element_blank() ,
            axis.text.y = element_text(size=8 ),
            axis.title.y = element_text(size=10),
            axis.line = element_line(colour = "black"),
            plot.title = element_text(size=10),
            legend.justification=c(0,0), legend.position=c(0,0.3),
            legend.title = element_text(size=8),
            legend.text = element_text(size=8),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black"),
            panel.background = element_blank())

pdf("HMR_size_age.pdf", width=6, height=3.5)
p
dev.off() 
