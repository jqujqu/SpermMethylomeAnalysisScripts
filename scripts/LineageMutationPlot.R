
# Change to proper path
ProjDir="~/sperm_methylome_evolution"

library(ggplot2)
MutDir= paste(ProjDir, "LineageMutation/regions_mutation", sep="/")
s <- "Region Lineage size
Human_Sperm.hmr.specific Human 6
Human_Sperm.hmr.specific Chimp 6
Chimp_Sperm.hmr.specific Human 6
Chimp_Sperm.hmr.specific Chimp 6
Gorilla_Sperm.hmr.specific Gorilla 4
Gorilla_Sperm.hmr.specific Rhesus 4
Rhesus_Sperm.hmr.specific Gorilla 4
Rhesus_Sperm.hmr.specific Rhesus 4
Mouse_Sperm.hmr.specific Mouse 60
Mouse_Sperm.hmr.specific Rat 60
Rat_Sperm.hmr.specific Mouse 60
Rat_Sperm.hmr.specific Rat 60
Primate_core.hmr.specific Catarrhini 8
Primate_core.hmr.specific Murinae 8
Rodent_core.hmr.specific Catarrhini 8
Rodent_core.hmr.specific Murinae 8"

d <- read.delim(textConnection(s), header=T,sep=" ", strip.white=TRUE, as.is=T)
d$Sub <- 0

for (i in 1:nrow(d)) {
region <- d$Region[i]
lineage <- d$Lineage[i]
f <- paste(MutDir,"/",region, ".txt", sep="")
x <- read.table(f, header=T, as.is=T); row.names(x) <- x$name
d$Sub[i] <- sum(x[lineage, 4:9])
}


dr <- read.table(paste(ProjDir, "LineageMutation/random_sample_region_0001/sampling_results.txt", sep="/"), header=T)

regions <- unique(d$Region)
p <- list()

for (i in seq(1,length(regions),2)) {
  region1 <- regions[i]; region2 <- regions[i+1]
  region1rows <- which(d$Region == region1 )
  region2rows <- which(d$Region == region2 )
  region1rate <- d$Sub[region1rows[1]]/(d$Sub[region1rows[2]])
  region2rate <- d$Sub[region2rows[1]]/(d$Sub[region2rows[2]])
  sp1 <- d$Lineage[region1rows[1]]
  sp2 <- d$Lineage[region1rows[2]]
  dat <- data.frame(ratio=rep(0,5000))
  set.seed(1234) # set seed to be reproducible
  for (j in 1:5000) {
    id <- sample(1:nrow(dr), size=d$size[region1rows[1]]);
    dat$ratio[j] <- sum(dr[id,sp1])/sum(dr[id,sp2])
  }
  m <- mean(dat$ratio)
  std <- sd(dat$ratio)
  ci <- qnorm(p=0.95, mean=m, sd=std) - m
  l <- max( abs( region1rate - m), abs(region2rate-m))*1.1

  # pvalues from distribution
  pval1 <- 1- pnorm(q=region1rate, mean=m, sd=std)
  pval2 <- pnorm(q=region2rate, mean=m, sd=std)

  # normal distribution fit contour
  dd <- data.frame(x=seq(m-l, m+l, length.out=1000), y=dnorm(x=seq(m-l, m+l, length.out=1000), mean=m, sd=std ))
  dd <- rbind(dd, dd[1,])
  dd$z <- dd$y
  dd$z[abs(dd$x-m) > ci] <- 0

  p[[(i+1)/2]] <- ggplot(dat) +
    geom_histogram(aes(x=ratio, y=..density..), breaks=seq(m-l, m+l, length.out=50), fill="gray") +
    layer(data = dd, stat="identity", position = "identity", mapping = aes(x=x, y=z),
          geom = "area", params=list(fill="#FFCC00", alpha=.6)) +
    stat_function(fun=dnorm, xlim=c(m-l, m+l),
              args=list(mean=mean(dat$ratio), sd=sd(dat$ratio))) +
    labs(title=paste(sp1,"/",sp2), y="density", x="") +
    geom_segment(mapping=aes(x = x2, y = 0, xend = x2, yend = height), 
                 data=data.frame(x2=region1rate, height=max(dd$y)*0.3), color="#E14B3B") +
    geom_point(mapping=aes(x=x2, y=height), data=data.frame(x2=region1rate, height=max(dd$y)*0.3), 
               shape=16, size=4, color="#E14B3B" ) +
    geom_text(mapping=aes(x=x2, y=y2, label=text2), vjust=0, hjust=0.5, colour="#E14B3B", 
              data=data.frame(x2=region1rate, y2= max(dd$y)*0.45, text2=paste(sp1)), size=2.5) +
    geom_text(mapping=aes(x=x2, y=y2, label=text2), vjust=1, hjust=0.5, colour="#E14B3B", 
              data=data.frame(x2=region1rate, y2= max(dd$y), text2=format.pval(pval1, digits=4)), size=2.5) +
    geom_segment(mapping=aes(x = x2, y = 0, xend = x2, yend = height), 
                 data=data.frame(x2=region2rate, height=max(dd$y)*0.3), color="#1C8B98") +
    geom_point(mapping=aes(x=x2, y=height), data=data.frame(x2=region2rate, height=max(dd$y)*0.3), 
               shape=16, size=4, color="#1C8B98" ) +
    geom_text(mapping=aes(x=x2, y=y2, label=text2), vjust=0, hjust=0.5, colour="#1C8B98", 
              data=data.frame(x2=region2rate, y2= max(dd$y)*0.45, text2=paste(sp2)), size=2.5) +
    geom_text(mapping=aes(x=x2, y=y2, label=text2), vjust=1, hjust=0.5, colour="#1C8B98", 
              data=data.frame(x2=region2rate, y2= max(dd$y), text2=format.pval(pval2, digits=4)), size=2.5) +
    theme_bw() +
    theme(plot.title = element_text(size=10),
      axis.text= element_text(size = 8),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"))
}

source(paste(ProjDir, "scripts/multiplot.R", sep="/"))

pdf("Ratio_distribution_0001.pdf", width=6.6, height=3, useDingbats=FALSE)
multiplot(p[[1]], p[[2]], p[[3]], p[[4]], cols=2)
dev.off()


