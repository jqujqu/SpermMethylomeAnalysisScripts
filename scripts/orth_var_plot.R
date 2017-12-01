dirs <- dir()
f <- character(0)
for (path in dirs ){
  files <- dir(path=path)
  files <- files[c(grep ("ORTH_", files), grep ("VAR_", files))]
  f <- c(f, paste(path,files, sep="/"))
}

mat <- data.frame(size= numeric(0), species=character(0),
                  region=character(0),  context = character(0),  stringsAsFactors = F)
for (i in f){
 x <- read.table(i, as.is=T)
 filename <- unlist(strsplit(i,"/"))[2]
 species <- unlist(strsplit(filename,"_"))[1]
 suffix <- unlist(strsplit(filename,"\\."))[3]
 region <- unlist(strsplit(suffix, "_"))[1]
 context <- unlist(strsplit(suffix, "_"))[2]
 mat <- rbind(mat, data.frame(size= x$V3-x$V2, species=rep(species, nrow(x)),
            region=rep(region, nrow(x)), context = rep(context, nrow(x)), stringsAsFactors = F))
}

mat <- mat[mat$size < quantile(mat$size, 0.995), ]
mat$species <- factor(mat$species, levels=c("Human","Chimp","Gorilla","Rhesus","Mouse", "Rat","Dog"))
mat$context <- factor(mat$context, levels = c("promoter","repeat","other"))
mat$region <- factor(mat$region, levels=c("ORTH", "VAR"))

library(ggplot2)
sp <- ggplot(mat, aes(x=context, y=size, fill=context, colour=region)) +
stat_boxplot(geom = "errorbar", position=position_dodge(1)) +
geom_boxplot(outlier.shape=NA, aes(colour=region), position=position_dodge(1))
p<- sp + facet_grid(. ~ species) +
  scale_fill_manual(values=c("#0071bc", "#29abe2", "#39b54a") )+
  scale_colour_manual(values=c("#000000", "#999999"))+
  theme_bw() +
  theme(axis.title.x=element_blank() ,
        axis.text.y = element_text(size=8 ),
        axis.title.y = element_text(size=10),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(size=10),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank())

pdf("OrthVar.pdf", width=6.6, height=2)
p
dev.off()

