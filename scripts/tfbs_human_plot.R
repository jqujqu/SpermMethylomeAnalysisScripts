
#birth
source("../scripts/volcano.R")
x <- read.table("Primate_lineage_specific.birth.filtered.TFcount", as.is=T)
x$V3<-x$V2/sum(x$V2)
labels = unlist(strsplit(unlist(x$V1), ".bed"))
y <- read.table("Human_Sperm.hmr.TFcount", as.is=T)
y$V2 <- y$V2 - x$V2
y$V3<-y$V2/sum(y$V2)
pdf("Lineage_to_human_merged_birth_TFBS_volcano.pdf", width=2.5, height=2.5, pointsize=8)
par(mar=c(3,3,1,1))
volcanoplot(x, y, orcutoff = 1, logpcutoff = 10,
            output = "Primate_lineage_specific.birth.filtered.TFcount.pval.tab",
            xlabel = "Enriched in Root-Human gain birth", 
            ylabel="Depleted in Root-Human gain birth", labels=labels)
dev.off()


#extension
x <- read.table("Primate_lineage_specific.extension.filtered.TFcount", as.is=T)
x$V3<-x$V2/sum(x$V2)
labels = unlist(strsplit(unlist(x$V1), ".bed"))
y <- read.table("Human_Sperm.hmr.TFcount", as.is=T)
y$V2 <- y$V2 - x$V2
y$V3<-y$V2/sum(y$V2)
pdf("Lineage_to_human_merged_extension_TFBS_volcano.pdf", width=2.5, height=2.5, pointsize=8)
par(mar=c(3,3,1,1))
volcanoplot(x, y,
            orcutoff = 1, logpcutoff = 10,
            output = "Primate_lineage_specific.extension.filtered.TFcount.pval.tab",
            xlabel = "Enriched in Root-Human gain extension", 
            ylabel = "Depleted in Root-Human gain extension", labels=labels)
dev.off()

# nonpromoter gain
x <- read.table("Primate_lineage_specific.gain.filtered.non_promoter.TFcount", as.is=T)
x$V3<-x$V2/sum(x$V2)
labels = unlist(strsplit(unlist(x$V1), ".bed"))
y <- read.table("Human_Sperm.hmr.TFcount", as.is=T)
y$V2 <- y$V2 - x$V2
y$V3<-y$V2/sum(y$V2)
pdf("Lineage_to_human_merged_gain_nonpromoter_TFBS_volcano.pdf", width=2.5, height=2.5, pointsize=8)
par(mar=c(3,3,1,1))
volcanoplot(x, y, orcutoff = 1, logpcutoff = 10,
            output = "Primate_lineage_specific.gain.filtered.non_promoter.TFcount.pval.tab",
            xlabel = "Enriched in Root-Human gain nonpromoter",
            ylabel="Depleted in Root-Human gain nonpromoter", labels=labels)
dev.off()


# promtoer HMR extension
x <- read.table("Primate_lineage_specific.extension.filtered.at_promoter.TFcount", as.is=T)
x$V3<-x$V2/sum(x$V2)
labels = unlist(strsplit(unlist(x$V1), ".bed"))
y <- read.table("Human_Sperm.hmr.TFcount", as.is=T)
y$V2 <- y$V2 - x$V2
y$V3<-y$V2/sum(y$V2)
pdf("Lineage_to_human_merged_extension_promoter_TFBS_volcano.pdf", width=2.5, height=2.5, pointsize=8)
par(mar=c(3,3,1,1))
volcanoplot(x, y,
            orcutoff = 1, logpcutoff = 5,
            output = "Primate_lineage_specific.extension.filtered.at_promoter.TFcount.pval.tab",
            xlabel = "Enriched in Root-Human extension at promoter",
            ylabel = "Depleted in Root-Human extension at promtoer", labels=labels)
dev.off()

