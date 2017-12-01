#birth
source("../scripts/volcano.R")
x <- read.table("Rodent_lineage_specific.birth.filtered.MouseTFcount", as.is=T)
x$V3<-x$V2/sum(x$V2)
labels = unlist(strsplit(unlist(x$V1), "-mouse.merged.bed.lift2hg19"))
y <- read.table("Mouse_Sperm_Hammoud2014.hmr.MouseTFcount", as.is=T)
y$V2 <- y$V2 - x$V2
y$V3<-y$V2/sum(y$V2)
pdf("Lineage_to_mouse_merged_birth_TFBS_volcano.pdf", width=2.5, height=2.5, pointsize=8)
par(mar=c(3,3,1,1))
volcanoplot(x, y, orcutoff = 1, logpcutoff = 10,
            output = "Rodent_lineage_specific.birth.filtered.TFcount.pval.tab",
            xlabel = "Enriched in Root-Mouse gain birth", 
            ylabel="Depleted in Root-Mouse gain birth", labels=labels)
dev.off()


#extension
x <- read.table("Rodent_lineage_specific.extension.filtered.MouseTFcount", as.is=T)
x$V3<-x$V2/sum(x$V2)
labels = unlist(strsplit(unlist(x$V1), "-mouse.merged.bed.lift2hg19"))
y <- read.table("Mouse_Sperm_Hammoud2014.hmr.MouseTFcount", as.is=T)
y$V2 <- y$V2 - x$V2
y$V3<-y$V2/sum(y$V2)
pdf("Lineage_to_mouse_merged_extension_TFBS_volcano.pdf", width=2.5, height=2.5, pointsize=8)
par(mar=c(3,3,1,1))
volcanoplot(x, y,
            orcutoff = 1, logpcutoff = 10,
            output = "Rodent_lineage_specific.extension.filtered.TFcount.pval.tab",
            xlabel = "Enriched in Root-Mouse gain extension", 
            ylabel = "Depleted in Root-Mouse gain extension", labels=labels)
dev.off()

# non-promoter HMR gain
x <- read.table("Rodent_lineage_specific.gain.filtered.non_promoter.MouseTFcount", as.is=T)
x$V3<-x$V2/sum(x$V2)
labels = unlist(strsplit(unlist(x$V1), "-mouse.merged.bed.lift2hg19"))
y <- read.table("Mouse_Sperm_Hammoud2014.hmr.MouseTFcount", as.is=T)
y$V2 <- y$V2 - x$V2
y$V3<-y$V2/sum(y$V2)
pdf("Lineage_to_mouse_merged_gain_nonpromoter_TFBS_volcano.pdf", width=2.5, height=2.5, pointsize=8)
par(mar=c(3,3,1,1))
volcanoplot(x, y, orcutoff = 1, logpcutoff = 10,
            output = "Rodent_lineage_specific.gain.filtered.non_promoter.TFcount.pval.tab",
            xlabel = "Enriched in Root-Mouse gain nonpromoter",
            ylabel="Depleted in Root-Mouse gain nonpromoter", labels=labels)
dev.off()


# promoter HMR extension
x <- read.table("Rodent_lineage_specific.extension.filtered.at_promoter.MouseTFcount", as.is=T)
x$V3<-x$V2/sum(x$V2)
labels = unlist(strsplit(unlist(x$V1), "-mouse.merged.bed.lift2hg19"))
y <- read.table("Mouse_Sperm_Hammoud2014.hmr.MouseTFcount", as.is=T)
y$V2 <- y$V2 - x$V2
y$V3<-y$V2/sum(y$V2)
pdf("Lineage_to_mouse_merged_extension_promoter_TFBS_volcano.pdf", width=2.5, height=2.5, pointsize=8)
par(mar=c(3,3,1,1))
volcanoplot(x, y,
            orcutoff = 1, logpcutoff = 2,
            output = "Rodent_lineage_specific.extension.filtered.at_promoter.TFcount.pval.tab",
            xlabel = "Enriched in Root-Mouse gain extension at promoter",
            ylabel = "Depleted in Root-Mouse gain extension at promoter", labels=labels)
dev.off()


