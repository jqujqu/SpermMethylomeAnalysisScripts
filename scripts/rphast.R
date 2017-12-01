library(rphast)
chroms <- c(1:22)
for (i in chroms) {
chrom <- paste("chr", i, sep="")
tree <- "(((((hg19,panTro4),gorGor3),rheMac3),(mm10,rn5)),canFam3)"
m <- read.msa(paste(chrom,".7sp.fasta", sep=""), format="FASTA")
mod <- phyloFit(m, tree=tree, subst.mod="UNREST")
write(mod$tree,file=paste(chrom,".7sp.UNREST_tree.nwk", sep=""))
}

