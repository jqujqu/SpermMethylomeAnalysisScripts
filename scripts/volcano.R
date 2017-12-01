
volcanoplot <- function(x, y, orcutoff=1.5, logpcutoff=10,
                        xlabel = "", ylabel="", labels="", output=""){
    pvals <- data.frame(V1=x$V1)
    for (i in 1:nrow(x) ){
        mat <- matrix(c(x$V2[i], y$V2[i], sum(x$V2[-i]), sum(y$V2[-i])), nrow=2, byrow=T); 
        m <- y$V2[i] + x$V2[i]
        n <- sum(x$V2[-i]) + sum(y$V2[-i])
        k <- sum(x$V2)
        q <- x$V2[i]
        pvals$hyper.log10.P.val[i] <- min(phyper(q, m, n, k, lower.tail = F, log.p = T)/log(10), 
                                          phyper(q, m, n, k, lower.tail = T, log.p = T)/log(10))
        pvals$Pval[i] <- fisher.test(mat)$p.value 
        pvals$OR[i] <- fisher.test(mat)$estimate
    }
    pvals$log10.p.adjust <- pvals$hyper.log10.P.val+log(nrow(x))/log(10)
    pvals$log10.p.adjust[which(pvals$log10.p.adjust >0)]  <- 0
    if (output!="") {
      pvals$FG <- x$V2
      pvals$BG <- y$V2
      names(pvals)[1] <- "TF"
      write.table(pvals, output, row.names=F, col.names=T,sep="\t", quote=F)
    }
    logOR <- log(pvals$OR)/log(2)
    logP <- -pvals$log10.p.adjust 
    indx <-  (logOR > log(orcutoff)/log(2) & logP > logpcutoff)
    indy <-  (logOR < -log(orcutoff)/log(2) & logP > logpcutoff)

    # plot 
    plot(logOR, logP, xlim=c(-4,4), pch="", xlab="log2 Odds Ratio", ylab="-log10 Adjusted P-value")
    points(logOR[!(indx | indy)], logP[!(indx | indy)], pch=19, col="gray")
    points(logOR[indx], logP[indx], col="blue", pch=19)
    points(logOR[indy], logP[indy], col="red", pch=19)
    abline(v=c(-log(orcutoff)/log(2), log(orcutoff)/log(2)), lty="dashed")
    abline(h=logpcutoff, lty="dashed")
    legend("topleft", fill=c("red","blue"), legend=c(ylabel, xlabel))
    chw <- par()$cxy[ 1 ]  ##  character width
    chh <- par()$cxy[ 2 ]  ##  character height
    if(length(indx) > 0) {
        text(x = logOR[indx] + chw*0.2, y = logP[indx] + chh*0,
            labels = labels[indx], adj = 0, cex=0.8)
    }
    if(length(indy) > 0) {
        text(x = logOR[indy] + chw*0.2, y = logP[indy] + chh*0,
            labels = labels[indy], adj = 2, cex=0.8)
    }
}


