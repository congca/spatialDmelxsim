load("combs_fraser.rda")
suppressPackageStartupMessages(library(SummarizedExperiment))
se
library(org.Dm.eg.db)
mcols(se)$symbol <- mapIds(org.Dm.eg.db, rownames(se), "SYMBOL", "ENSEMBL")
rownames(se) <- mcols(se)$symbol
# a1 = D simulans
# a2 = D melanogaster
assay(se, "total") <- assay(se, "a1") + assay(se, "a2") 
assay(se, "ratio") <- assay(se, "a1") / assay(se, "total")
plotGene <- function(gene) {
  x <- se$slice
  # deal with diff number of slices
  x[se$rep == 2] <- x[se$rep == 2] * 27/26
  x[se$rep == 4] <- x[se$rep == 4] * 27/25
  y <- assay(se, "ratio")[gene,]
  plot(x, y, xlab="slice", ylab="ratio", ylim=c(0,1), main=gene)
  lw <- loess(y ~ x, data=data.frame(x,y))
  lines(sort(lw$x), lw$fitted[order(lw$x)], col="red", lwd=2)
  abline(h=0.5, col="grey")
}
png(file="fig/slam.png")
plotGene("slam")
dev.off()
png(file="fig/uif.png")
plotGene("uif")
dev.off()
