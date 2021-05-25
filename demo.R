load("combs_fraser.rda")
suppressPackageStartupMessages(library(SummarizedExperiment))
se
library(org.Dm.eg.db)
mcols(se)$symbol <- mapIds(org.Dm.eg.db, rownames(se), "SYMBOL", "ENSEMBL")
rownames(se) <- mcols(se)$symbol
assay(se, "total") <- assay(se, "a1") + assay(se, "a2") 
assay(se, "ratio") <- assay(se, "a1") / assay(se, "total")
plotGene <- function(gene) {
  x <- se$slice
  y <- assay(se, "ratio")[gene,]
  plot(x, y, xlab="slice", ylab="ratio", ylim=c(0,1), main=gene)
  lw <- loess(y ~ x, data=data.frame(x,y))
  lines(sort(unique(x)), lw$fitted[order(unique(x))], col="red", lwd=2)
  abline(h=0.5, col="grey")
}
png(file="fig/slam.png")
plotGene("slam")
dev.off()
png(file="fig/uif.png")
plotGene("uif")
dev.off()
