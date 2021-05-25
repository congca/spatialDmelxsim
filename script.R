library(GEOquery)
g <- getGEO("GSE102233")
titles <- c(pData(g[[2]])$title,pData(g[[3]])$title)
slice <- as.numeric(sub(".*sl(.*?)$","\\1",titles))
table(slice)
strain <- substr(titles, 1, 7)
table(strain)
rep <- as.numeric(sub(".*rep(.*?)_.*","\\1",titles))
rep <- ifelse(strain == "simXmel", rep + 3, rep)
table(rep, strain)
table(rep, slice)
length(titles)

## for (i in 2:3) {
##   for (j in seq_along(pData(g[[i]])$title)) {
##     download.file(pData(g[[i]])$supplementary_file_1[j],
##                   destfile=paste0("ase_files/",i,"_",j,".txt.gz"), method="wget")
##   }
## }


n <- length(titles)
g <- 13620
a1 <- matrix(nrow=g, ncol=n) # effect / alt allele (D sim)
a2 <- matrix(nrow=g, ncol=n) # non-effect / ref allele (D mel)
i <- 1
j <- 1
x <- read.delim(paste0("ase_files/",i+1,"_",j,".txt.gz"))
rownames(a1) <- rownames(a2) <- x$FEATURE
for (i in 1:2) {
  for (j in 1:(c(52,80)[i])) {
    cat(j)
    x <- read.delim(paste0("ase_files/",i+1,"_",j,".txt.gz"))
    jj <- j + (i-1)*52
    a1[,jj] <- x$ALT_COUNTS
    a2[,jj] <- x$REFERENCE_COUNTS
  }
}

range <- strsplit(x[,4], "-")
suppressPackageStartupMessages(library(SummarizedExperiment))
genes <- GRanges(sub("dmel_","",x[,2]), IRanges(
                          as.numeric(sapply(range, `[`, 1)),
                          as.numeric(sapply(range, `[`, 2))),
                 strand=x[,3], gene_id=x[,1])
genome(genes) <- "dmr5.57"
names(genes) <- genes$gene_id
coldata <- data.frame(strain=factor(strain), slice, rep)
se <- SummarizedExperiment(assays=list(a1=a1, a2=a2), rowRanges=genes, colData=coldata)
colnames(se) <- paste0("r", se$rep, "s", se$slice)
metadata(se) <- list(a1="alt / D simulans",
                     a2="ref / D melanogaster",
                     author="Combs PA, Fraser HB",
                     title="Spatially varying cis-regulatory divergence in Drosophila embryos elucidates cis-regulatory logic",
                     journal="PLOS Genetics",
                     year="2018",
                     additional="14(11):e1007631",
                     doi="10.1371/journal.pgen.1007631")

# FWIW ~99% of the genes are in dm6 and have same coordinates
#library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
#g <- genes(TxDb.Dmelanogaster.UCSC.dm6.ensGene)


