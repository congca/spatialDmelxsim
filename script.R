library(GEOquery)
g <- getGEO("GSE102233")
i <- 3
#pData(g[[i]])[1,]
pData(g[[i]])$title
pData(g[[i]])$"hybrid:ch1"
pData(g[[i]])$supplementary_file_1[1]

titles <- c(pData(g[[2]])$title,pData(g[[3]])$title)
table(sub(".*(sl.*)$","\\1",titles))
strain <- substr(titles, 1, 7)
table(strain)
length(titles)

## for (i in 2:3) {
##   for (j in seq_along(pData(g[[i]])$title)) {
##     download.file(pData(g[[i]])$supplementary_file_1[j],
##                   destfile=paste0("ase_files/",i,"_",j,".txt.gz"), method="wget")
##   }
## }

means <- list(numeric(52), numeric(80))
for (i in 1:2) {
  for (j in 1:(c(52,80)[i])) {
    cat(j)
    x <- read.delim(paste0("ase_files/",i+1,"_",j,".txt.gz"))
    head(x[,1:4]) # these are dm6
    not_xy <- !x$CHROMOSOME %in% c("dmel_X","dmel_Y")
    standard <- nchar(x$CHROMOSOME <= 7)
    table(not_xy & standard)
    y <- x[not_xy & standard,c(1,5:6)]
    y <- y[complete.cases(y),]
    y <- y[rowSums(as.matrix(y[,2:3])) >= 10,]
    pc <- 1
    ratio <- (y[,3] + pc)/(rowSums(y[,2:3]) + 2*pc)
    means[[i]][j] <- mean(ratio)
  }
}
boxplot(means)    
mean_ratio <- unlist(means)
boxplot(mean_ratio ~ strain, ylim=c(0, 1))
