library(GEOquery)
g <- getGEO("GSE102233")
i <- 2
#pData(g[[i]])[1,]
pData(g[[i]])$title
pData(g[[i]])$"hybrid:ch1"
pData(g[[i]])$supplementary_file_1[1]

titles <- sort(c(pData(g[[2]])$title,pData(g[[3]])$title))
table(sub(".*(sl.*)$","\\1",titles))

download.file(pData(g[[i]])$supplementary_file_1[1], "ase.txt.gz")
x <- read.delim("ase.txt.gz")
head(x[,1:4]) # these are dm6
y <- x[,c(1,5:6)]
y <- y[complete.cases(y),]
y <- y[rowSums(as.matrix(y[,2:3])) > 10,]
hist((y[,3] + 1) / (y[,2] + y[,3] + 2), breaks=100)
