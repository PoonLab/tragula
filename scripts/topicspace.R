require(wordspace)
require(Matrix)
require(irlba)

setwd("~/git/tragula")

# co-occurrence matrix
ccm <- read.csv("test.csv", header=F)
sparse <- sparseMatrix(i=ccm[,1], j=ccm[,2], x=1, index1=FALSE)
#sparse <- sparseMatrix(i=ccm[,1], j=ccm[,2], x=ccm[,3], index1=FALSE)

index <- read.csv("index.txt")
index <- index[!is.na(index$index), ]

p1 <- prcomp_irlba(t(sparse), n=6)
plot(p1$x)

plot(p1$x, type='n')
text(p1$x, labels=index$word, cex=0.5)

plot(p1$x[,3:4], type='n')
text(p1$x[,3:4], labels=index$word, cex=0.5)

plot(p1$x[,5:6], type='n')
text(p1$x[,5:6], labels=index$word, cex=0.5)

#dm <- wordspace::dist.matrix(t(sparse), method="euclidean")
#p2 <- prcomp(dm)
#plot(p2$x)

#require(uwot)
#um <- umap(sparse)
