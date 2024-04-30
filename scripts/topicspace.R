require(Matrix)
require(irlba)

setwd("~/git/tragula")

# co-occurrence matrix
ccm <- read.csv("results/cooccur.csv", header=F)

# presence/absence seems to work better
sparse <- sparseMatrix(i=ccm[,1], j=ccm[,2], x=1, index1=FALSE)
#sparse <- sparseMatrix(i=ccm[,1], j=ccm[,2], x=ccm[,3], index1=FALSE)

index <- read.csv("results/index.csv")
index <- index[!is.na(index$index), ]

p1 <- prcomp_irlba(t(sparse), n=6)
plot(p1$x)

plot(p1$x, type='n')
text(p1$x, labels=index$word, cex=0.5)

plot(p1$x[,3:4], type='n')
text(p1$x[,3:4], labels=index$word, cex=0.5)

plot(p1$x[,5:6], type='n')
text(p1$x[,5:6], labels=index$word, cex=0.5)
