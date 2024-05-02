require(Matrix)
require(irlba)
require(transport)
require(jsonlite)

setwd("~/git/tragula")

# co-occurrence matrix
ccm <- read.csv("results/cooccur.csv", header=F)

# presence/absence seems to work better
sparse <- sparseMatrix(i=ccm[,1], j=ccm[,2], x=1, index1=FALSE)
#sparse <- sparseMatrix(i=ccm[,1], j=ccm[,2], x=ccm[,3], index1=FALSE)

# load global index of all words
index <- read.csv("results/index.csv")
index <- index[!is.na(index$index), ]

# PCA projection of word co-occurrence matrix
p1 <- prcomp_irlba(t(sparse), n=8)

# load author-specific word counts
by.author <- read_json("results/by_author.json", simplifyVector = TRUE)
by.author <- sapply(by.author, unlist)  # more convenient named vectors

# map each author's word counts to global index and convert to a 
# weighted point pattern (discrete probability distribution over a 
# finite number of points in a continuous space)
objs <- lapply(by.author, function(wordcounts) {
  idx <- match(index$word, names(wordcounts))
  mass <- wordcounts[idx]
  mass[is.na(mass)] <- 0
  mass <- mass / sum(mass)
  wpp(p1$x, mass=mass)
})
n <- length(objs)

# finally, calculate the pairwise Wasserstein distance matrix
wdist <- matrix(0, nrow=n, ncol=n)
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    wdist[i,j] <- wdist[j,i] <- wasserstein(objs[[i]], objs[[j]], prob=TRUE)
  }
}

mds <- cmdscale(wdist, k=2)

par(mar=rep(2,4))
plot(mds, type='n', bty='n', xaxt='n', yaxt='n', xlab=NA, ylab=NA)
text(mds, labels=names(by.author), cex=0.7, xpd=NA)
