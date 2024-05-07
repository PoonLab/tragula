# python scripts/analyze.py --counts results/by_author.json --matrix results/cooccur.csv --index results/index.csv data

require(Matrix)
require(irlba)
require(transport)
require(jsonlite)
require(igraph)
require(parallel)

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
res <- mclapply(1:(n-1), function(i) {
  sapply((i+1):n, function(j) {
    wasserstein(objs[[i]], objs[[j]], prob=TRUE)  
  })
}, mc.cores = 4)  # careful of RAM usage!

for (i in 1:(n-1)) {
  #print(i)
  for (j in (i+1):n) {
    wdist[i,j] <- wdist[j,i] <- wasserstein(objs[[i]], objs[[j]], prob=TRUE)
  }
}

# write matrix out to file
rownames(wdist) <- names(objs)
colnames(wdist) <- names(objs)
write.csv(wdist, "results/wdist.csv")

# project distance matrix into 2/3 dimensions
mds <- cmdscale(wdist, k=4)
par(mar=rep(2,4))
plot(mds, type='n', bty='n', xaxt='n', yaxt='n', xlab=NA, ylab=NA)
text(mds, labels=names(by.author), cex=0.7, xpd=NA)

# bioinformatics/sequence analysis is the 3rd dimension...
plot(mds[,3:4], type='n', bty='n', xaxt='n', yaxt='n', xlab=NA, ylab=NA)
text(mds[,3:4], labels=names(by.author), cex=0.7, xpd=NA)

# try network visualization instead
hist(wdist[upper.tri(wdist)])
#cutoff <- quantile(wdist[upper.tri(wdist)], 0.2)

# everyone should be connected to at least one other 
cutoff <- max(apply(wdist, 1, function(x) min(x[x>0]))) * 1.01

adj.mat <- wdist < cutoff
g <- graph_from_adjacency_matrix(adj.mat, mode="undirected", diag=F)
plot(g)
