# python scripts/analyze.py --counts results/by_author.json --matrix results/cooccur.csv --index results/index.csv data

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
if (FALSE) {
  # this doesn't work very well!
  p1 <- prcomp_irlba(t(sparse), n=8)
  par(mar=c(5,5,1,1))
  plot(p1$x[,3:4], type='n', bty='n')
  text(p1$x[1:1000,3:4], labels=index$word[1:1000], cex=0.5)  
}

# try UMAP instead
require(wordspace)
#d1 <- wordspace::dist.matrix(t(sparse))
d1 <- wordspace::dist.matrix(t(sparse)[1:1000,], as.dist=TRUE)

require(uwot)
u1 <- uwot::umap(d1, n_components = 3)


par(mar=rep(0, 4))
plot(u1, type='n', bty='n', xaxt='n', yaxt='n', xlab=NA, ylab=NA)
text(u1, labels=index$word[1:1000], cex=0.5)


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


# finally, calculate the pairwise Wasserstein distance matrix
n <- length(objs)
mc.cores <- 10  # careful of RAM usage! about 2GB per core
if (require(parallel, quietly=TRUE)) {
  res <- mclapply(0:(n*n-1), function(k) {
    i <- k %/% n + 1
    j <- k %% n + 1
    if (i < j) { 
      wasserstein(objs[[i]], objs[[j]], prob=TRUE) 
    } else {
      0
    }
  }, mc.cores = mc.cores)
  wdist <- matrix(unlist(res), nrow=n, ncol=n, byrow=T)
  # reflect upper triagonal
  ix <- lower.tri(wdist, diag=FALSE)
  wdist[ix] <- t(wdist)[ix]  
} else {
  # single-threaded version
  wdist <- matrix(0, nrow=n, ncol=n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      wdist[i,j] <- wdist[j,i] <- wasserstein(objs[[i]], objs[[j]], prob=TRUE)
    }
  }
}

# write matrix out to file
rownames(wdist) <- names(objs)
colnames(wdist) <- names(objs)
write.csv(wdist, "results/wdist.csv")
