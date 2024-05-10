# python scripts/analyze.py --counts results/by_author.json --matrix results/cooccur.csv --index results/index.csv data

require(Matrix)
require(irlba)
require(transport)
require(jsonlite)
require(wordspace)
require(uwot)

#setwd("~/git/tragula")

#' topicspace
#' Map the word-context co-occurrence matrix to a lower-dimensional space
#' @param index.path:  character, path to CSV file with word index
#' @param cooccur.path:  character, path to CSV file with co-occurrence counts
#' @param max.words:  integer, maximum number of words to analyze
#' @param n.comp:  integer, number of components for UMAP
#' @param ...:  other arguments to pass to uwot::umap()
#' @return matrix, UMAP coordinates for words
topicspace <- function(index.path, cooccur.path, 
                       max.words=5000, n.comp=3, ...) {
  # load global index of all words
  index <- read.csv(index.path)
  # check file integrity
  stopifnot(
    is.element(c("word", "count"), names(index)),
    index$index==(0:(nrow(index)-1)), 
    diff(index$count) <= 0
    )
  index <- index[1:max.words, ]
  
  # load co-occurrence as a sparse matrix (doc #, word #, count)
  ccm <- read.csv(cooccur.path, header=F)
  stopifnot(ncol(ccm)==3, apply(ccm, 2, class)=="integer")
  smx <- sparseMatrix(i=ccm[,1], j=ccm[,2], x=1, index1=FALSE)
  
  # generate cosine distance matrix
  d1 <- wordspace::dist.matrix(t(sparse)[1:max.words,], as.dist=TRUE)
  u1 <- uwot::umap(d1, n_components=n.comp, ...)  # run UMAP
  row.names(u1) <- index$word
  u1
}



# sort by frequency in descending order
#index <- index[order(index$count, decreasing = TRUE), ]

par(mar=rep(0, 4))
plot(u1, type='n', bty='n', xaxt='n', yaxt='n', xlab=NA, ylab=NA)
text(u1[1:5000,], labels=index$word[1:5000], cex=0.5)


# load author-specific word counts
by.author <- read_json("results/by_author.json", simplifyVector = TRUE)
by.author <- sapply(by.author, unlist)  # more convenient named vectors

if (FALSE) {
  # top two words are negatively correlated
  y <- sapply(by.author, function(a) 
    ifelse(is.element("patient", names(a)), a[["patient"]]/length(a), NA))
  x <- sapply(by.author, function(a) 
    ifelse(is.element("cell", names(a)), a[["cell"]]/length(a), NA))
  
  plot(x,y, type='n', log='xy', bty='n', xlab="f(Cell)", ylab="f(Patient)")
  text(x, y, labels=surname, cex=0.5, xpd=NA)
  cor.test(x, y)  
}


fingerprint <- function(author) {
  # a useful way of viewing each author's word frequency in UMAP space
  idx <- match(index$word[1:5000], names(by.author[[author]]))
  plot(u1, pch=19, cex=0.1, col='grey')
  points(u1, cex=sqrt(by.author[[author]][idx])/2, pch=19, 
         col=rgb(0,0,0,0.2))  
}


# map each author's word counts to global index and convert to a 
# weighted point pattern (discrete probability distribution over a 
# finite number of points in a continuous space)
objs <- lapply(by.author, function(wordcounts) {
  idx <- match(index$word[1:5000], names(wordcounts))
  mass <- wordcounts[idx]
  mass[is.na(mass)] <- 0
  mass <- mass / sum(mass)
  wpp(u1, mass=mass)
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
