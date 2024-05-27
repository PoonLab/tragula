require(Matrix)
require(transport)
require(jsonlite)
require(wordspace)
require(uwot)


#' topicspace
#' Map the word-context co-occurrence matrix to a lower-dimensional space
#' using uniform manifold approximation and projection (UMAP)
#' 
#' @param index.path:  character, path to CSV file with word index
#' @param cooccur.path:  character, path to CSV file with co-occurrence counts
#' @param max.words:  integer, maximum number of words to analyze
#' @param n.comp:  integer, number of components for UMAP
#' @param ...:  other arguments to pass to uwot::umap()
#' @return S3 object of class 'topicspace"
topicspace <- function(index.path, cooccur.path, author.path,
                       max.words=5000, n.comp=2, ...) {
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
  smx <- Matrix::sparseMatrix(i=ccm[,1], j=ccm[,2], x=1, index1=FALSE)
  
  # generate cosine distance matrix
  d1 <- wordspace::dist.matrix(t(smx)[1:max.words,], as.dist=TRUE)
  u1 <- uwot::umap(d1, n_components=n.comp, ...)  # run UMAP
  row.names(u1) <- index$word
  
  # load word counts by author
  by.author <- jsonlite::read_json(author.path, simplifyVector = TRUE)
  by.author <- sapply(by.author, unlist)
  
  obj <- list(um=u1, index=index, by.author=by.author)
  class(obj) <- c("topicspace", "list")
  obj
}

#' Generic plot function for topicspace S3 class
#' 
#' Plot the first two components of the distribution of words in the UMAP 
#' representation.  Optionally visualize the word frequency of a specific 
#' author instead of labeling words.
#' 
#' @param obj:  S3 object of class 'topicspace'
#' @param author:  character, if specified, then overlay the word frequency 
#'                 distribution for a specific author
#' @param limit:  integer, number of words to display (in decreasing order
#'                of global word frequency); defaults to 100
#' @param pt.cex:  numeric, character expansion factor for points, not 
#'                 including author-specific point scaling
#' @param text.cex:  numeric, character expansion factor for text
#' @param scale:  numeric, scaling factor for author-specific points
#' @param ...:  other options passed to the initial call to plot()
plot.topicspace <- function(obj, i=1, j=2, author=NA, limit=100, text.cex=0.5, 
                            pt.cex=0.1, scale=3, ...) {
  x <- obj$um[,i]
  y <- obj$um[,j]
  plot(x, y, type='n', bty='n', xaxt='n', yaxt='n', xlab=NA, ylab=NA,
       mar=rep(0,4), ...)
  points(x, y, pch=19, cex=pt.cex, col='grey')
  if (is.na(author)) {
    text(x[1:limit], y[1:limit], labels=obj$index$word[1:limit], cex=text.cex)
  } else {
    stopifnot(is.element(author, names(obj$by.author)))
    # a useful way of viewing each author's word frequency in UMAP space
    counts <- obj$by.author[[author]]
    idx <- match(obj$index$word, names(counts))
    # use relative frequencies so we do not over-emphasize authors who 
    # write more
    freq <- 100 * counts[idx] / sum(counts)
    points(x, y, cex=scale*sqrt(freq), pch=19, col=rgb(0,0,0,0.2))
  }
}

summary.topicspace <- function(obj) {
  cat("topicspace object\n")
  cat("  ", length(obj$by.author), "authors,", nrow(obj$index), "words,", 
      ncol(obj$um), "UMAP components\n")
  cat("  ", "Top words:", paste(head(obj$index$word), collapse=", "), 
      "\n")
  totals <- sapply(obj$by.author, length)
  cat("  ", "Words per author:", 
      paste(round(mean(totals), 2), " (range ", min(totals), "-", max(totals), ")", 
            sep=""))
}

print.topicspace <- function(obj) {
  summary(obj)
}


#' get.dist
#' 
#' Map each author's word counts to a topicspace and convert to 
#' weighted point patterns (discrete probability distribution over a 
#' finite number of points in a continuous space).  Calculate a 
#' pairwise matrix using Wasserstein (earth mover's) distances.
#' 
#' @param obj:  S3 object of class 'topicspace'
#' @param mc.cores:  integer, number of cores if using parallel
#' @return  dist object
get.dist <- function(obj, mc.cores=1) {
  # calculate weighted point patterns
  wpps <- lapply(obj$by.author, function(counts) {
    idx <- match(obj$index$word, names(counts))
    mass <- counts[idx]
    mass[is.na(mass)] <- 0  # handle missing entries
    mass <- mass / sum(mass)  # normalize
    transport::wpp(obj$um, mass=mass)
  })
  
  # calculate the pairwise Wasserstein distance matrix
  n <- length(wpps)
  if (require(parallel, quietly=TRUE)) {
    res <- mclapply(0:(n*n-1), function(k) {
      i <- k %/% n + 1
      j <- k %% n + 1
      if (i < j) { 
        transport::wasserstein(wpps[[i]], wpps[[j]], prob=TRUE) 
      } else {
        0
      }
    }, mc.cores = mc.cores)
    
    wdist <- matrix(unlist(res), nrow=n, ncol=n, byrow=T)
    
    # reflect upper triangular portion of matrix
    ix <- lower.tri(wdist, diag=FALSE)
    wdist[ix] <- t(wdist)[ix]  
  } else {
    # single-threaded version
    wdist <- matrix(0, nrow=n, ncol=n)
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        wdist[j, i] <- transport::wasserstein(wpps[[i]], wpps[[j]], prob=TRUE)
        wdist[i, j] <- wdist[j, i]  # symmetric matrix
      }
    }
  }
  
  rownames(wdist) <- names(wpps)
  colnames(wdist) <- names(wpps)  
  wdist
}

