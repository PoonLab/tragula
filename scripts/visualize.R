require(igraph)


#' Use multi-dimensional scaling to visualize author distance matrix
#' @param wdist:  matrix or data.frame, pairwise distance matrix for authors
#' @param 
plot.wdist <- function(wdist, k=2) {
  mds <- cmdscale(wdist, k=k)
  par(mar=rep(1,4))
  plot(mds, type='n', bty='n', xaxt='n', yaxt='n', xlab=NA, ylab=NA)
  text(mds, labels=row.names(wdist), cex=0.7, xpd=NA)  
}


#' Make k-nearest neighbour graph from a distance matrix.
#' @param dmx:  matrix, distance matrix
#' @param k:  integer, number of nearest neighbors
#' @param cutoff:  numeric, if set, then nearest neighbors are excluded
#'                 if their distance exceeds this value.
#' @param names:  character, override input `dmx` matrix row.names
#' @param undirected:  bool, if TRUE, then build an undirected graph
#'                     from an asymmetric adjacency matrix
knn <- function(dmx, k=3, cutoff=NA, names=NA, undirected=TRUE) {
  n <- nrow(dmx)
  # handle default values
  if (is.na(names)) {
    names <- row.names(dmx)
  }
  if (is.na(cutoff)) {
    cutoff <- max(apply(dmx, 1, function(x) min(x[x>0])))
  }
  
  adj <- matrix(0, nrow=n, ncol=n, dimnames=list(names, names))
  for (i in 1:n) {
    row <- as.numeric(dmx[i,])
    ranks <- order(row)
    # exclude self and find nearest neighbors
    nn <- ranks[-(ranks==i)][1:k]
    nn <- nn[row[nn] <= cutoff]  # exclude edges that are too long
    adj[i,nn] <- 1
  }
  mode <- ifelse(undirected, "undirected", "directed")
  igraph::graph_from_adjacency_matrix(adj, diag=F, mode=mode)
}


#' Write graph to GraphViz DOT file
#' @param g:  'igraph' class object, graph to export
#' @param fn:  character, filename (path) to write DOT file
write.dot <- function(g, fn) {
  conn <- file(fn, "w")  # open connection to file
  cat("graph {\n", file=conn)
  cat("\toutputorder=edgesfirst;\n", file=conn)
  cat("\tnode [shape=rect, style=filled, fillcolor=white, margin=0.05, height=0];\n", file=conn)
  cat("\tedge [len=1.1];\n", file=conn)
  
  # write edge list
  el <- igraph::as_edgelist(g)
  for (i in 1:nrow(el)) {
    cat(paste("\t\"", el[i,1], "\"--\"", el[i,2], "\";\n", sep=""), file=conn)
  }
  cat("}\n", file=conn)
  close(conn)
}



