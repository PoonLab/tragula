require(igraph)
require(jsonlite)

#' Use multi-dimensional scaling to visualize author distance matrix
#' @param wdist:  S3 object of class `wdist`, pairwise distance matrix for authors
#' @param i:  integer, select component for x-axis (default 1)
#' @param j:  integer, select component for y-axis (default 2)
#' @param type:  character, 'u' for UMAP, 'm' for multidimensional scaling
#' @param k:  integer, number of components for dimensionality reduction
#' @param labels:  character, optionally specify a custom set of labels
#' @param col:  colour for text
plot.wdist <- function(wdist, i=1, j=2, type='u', k=2, labels=NA, col=1, 
                       ...) {
  if (all(is.na(labels))) {
    labels <- attr(wdist, "Labels")
  }
  if (type=='m') {
    proj <- cmdscale(wdist, k=k)
  } else if (type=='u') {
    proj <- umap(wdist, n_components=k)
  }
  x <- proj[,i]
  y <- proj[,j]
  par(mar=rep(1,4))
  plot(x, y, type='n', bty='n', xaxt='n', yaxt='n', xlab=NA, ylab=NA, ...)
  text(x, y, labels=labels, cex=0.7, xpd=NA, col=col)
}


#' Make k-nearest neighbour graph from a distance matrix.
#' @param wdist:  S3 object of class 'wdist', Wasserstein distance matrix
#' @param k:  integer, number of nearest neighbors
#' @param cutoff:  numeric, if set, then nearest neighbors are excluded
#'                 if their distance exceeds this value.
#' @param names:  character, override input `dmx` matrix row.names
#' @param undirected:  bool, if TRUE, then build an undirected graph
#'                     from an asymmetric adjacency matrix
make.knn <- function(wdist, k=3, cutoff=NA, names=NA, undirected=TRUE) {
  dmx <- as.matrix(wdist)
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
  g <- igraph::graph_from_adjacency_matrix(adj, diag=F, mode=mode)
  if (mode=="directed") {
    idx <- which(as.logical(adj))
    E(g)$weight <- dmx[idx]    
  } else {
    eids <- sapply(as_ids(E(g)), function(x) strsplit(x, "\\|")[[1]])
    parent <- match(eids[1,], as_ids(V(g)))
    child <- match(eids[2,], as_ids(V(g)))
    E(g)$weight <- sapply(1:ncol(eids), function(i) dmx[parent[i], child[i]])
  }
  g
}


#' Write graph to GraphViz DOT file
#' @param g:  'igraph' class object, graph to export
#' @param fn:  character, filename (path) to write DOT file
write.dot <- function(g, fn, labels=NA, groups=NA, pal=NA) {
  if (all(is.na(pal)) & all(!is.na(groups))) {
    # default palette
    pal <- hcl.colors(n=length(unique(groups)), palette="Set3")
  }
  
  conn <- file(fn, "w")  # open connection to file
  cat("graph {\n", file=conn)
  cat("\toutputorder=edgesfirst;\n", file=conn)
  cat("\tnode [shape=rect, style=filled, margin=0.05, height=0];\n", file=conn)
  cat("\tedge [len=1.1];\n", file=conn)
  
  # write node list
  nodes <- as_ids(V(g))
  if (all(is.na(labels))) {
    labels <- nodes  # by default use node names as labels
  }
  for (i in 1:length(nodes)) {
    cat(paste("\t\"", nodes[i], "\" [label=\"", labels[i], "\"", sep=""), file=conn)
    if (any(!is.na(groups))) {
      cat(", fillcolor=\"", pal[groups[i]], "\"", file=conn)
    }
    cat("];\n", file=conn)  # EOL
  }
  
  # write edge list
  el <- igraph::as_edgelist(g)
  for (i in 1:nrow(el)) {
    cat(paste("\t\"", el[i,1], "\"--\"", el[i,2], "\";\n", sep=""), file=conn)
  }
  cat("}\n", file=conn)
  close(conn)
}


#' Generate a JSON object to pass to JavaScript layer via r2d3
#' @param g:  'igraph' class object
#' @param labels:  character, optionally use custom node labels
export.json <- function(g, labels=NA, groups=NA) {
  nodes <- data.frame(id=as_ids(V(g)))
  if (all(is.na(labels))) {
    nodes$label <- nodes$id
  } else {
    nodes$label <- labels
  }
  if (all(is.na(groups))) {
    nodes$group <- 1
  } else {
    nodes$group <- as.integer(as.factor(groups))
  }
  edges <- as.data.frame(igraph::as_edgelist(g))
  names(edges) <- c('source', 'target')
  edges$distance <- E(g)$weight
  jsonlite::toJSON(list(nodes=nodes, edges=edges))
}

