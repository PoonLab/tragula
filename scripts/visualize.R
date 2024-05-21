require(igraph)
require(jsonlite)

# load data
wdist <- read.csv("results/wdist.csv", row.names=1)
by.author <- read_json("results/by_author.json", simplifyVector = TRUE)
by.author <- sapply(by.author, unlist)  # more convenient named vectors

# project distance matrix into 2/3 dimensions
mds <- cmdscale(wdist, k=6)
par(mar=rep(2,4))
plot(mds, type='n', bty='n', xaxt='n', yaxt='n', xlab=NA, ylab=NA)
text(mds, labels=names(by.author), cex=0.7, xpd=NA)

# bioinformatics/sequence analysis is the 3rd dimension...
plot(mds[,3:4], type='n', bty='n', xaxt='n', yaxt='n', xlab=NA, ylab=NA)
text(mds[,3:4], labels=names(by.author), cex=0.7, xpd=NA)

plot(mds[,5:6], type='n', bty='n', xaxt='n', yaxt='n', xlab=NA, ylab=NA)
text(mds[,5:6], labels=names(by.author), cex=0.7, xpd=NA)

# try network visualization instead
par(mar=c(5,5,1,1))
hist(wdist[upper.tri(wdist)], breaks=50)
#cutoff <- quantile(wdist[upper.tri(wdist)], 0.2)

# everyone should be connected to at least one other 
cutoff <- max(apply(wdist, 1, function(x) min(x[x>0]))) * 1.01
adj.mat <- wdist < cutoff
g <- graph_from_adjacency_matrix(adj.mat, mode="undirected", diag=F)
plot(g, vertex.shape="none", edge.width=2)

lc <- cluster_leiden(g, resolution_parameter=0.8)
sizes(lc)
plot(lc, g, vertex.shape="none", edge.width=2)


# k-nearest neighbour graph



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
  graph_from_adjacency_matrix(adj, diag=F, mode=mode)
}


#' Write graph to GraphViz DOT file
#' @param g:  'igraph' class object, graph to export
#' @param fn:  character, filename (path) to write DOT file
#' 
write.dot <- function(g, fn) {
  conn <- file(fn, "w")  # open connection to file
  cat("graph {\n", file=conn)
  cat("\toutputorder=edgesfirst;\n", file=conn)
  cat("\tnode [shape=rect, style=filled, fillcolor=white, margin=0.05, height=0];\n", file=conn)
  cat("\tedge [len=1.1];\n", file=conn)
  
  # write edge list
  el <- as_edgelist(g)
  for (i in 1:nrow(el)) {
    cat(paste("\t\"", el[i,1], "\"--\"", el[i,2], "\";\n", sep=""), file=conn)
  }
  cat("}\n", file=conn)
  close(conn)
}
