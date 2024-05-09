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

make_ego_graph(g, order=3)

# k-nearest neighbour graph
cutoff <- max(apply(wdist, 1, function(x) min(x[x>0])))*1.01
n <- nrow(wdist)
nm <- row.names(wdist)
mx <- matrix(0, nrow=n, ncol=n, dimnames=list(nm, nm))
max.deg <- 3
for (i in 1:nrow(wdist)) {
  row <- as.numeric(wdist[i,])
  temp <- order(row)
  idx <- temp[-(temp==i)][1:max.deg]
  j <- idx[row[idx] < cutoff]
  mx[i,j] <- 1
}

g <- graph_from_adjacency_matrix(mx, diag=F, mode="undirected")

plot(g, vertex.shape="none", label=row.names(wdist),
     edge.width=2, edge.arrow.mode='-',
     layout=layout_with_kk(g))

write.dot <- function(g, fn) {
  conn <- file(fn, "w")
  cat("graph {\n", file=conn)
  cat("\tnode [shape=box];\n", file=conn)
  
  # write edge list
  el <- as_edgelist(g)
  for (i in 1:nrow(el)) {
    cat(paste("\t\"", el[i,1], "\"--\"", el[i,2], "\";\n", sep=""), file=conn)
  }
  cat("}\n", file=conn)
  close(conn)
}
