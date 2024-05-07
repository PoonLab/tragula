require(igraph)

wdist <- read.csv("results/wdist.csv")

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

