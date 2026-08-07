
setwd("~/git/tragula")

index <- read.csv("results/palm_index.csv")
nrow(index)  # number of words
sum(index$count < 10)

# fit Zipf-Mandelbrot distribution using least-squares
zmdist <- function(rank, alpha, beta, C) {
  C/(rank+beta)^alpha
}
objfunc <- function(par) {
  pred <- zmdist(1:nrow(index), par[1], par[2], par[3])
  sum((log(pred)-log(index$count))^2)
}
res <- optim(c(1,10,5000), objfunc)

pdf("doc/wordfreq.pdf", width=4, height=4)
par(mar=c(5,5,1,1), family='Palatino')
plot(index$count, log='xy', type='s', bty='n', yaxt='n',
     ylab='Word frequency', xlab='Log rank', cex.axis=0.8)
lines(1:nrow(index), 
      zmdist(1:nrow(index), alpha=res$par[1], beta=res$par[2], C=res$par[3]),
      col='red')
x <- c(2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000)
axis(side=2, at=x, label=x, las=2, cex.axis=0.8)
abline(h=c(10, 100, 1000), lty=2, col=rgb(0,0,0,0.5))
dev.off()

# limiting analysis to top 5,000 words still covers 86% of corpus
csum <- cumsum(index$count) / sum(index$count)
plot(csum)
csum[5000]


# examine the two most frequent words, "cell" and "patient"
require(jsonlite)

by.author <- read_json("results/palm_author.json", simplifyVector = TRUE)
cell <- unlist(sapply(by.author, function(x) 
  ifelse(is.null(x[['cell']]), 0, x[['cell']] / sum(unlist(x)))
  ))
patient <- unlist(sapply(by.author, function(x) 
  ifelse(is.null(x[['patient']]), 0, x[['patient']] / sum(unlist(x)))
))

# not a significant correlation
cor.test(cell, patient, method='spearman')
layout(1)
par(mar=c(5,5,1,1))
plot(cell, patient)

# what about among documents?
coo <- read.csv("results/palm_cooccur.csv", header=F)
names(coo) <- c("document", "term", "count")
length(unique(coo$document))  # number of documents

cell.patient <- sapply(0:max(coo$document), function(i) {
  c(sum(coo$term[coo$document==i] == 0), sum(coo$term[coo$document==i] == 1))
})
foo <- matrix(cell.patient, byrow=T, ncol=2)
fisher.test(table(foo[,1], foo[,2]))

patient.cancer <- sapply(0:max(coo$document), function(i) {
  c(sum(coo$term[coo$document==i] == 1), sum(coo$term[coo$document==i] == 2))
})
foo <- matrix(patient.cancer, byrow=T, ncol=2)
fisher.test(table(foo[,1], foo[,2]))


require(wordcloud, quietly=TRUE)
set.seed(4)
pdf("doc/wordcloud.pdf", height=5, width=5)
wordcloud(index$word[1:80], index$count[1:80], 
          colors=hcl.colors(20, 'Berlin'), random.color=TRUE)
dev.off()

# sparse document-term matrix
ccm <- read.csv("results/palm_cooccur.csv", header=F)
names(ccm) <- c('document', 'term', 'count')
smx <- Matrix::sparseMatrix(i=ccm$document, j=ccm$term, x=1, index1=FALSE)
checksums <- sapply(1:nrow(smx), function(i) {
  sum(which(smx[i,]==1))
})
dups <- which(duplicated(checksums))
length(dups)  # number of duplicate documents

source("scripts/topicspace.R")
#ts <- topicspace(
#  index.path = "results/palm_index.csv",
#  cooccur.path = "results/palm_cooccur.csv",
#  author.path = "results/palm_author.json",
#  max.words=5000, n.comp=3)
#save(ts, file="results/palm_ts_mw5000_nc3.RData")
load("results/palm_ts_mw5000_nc3.RData")


png("doc/topicspace.png", width=6*300, height=4.5*300, res=300)
set.seed(1)
par(family='Palatino', mar=rep(0, 4))
layout(matrix(c(1, 1, 1, 2, 1, 1, 1, 3, 1, 1, 1, 4), byrow=T, nrow=3))
plot(ts, idx=sample(1:1000, 100), text.cex=0.9, xpd=NA,
     text.col=hcl.colors(6, 'Berlin')[2:5])  #gg.rainbow(3, l=50))
title("A", cex.main=2.2, line=-1.6, adj=0, font.main=1)
plot(ts, author=39)
title("B", cex.main=2.2, line=-1.6, adj=0, font.main=1)
plot(ts, author=20)
plot(ts, author=77)
dev.off()


## this takes about a minute with parallel processing
wdist <- get.dist(ts, p=2, mc.cores=24)
#save(wdist, file="results/palm_wdist_L2.RData")
#load("results/palm_wdist.RData")  # L1-norm
load("results/palm_wdist_L2.RData")  # L1-norm

# time complexity (24 cores)
t0 <- Sys.time()
. <- get.dist(ts, p=2, mc.cores=24)
print(Sys.time() - t0)
x <- c(200, 200, 500, 1000, 2000, 2000, 5000, 5000)  # varying max words
y <- c(1.981, 1.963, 3.862, 9.197, 25.49, 25.437, 65.19906, 66.413) # seconds
par(mar=c(5,5,1,1))
plot(unique(x), sapply(split(y, x), mean), type='b')

# 1 core
t0 <- Sys.time()
. <- get.dist(ts, p=2, mc.cores=1)
print(Sys.time() - t0)
x <- c(100, 200, 300, 400, 500, 5000)
y <- c(5.757, 12.161, 19.748, 28.938, 37.782, 429.626)
plot(x, y, type='b')


wmx <- as.matrix(wdist)
i1 <- which(grepl("Poon", colnames(wmx)))
i2 <- which(grepl("Quin", colnames(wmx)))
i3 <- which(grepl("Garcia", colnames(wmx)))
hist(wmx[upper.tri(wmx)], freq=F, main=NA)
par(mfrow=c(1,1), mar=c(5,5,1,1))
plot(density(wmx[upper.tri(wmx)]), bty='n', main=NA,
     xlab="Wasserstein distance")
abline(v=wmx[i1, i2], col='red')  # should be more similar
abline(v=wmx[i1, i3])
abline(v=wmx[i2, i3])



# visualize as a network
source("scripts/visualize.R")

fac <- read.csv("data/palm.tsv", sep='\t')
fac$filename <- paste(fac$lastname, fac$forename, sep="_")
fac$filename <- gsub(" ", "_", fac$filename)
idx <- match(attr(wdist, "Labels"), fac$filename)
fac <- fac[idx, ]  # sort in order of graph data
labels <- fac$initials

#idx <- match(fac$filename, names(by.author))
#fac$n.abs <- sapply(idx, function(i) length(by.author[[i]]))
fac$n.abs <- sapply(fac$filename, function(fn) {
  json <- read_json(paste("data/palm/", fn, '.json', sep=''))
  length(json)
})

require(igraph)

g <- make.knn(wdist, k=3)
set.seed(2)
cl <- cluster_louvain(g)
table(cl$membership)


pdf("doc/knn.pdf", width=5, height=5)
par(mar=c(0,0,0,0), family='Palatino')
layout <- layout_with_fr(g)
plot(g, layout=layout,
     mark.groups=cl, mark.col=hcl.colors(length(cl), "Pastel 1"),
     vertex.shape="none", vertex.label=tolower(labels), vertex.label.color="black",
     vertex.label.cex=1, vertex.size=10,
     edge.width=2, edge.arrow.mode='-')
text(x=0.1, y=0.75, label='1', cex=1.2)
text(x=0.9, y=0.5, label='2', cex=1.2)
text(x=0.4, y=0.3, label='3', cex=1.2)
text(x=-0.1, y=0.2, label='4', cex=1.2)
text(x=-0.65, y=0.4, label='5', cex=1.2)
text(x=0.6, y=-0.2, label='6', cex=1.2)
text(x=-0.1, y=-0.8, label='7', cex=1.2)
text(x=-0.2, y=-0.1, label='8', cex=1.2)
dev.off()

# guide for manual labeling of clusters
cl$names[cl$membership==7]


# sensitivity analysis
num.cl <- sapply(1:1000, function(i) {
  max(cluster_louvain(g)$membership)
})
table(num.cl)


# characteristic words



# https://www.schulich.uwo.ca/pathol//research/areas_of_research/index.html
areas <- c('Bioinformatics and data science\n', 'Cancer biology\n', 
            'Cardiovascular, respiratory health\nand metabolic diseases', 
            'Computational pathology\n', 'Epigenetics\n',
            'Infection, immunity, and\ninflammation',
            'Interdisciplinary health and\neducation',
            'One Health\n', 'Oral biology and medicine\n',
            'Pathobiology of neurologic\ndiseases',
            'Regenerative and transplantation\nmedicine',
            'Test utilization, optimization\nand quality assurance')
abbrv <- c('BD', 'CB', 'CRM', 'CP', 'E', 'III', 'IHE', 'OH', 'OBM', 'PND', 
           'RTM', 'TOQ')
pdf("doc/areas.pdf", width=9, height=12)
par(mfrow=c(4,3), mar=c(2,0,0,0), xpd=NA, family='Palatino')
for (i in 5:16) {
  plot(g, layout=layout, vertex.label='',
       vertex.color=ifelse(is.na(fac[,i]), 'white', 'black'))
  title(main=paste(abbrv[i-4], "   "), line=-3, cex.main=2.7, adj=1, font.main=1)
}
par(xpd=FALSE)
dev.off()


# test assortativity by research area (Table 2)
for (i in 5:16) {
  obs.r <- assortativity(g, !is.na(fac[,i]))  # NA is FALSE
  null.r <- sapply(1:1e6, function(j) {
    p <- sample(fac[,i], nrow(fac), replace=FALSE)
    assortativity(g, !is.na(p))
  })
  print(paste(areas[i-4], obs.r, sum(null.r > obs.r)/1e6))
}




table(replicate(1000, max(cluster_louvain(g)$membership)))

for (i in 1:1000) {
  set.seed(i)
  cl8 <- cluster_louvain(g)
  if (max(cl8$membership)==8) {
    temp <- split(cl8$names, cl8$membership)
    print(temp[[8]])
  }
}
set.seed(i)
cl8 <- cluster_louvain(g)
split(cl8$names, cl8$membership)

# word-to-word cosine distances
n <- length(by.author)
w2w <- matrix(0, nrow=n, ncol=n)
for (i in 1:(n-1)) {
  ai <- by.author[[i]]
  ai2 <- sum(as.integer(ai)*as.integer(ai))
  for (j in 2:n) {
    aj <- by.author[[j]]
    aj2 <- sum(as.integer(aj)*as.integer(aj))
    common <- intersect(names(ai), names(aj))
    aij <- sum(as.integer(ai[common]) * as.integer(aj[common]))
    w2w[i,j] <- w2w[j,i] <- 1-exp(log(aij) - 0.5*log(ai2) - 0.5*log(aj2)) # sqrt(ai2*aj2)
  }
}
wmx <- as.matrix(wdist)  # in case this hasn't been run yet
row.names(w2w) <- row.names(wmx)

par(mfrow=c(1,1), mar=c(5,5,1,1), xpd=F)
plot(density(w2w[upper.tri(w2w)]))
abline(v=w2w[i1, i2])
abline(v=w2w[i1, i3])
abline(v=w2w[i2, i3])

g.w2w <- make.knn(w2w)
set.seed(3)
cl.w2w <- cluster_louvain(g.w2w)
pdf("doc/word2word.pdf", width=5, height=5)
par(mfrow=c(1,1), mar=c(0,0,0,0))
plot(g.w2w, layout=layout_with_fr, 
     mark.groups=cl.w2w, #mark.col=hcl.colors(length(cl.w2w), "Pastel 1"),
     vertex.shape="none", vertex.label=tolower(labels), vertex.label.color="black",
     vertex.label.cex=1, vertex.size=10,
     edge.width=2, edge.arrow.mode='-')
dev.off()



# compare network statistics
summary(degree(g.w2w))
summary(degree(g))

plot(NA, xlim=c(1, 19), ylim=c(-0.5, 0.5), bty='n', xaxt='n',
     xlab="Degree size", yaxt='n', ylab='')
abline(h=0)
tab <- table(factor(degree(g.w2w), levels=1:19))
segments(x0=1:19, y0=0, y1=as.integer(tab)/sum(tab), lwd=5, col='salmon', lend=1)

plot(density(closeness(g.w2w)))

tab <- table(factor(degree(g), levels=1:9))
segments(x0=1:9, y0=0, y1=-as.integer(tab)/sum(tab), lwd=5, 
         col='dodgerblue', lend=1)

x <- closeness(g)
y <- closeness(g.w2w)
plot(x, y, log='xy')
text(x, y, label=cl.w2w$membership)


require(ggfree)
slopegraph(x, y, log='y', color=T)
drop <- which(cl$membership %in% c(1,2,3))
slopegraph(closeness(g.w2w)[-drop], closeness(g)[-drop], color=T)
wilcox.test(closeness(g.w2w)[-drop], closeness(g)[-drop], paired=T)
wilcox.test(closeness(g.w2w), closeness(g), paired=T)

# outcomes are not normally distributed
#require(PairedData)
#Var.test(closeness(g), closeness(g.w2w))

# test for difference in variation by bootstrap resampling
null.dist <- sapply(1:1e5, function(.) {
  idx <- sample(c(FALSE, TRUE), size=length(x), replace=T)
  var(ifelse(idx, x, y)) - var(ifelse(idx, y, x))
})
var(x)-var(y)


pdf("doc/correlation.pdf", width=5, height=5)
par(mar=c(5,5,1,1), family='Palatino')
plot(wmx[1,-1], w2w[1,-1], type='n', xlim=range(wmx), ylim=c(0.2, 1),
     xlab="Wasserstein distance",
     ylab="Word-to-word cosine distance", bty='n', las=1)
for (i in 1:n) {
  lines(smooth.spline(wmx[i,-i], w2w[i,-i], spar=1.), col=rgb(0,0,0,0.1))
}
points(wmx[1,-1], w2w[1,-1], pch=19, cex=0.5, col='darkorange')
lines(smooth.spline(wmx[1,-1], w2w[1,-1], spar=1.0), lwd=2, col='darkorange')
points(wmx[2,-2], w2w[2,-2], pch=19, cex=0.5, col='firebrick')
lines(smooth.spline(wmx[2,-2], w2w[2,-2], spar=1.0), col='firebrick', lwd=2)
points(wmx[3,-3], w2w[3,-3], pch=19, cex=0.5, col='cadetblue')
lines(smooth.spline(wmx[3,-3], w2w[3,-3], spar=1.0), col='cadetblue', lwd=2)
dev.off()

cor.test(as.numeric(wmx[upper.tri(wmx)]), as.numeric(w2w[upper.tri(w2w)]), method='spearman')
require(DescTools)
SpearmanRho(as.numeric(wmx[upper.tri(wmx)]), as.numeric(w2w[upper.tri(w2w)]), 
            conf.level = 0.95)


#plot(g, vertex.shape="circle", vertex.label=NA, 
#     vertex.label.cex=1, vertex.size=sqrt(fac$n.abs), 
#     edge.width=2, edge.arrow.mode='-')


# characterize each cluster
our.words <- function(i, limit=500, send=10) {
  members <- cl$names[cl$membership==i]
  df <- sapply(members, function(nm) {
    sapply(index$word, function(w) {
      count <- by.author[[nm]][[w]]
      ifelse (is.null(count), 0, count)
    })
  })
  counts <- apply(df, 1, sum)
  #counts <- sort(counts, decreasing = TRUE)
  denom <- index$count/sum(index$count)
  resid <- (counts/sum(counts) - denom) / denom
  names(sort(resid[1:limit], decreasing = TRUE)[1:send])  
}

# for Table 2
sapply(1:max(cl$membership), our.words)
table(cl$membership)

############ SCRAP

c1 <- by.author[[n1[1]]]
for (i in 2:length(n1)) {
  
}

(cl$names[cl$membership==1], function(x) {
  temp <- by.author[[x]]
  counts <- as.integer(temp)
  names(counts) <- names(temp)
  
})



# remove duplicate columns from term-document matrix
index <- read.csv("results/palm_index.csv")
sparse <- read.csv("results/palm_cooccur.csv", header=F)
names(sparse) <- c('doc.idx', 'word.idx', 'count')
n.words <- 5000

foo <- sapply(split(sparse[,2:3], sparse$doc.idx), function(x) {
  
})
bar <- unique(foo)

top5 <- sparse[sparse$word.idx < n.words, ]  # most frequent words
totals <- sapply(split(top5$count, top5$doc.idx), sum)
top10 <- unique(top5$doc.idx)[order(totals, decreasing = TRUE)[1:n.docs]]

temp <- sparse[sparse$word.idx < n.words & sparse$doc.idx %in% top10, ]
temp$doc.idx <- as.integer(as.factor(temp$doc.idx))  # renumber
tdm <- matrix(0, nrow=n.words, ncol=n.docs)
rownames(tdm) <- index$word[1:n.words]
colnames(tdm) <- paste("D", 1:n.docs, sep="")
for (row in 1:nrow(temp)) {
  tdm[temp$word.idx[row]+1, temp$doc.idx[row]] <- temp$count[row]
}
foo <- unique(tdm, margin=2)
