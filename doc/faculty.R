## To generate input files:
## python scripts/analyze.py data/acb data/biochem data/epibio data/mbp \ 
##   data/mimm data/palm data/ppharm --counts results/b7_counts.json \
##   --matrix results/b7_cooccur.csv --index results/b7_index.csv

setwd("~/git/tragula")

files <- c("data/acb.tsv", "data/biochem.tsv", "data/epibio.tsv", 
           "data/mbp.tsv", "data/mimm.tsv", "data/palm.tsv", 
           "data/ppharm.tsv")

fac <- do.call(rbind, lapply(files, function(f) {
  dept <- gsub("\\.tsv", "", basename(f))
  tsv <- read.csv(f, sep='\t')
  tsv <- tsv[,1:4]
  names(tsv) <- c("lastname", "forename", "initials", "query")
  tsv$affil <- dept
  tsv
}))

# make some more convenient labels
fac$labels <- paste(substr(fac$forename, 1, 1),
                    substr(fac$lastname, 1, 4), sep="")
#paste(fac$initials, substr(fac$lastname, 2, 2), sep="")
fac$filename <- gsub(" ", "_", paste(fac$lastname, fac$forename, sep="_"))

# UPDATE: manually added centre affiliations
fac <- read.csv("data/faculty-centres.csv")

require(jsonlite)
fac$n.abs <- sapply(1:nrow(fac), function(i) {
  fn <- fac$filename[i]
  af <- fac$affil[i]
  fn <- paste("data/", af, "/", fn, '.json', sep='')
  if (file.exists(fn)) {
    json <- read_json(fn)
    return(length(json))
  }
  NA
})

# Median abstracts per member
sapply(split(fac$n.abs, fac$affil), median)

sapply(split(fac$n.abs, fac$affil), sum)


# generate word cloud (Figure S2)
index <- read.csv("results/b7_index.csv")

require(wordcloud, quietly = TRUE)
set.seed(70)
pdf("doc/faculty-cloud.pdf", width=6, height=6)
wordcloud(index$word[1:200], index$count[1:200], 
          colors=hcl.colors(20, 'Berlin'), random.color=TRUE)
dev.off()

source("scripts/topicspace.R")
source("scripts/visualize.R")

## warning, screening for duplicates involves a sparse-to-dense conversion 
## of a matrix - in this case, it consumes 6.4 GB of RAM.
set.seed(1)
ts <- topicspace(
  index.path = "results/b7_index.csv",
  cooccur.path = "results/b7_cooccur.csv",
  author.path = "results/b7_counts.json",
  max.words=10000, n.comp=3)
#save(ts, file="results/b7_ts_mw1e4_n3.RData")
load("results/b7_ts_mw1e4_n3.RData")

idx <- match(names(ts$by.author), fac$filename)
fac <- fac[idx, ]


t0 <- Sys.time()  # timing (about half an hour)
wdist <- get.dist(ts, p=2, mc.cores=23)
print(Sys.time() - t0)
#save(wdist, file="results/b7_wdist_L2.RData")
#load("results/b7_wdist_10E4.RData")
load("results/b7_wdist_L2.RData")

set.seed(1)
g <- make.knn(wdist, k=3) #, cutoff=0.75)
layout <- layout_with_kk(g)


# plot the k-nearest neighbors (k-NN) graph
shortened <- list(acb='a', biochem='b', epibio='e', mbp='mb', mimm='mi', 
                  palm='pl', ppharm='ph')
labels <- fac$label
groups <- sapply(fac$affil, function(x) shortened[[x]])
pal <- hcl.colors(n=7, palette="Dark3")
colors <- pal[as.integer(as.factor(groups))]

pdf("doc/faculty.pdf", width=8, height=8)
par(mar=rep(0, 4), family='Palatino', mfrow=c(1,1))
plot(g, layout=layout,
     vertex.shape="circle", vertex.label=groups, vertex.label.cex=0.7, 
     vertex.size=sqrt(fac$n.abs)/2+3, 
     vertex.color=colors, vertex.frame.color='white', vertex.label.color='white',
     edge.width=2, edge.arrow.mode='-')
dev.off()


# show k-NN graphs for other values of k
graphs <- list('2'=NA, '3'=g, '4'=NA, '5'=NA, '6'=NA)

pdf("doc/other-k.pdf", width=8, height=8)
par(mar=rep(0, 4), family='Palatino', mfrow=c(2,2))

for (k in c(2, 4, 5, 6)) {
  g <- make.knn(wdist, k=k)
  graphs[[as.character(k)]] <- g
  layout <- layout_with_kk(g)
  plot(g, layout=layout,
       vertex.shape="circle", vertex.label=groups, vertex.label.cex=0.5, 
       vertex.size=sqrt(fac$n.abs)/2+3, 
       vertex.color=colors, vertex.frame.color='white', vertex.label.color='white',
       edge.width=2, edge.arrow.mode='-')
  title(main=paste(" k=",k), font.main=1, adj=0, line=-1)
}
abline(h=par('usr')[4], xpd=NA, lwd=0.5)
abline(v=par('usr')[1], xpd=NA, lwd=0.5)
dev.off()


# generate data for Supplementary Table (other-k)
types <- as.integer(as.factor(groups))
for (k in 1:6) {
  g <- make.knn(wdist, k=k)
  asrt <- assortativity_nominal(g, types, directed=FALSE)
  cl <- cluster_louvain(g)
  ari <- adjustedRandIndex(fac$affil, cl$membership)
  t0 <- transitivity(g)
  temp <- sapply(unique(fac$affil), function(a) {
    transitivity(g, type='localundirected', vids=fac$filename[fac$affil==a])
  })
  t1 <- mean(sapply(temp, function(x) mean(x[!is.nan(x)])))
  print(paste(round(c(asrt, ari, t0, t1, (t1-t0)/t0), 3)))
}

# evaluate assortativity coefficient for departments
g <- graphs[['3']]
types <- as.integer(as.factor(groups))
assortativity_nominal(g, types, directed=FALSE)
null.asrt <- sapply(1:1e5, function(i) {
  permut <- sample(types, length(types))
  assortativity_nominal(g, permut, directed=FALSE)  
})


# data for Table 3
for (group in unique(groups)) {
  types <- (groups==group)+1
  obs.asrt <- assortativity_nominal(g, types, directed=FALSE)
  null.asrt <- sapply(1:1e5, function(.) {
    permut <- sample(types, length(types))
    assortativity_nominal(g, permut, directed=FALSE)
  })
  print(paste(group, obs.asrt, sum(null.asrt > obs.asrt)))
}


# examine clustering coefficients for departments in faculty graph
depts <- c('acb', 'biochem', 'epibio', 'mbp', 'mimm', 'ppharm', 'palm')
close.cent <- sapply(depts, function(a) { 
  closeness(g, v=fac$filename[fac$affil==a], normalized=T) 
  })

sapply(close.cent, mean)

sapply(depts, function(a) {
  edge_density(subgraph(g, fac$filename[fac$affil==a]))
})

# path lengths from a member of cluster to every other member (n choose 2)
d <- distances(g, algorithm='unweighted')
within.d <- lapply(depts, function(a) {
    i <- which(fac$affil==a)
    dd <- d[i,i]
    dd[upper.tri(dd)]
    })
names(within.d) <- depts
sapply(within.d, mean)

null.d <- sapply(1:1e5, function(.) {
  # random permutation of dept labels in graph
  prm <- sample(fac$affil, nrow(fac), replace=FALSE)
  sapply(depts, function(a) {
    i <- which(prm==a)
    dd <- d[i,i]
    mean(dd[upper.tri(dd)])  # mean path length
  })
})
apply(null.d, 1, summary)

sapply(1:7, function(i) {
  sum(null.d[i,] < mean(within.d[[depts[i]]]))
})


pdf("doc/graphstats.pdf", width=8, height=4)
par(mar=c(5,4,1,1), mfrow=c(1,2))
boxplot(close.cent, col=hcl.colors(7, 'Dark3'), frame=F, 
        ylab='Closeness centrality', las=2, cex.axis=0.8)
boxplot(within.d, names=depts, col=hcl.colors(7, 'Dark3'), 
        ylab="Within cluster distance", las=2, frame=F, cex.axis=0.8)
dev.off()




### LOUVAIN CLUSTERING ###

# random restarts to find best clustering
clusterings <- sapply(1:1000, function(i) {
  set.seed(i)
  cluster_louvain(g)
})
x <- sapply(clusterings, function(cl) max(cl$modularity))
cl <- clusterings[[which.max(x)]]

set.seed(2)
cl <- cluster_louvain(g)

# examine clusters at faculty level
pdf("doc/faculty-clusters.pdf", width=8, height=8)
par(mar=c(0,0,0,0), mfrow=c(1,1))
plot(g, layout=layout, mark.groups = cl, 
     vertex.label=groups, #fac$initials, 
     vertex.shape="circle", vertex.label.cex=0.7,
     vertex.size=sqrt(fac$n.abs)/2+3, vertex.color='white',
     edge.width=2, edge.arrow.mode='-')
dev.off()

# concordance between academic departments and clusters?
table(fac$affil, cl$membership)

# stability of clustering
table(sapply(clusterings, function(cl) max(cl$membership)))



require(mclust)
obs.ARI <- adjustedRandIndex(fac$affil, cl$membership)  
null.ARI <- sapply(1:1e5, function(i) {
  p <- sample(cl$membership, nrow(fac), replace=FALSE)
  adjustedRandIndex(fac$affil, p)
})
table(null.ARI < obs.ARI)
obs.ARI
summary(null.ARI)


# characterize each cluster
by.author <- read_json("results/b7_counts.json")
index <- read.csv("results/b7_index.csv")
our.words <- function(i, limit=500, send=10) {
  members <- cl$names[cl$membership==i]
  df <- sapply(members, function(nm) {
    sapply(index$word, function(w) {
      count <- by.author[[nm]][[w]]
      ifelse (is.null(count), 0, count)
    })
  })
  counts <- apply(df, 1, sum)
  denom <- index$count/sum(index$count)
  resid <- (counts/sum(counts) - denom) / denom
  names(sort(resid[1:limit], decreasing = TRUE)[1:send])  
}
sapply(1:max(cl$membership), our.words)

table(cl$membership, fac$affil)



# show centre affiliations (omit CQuInS, only one member)
centres <- c('AFRI', 'BJI', 'BMI', 'CERI', 'CSFM', 'CTCR', 'WIN', 'WI3C')

pdf("doc/centres.pdf", width=8, height=8)
par(mfrow=c(3, 3), mar=rep(0,4), family="Palatino")
for (centre in centres) {
  o <- order(fac[[centre]])
  p <- rank(fac[[centre]], ties.method='first')
  g2 <- permute(g, p)
  l2 <- layout[o, ]
  fac2 <- fac[o, ]
  plot(g2, layout=l2, 
       vertex.shape="circle", vertex.label=NA,
       vertex.size=sqrt(fac2$n.abs)/2+3, 
       vertex.color=ifelse(fac2[[centre]], 'black', 'white'), 
       vertex.frame.color=ifelse(fac2[[centre]], 'white', 'grey50'),
       edge.width=2, edge.arrow.mode='-')  
  text(x=par('usr')[2]*0.75, y=par('usr')[3]*0.75, label=centre, 
       cex=2, adj=1, xpd=NA)
}
dev.off()


# for Table 3
result <- data.frame(centre=centres, obs=numeric(8),
                     min=numeric(8), max=numeric(8), pval=numeric(8))
for(i in 1:8) {
  centre <- centres[i]
  types <- as.integer(fac[[centre]])+1L  # next function requires 1-index
  obs.asrt <- assortativity_nominal(g, types, directed=FALSE)
  null.asrt <- sapply(1:1e5, function(i) {
    permut <- sample(types, length(types))
    assortativity_nominal(g, permut, directed=FALSE)  
  })
  result[i, ] <- list(centre, obs.asrt, min(null.asrt), max(null.asrt), sum(null.asrt > obs.asrt)/1e5)
}

#require(knitr)
#kable(result, format="latex", booktabs=T)

tab <- table(fac$WIN, cl$membership)
chisq.test(tab)$stdres
2*(1-pnorm(7.4261539))
2*(1-pnorm(5.3741853))

# drill down on CSFM membership
members <- cl$names[fac$CSFM]
df <- sapply(members, function(nm) {
  sapply(index$word, function(w) {
    count <- by.author[[nm]][[w]]
    ifelse (is.null(count), 0, count)
  })
})
counts <- apply(df, 1, sum)
denom <- index$count/sum(index$count)
resid <- (counts/sum(counts) - denom) / denom
names(sort(resid[1:500], decreasing = TRUE)[1:10])  

members <- cl$names[
  fac$affil=='epibio' & 
    fac$initials %in% c('SA', 'NR', 'KKA', 'SS', 'PW', 'AD')]
    #fac$initials %in% c('YHC', 'OEG', 'EK', 'IK', 'GZ', 'JJK')]
members <- members[-5]
