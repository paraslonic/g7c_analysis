setwd("/data7a/bio/operonTravel/salmonella/tmp/")

library("seqinr")
library("stringr")
library("ape")
library("phangorn")
library("gplots")
library("abind")

### CORE GENES
dist.list = list()
FILES = list.files("coreogaligns/",pattern = ".*fasta")
FILES.count = length(FILES)
FILES.count = 100
for(f in 1:FILES.count)
{
  if(f%%100==0) print(f)
  tryCatch({
    name = sub(".fasta", "", FILES[f])
    #print(name)
    A = read.dna(paste0("coreogaligns/",FILES[f]), format = "fasta")
    D = dist.dna(A, model="K80",pairwise.deletion=TRUE)
    #D = dist.ml(A,model="Blosum62")
    D = as.matrix(D)
    if(max(D) == 0) { next }
    colnames(D) = sub("\\|.+","",colnames(D))
    rownames(D) = sub("\\|.+","",rownames(D))
    .order = order(rownames(D))
    D = D[.order,.order]
    #heatmap.2(as.matrix(D),  trace ="none", cexRow = 0.8, margins = c(6,12), cex = 0.8)
    dist.list[[name]] = as.matrix(D)}, error=function(e){cat("ERROR in",name, "\n")})
}
save(dist.list, file="dist.list.rdata")

dists.3d <- abind(dist.list, along=3)
dists.mean <- apply(dists.3d, c(1,2), mean)

save(dists.mean, file="../core.dists.mean.nuc.rdata")

hist((unlist(dists.mean)), col ="dodgerblue", breaks = 50)
hist(log(unlist(dists.mean)), col ="dodgerblue")

#heatmap.2(as.matrix(dists.mean),trace="none")

##### OPERON GENES

operon.dist.list = list()
FILES = list.files("operon_genes_aligned/",pattern = ".*fasta")
FILES.count = length(FILES)
#FILES.count = 100
for(f in 1:FILES.count)
{
  if(f%%100==0) print(f)
  tryCatch({
    name = sub(".fasta", "", FILES[f])
    #print(name)
    A = read.dna(paste0("operon_genes_aligned/",FILES[f]), format = "fasta")
    if(length(labels(A)) > 68) { next }
    D = dist.dna(A, model="K80",pairwise.deletion=TRUE)
    #D = dist.ml(A,model="Blosum62")
    D = as.matrix(D)
    if(max(D) == 0) { next }
    colnames(D) = sub("\\|.+","",colnames(D))
    rownames(D) = sub("\\|.+","",rownames(D))
    .order = order(rownames(D))
    D = D[.order,.order]
    #heatmap.2(as.matrix(D),  trace ="none", cexRow = 0.8, margins = c(6,12), cex = 0.8)
    operon.dist.list[[name]] = as.matrix(D)}, error=function(e){cat("ERROR in",name, "\n")})
}
save(operon.dist.list, file="operon.dist.list.rdata")

#### CORRELATION

og.dists = operon.dist.list[[25]]

.selection = which(colnames(dists.mean) %in% colnames(og.dists))
.dists.mean = dists.mean[.selection, .selection]

all(colnames(.dists.mean) == colnames(og.dists))

plot((unlist(.dists.mean)), (unlist(og.dists)), pch = 16, cex = 1.2, col = rgb(0.3,0.3,0.3,0.3),
     xlab="core genes", ylab = "og", main = "DNA sequence distance", cex.lab=1.4)

plot((unlist(.dists.mean)), (unlist(og.dists)), pch = 16, cex = 1.2, col = rgb(0.3,0.3,0.3,0.3),
     xlab="core genes", ylab = "og", main = "DNA sequence distance", cex.lab=1.4, xlim = c(0, 0.02), ylim = c(0,0.02))


dists.mean.nuc = dists.mean
cor.list = list()
for(q in names(dist.list.core)){
  tryCatch({
  print(q)
  cor.list[[q]] = cor(vec(dist.list.core[[q]]), vec(dists.mean.nuc), method = "spearman")
  #plot(vec(q), vec(dists.mean))
  }, error=function(e){cat("ERROR in",name, "\n")})
}
hist(unlist(cor.list), col = "orange")
save(cor.list, file="corlist.rdata")

# compare with core -------------------------------------------------------

rownames(core.dist) = colnames(core.dist)
.selection = which(colnames(core.dist) %in% colnames(dists.mean))
core.dist = core.dist[.selection,.selection]
.order = order(colnames(core.dist))
core.dist = core.dist[.order,.order]

.selection = which(colnames(dists.mean) %in% colnames(core.dist))
dists.mean = dists.mean[.selection,.selection]


plot((unlist(dists.mean)), (unlist(core.dist)), pch = 16, cex = 1.2, col = rgb(0.3,0.3,0.3,0.3),
     xlab="pdu operon", ylab = "core genes", main = "DNA sequence distance", cex.lab=1.4)

heatmap.2(as.matrix(core.dist),  trace ="none", cexRow = 0.8, margins = c(6,12), cex = 0.8)

heatmap.2(as.matrix(dists.mean),  trace ="none", cexRow = 0.8, margins = c(6,12), cex = 0.8)

####
