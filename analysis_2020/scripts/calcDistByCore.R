
library("stringr")
library("ape")
library("phangorn")
library("gplots")
library("abind")

### CORE GENES
dist.list = list()
aligned_dir <- "Results/ortho/coreogs_aligned/"
out_dir <- "Results/rdata/"

FILES = list.files(aligned_dir,pattern = ".*fasta")
FILES.count = length(FILES)


for(f in 1:FILES.count)
{
  if(f%%100==0) print(f)
  tryCatch({
    name = sub(".fasta", "", FILES[f])
    A = read.dna(paste0(aligned_dir,FILES[f]), format = "fasta")
    D = dist.dna(A, model="K80",pairwise.deletion=TRUE)
    D = as.matrix(D)
    if(max(D) == 0) { next }
    colnames(D) = sub("\\|.+","",colnames(D))
    rownames(D) = sub("\\|.+","",rownames(D))
    .order = order(rownames(D))
    D = D[.order,.order]
    #heatmap.2(as.matrix(D),  trace ="none", cexRow = 0.8, margins = c(6,12), cex = 0.8)
    dist.list[[name]] = as.matrix(D)},
  error=function(e){cat("ERROR in",name, "\n")})
}

dists.3d <- abind(dist.list, along=3)
dists.mean <- apply(dists.3d, c(1,2), mean)

system(paste("mkdir",out_dir))
save(dist.list, dists.mean, file=paste0(out_dir,"core.dists.rdata"))

