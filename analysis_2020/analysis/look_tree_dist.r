library("phangorn")

setwd("/data11/bio/PROJECTS/G7C/analysis_2020/")
tree.dir <- "Results/ortho/coreogs_aligned/"
files <- list.files(tree.dir,"*.fasta.treefile")
tree.list <- lapply(files, function(f){
  t <- read.tree(paste0(tree.dir,f))
  t$tip.label <- gsub("_.+","",t$tip.label)
  return(t)
  })

names(tree.list) <- gsub(".fasta.treefile","",files)

class(tree.list)<-"multiPhylo"
M <- path.dist(tree.list)
M <- as.matrix(M)

ot <- read.delim("Results/ortho/ortho_table.txt")
ognames <- data.frame(og=gsub(":","",ot$id), name=ot$product)

qu <- ognames[match(names(tree.list),ognames$og),]
all(names(tree.list) == qu$og)
rownames(M) <- qu$name
colnames(M) <- qu$name
pheatmap(M, trace = "none", fontsize_row = 6,fontsize_col = 6)

pheatmap(Muu)
