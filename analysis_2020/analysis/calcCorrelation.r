library("ape")
library("abind")
library("gplots")
library("RColorBrewer")
library("phytool")
library("phangorn")

out_dir <- "Results/rdata/"
load(paste0(out_dir,"core.dists.rdata"))

M.list = dist.list

sortmat = function(x) {x = as.matrix(x); x = x[order(rownames(x)),]; x = x[,order(colnames(x))]; as.data.frame((x)) }
M.list = lapply(M.list,sortmat)


n = length(M.list)
C  = matrix(0, ncol = n, nrow = n)
Csp  = matrix(0, ncol = n, nrow = n)

for(f in 1:n){
  for(g in 1:n){
    Csp[f,g] = cor(unlist(M.list[[f]]),unlist(M.list[[g]]), method="spearman")   
    C[f,g] = cor(unlist(M.list[[f]]),unlist(M.list[[g]]))   
  }
}

C[is.na(C)] = 1
Csp[is.na(Csp)] = 1
View(C)

ot <- read.delim("Results/ortho/ortho_table.txt")
ognames <- data.frame(og=gsub(":","",ot$id), name=ot$product)

qu <- ognames[match(names(M.list),ognames$og),]


rownames(C) = qu$name
colnames(C) = rownames(C)
rownames(Csp) = rownames(C)
colnames(Csp) = rownames(C)


pdf("chooseColors.pdf")

myPalette <-  colorRampPalette(c("#e34a33", "#e34a33","#e34a33", "#ffeda0", "#3182bd"))(n = 200)
heatmap.2(Csp, trace="none",  hclustfun = function(x) hclust(x, method = 'ward.D'), 
          col =myPalette, cexRow = 0.5, cexCol = 0.5, dendrogram = "row",
          margins = c(14,14))


myPalette <-  colorRampPalette(c("white","white","gray95",  "#3182bd"))(n = 200)
heatmap.2(Csp, trace="none",  hclustfun = function(x) hclust(x, method = 'ward.D'), 
          col =myPalette, cexRow = 0.5, cexCol = 0.5, dendrogram = "row")

dev.off()

pdf("correlation.pdf")
heatmap.2(C, trace="none", main = "pearson", hclustfun = function(x) hclust(x,method = 'ward.D'))
heatmap.2(Csp, trace="none", main = "spearman", hclustfun = function(x) hclust(x, method = 'ward.D'), col = myPalette)
dev.off()


colorRampPalette(colors = c("white","blue"))

hm <- heatmap.2(Csp, trace = "none", hclustfun = function(x) hclust(x,method = 'ward.D'))
hc <- as.hclust( hm$rowDendrogram )
plot(hc)
rect.hclust(hc, k=3)
groups = cutree(hc, k = 3)
s = silhouette(groups, Csp)
mean(s)

g1 = names(groups[groups == 1])
g2 = names(groups[groups == 2])
g3 = names(groups[groups == 3])

names(M.list) = gsub("]","",names(M.list))
rename = function(d)
{
  colnames(d) = gsub("\\s.+","", colnames(d))
  colnames(d)[grepl("NC_015933", colnames(d))] = "G7C"
  rownames(d) = colnames(d)
  return(d)
}


G1 = abind(M.list[g1], along = 3)
G2 = abind(M.list[g2], along = 3)
G3 = abind(M.list[g3], along = 3)
G1 = rename(G1)
G2 = rename(G2)
G3 = rename(G3)

dG1 <- apply(G1, c(1,2), mean)
dG2 <- apply(G2, c(1,2), mean)
dG3 <- apply(G3, c(1,2), mean)

colnames(dG1)
colnames(dG2)

pdf("three_clusters.pdf")
heatmap.2(as.matrix(dG1),trace="none", main = "group 1", margins = c(10,10))
heatmap.2(as.matrix(dG3),trace="none", main = "group 2", margins = c(10,10))
heatmap.2(as.matrix(dG2),trace="none", main = "group 3", margins = c(10,10))
dev.off()

write.tree(midpoint(bionj(as.matrix(dG1))), file = "g1.nwk")
write.tree(midpoint(bionj(as.matrix(dG2))), file="g2.nwk")
write.tree(midpoint(bionj(as.matrix(dG3))), file="g3.nwk")




pdf("bygroup.pdf")

for(i in 1:n){
  if(max(M.list[[i]]) == 0) { next }
  tryCatch({
    heatmap.2(as.matrix(M.list[[i]]),trace="none", margins = c(12,12), main = core_gp$V1[i])
  }, error=function(e){ print(e)}) 
  print(core_gp$V1[i])
  }

dev.off()

