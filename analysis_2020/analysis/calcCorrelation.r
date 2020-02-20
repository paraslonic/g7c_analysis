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
rownames(C) = names(M.list)
colnames(C) = rownames(C)
rownames(Csp) = rownames(C)
colnames(Csp) = rownames(C)


# rownames(Csp) = gsub("]","",rownames(Csp))
# colnames(Csp) = gsub("]","",colnames(Csp))

pdf("chooseColors.pdf")

myPalette <-  colorRampPalette(c("#e34a33", "#e34a33","#e34a33", "#ffeda0", "#3182bd"))(n = 200)
heatmap.2(Csp, trace="none",  hclustfun = function(x) hclust(x, method = 'ward.D'), 
          col =myPalette, cexRow = 0.5, cexCol = 0.5, dendrogram = "row")


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


# *** ---------------------------------------------------------------------

getDist1234 = function(gene)
{
  d = M.list[[gene]] 
  colnames(d) = gsub("\\s.+","", colnames(d))
  colnames(d)[grepl("NC_015933", colnames(d))] = "G7C"
  rownames(d) = colnames(d)
  d1 = mean(unlist(d["10GB",c("G8C","G7C", "Alt63")]))
  d2 = mean(unlist(d["10GB",c("N4G2", "N4Sz33")]))
  d3 = mean(unlist(d["St10y",c("G8C","G7C", "Alt63")]))
  d4 = mean(unlist(d["St10y",c("N4G2", "N4Sz33")]))
  return(c(d1,d2,d3,d4))  
}

getDist1234. = function(gene)
{
  d = M.list[[gene]] 
  colnames(d) = gsub("\\s.+","", colnames(d))
  colnames(d)[grepl("NC_015933", colnames(d))] = "G7C"
  rownames(d) = colnames(d)
  d1 = mean(unlist(d["10GB",c("G8C")]))
  d2 = mean(unlist(d["10GB",c("N4Sz33")]))
  d3 = mean(unlist(d["St10y",c("G8C")]))
  d4 = mean(unlist(d["St10y",c("N4Sz33")]))
  return(c(d1,d2,d3,d4))  
}


L = list()
for (g in as.character(core_gp$V1)){
  print(g)
  L[[g]] = getDist1234.(g)
}

dist.tab = do.call(rbind, L)

write.table(dist.tab, "dist1234_krai.txt", row.names=FALSE, quote=FALSE)

### 

getDistAll.g8c = function(gene)
{
  d = M.list[[gene]] 
  colnames(d) = gsub("\\s.+","", colnames(d))
  colnames(d)[grepl("NC_015933", colnames(d))] = "G7C"
  rownames(d) = colnames(d)
  res = sapply(rownames(d), function(x){
    d1 = mean(unlist(d[x,c("G8C")]))
    return(d1)
    })
  return(res)  
}

getDistAll.g7c = function(gene)
{
  d = M.list[[gene]] 
  colnames(d) = gsub("\\s.+","", colnames(d))
  colnames(d)[grepl("NC_015933", colnames(d))] = "G7C"
  rownames(d) = colnames(d)
  res = sapply(rownames(d), function(x){
    d1 = mean(unlist(d[x,c("G7C")]))
    return(d1)
  })
  return(res)  
}


getDistAll.n4 = function(gene)
{
  d = M.list[[gene]] 
  colnames(d) = gsub("\\s.+","", colnames(d))
  colnames(d)[grepl("NC_015933", colnames(d))] = "G7C"
  rownames(d) = colnames(d)
  res = sapply(rownames(d), function(x){
    d1 = mean(unlist(d[x,c("N4Sz33")]))
    return(d1)
  })
  return(res)  
}

L = list()
for (g in as.character(core_gp$V1)){
  print(g)
  L[[g]] = getDistAll.g7c(g)
}

dist.tab = do.call(rbind, L)
write.table(dist.tab, "distAll_g7c.txt", row.names=FALSE, quote=FALSE)

L = list()
for (g in as.character(core_gp$V1)){
  print(g)
  L[[g]] = getDistAll.n4(g)
}

dist.tab = do.call(rbind, L)
write.table(dist.tab, "distAll_n4.txt", row.names=FALSE, quote=FALSE)

dist.tab = t(dist.tab)
plot(log(dist.tab[4,]),type='l',ylim = c(0, -12))
for (i in 1:nrow(dist.tab)){
  lines(log(0.001+dist.tab[i,]), col = 1+i, lwd = 2)  
}
legend("top",row.names(dist.tab), col=1+1:nrow(dist.tab), lwd = 1)

