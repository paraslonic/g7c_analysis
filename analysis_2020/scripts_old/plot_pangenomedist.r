library(vegan)
library("gplots")

gene.table = read.delim("ortho_table.txt", check.names = F)
gene.table = gene.table[,-c(1,2,match("strains", colnames(gene.table)))]
gene.bool = gene.table
gene.bool[gene.bool > 1] = 1

### plot pangenome distances

pg.dist = as.matrix(vegdist(t(gene.bool), method="bray"))
pdf("../PangenomeProfileDistances.pdf")
heatmap.2(as.matrix(pg.dist),trace="none",margin=c(11,11), main = "Pangenome profile distance",
          dendrogram="row", keysize=1.5, cexRow = 0.3, cexCol = 0.3)
dev.off()
