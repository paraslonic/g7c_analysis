library(ape)
library(phangorn) 

strains_num <-  read.delim("../tmp/strains_num", head = F)
colnames(strains_num) <- c("strain","num")

dnml.tree <- read.tree("../tmp/dnaml.tree")
dnml.tree$tip.label <-  strains_num$strain[as.integer(dnml.tree$tip.label)]
write.tree(dnml.tree, file = "../CoreGenesMLTree.nwk");

source("ecoliGroups.r")

tree <- dnml.tree
# remove redundant (too similar) strains
tree$tip.label = gsub("__","_", tree$tip.label)
tree = drop.tip(tree, redundant)

strains.long = tree$tip.label

# colors
rce.colors = "black"
litc.colors = "red"
path.colors = "pink"
nonpath.colors = "darkolivegreen3"

count = length(tree$tip.label)
colors = rep("gray100", count)
colors[grep("RCE", tree$tip.label)] = rce.colors
colors[tree$tip.label %in% pathogenic ] = path.colors
colors[tree$tip.label %in% healthy ] = nonpath.colors
colors[tree$tip.label %in% crohns.ref ] = litc.colors

tree$tip.label = sub("Escherichia_coli_", "", tree$tip.label)
tree$tip.label = sub("_uid.+", "", tree$tip.label,perl = T)


## add plasmid labels
pjj = strains.long %in% have.jj
plf = strains.long %in% have.pLF82
tree$tip.label[pjj] = paste(tree$tip.label[pjj], "(pJJ)")
tree$tip.label[plf] = paste(tree$tip.label[plf], "(pLF82)")

### P L O T   T R E E 

pdf("../Figure3.pdf")
par(mar = c(1,1,1,1))
rech = 0.4
plot(midpoint(tree), tip.color = "gray30", cex = 0.8,  
     label.offset = 0.0025,edge.lty = 1, edge.width = 1.4)

## add rectangles
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
tip <- 1:lastPP$Ntip
XX <- lastPP$xx[tip]
YY <- lastPP$yy[tip]
rect(XX + 0.0005,YY-rech,XX + 0.002, YY+rech, lwd = 1, col = colors, border = "gray")

## scale bar
add.scale.bar(x = 0.0, y = 30, cex = 0.6,  col = "gray40",lcol="gray40")
legend(-0.002,18, c("Crohn","CrohnLit","Pathogenic","Non pathogenic","Other"), fill=  c(rce.colors, litc.colors, path.colors, nonpath.colors, "gray90"), cex = 0.8,  bty = "n", lwd = 0, lty = 0, col="white", 
       y.intersp=1.2,  x.intersp = -0.5)

## add plasmid symbols
pjj = strains.long %in% have.jj
plf = strains.long %in% have.pLF82
dev.off()


