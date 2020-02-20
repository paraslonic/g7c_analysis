library("ape")
library("gplots")

setwd("/data5/bio/runs-danio/letarov_fag/annotationFAG/allffn")

core_gp = read.delim("gp_list",head = FALSE)

for(f in core_gp$V1){
  data <- read.dna(paste0("gp_ffn_aligned/",f,".ffn"), format = "fasta")
  M = dist.dna(data)
  M.list[[f]] = M
}

n = length(core_gp$V1)

pdf("bygroup_all.pdf")

for(i in 1:n){
  if( length(M.list[[i]]) == 0 | max(M.list[[i]]) == 0 | is.nan(max(M.list[[i]]))) { next }
  tryCatch({
    heatmap.2(as.matrix(M.list[[i]]),trace="none", margins = c(12,12), main = core_gp$V1[i])
  }, error=function(e){ print(e)}) 
  print(core_gp$V1[i])
  }

dev.off()
