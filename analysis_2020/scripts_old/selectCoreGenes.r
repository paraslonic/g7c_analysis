library(data.table)
length_maxvar <- 0.4

ogtab <- fread("ortho_table.txt")
strain.count <- ncol(ogtab) - 3

og.core <- subset(ogtab, ogtab$strains == strain.count)
if(nrow(og.core) < 1) stop ("STOP: core genes less then 1")
print(paste("total core genes:",nrow(og.core)))

og.core <- subset(og.core, rowSums(og.core[,-c(1,2,ncol(og.core)), with=FALSE]) == strain.count)
print(paste("one copy core genes:",nrow(og.core)))

## add length and filter
ognames <- fread("ortho_table_names.txt")

get_length <-  function(x) {
  if(x == "") return(NA)
  parts <- strsplit(x, "\\|")[[1]]; 
  if(length(parts) != 6) return(NA)
  end <- as.integer(parts[6])
  start <- as.integer(parts[5])
  return (abs(end-start))
}

ognames = subset(ognames, ognames$id %in% og.core$id)
oglen <- apply(ognames[,-c(1,2),with=FALSE], c(1,2), get_length)

og.maxlen <- apply(oglen,1, max, na.rm = TRUE)
og.minlen <- apply(oglen,1, min, na.rm = TRUE)
og.medlen <- apply(oglen,1, median, na.rm = TRUE)

og.core <- subset(og.core, 
                  (og.medlen-og.minlen)/og.medlen < length_maxvar &
                  (og.maxlen - og.medlen)/og.medlen < length_maxvar)

print(paste("one copy good core genes",nrow(og.core)))

write.table(og.core$id,"../tmp/coreog", quote = FALSE, row.names=FALSE, col.names=FALSE)

