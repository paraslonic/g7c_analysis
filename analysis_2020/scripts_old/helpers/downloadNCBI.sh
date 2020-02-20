rsync -av rsync://ftp.ncbi.nlm.nih.gov/genomes/Bacteria --include "*/" --include "Bacteria/$1*/*.fna" --exclude=* .
rsync -av rsync://ftp.ncbi.nlm.nih.gov/genomes/Bacteria_DRAFT --include "*/" --include "Bacteria_DRAFT/$1*/*.fna.tgz" --exclude=* .

find Bacteria -empty -type d -delete
find Bacteria_DRAFT -empty -type d -delete
find Bacteria -mindepth 1 -type d | xargs -L 1 -I '{}' find {} -name "*.fna" | while read i ; do cat "$i" >> `dirname "$i"`.fasta ; done
find Bacteria_DRAFT -mindepth 1 -type d | xargs -L 1 -I '{}' find {} -name "*.fna.tgz" | while read i ; do tar zxvf "$i" -C `dirname "$i"`; done
find Bacteria_DRAFT -mindepth 1 -type d | xargs -L 1 -I '{}' find {} -name "*.fna" | while read i ; do cat "$i" >> `dirname "$i"`.fasta ; done
 
