mkdir -p ../tmp/operon_genes
#perl splitToOg.pl ../OrthologousGroups.txt ../tmp/allgenes.fasta ../tmp/operon_genes spi_1.og

echo "aligning core genes "
mkdir -p ../tmp/operon_genes_aligned
for f in ../tmp/operon_genes/*.fasta
do
	name=$(basename $f .fasta)
	echo "./muscle -in ../tmp/operon_genes/$name.fasta -out ../tmp/operon_genes_aligned/$name.fasta -quiet" | qsub -cwd
	printf "."
done

