#mkdir -p ../tmp
#mkdir -p ../faa
#mkdir -p ../fna
#
#echo "converting genbank files "
#
#for f in ../gb/*.gbf
#do
#	name=$(basename $f .gbf)
##	perl gb2faa.pl $f > ../faa/$name.fasta
#	perl gb2fna.pl $f > ../fna/$name.fasta
#	printf "."
#done
#echo 
#
#echo "selecting core genes"
#
#sed -i 's/://' ../OrthologousGroups.txt
#perl namedGroups2table.pl ../OrthologousGroups.txt
#
##Rscript plot_pangenomedist.r
#
#Rscript selectCoreGenes.r
#
mkdir -p ../tmp
#cat ../fna/*.fasta > ../tmp/allgenes.fasta

rm -r -f ../tmp/ogfasta
mkdir -p ../tmp/ogfasta 
perl splitToOg.pl ../OrthologousGroups.txt ../tmp/allgenes.fasta  ../tmp/ogfasta og_of_50

echo "aligning core genes "
rm -r -f ../tmp/ogaligns
mkdir -p ../tmp/ogaligns
for f in ../tmp/ogfasta/*.fasta
do
	name=$(basename $f .fasta)
	echo "./muscle -in ../tmp/ogfasta/$name.fasta -out ../tmp/ogaligns/$name.fasta -quiet" | qsub -cwd
	printf "."
done
