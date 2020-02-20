mkdir -p ../tmp
mkdir -p ../faa
mkdir -p ../fna

echo "converting genbank files "

for f in ../gb/*.gbf
do
	name=$(basename $f .gbf)
#	perl gb2faa.pl $f > ../faa/$name.fasta
#	perl gb2fna.pl $f > ../fna/$name.fasta
	printf "."
done
echo 

echo "selecting core genes"

sed -i 's/://' ../OrthologousGroups.txt
perl namedGroups2table.pl ../OrthologousGroups.txt

Rscript plot_pangenomedist.r

Rscript selectCoreGenes.r

cat ../fna/*.fasta > ../tmp/allgenes.fasta

rm -r -f ../tmp/coreogfasta
mkdir -p ../tmp/coreogfasta 
perl splitToOg.pl ../OrthologousGroups.txt ../tmp/allgenes.fasta ../tmp/coreogfasta ../tmp/coreog

echo "aligning core genes "
rm -r -f ../tmp/coreogaligns
mkdir -p ../tmp/coreogaligns
for f in ../tmp/coreogfasta/*.fasta
do
	name=$(basename $f .fasta)
	./muscle -in ../tmp/coreogfasta/$name.fasta -out ../tmp/coreogaligns/$name.fasta -quiet
	printf "."
done

rm -f ../tmp/coreogaligned.fasta

perl concatenate_core.pl
perl renamer_simple.pl ../tmp/coreogaligned.fasta fasta ../tmp/coreogaligned_renamed.fasta > ../tmp/strains_num

echo "building tree (takes about 8hr)"
fdnaml -sequence ../tmp/coreogaligned_renamed.fasta -outfile ../tmp/dnaml.out -outtreefile ../tmp/dnaml.tree -progress -treeprint -intreefile ''

Rscript plot_tree.r

