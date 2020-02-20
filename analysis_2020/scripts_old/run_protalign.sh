cat ../faa/*.fasta > ../tmp/allgenes_aa.fasta

rm -r -f ../tmp/coreogfasta_aa
mkdir -p ../tmp/coreogfasta_aa
perl splitToOg.pl ../OrthologousGroups.txt ../tmp/allgenes_aa.fasta ../tmp/coreogfasta_aa ../tmp/coreog

echo "aligning core genes "
rm -r -f ../tmp/coreogaligns_aa
mkdir -p ../tmp/coreogaligns_aa
for f in ../tmp/coreogfasta_aa/*.fasta
do
	name=$(basename $f .fasta)
	./muscle -in ../tmp/coreogfasta_aa/$name.fasta -out ../tmp/coreogaligns_aa/$name.fasta -quiet
	sed -i 's/|.\+//g' ../tmp/coreogaligns_aa/$name.fasta 			# remove protein id
	printf "."
done

