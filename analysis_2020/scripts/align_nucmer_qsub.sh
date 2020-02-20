mkdir -p nucmer

for f in ../genomes/*.fasta
do 
	name=$(basename $f .fasta)
	echo $name
	echo "nucmer -p $name  ../ref/reference.fasta $f; show-coords -d -T -H $name.delta -L 500 > nucmer/$name.coords; rm $name.delta" | qsub -cwd 
done



## 2do
# select reference automaticly if no refence is defined. Largerst and finilized as cirteria
