threads=$1
orthofinder -t $threads -a $threads -og -f faa
mkdir -p Results
mkdir -p Results/ortho

find faa -name 'Orthogroups.txt' -exec cp {} Results/ortho \; 
find faa -name 'SingleCopyOrthogroups.txt' -exec cp {} Results/ortho \; 


