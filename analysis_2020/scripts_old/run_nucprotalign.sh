for i in ../ogaligned/*.fasta
do
        name=$(basename $i .fasta)
        grep ">" $i | sed 's/>//' > order
        ./reorder.py ../ogfastanuc/unal/$name.fasta order > ../ogfastanuc/$name.fasta
        perl pal2nal.v14/pal2nal.pl $i  ../ogfastanuc/$name.fasta -output fasta > ../ogalignednuc/$name.fasta
#       echo "muscle -in ../ogfastanuc/$name.fasta -out ../ogaligned/$f.fasta" | qsub -N "mus" -cwd 

done
