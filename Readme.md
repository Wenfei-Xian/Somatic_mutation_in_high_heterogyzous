# Classification the position of short reads alignment

## step1 get the position of diploid assembly in haploid assemby
```
minimap2 -ax asm20 -t 128 ../01.scaffolding/Quercus_robur.fasta Verkko.fasta | samtools sort -@24 -o Quercus_robur.fasta.verkko.sorted.bam -
python3 /ebio/abt6/wxian/script/HQSNV/step1.contigs.refspan.orientation.py --input Quercus_robur.fasta.verkko.sorted.bam --output Quercus_robur.fasta.verkko.sorted.bam.alignment.txt --fasta Verkko.fasta
```
## step2 get the approximate position of short reads in haploid assemby
```
grep ">" Verkko.fasta | sed "s/>//" > Verkko.fasta.ID
for i in `cat Verkko.fasta.ID `;do echo "perl /ebio/abt6/wxian/script/HQSNV/step2.allowed.goodreads.coordination.pl Quercus_robur.fasta.verkko.sorted.bam.alignment.txt upper.sorted.mdup.bam $i > tmp.$i.out" ;done > tmp
cat tmp | parallel -j 64
cat tmp.*.out > upper.reads.coordination.txt
rm tmp tmp.*.out
```
