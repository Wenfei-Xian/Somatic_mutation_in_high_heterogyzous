# 1) Identification misalignment of short reads in haploid assembly

## step1 get the position of diploid assembly in haploid assembly
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
## step3 check the position of short in haploid assembly located in the suitable position
```
perl /ebio/abt6/wxian/script/HQSNV/split.pl upper.reads.coordination.txt
for i in {1..12};do echo "python3 /ebio/abt6/wxian/script/HQSNV/step3.contigsfilter.perfectalign.py --ranges upper.reads.coordination.txt.Chr$i.ID --input_bam  ../01.hifiasm_0.24.0/09.bwa_alignment_v1/upper.sorted.mdup.bam --output_bam hap.Chr$i.upper.sorted.mdup.bam --output_txt hap.Chr$i.upper.sorted.mdup.bam.out.txt --target_chr Chr$i" ;done > tmp1
cat tmp1 | parallel -j 24
```
