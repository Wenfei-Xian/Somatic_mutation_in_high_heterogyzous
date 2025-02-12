# 1) Identification misalignment of short reads in haploid assembly

## step1 extract the flanking 2 kb sequence of each reads in the diploid assembly
```
for i in `cat Verkko.fasta.ID`;do echo "python3 extend.readsto2kb.indip.py -b lower.sorted.mdup.bam -r Verkko.fasta --chrom $i | /tmp/global2/wxian/software/pigz-master/pigz -p 2 > $i.2000bp.fasta.gz";done > tmp
cat tmp | parallel -j 100
```
## step2 align the extended sequence to haploid assembly
```
for i in `ls | grep "contig" | grep "gz$"`;do echo "/tmp/global2/wxian/software/minimap2-2.24_x64-linux/minimap2 -c --eqx -t 8 Quercus_robur.fasta $i -o Quercus_robur.fasta.$i.paf" > minimap2.$i.sh ;done
```
## step3 get the expected position of short reads in the haploid assembly
```
for i in `cat Verkko.fasta.ID`;do echo "perl best.paf.pl map_extended_sequences/Quercus_robur.fasta.$i.2000bp.fasta.gz.paf > tmp.$i.out" ;done > tmp
cat tmp | parallel -j 168
cat tmp.*.out > lower.reads.coordination.txt
rm tmp tmp.*.out
```
