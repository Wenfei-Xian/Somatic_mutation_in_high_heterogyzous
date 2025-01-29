# Classification the position of short reads alignment

## step1 get the position of diploid assembly on haploid assemby
```
minimap2 -ax asm20 -t 128 ../01.scaffolding/Quercus_robur.fasta Verkko.fasta | samtools sort -@24 -o Quercus_robur.fasta.verkko.sorted.bam -
python3 /ebio/abt6/wxian/script/HQSNV/step1.contigs.refspan.orientation.py --input /ebio/abt6_projects/Oak_somatic_mutation/data/01.hifiasm_0.24.0/14.graph/Quercus_robur.fasta.verkko.sorted.bam --output Quercus_robur.fasta.verkko.sorted.bam.alignment.txt --fasta Verkko.fasta
```
