#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pysam
import argparse
import os
import sys

def main():
    parser = argparse.ArgumentParser(
        description="Extract reference subsequences around properly paired, non-clipped reads mapped to a specific chromosome."
    )
    parser.add_argument("-b", "--bam", required=True, help="Input BAM file")
    parser.add_argument("-r", "--ref", required=True, help="Reference genome (FASTA)")
    parser.add_argument("-o", "--out", default="-", help="Output FASTA file (default: stdout)")
    parser.add_argument("--chrom", default="contig-0000001", help="Chromosome name to extract (default: Chr2)")
    parser.add_argument("--flank", type=int, default=2000, help="Flanking size to extend (default: 2000bp)")
    parser.add_argument("--exclude-clipped", action="store_true", help="Exclude soft/hard clipped reads")

    args = parser.parse_args()

    # 检查 BAM 索引
    if not os.path.exists(args.bam + ".bai"):
        raise FileNotFoundError(f"BAM index file ({args.bam}.bai) not found! Please run 'samtools index {args.bam}'.")

    with pysam.AlignmentFile(args.bam, "rb", threads=4) as bamfile, \
         pysam.FastaFile(args.ref) as ref_fa, \
         (open(args.out, "w") if args.out != "-" else sys.stdout) as outfile:

        # 获取染色体长度
        try:
            chrom_length = ref_fa.get_reference_length(args.chrom)
        except ValueError:
            raise ValueError(f"Chromosome '{args.chrom}' not found in reference {args.ref}.")

        # 遍历 BAM
        for read in bamfile.fetch(args.chrom):

            # 筛选条件
            if not read.is_paired or read.is_unmapped or read.mate_is_unmapped or not read.is_proper_pair:
                continue

            if read.is_secondary or read.is_supplementary:
                continue

            # 可选：排除 clipped reads
            if args.exclude_clipped and any(op in {4, 5} for op, length in (read.cigartuples or [])):
                continue

            # 获取对齐区域
            ref_start = read.reference_start
            ref_end = read.reference_end

            # 计算 flanking 区域
            flank = args.flank
            start_pos = max(0, ref_start - flank)
            end_pos = min(chrom_length, ref_end + flank)

            if start_pos == 0:
                continue

            if end_pos == chrom_length:
                continue

            # 提取参考序列
            seq = ref_fa.fetch(args.chrom, start_pos, end_pos)

            # 读的方向信息
            pair_label = "1" if read.is_read1 else "2"

            # 生成 FASTA ID
            #fasta_id = f"{read.query_name}_{pair_label}_{args.chrom}_{start_pos}_{end_pos}"
            fasta_id = f"{read.query_name}_{pair_label}"
            outfile.write(f">{fasta_id}\n{seq}\n")

    print(f"Done. Sequences saved to {args.out if args.out != '-' else 'stdout'}", file=sys.stderr)

if __name__ == "__main__":
    main()

