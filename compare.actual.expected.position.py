import pysam
import csv
import argparse

def load_acceptable_ranges(file_path, target_chr):
    acceptable_ranges = {}
    with open(file_path, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            chr_ = row[2]  # 染色体
            if chr_ != target_chr:  # 只处理目标染色体
                continue
            read_id = row[0]  # Read ID
            pair_label = row[1]  # 1: first in pair, 2: second in pair
            start = int(row[3])  # 起始位置
            end = int(row[4])  # 终止位置
            
            key = (read_id, pair_label)  # 使用 (read_id, pair_label) 作为键
            if key not in acceptable_ranges:
                acceptable_ranges[key] = []
            acceptable_ranges[key].append((chr_, start, end))
    return acceptable_ranges

def is_in_acceptable_range(read, acceptable_ranges):
    pair_label = "1" if read.is_read1 else "2"  # 读取 read 的方向信息
    key = (read.query_name, pair_label)  # 组合 key

    if key not in acceptable_ranges:  # 判断 key 是否在可接受范围内
        return False

    for chr_, start, end in acceptable_ranges[key]:  # 遍历可接受区间
        if read.reference_name == chr_ and start <= read.reference_start < end:
            return True
    return False

def process_bam(input_bam, output_bam, output_txt, acceptable_ranges, target_chr):
    input_bam_file = pysam.AlignmentFile(input_bam, "rb")
    output_bam_file = pysam.AlignmentFile(output_bam, "wb", template=input_bam_file)

    with open(output_txt, 'w') as txt_file:
        for read in input_bam_file.fetch(contig=target_chr):
            pair_label = "1" if read.is_read1 else "2"  # 读取方向信息
            key = (read.query_name, pair_label)  # 组合 key
            
            expected_positions = acceptable_ranges.get(key, [])  # 获取 read 期望的范围
            
            # 如果没有期望范围，显示 "None"，否则格式化期望范围
            expected_str = ";".join([f"{chr_}:{start}-{end}" for chr_, start, end in expected_positions]) if expected_positions else "None"

            if is_in_acceptable_range(read, acceptable_ranges) :
                output_bam_file.write(read)  # 只写入符合范围且剪切不过多的 reads
            else:
                txt_file.write(f"{read.query_name}\t{pair_label}\t{read.reference_name}\t{read.reference_start}\t{read.reference_end}\t{expected_str}\n")
    
    input_bam_file.close()
    output_bam_file.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter reads based on acceptable ranges and exclude reads with clipping > 1bp.")
    parser.add_argument("--ranges", required=True, help="Path to the acceptable ranges file.")
    parser.add_argument("--input_bam", required=True, help="Path to the input BAM file.")
    parser.add_argument("--output_bam", required=True, help="Path to the output BAM file for filtered reads.")
    parser.add_argument("--output_txt", required=True, help="Path to the output text file for unmatched reads.")
    parser.add_argument("--target_chr", required=True, help="Specify a target chromosome to filter reads.")

    args = parser.parse_args()

    acceptable_ranges = load_acceptable_ranges(args.ranges, args.target_chr)

    process_bam(args.input_bam, args.output_bam, args.output_txt, acceptable_ranges, args.target_chr)

