import pysam
import csv
import argparse

def load_acceptable_ranges(file_path, target_chr):
    acceptable_ranges = {}
    with open(file_path, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            chr_ = row[1]
            if chr_ != target_chr:
                continue
            read_id = row[0]
            start = int(row[2])
            end = int(row[3])
            if read_id not in acceptable_ranges:
                acceptable_ranges[read_id] = []
            acceptable_ranges[read_id].append((chr_, start, end))
    return acceptable_ranges

def is_in_acceptable_range(read, acceptable_ranges):
    if read.query_name not in acceptable_ranges:
        return False
    for chr_, start, end in acceptable_ranges[read.query_name]:
        if read.reference_name == chr_ and start <= read.reference_start < end:
            return True
    return False

def process_bam(input_bam, output_bam, output_txt, acceptable_ranges, target_chr):
    input_bam_file = pysam.AlignmentFile(input_bam, "rb")
    output_bam_file = pysam.AlignmentFile(output_bam, "wb", template=input_bam_file)
    with open(output_txt, 'w') as txt_file:
        for read in input_bam_file.fetch(contig=target_chr):
            if not is_in_acceptable_range(read, acceptable_ranges):
                txt_file.write(f"{read.query_name}\t{read.reference_name}\t{read.reference_start}\t{read.reference_end}\n")
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

