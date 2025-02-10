import pysam
import argparse

def load_contig_lengths(fasta_file):
    """
    Load contig lengths from a FASTA file.

    Args:
        fasta_file (str): Path to the FASTA file.

    Returns:
        dict: A dictionary mapping contig names to their lengths.
    """
    from Bio import SeqIO
    contig_lengths = {}
    with open(fasta_file, "r") as fasta:
        for record in SeqIO.parse(fasta, "fasta"):
            contig_lengths[record.id] = len(record.seq)
    return contig_lengths

def process_bam(bam_file, fasta_file, output_file):
    """
    Process a BAM file to extract contig alignments and their positions on the reference genome,
    handling clipping information and outputting strand direction.

    Args:
        bam_file (str): Path to the input BAM file.
        fasta_file (str): Path to the contigs FASTA file.
        output_file (str): Path to the output file.
    """
    # Load contig lengths from the FASTA file
    contig_lengths = load_contig_lengths(fasta_file)

    # Open the BAM file
    bam = pysam.AlignmentFile(bam_file, "rb")

    # Open output file
    with open(output_file, "w") as out:
        # Write header
        out.write("Contig_Name\tContig_Direction\tContig_Start\tContig_End\tReference_Chromosome\tReference_Start_Position\tReference_End_Position\tContig_Length\tPositive_Strand_Contig_Name\tPositive_Strand_Contig_Start\tPositive_Strand_Contig_End\n")

        # Iterate over each alignment in the BAM file
        for read in bam:
            # Skip unmapped reads or reads with MAPQ not equal to 60
            if read.is_unmapped or read.mapping_quality != 60:
                continue

            # Get the full contig length from the FASTA file
            contig_name = read.query_name
            contig_length = contig_lengths.get(contig_name, None)

            if contig_length is None:
                raise ValueError(f"Contig {contig_name} not found in the provided FASTA file.")

            # Initialize contig start and end
            contig_start = 1
            contig_end = contig_length

            # Check for clipping in CIGAR
            if read.cigartuples:
                # Adjust for clipping at the start and end
                if read.cigartuples[0][0] in [4, 5]:  # Soft clip (4) or Hard clip (5)
                    contig_start = read.cigartuples[0][1] + 1
                if read.cigartuples[-1][0] in [4, 5]:  # Soft clip (4) or Hard clip (5)
                    contig_end = contig_length - read.cigartuples[-1][1]

            # Extract additional information
            contig_direction = "+" if not read.is_reverse else "-"
            reference_name = bam.get_reference_name(read.reference_id)
            reference_start_position = read.reference_start + 1  # Convert to 1-based
            reference_end_position = read.reference_end  # End position is already 1-based

            # Convert to positive strand coordinates if necessary
            if contig_direction == "-":
                positive_contig_start = contig_length - contig_end + 1
                positive_contig_end = contig_length - contig_start + 1
            else:
                positive_contig_start = contig_start
                positive_contig_end = contig_end

            # Write to output file
            out.write(f"{contig_name}\t{contig_direction}\t{contig_start}\t{contig_end}\t{reference_name}\t{reference_start_position}\t{reference_end_position}\t{contig_length}\t{contig_name}\t{positive_contig_start}\t{positive_contig_end}\n")

    # Close BAM file
    bam.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process a BAM file to extract alignment positions.")
    parser.add_argument("--input", required=True, help="Path to the input BAM file.")
    parser.add_argument("--fasta", required=True, help="Path to the contigs FASTA file.")
    parser.add_argument("--output", required=True, help="Path to the output file.")

    args = parser.parse_args()

    process_bam(args.input, args.fasta, args.output)
    print(f"Processed BAM file and saved results to {args.output}")

