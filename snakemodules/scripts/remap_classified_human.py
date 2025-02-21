from Bio import SeqIO
import argparse
import os
import logging

def extract_human_read_ids(classification_file):
    human_reads = set()
    lines_processed = 0
    with open(classification_file, 'r') as f:
        for line in f:
            lines_processed += 1
            parts = line.strip().split()
            if len(parts) < 3:
                logging.warning(f"Line {lines_processed} in {classification_file} skipped (insufficient columns): {line.strip()}")
                continue
            # First column: classification (C/U)
            classification = parts[0]
            # Second column: read ID
            read_id = parts[1]
            # Third column: taxonomy ID (human is '9606')
            taxonomy = parts[2]
            if classification == 'C' and taxonomy == '9606':
                human_reads.add(read_id)
    logging.info(f"Extracted {len(human_reads)} human read IDs from {classification_file} (processed {lines_processed} lines)")
    return human_reads

def extract_reads_from_fastq(fastq_file, human_read_ids, output_fastq):
    total_records = 0
    human_records = 0
    with open(output_fastq, 'w') as out_fq:
        for record in SeqIO.parse(fastq_file, "fastq"):
            total_records += 1
            if record.id in human_read_ids:
                SeqIO.write(record, out_fq, "fastq")
                human_records += 1
    logging.info(f"Processed {total_records} records from {fastq_file}; extracted {human_records} human reads into {output_fastq}")

def main(single_output, paired_output, single_fq, paired_fq1, paired_fq2, output_prefix):
    # Extract human read IDs from the classification files
    logging.info("Starting extraction of human read IDs from classified single reads...")
    single_human_reads = extract_human_read_ids(single_output)
    logging.info("Starting extraction of human read IDs from classified paired reads...")
    paired_human_reads = extract_human_read_ids(paired_output)

    # Combine the read IDs from single and paired outputs
    all_human_reads = single_human_reads.union(paired_human_reads)
    logging.info(f"Total unique human read IDs combined: {len(all_human_reads)}")

    # Construct output file paths using os.path.join
    single_out = os.path.join(output_prefix, 'human_classified_single.fastq')
    paired_out1 = os.path.join(output_prefix, 'human_classified_paired_1.fastq')
    paired_out2 = os.path.join(output_prefix, 'human_classified_paired_2.fastq')

    # Extract human reads from the FASTQ files
    logging.info("Extracting human reads from FASTQ files...")
    extract_reads_from_fastq(single_fq, all_human_reads, single_out)
    extract_reads_from_fastq(paired_fq1, all_human_reads, paired_out1)
    extract_reads_from_fastq(paired_fq2, all_human_reads, paired_out2)

    logging.info("Extraction complete.")
    logging.info(f"Outputs saved to:\n  {single_out}\n  {paired_out1}\n  {paired_out2}")

if __name__ == "__main__":
    # Set up logging to include the time, log level, and message
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    
    parser = argparse.ArgumentParser(description="Extract human reads based on classification and taxonomy.")
    parser.add_argument("--classified_single", required=True, help="Path to the single_output.txt file")
    parser.add_argument("--classified_paired", required=True, help="Path to the paired_output.txt file")
    parser.add_argument("--single_fq", required=True, help="Path to the single.fq file")
    parser.add_argument("--paired_fq1", required=True, help="Path to the paired.1.fq file")
    parser.add_argument("--paired_fq2", required=True, help="Path to the paired.2.fq file")
    parser.add_argument("--output_prefix", required=True, help="Directory for the output FASTQ files")
    
    args = parser.parse_args()
    main(args.classified_single, args.classified_paired, args.single_fq, args.paired_fq1, args.paired_fq2, args.output_prefix)

