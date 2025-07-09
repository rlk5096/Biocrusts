#!/usr/bin/env python3

import sys
from Bio import SeqIO
import os

def split_fasta(input_fasta, output_dir, sequences_per_file):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    basename = os.path.basename(input_fasta).replace(".fa", "").replace(".fasta", "").replace(".fna", "")
    record_iter = SeqIO.parse(input_fasta, "fasta")

    file_count = 1
    records = []

    for i, record in enumerate(record_iter, 1):
        records.append(record)
        if i % sequences_per_file == 0:
            out_file = os.path.join(output_dir, f"{basename}{file_count}.fa")
            SeqIO.write(records, out_file, "fasta")
            records = []
            file_count += 1

    # Write remaining sequences
    if records:
        out_file = os.path.join(output_dir, f"{basename}{file_count}.fa")
        SeqIO.write(records, out_file, "fasta")

    print(f"Split complete. {file_count} files written to: {output_dir}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: split_fasta.py <input_fasta> <output_dir> <sequences_per_file>")
        sys.exit(1)

    input_fasta = sys.argv[1]
    output_dir = sys.argv[2]
    sequences_per_file = int(sys.argv[3])

    split_fasta(input_fasta, output_dir, sequences_per_file)

