"""
This script reads in a fasta file, computes statistics on each sequence
contained within the file, then outputs the stats in a tab
delimited format.

~~~~~~~~~~~~~~~~~~
USAGE
~~~~~~~~~~~~~~~~~~

python compute_fasta_stats.py input.fa path/to/output.txt

"""

import argparse
from collections import namedtuple
from typing import Generator


def read_seqs_from_fasta(fasta_file: str) -> Generator[tuple[str, str], None, None]:
    """
    This function takes a path to a fasta file
    and returns a generator containing the ID and sequence
    for every sequence in the file. ID refers to the first word
    after the '>'. Only fasta files containing DNA
    will be considered. It skips sequences with invalid characters.
    It also skips empty sequences.

    Args:
        fasta_file : path to fasta file
    Returns:
        Iterator: The seqID and its corresponding sequence

    """

    with open(fasta_file, "r") as fasta_reader:
        line = fasta_reader.readline().strip()
        if line == "":
            raise ValueError("File is empty")
        if not line.startswith(">"):
            raise ValueError("Fasta files should start with '>'")

        while line:
            seqid = line.split()[0][1:]
            sequence = ""
            line = fasta_reader.readline().strip()
            # get the sequence
            while not line.startswith(">") and line:
                sequence += line.upper()
                line = fasta_reader.readline().strip()

            if len(sequence) == 0:
                print(f"Skipping {seqid} as there it has no corresponding sequence")
            elif all(base in "ATCGN" for base in sequence):
                yield seqid, sequence
            else:
                print(f"Skipping sequence {seqid} as an unknown character was detected")


# Object to store stats values
seq_stats = namedtuple("seq_stats", ["length", "gc_frac", "cpgs", "max_hp"])


def calculate_seq_stats(sequence: str) -> seq_stats:
    """
    This function calculates sequence length,
    fraction of bases that are GC, the CpG count,
    and the length of the longest homopolymer
    for a given sequence. Returns stats in a named tuple
    seq_stats object

    Args:
        sequence: string containing DNA characters
    Returns:
        seq_stats: named tuple containing length, fraction of
        gc bases, # of cpgs and max homopolymer length for the
        sequence

    """
    gc_count, cpg_count, homopolymer_length, max_hp_length = 0, 0, 0, 0

    # Loop through bases in sequence. Keep track of prev base
    # in order to detect cpgs and to check for homopolymers
    prev_base = ""
    for base in sequence:
        if base in "CG":
            gc_count += 1
        if base == "G" and prev_base == "C":
            cpg_count += 1
        if base == prev_base:
            homopolymer_length += 1
        else:
            homopolymer_length = 1
        max_hp_length = max(homopolymer_length, max_hp_length)
        prev_base = base

    # calculate GC content if len sequence greater than 0
    if sequence:
        frac_gc = gc_count / len(sequence)
    else:
        frac_gc = 0
    stats = seq_stats(len(sequence), frac_gc, cpg_count, max_hp_length)
    return stats


def compute_stats(fasta_file: str, output_file: str) -> None:
    """
    This function processes a fasta file and writes statistics
    for each sequence to an output file.

    Args:
        fasta_file: path to fasta file
        output_file: path to desired output file

    """
    header = ["seqid", "length", "frac_gc", "cpgs", "longest_homopolymer"]
    with open(output_file, "w") as output_writer:
        output_writer.write("\t".join(header))
        output_writer.write("\n")
        for seqid, sequence in read_seqs_from_fasta(fasta_file):
            stats = calculate_seq_stats(sequence)
            output_writer.write(
                f"{seqid}\t{stats.length}\t{stats.gc_frac:.2f}\t{stats.cpgs}\t{stats.max_hp}\n"
            )


def main():
    parser = argparse.ArgumentParser(
        description="Computes statistics from a fasta file"
    )
    parser.add_argument(
        "-i",
        "--input_fasta",
        type=str,
        help="Path to the fasta file that stats will be computed on",
    )
    parser.add_argument(
        "-o",
        "--output_file",
        type=str,
        help="Path to the output file where the stats will be saved",
        default="output_stats.txt",
    )

    args = parser.parse_args()
    print(
        f"Calculating stats from {args.input_fasta} and saving them to {args.output_file}"
    )
    compute_stats(args.input_fasta, args.output_file)


if __name__ == "__main__":
    main()
