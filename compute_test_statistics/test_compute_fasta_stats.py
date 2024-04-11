"""
This file contains unit tests for the compute_fasta_stats.py script.

"""
import os
from tempfile import NamedTemporaryFile as NamedTemp

import pytest
from compute_fasta_stats import (
    calculate_seq_stats,
    compute_stats,
    read_seqs_from_fasta,
    seq_stats,
)

test_fasta = "example.fa"


def create_tempfasta(sequence: str, seqid: str = ">testid") -> str:
    """
    This function creates a temporary fasta file based on the sequence.

    Args:
        sequence (str): sequence to be written to the FASTA file.
        seqid(str): value of header for the seqid

    Returns:
        str: path to created FASTA file.
    """

    with NamedTemp(mode="w", suffix=".fasta", delete=False) as temp_fasta:
        temp_fasta.write(f"{seqid}\n{sequence}")
        temp_fasta.flush()
        temp_fasta_path = temp_fasta.name

    return temp_fasta_path


def test_emptyfile():
    """
    Tests whether empty file is handled properly
    """
    testfa_path = create_tempfasta("", "")
    with pytest.raises(ValueError):
        for _ in read_seqs_from_fasta(testfa_path):
            pass
    os.remove(testfa_path)


def test_invalid_fasta_detected():
    """
    Tests whether fasta file not starting
    with a '>' is detected
    """
    testfa_path = create_tempfasta("CTCTC", "invalid")

    with pytest.raises(ValueError):
        for _ in read_seqs_from_fasta(testfa_path):
            pass
    os.remove(testfa_path)


def test_validseq():
    """
    Tests whether sequence without valid
    characters are skipped
    """
    testfa_path = create_tempfasta("CTCMMTC")

    sequences = list(read_seqs_from_fasta(testfa_path))

    assert (
        len(sequences) == 0
    ), "Sequence containing invalid characters should be skipped"
    os.remove(testfa_path)


def test_proper_case():
    """
    Tests that different cases in sequence is properly converted to upper
    """
    testfa_path = create_tempfasta("aaTacGG")

    _, sequences = list(read_seqs_from_fasta(testfa_path))[0]

    assert sequences == "AATACGG", "Cases should have been converted to upper"
    os.remove(testfa_path)


def test_empty_sequence_properly():
    """
    Tests whether an empty sequence
    is properly skipped
    """
    testfa_path = create_tempfasta("")

    sequences = list(read_seqs_from_fasta(testfa_path))

    assert (
        len(sequences) == 0
    ), "Expected no sequences to be yielded for an empty sequence"

    os.remove(testfa_path)


def test_allN_seqs():
    """
    Test situation where all Ns
    """
    testfa_path = create_tempfasta("NNNN")
    _, sequence = list(read_seqs_from_fasta(testfa_path))[0]
    stats = calculate_seq_stats(sequence)

    assert (
        stats == seq_stats(4, 0, 0, 4)
    ), "Sequence with all Ns should have length and max homoplymer of 4 and gc frac/cpg count of 0"
    os.remove(testfa_path)


def test_no_gc_nor_cpg():
    """
    Test situation where there's no Gs or Cs
    """
    testfa_path = create_tempfasta("ATATATAT")
    _, sequence = list(read_seqs_from_fasta(testfa_path))[0]
    stats = calculate_seq_stats(sequence)

    assert stats.gc_frac == 0, "Sequence without Gs or Cs should have 0 GC fraction"
    assert stats.cpgs == 0, "Sequence without Gs or Cs should have no cpgs"
    os.remove(testfa_path)


def test_homopolymer():
    """
    Tests situation where there's homoplymers at end of sequence
    """
    testfa_path = create_tempfasta("ATATATAAAAA")
    _, sequence = list(read_seqs_from_fasta(testfa_path))[0]
    stats = calculate_seq_stats(sequence)

    assert stats.max_hp == 5, "Incorrect calculation of max homopolymer length"
    os.remove(testfa_path)


def test_multiple_sequences():
    """
    Tests whether multiple sequences are handled properly.
    """
    correct_stats = {
        "chr1": seq_stats(331, pytest.approx(0.42, 0.1), 5, 5),
        "chr2": seq_stats(350, pytest.approx(0.43, 0.1), 2, 5),
        "chr3": seq_stats(463, pytest.approx(0.4, 0.1), 7, 10),
        "chr4": seq_stats(318, pytest.approx(0.35, 0.1), 4, 6),
        "chr5": seq_stats(326, pytest.approx(0.42, 0.1), 2, 9),
    }

    for seqid, sequence in read_seqs_from_fasta(test_fasta):
        stats = calculate_seq_stats(sequence)
        assert stats == correct_stats[seqid], f"Stats do not match for {seqid}"


def test_output_written():
    """
    Test whether output is written as expected properly
    """
    output_path = "/tmp/results.txt"
    compute_stats(test_fasta, output_path)
    assert os.path.exists(output_path), f"File '{output_path}' was not created"
    os.remove(output_path)


def test_proper_output_data_written():
    """
    Tests whether the saved output file contains the expected text
    """
    testfa_path = create_tempfasta("GATAA")
    output_path = "/tmp/results.txt"

    compute_stats(testfa_path, output_path)

    # open file and check for header and the data is present
    with open(output_path, "r") as output_file:
        header, data = output_file.readlines()

        expected_header = "seqid\tlength\tfrac_gc\tcpgs\tlongest_homopolymer\n"
        expected_data = "testid\t5\t0.20\t0\t2\n"

        assert expected_header == header, "The header in the file is not as expected"

        assert (
            expected_data == data
        ), "The data in the file does not match expected data"

    # Cleanup: Remove temporary files
    os.remove(testfa_path)
    os.remove(output_path)
