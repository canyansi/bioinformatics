# Sequence Analysis Tools

A collection of Python tools for analyzing DNA sequences, including FASTA file statistics computation and CNV (Copy Number Variation) analysis.

## Compute Test Statistics

### Components
- `compute_fasta_stats.py`: Main script for processing FASTA files
- `test_compute_fasta_stats.py`: Test suite
- `example.fa`: Example FASTA file

### Description
The `compute_fasta_stats.py` script reads a FASTA file, computes various statistics for each DNA sequence contained within the file, and outputs the results in a tab-delimited format. The statistics calculated include sequence length, GC content, the number of CpG sites, and the length of the longest homopolymer.

### Input
- **Input File**: The script accepts a FASTA file as input. The file should contain DNA sequences formatted according to the FASTA standard, where each sequence starts with a header line beginning with `>`, followed by lines of sequence data.

### Output
- The results of the analysis are saved in a specified output file, typically including:
  - A tab-delimited file containing the following columns:
    - `seqid`: The identifier for each sequence.
    - `length`: The length of the sequence.
    - `frac_gc`: The fraction of bases that are GC.
    - `cpgs`: The number of CpG sites in the sequence.
    - `longest_homopolymer`: The length of the longest homopolymer in the sequence.

### Usage
To run the FASTA statistics script, use the following command:
```bash
python compute_test_statistics/compute_fasta_stats.py -i <input_fasta> -o <output_file>
```

## CNV Analysis

Tools for analyzing Copy Number Variations in genomic data.

### Components
- `CNV_identification.py`: Main script for CNV detection

### Description
The `CNV_identification.py` script is designed to identify copy number variations (CNVs) from genomic data. It utilizes normalization techniques to adjust for sample variability and probe variability, followed by statistical methods to estimate copy numbers. The script detects regions where the number of copies of a particular gene or genomic region differs from the expected number, categorizing these variations as deletions or duplications.

### Input
- **Input File**: The script accepts genomic data in a CSV format. For example, the `cnsl_data.csv` file contains:
  - **Columns**: 
    - `ethnicity`: The ethnicity of each sample.
    - Other columns represent probe data with copy number values.
  - **Format**: The first column is the index (sample identifiers), followed by the `ethnicity` column and then the probe data columns containing numeric values representing the copy number for each probe.

### Output
- The results of the CNV analysis are saved in the specified output directory, typically including:
  - A summary report of detected CNVs, indicating the number of different CNV/breakpoint combinations detected.
  - Detailed information about the genomic regions affected, including the counts of CNVs per ethnicity and the proportions of each CNV type.

### Usage
To run the CNV identification script, use the following command:
```bash
python cnv_analysis/CNV_identification.py -i <input_file> -o <output_directory> [options]
```

## Installation

Create conda environment:

    conda env create -f environment.yml
    conda activate sequence_analysis

## Testing

    cd compute_test_statistics
    pytest test_compute_fasta_stats.py

## License

MIT License

