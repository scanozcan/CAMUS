#!/usr/bin/env python3
"""
Simple sgRNA counting from demuxed FASTQ files.
"""

import pandas as pd
from pathlib import Path
from typing import Dict
from collections import defaultdict


def count_sgrnas_from_fastq(
    fastq_file: str,
    library_file: str,
    output_file: str
) -> pd.DataFrame:
    """
    Count sgRNA occurrences from a demuxed FASTQ file.

    Args:
        fastq_file: Path to demuxed FASTQ file (R1)
        library_file: Path to library file with sgRNA sequences
        output_file: Path to save count output

    Returns:
        DataFrame with sgRNA counts
    """
    # Load library
    library_df = pd.read_csv(library_file, sep='\t')

    # Create hash table for fast lookup
    sgrna_to_gene = {}
    for _, row in library_df.iterrows():
        sgrna_to_gene[row['sgRNA_sequence']] = (row['sgRNA'], row['Gene'], row['in_library'])

    # Count sgRNAs
    counts = defaultdict(int)
    total_reads = 0
    matched_reads = 0

    with open(fastq_file, 'r') as f:
        while True:
            # Read FASTQ record (4 lines)
            header = f.readline()
            if not header:
                break
            seq = f.readline().strip()
            plus = f.readline()
            qual = f.readline()

            total_reads += 1

            # Check if sequence matches any sgRNA (first 20bp)
            seq_20bp = seq[:20]
            if seq_20bp in sgrna_to_gene:
                sgrna_name, gene, in_lib = sgrna_to_gene[seq_20bp]
                counts[sgrna_name] += 1
                matched_reads += 1

    # Create count DataFrame
    count_data = []
    for _, row in library_df.iterrows():
        count_data.append({
            'sgRNA': row['sgRNA'],
            'Gene': row['Gene'],
            'in_library': row['in_library'],
            'count': counts.get(row['sgRNA'], 0)
        })

    count_df = pd.DataFrame(count_data)

    # Save
    count_df.to_csv(output_file, sep='\t', index=False)

    print(f"  Total reads: {total_reads:,}")
    print(f"  Matched reads: {matched_reads:,} ({100*matched_reads/total_reads:.1f}%)")
    print(f"  Saved to: {output_file}")

    return count_df


def combine_counts_to_matrix(
    count_files: Dict[str, str],
    library_file: str,
    output_file: str
) -> pd.DataFrame:
    """
    Combine individual count files into a count matrix.

    Args:
        count_files: Dictionary mapping sample names to count file paths
        library_file: Path to library file
        output_file: Path to save count matrix

    Returns:
        Combined count matrix DataFrame
    """
    # Load library as template
    library_df = pd.read_csv(library_file, sep='\t')

    # Start with library columns
    matrix_df = library_df[['sgRNA', 'Gene', 'in_library']].copy()

    # Add count column for each sample
    for sample_name, count_file in count_files.items():
        count_df = pd.read_csv(count_file, sep='\t')
        matrix_df[sample_name] = count_df['count'].values

    # Save
    matrix_df.to_csv(output_file, sep='\t', index=False)

    print(f"  Combined {len(count_files)} samples")
    print(f"  Saved to: {output_file}")

    return matrix_df
