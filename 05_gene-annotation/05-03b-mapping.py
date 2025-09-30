#!/usr/bin/python

import pandas as pd
import argparse

def main(contig2acc_file, gene_map_file, output_file):
    # Define column names for the DIAMOND output file (12-column tab-separated)
    col_names = [
        "query", "accession", "pident", "length", "mismatch", "gapopen", 
        "qstart", "qend", "sstart", "send", "evalue", "bitscore"
    ]
    # Step 1: Read contig2acc (DIAMOND result)
    contig2acc = pd.read_csv(contig2acc_file, sep='\t', header=None, names=col_names)

    # Step 2: Read gene mapping file (2-column, tab-separated, no headers)
    gene_df = pd.read_csv(gene_map_file, sep='\t', header=None, names=['accession', 'gene'])

    # Step 3: Build the gene map dictionary
    geneMap = pd.Series(gene_df['gene'].values, index=gene_df['accession']).to_dict()

    # Step 4: Map accessions to gene names
    contig2acc['gene'] = contig2acc['accession'].map(geneMap)

    # Step 5: Drop rows without a gene match
    contig2acc = contig2acc.dropna(subset=['gene'])

    # Step 6: Save to output file
    contig2acc.to_csv(output_file, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Map genes to DIAMOND hits and filter unmatched entries.")
    parser.add_argument("contig2acc_file", help="Tab-separated DIAMOND output file (12 columns, no header)")
    parser.add_argument("gene_map_file", help="2-column tab-separated file with accession-to-gene mappings")
    parser.add_argument("output_file", help="Path to write filtered output CSV")

    args = parser.parse_args()
    main(args.contig2acc_file, args.gene_map_file, args.output_file)