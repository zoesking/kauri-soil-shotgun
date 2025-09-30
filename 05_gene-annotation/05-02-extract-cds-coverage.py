#!/usr/bin/env python3

import argparse
import pandas as pd
import os

def load_dataframe(file_path):
    ext = os.path.splitext(file_path)[1]
    if ext == '.csv':
        return pd.read_csv(file_path)
    elif ext in ['.tsv', '.txt']:
        return pd.read_csv(file_path, sep='\t')
    else:
        raise ValueError(f"Unsupported file type: {ext}. Use .csv, .tsv, or .txt.")

def extract_coverage(coverage_df, annotation_df, coverage_id_col, annotation_id_col):
    # Ensure both ID columns exist
    if coverage_id_col not in coverage_df.columns:
        raise ValueError(f"Column '{coverage_id_col}' not found in coverage file.")
    if annotation_id_col not in annotation_df.columns:
        raise ValueError(f"Column '{annotation_id_col}' not found in annotation file.")

    annotated_ids = annotation_df[annotation_id_col].unique()
    filtered_df = coverage_df[coverage_df[coverage_id_col].isin(annotated_ids)]
    return filtered_df

def main():
    parser = argparse.ArgumentParser(description="Extract coverage for annotated genes.")
    parser.add_argument("-c", "--coverage", required=True, help="Path to coverage file (.txt, .tsv, .csv)")
    parser.add_argument("-a", "--annotation", required=True, help="Path to annotation file (.csv, .tsv)")
    parser.add_argument("--coverage_id", default="Contig", help="Column in coverage file with gene IDs (default: Contig)")
    parser.add_argument("--annotation_id", required=True, help="Column in annotation file with gene IDs (e.g., 'query' or '#query')")
    parser.add_argument("-o", "--output", required=True, help="Output file for filtered coverage table (.csv or .tsv)")

    args = parser.parse_args()

    # Load data
    coverage_df = load_dataframe(args.coverage)
    annotation_df = load_dataframe(args.annotation)

    # Filter coverage data for annotated genes
    filtered_df = extract_coverage(coverage_df, annotation_df, args.coverage_id, args.annotation_id)

    # Save output
    if args.output.endswith(".csv"):
        filtered_df.to_csv(args.output, index=False)
    elif args.output.endswith((".tsv", ".txt")):
        filtered_df.to_csv(args.output, sep='\t', index=False)
    else:
        raise ValueError("Output file must end with .csv, .tsv, or .txt")

if __name__ == "__main__":
    main()