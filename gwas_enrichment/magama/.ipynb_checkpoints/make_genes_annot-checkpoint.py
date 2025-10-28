#!/usr/bin/env python3

"""
make_genes_annot.py

Convert SNP–cCRE–gene intersection file into MAGMA .genes.annot format
using real gene coordinates (only autosomes, chr1–22).

Inputs:
  - --intersect: Output from bedtools intersect (snp–cCRE–gene file)
  - --gene-pos: Gene position file with columns: gene_name, chrom, start, end

Output:
  - .genes.annot file with: gene_id, CHR:FROM:TO, SNP IDs
"""

import argparse
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(description="Generate MAGMA .genes.annot using real gene coordinates.")
    parser.add_argument("--intersect", required=True,
                        help="Input file: SNP–cCRE–gene file from bedtools intersect")
    parser.add_argument("--gene-pos", required=True,
                        help="Gene position file: 4 columns [gene_name, chrom, start, end]")
    parser.add_argument("--out", required=True,
                        help="Output MAGMA .genes.annot file")
    return parser.parse_args()

def main():
    args = parse_args()

    # Load intersect file
    print(f"Reading intersect file: {args.intersect}")
    df = pd.read_csv(args.intersect, sep="\t", header=None,
                     names=["snp_chr", "snp_start", "snp_end", "snp_id",
                            "ccre_chr", "ccre_start", "ccre_end", "gene"])

    # Load gene position file
    print(f"Reading gene position file: {args.gene_pos}")
    gene_pos = pd.read_csv(args.gene_pos, sep="\t", header=None,
                           names=["gene", "chr", "start", "end"])

    # Remove "chr" prefix if present
    gene_pos["chr"] = gene_pos["chr"].astype(str).str.replace("^chr", "", regex=True)

    # Keep only autosomes (1–22)
    gene_pos = gene_pos[gene_pos["chr"].isin([str(c) for c in range(1, 23)])].copy()

    # Construct CHR:FROM:TO
    gene_pos["coord"] = gene_pos["chr"] + ":" + gene_pos["start"].astype(str) + ":" + gene_pos["end"].astype(str)
    gene_pos = gene_pos.set_index("gene")

    # Group SNPs by gene
    print("Aggregating SNPs by gene...")
    gene2snps = (
        df.groupby("gene")["snp_id"]
        .apply(lambda x: sorted(set(x)))
        .reset_index()
    )

    print(f"Writing MAGMA .genes.annot to {args.out}")
    with open(args.out, "w") as fout:
        for _, row in gene2snps.iterrows():
            gene = row["gene"]
            snps = row["snp_id"]
            if gene in gene_pos.index:
                coord = gene_pos.loc[gene, "coord"]
                fout.write(f"{gene}\t{coord}\t{' '.join(snps)}\n")
            else:
                # skip genes without coordinate info or on non-autosomes
                continue

    print("Done.")

if __name__ == "__main__":
    main()
