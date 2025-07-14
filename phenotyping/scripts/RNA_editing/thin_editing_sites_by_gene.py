#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, fcluster
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description="Reduce RNA editing site redundancy by clustering correlated sites per gene")
    parser.add_argument("--edit_matrix", required=True, help="Input edit matrix (site x sample, tab-separated)")
    parser.add_argument("--site_to_gene", required=True, help="TSV with columns: site<TAB>gene_id")
    parser.add_argument("--output_matrix", required=True, help="Output reduced matrix (cluster-averaged edit levels)")
    parser.add_argument("--output_clusters", required=True, help="Output TSV with columns: site, gene_id, cluster_id")
    parser.add_argument("--correlation_threshold", type=float, default=0.9, help="Min Pearson r to group sites into a cluster (default: 0.9)")
    return parser.parse_args()

def cluster_sites(edit_df: pd.DataFrame, threshold: float):
    """Cluster sites based on correlation threshold and return a mapping: cluster_id → list of site_ids"""
    if edit_df.shape[0] <= 1:
        return {"Cluster_1": edit_df.index.tolist()}

    # Calculate pairwise correlation matrix
    corr = edit_df.T.corr().fillna(0)

    # Clip correlation values to avoid floating-point precision issues
    corr = corr.clip(lower=-1.0, upper=1.0)

    # Convert to distance matrix
    dist = 1 - corr
    dist[dist < 0] = 0  # Ensure no negative distances

    dist_matrix = squareform(dist.values, checks=False)

    # Hierarchical clustering
    Z = linkage(dist_matrix, method="average")
    cluster_labels = fcluster(Z, t=1 - threshold, criterion="distance")

    # Build cluster mapping
    clusters = defaultdict(list)
    for site, label in zip(edit_df.index, cluster_labels):
        clusters[f"Cluster_{label}"].append(site)
    return clusters


def main():
    args = parse_args()

    # Load matrix and mapping
    matrix = pd.read_csv(args.edit_matrix, sep="\t", index_col=0)
    mapping = pd.read_csv(args.site_to_gene, sep="\t", header=None, names=["site", "gene"])

    print(f"Loaded matrix with {matrix.shape[0]} sites and {matrix.shape[1]} samples")
    print(f"Loaded site-to-gene map with {mapping.shape[0]} entries")

    # Group by gene and reduce
    output_rows = []
    cluster_mappings = []

    for gene, group in mapping.groupby("gene"):
        site_ids = group["site"].tolist()
        submatrix = matrix.loc[matrix.index.intersection(site_ids)]

        if submatrix.empty or submatrix.shape[0] < 1:
            continue

        # Convert ratios to float edit levels
        float_matrix = submatrix.applymap(lambda val: float(val.split("/")[0]) / float(val.split("/")[1]) if "/" in val and val != "0/0" and float(val.split("/")[1]) > 0 else np.nan)

        clusters = cluster_sites(float_matrix, args.correlation_threshold)
        for cluster_id, sites in clusters.items():
            cluster_mean = float_matrix.loc[sites].mean(axis=0, skipna=True)
            cluster_name = f"{gene}__{cluster_id}"
            output_rows.append(pd.DataFrame([cluster_mean], index=[cluster_name]))

            for site in sites:
                cluster_mappings.append((site, gene, cluster_name))

    # Save reduced matrix
    reduced_matrix = pd.concat(output_rows)
    reduced_matrix.index.name = "site"
    reduced_matrix.to_csv(args.output_matrix, sep="\t", na_rep="NA")

    # Save cluster mapping
    cluster_df = pd.DataFrame(cluster_mappings, columns=["site", "gene", "cluster_id"])
    cluster_df.to_csv(args.output_clusters, sep="\t", index=False)

    print(f"Wrote reduced matrix with {reduced_matrix.shape[0]} clusters to {args.output_matrix}")
    print(f"Wrote cluster mapping for {len(cluster_df)} sites to {args.output_clusters}")

if __name__ == "__main__":
    main()
