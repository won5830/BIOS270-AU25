
import sys
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

def parse_fasta_headers(faa_path):
    id_to_name = {}
    with open(faa_path) as f:
        for line in f:
            if line.startswith(">"):
                header = line[1:].strip()
                parts = header.split(None, 1)
                prot_id = parts[0]
                prot_name = parts[1] if len(parts) > 1 else ""
                id_to_name[prot_id] = prot_name
    return id_to_name


def main(faa_path, cluster_tsv, out_prefix):
    prot_name_map = parse_fasta_headers(faa_path)

    df = pd.read_csv(
        cluster_tsv,
        sep="\t",
        header=None,
        names=["cluster_id", "protein_id"]
    )

    df = df.drop_duplicates(subset=["cluster_id", "protein_id"])

    cluster_sizes = df.groupby("cluster_id").size().rename("copy_number")
    df = df.join(cluster_sizes, on="cluster_id")

    df["protein_name"] = df["protein_id"].map(prot_name_map).fillna("")

    summary_cols = ["protein_id", "protein_name", "copy_number", "cluster_id"]
    summary = df[summary_cols].sort_values("copy_number", ascending=False)
    summary_tsv = f"{out_prefix}_paralog_summary.tsv"
    summary.to_csv(summary_tsv, sep="\t", index=False)
    print(f"Wrote summary table: {summary_tsv}")

    paralog_df = summary[summary["copy_number"] > 1]

    if not paralog_df.empty:
        top = paralog_df.head(20).copy()
        labels = top["protein_name"].replace("", pd.NA).fillna(top["protein_id"])

        plt.figure(figsize=(10, 6))
        plt.barh(labels, top["copy_number"])
        plt.xlabel("Copy number")
        plt.ylabel("Protein (top paralogs)")
        plt.gca().invert_yaxis()
        plt.tight_layout()
        out_png = f"{out_prefix}_top_paralogs.png"
        plt.savefig(out_png, dpi=300)
        plt.close()
        print(f"Wrote bar plot of top paralogs: {out_png}")

    plt.figure(figsize=(6, 4))
    summary["copy_number"].plot.hist(
        bins=range(1, summary["copy_number"].max() + 2),
        align="left"
    )
    plt.xlabel("Copy number per cluster")
    plt.ylabel("Number of proteins")
    plt.tight_layout()
    hist_png = f"{out_prefix}_copy_number_hist.png"
    plt.savefig(hist_png, dpi=300)
    plt.close()
    print(f"Wrote copy-number histogram: {hist_png}")

if __name__ == "__main__":
    if len(sys.argv) == 3:
        faa_path = sys.argv[1]
        cluster_tsv = sys.argv[2]
        p = Path(cluster_tsv)
        out_prefix = p.stem                      

    elif len(sys.argv) == 4:
        faa_path = sys.argv[1]
        cluster_tsv = sys.argv[2]
        out_prefix = Path(sys.argv[3]).stem

    else:
        print("Usage: python summarize_paralogs.py assembly.faa cluster.tsv [out_prefix]")
        sys.exit(1)

    print("FASTA path:", faa_path)
    print("Cluster path:", cluster_tsv)
    print("Out prefix:", out_prefix)

    main(faa_path, cluster_tsv, out_prefix)
