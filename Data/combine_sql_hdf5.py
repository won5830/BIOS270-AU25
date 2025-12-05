#!/usr/bin/env python

import argparse
import numpy as np
import h5py

from query_bacteria_db import BacteriaDatabase


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--database_path",
        type=str,
        default="bacteria.db",
        help="Path to bacteria.db",
    )
    parser.add_argument(
        "--h5_path",
        type=str,
        default="/farmshare/home/classes/bios/270/data/processed_bacteria_data/protein_embeddings.h5",
        help="Path to protein_embeddings.h5",
    )
    parser.add_argument(
        "--record_id",
        type=str,
        required=True,
        help="record_id to query from gff table",
    )
    parser.add_argument(
        "--metric",
        type=str,
        choices=["mean", "mean_mid"],
        required=True,
        help="Which embedding metric to use (mean or mean_mid)",
    )
    parser.add_argument(
        "--output_path",
        type=str,
        default="embeddings.npy",
        help="Where to save resulting (N, D) matrix as .npy",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    db = BacteriaDatabase(args.database_path)
    protein_ids = db.get_protein_ids_from_record_id(args.record_id)
    protein_ids = [pid for pid in protein_ids if pid is not None]

    if not protein_ids:
        raise ValueError(f"No protein_ids found for record_id={args.record_id}")

    with h5py.File(args.h5_path, "r") as f:
        # protein_ids: (num_proteins,)
        h5_protein_ids = f["protein_ids"][:]

        # metric 이름 → 실제 HDF5 dataset 이름 매핑
        if args.metric == "mean":
            ds_name = "mean_embeddings"
        elif args.metric == "mean_mid":
            ds_name = "mean_mid_embeddings"
        else:
            raise ValueError(f"Unknown metric: {args.metric}")

        if ds_name not in f:
            raise KeyError(f"Dataset '{ds_name}' not found in HDF5. Available: {list(f.keys())}")

        emb_ds = f[ds_name]   # (num_proteins, 164) 
        if isinstance(h5_protein_ids[0], bytes):
            h5_protein_ids = [pid.decode("utf-8") for pid in h5_protein_ids]
        else:
            h5_protein_ids = [str(pid) for pid in h5_protein_ids]

        id2idx = {pid: i for i, pid in enumerate(h5_protein_ids)}

        indices = []
        for pid in protein_ids:
            if pid in id2idx:
                indices.append(id2idx[pid])

        if not indices:
            raise ValueError(
                f"No protein_ids for record_id={args.record_id} were found in HDF5 file"
            )
        indices = sorted(set(indices))

        embeddings = emb_ds[indices, :]        # shape (N, D)

    np.save(args.output_path, embeddings)
    print(
        f"Saved embeddings for record_id={args.record_id} "
        f"using metric={args.metric} with shape {embeddings.shape} "
        f"to {args.output_path}"
    )


if __name__ == "__main__":
    main()
