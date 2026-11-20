import gzip
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import anndata as ad
from scipy.sparse import csr_matrix


# =========================
# 1. Input / Output
# =========================

fragment_file = "/data/Fanyt/ATAC_0209_output/ATAC_sorted.bed.gz"

outdir = "/data/Fanyt/ATAC_0209_output/QC"

os.makedirs(outdir, exist_ok=True)


# =========================
# 2. Read fragments
# =========================

print("Reading fragments...")

fragments = pd.read_csv(
    fragment_file,
    sep="\t",
    header=None,
    compression="gzip"
)

# Standard 10x fragments format:
# chr, start, end, barcode, count

fragments.columns = [
    "chrom",
    "start",
    "end",
    "barcode",
    "count"
]

fragments["fragment_length"] = (
    fragments["end"] - fragments["start"]
)

print(f"Total fragment records: {len(fragments):,}")
print(f"Total barcodes: {fragments['barcode'].nunique():,}")


# =========================
# 3. Fragment size distribution
# =========================

plt.figure(figsize=(5, 4))

sns.histplot(
    fragments["fragment_length"],
    bins=100,
    binrange=(0, 700),
    color="steelblue",
    edgecolor=None
)

plt.xlim(0, 700)

plt.xlabel("Fragment size (bp)")
plt.ylabel("Count")

plt.title("Fragment size distribution")

plt.tight_layout()

plt.savefig(
    os.path.join(outdir, "fragment_size_distribution.png"),
    dpi=300,
    bbox_inches="tight"
)

plt.close()


# =========================
# 4. Barcode QC
# =========================

barcode_qc = (
    fragments.groupby("barcode")
    .agg(
        n_fragment=("count", "sum"),
        n_fragment_records=("fragment_length", "count"),
        mean_fragment_length=("fragment_length", "mean"),
        median_fragment_length=("fragment_length", "median")
    )
)

barcode_qc = barcode_qc.sort_values(
    "n_fragment",
    ascending=False
)

print(barcode_qc.head())


# =========================
# 5. Fragment number distribution
# =========================

plt.figure(figsize=(5, 4))

plt.hist(
    barcode_qc["n_fragment"],
    bins=100,
    color="steelblue",
    edgecolor="none"
)

plt.xlabel("Number of fragments")
plt.ylabel("Number of barcodes")

plt.title("Fragments per barcode")

plt.tight_layout()

plt.savefig(
    os.path.join(outdir, "barcode_fragment_distribution.png"),
    dpi=300,
    bbox_inches="tight"
)

plt.savefig(
    os.path.join(outdir, "barcode_fragment_distribution.pdf"),
    bbox_inches="tight"
)

plt.close()


# =========================
# 6. Log-scale fragment distribution
# =========================

plt.figure(figsize=(5, 4))

plt.hist(
    np.log10(barcode_qc["n_fragment"] + 1),
    bins=100,
    color="steelblue",
    edgecolor="none"
)

plt.xlabel("log10(Number of fragments + 1)")
plt.ylabel("Number of barcodes")

plt.title("Barcode fragment distribution")

plt.tight_layout()

plt.savefig(
    os.path.join(outdir, "barcode_fragment_distribution_log.png"),
    dpi=300,
    bbox_inches="tight"
)

plt.close()


# =========================
# 7. Construct AnnData
# =========================

print("Creating AnnData...")

barcodes = barcode_qc.index.tolist()

X = csr_matrix(
    barcode_qc["n_fragment"].values.reshape(-1, 1)
)

adata = ad.AnnData(
    X=X,
    obs=barcode_qc.copy(),
    var=pd.DataFrame(
        index=["total_fragments"]
    )
)

adata.obs_names = barcodes

adata.obs["n_fragment"] = (
    adata.obs["n_fragment"].astype(int)
)

adata.obs["n_fragment_records"] = (
    adata.obs["n_fragment_records"].astype(int)
)

adata.var_names = ["total_fragments"]


# =========================
# 8. Store metadata
# =========================

adata.uns["ATAC_QC"] = {
    "source_file": fragment_file,
    "total_fragment_records": int(len(fragments)),
    "total_barcodes": int(fragments["barcode"].nunique()),
    "genome": "mm10"
}


# =========================
# 9. Save H5AD
# =========================

h5ad_file = os.path.join(
    outdir,
    "ATAC_QC.h5ad"
)

adata.write_h5ad(h5ad_file)

print(f"H5AD saved: {h5ad_file}")


# =========================
# 10. QC summary
# =========================

qc_summary = pd.DataFrame({
    "Metric": [
        "Total fragment records",
        "Total barcodes",
        "Median fragments per barcode",
        "Mean fragments per barcode",
        "Median fragment length",
        "Mean fragment length"
    ],
    "Value": [
        len(fragments),
        fragments["barcode"].nunique(),
        barcode_qc["n_fragment"].median(),
        barcode_qc["n_fragment"].mean(),
        fragments["fragment_length"].median(),
        fragments["fragment_length"].mean()
    ]
})

qc_summary.to_csv(
    os.path.join(outdir, "QC_summary.csv"),
    index=False
)

print(qc_summary)

print("QC completed!")