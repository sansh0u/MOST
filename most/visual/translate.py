import gzip
import os
import sys
import warnings
import subprocess
import tempfile
import shutil
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import snapatac2 as snap
from scipy.stats import gaussian_kde
# ============================================================
# Utility functions
# ============================================================
def open_fragment_file(filename):
    if str(filename).endswith(".gz"):
        return gzip.open(filename, "rt")
    return open(filename, "rt")
def normalize_column_name(name):
    return (
        str(name)
        .strip()
        .lower()
        .replace(" ", "")
        .replace("-", "")
        .replace("_", "")
    )
def find_column(df, candidates):
    normalized = {
        normalize_column_name(col): col
        for col in df.columns
    }
    for candidate in candidates:
        candidate = normalize_column_name(candidate)
        if candidate in normalized:
            return normalized[candidate]
    return None
def read_chromap_summary(filename):
    """
    Automatically detect CSV or TSV.
    """
    summary = pd.read_csv(
        filename,
        sep=None,
        engine="python"
    )
    summary.columns = [
        str(col).strip()
        for col in summary.columns
    ]
    print("\nChromap summary columns:")
    print(summary.columns.tolist())
    return summary
def save_figure(filename):
    plt.tight_layout()
    plt.savefig(
        filename,
        dpi=300,
        bbox_inches="tight"
    )
    plt.close()
# ============================================================
# Fragment sampling
# ============================================================
def sample_fragments(
    fragment_file,
    max_sample=2000000,
    max_fragment_size=700
):
    """
    Read fragment BED and use reservoir sampling.
    Expected BED format:
        chromosome
        start
        end
        barcode
    Additional columns are allowed.
    Returns:
        fragment_sizes
        barcodes
        total_fragments
    """
    print("\nReading fragments for QC sampling...")
    rng = np.random.default_rng(1234)
    fragment_sizes = []
    barcodes = []
    total_fragments = 0
    with open_fragment_file(fragment_file) as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 4:
                continue
            try:
                start = int(fields[1])
                end = int(fields[2])
                barcode = fields[3]
                fragment_size = end - start
            except Exception:
                continue
            total_fragments += 1
            if len(fragment_sizes) < max_sample:
                fragment_sizes.append(fragment_size)
                barcodes.append(barcode)
            else:
                replacement_index = rng.integers(
                    0,
                    total_fragments
                )
                if replacement_index < max_sample:
                    fragment_sizes[
                        replacement_index
                    ] = fragment_size
                    barcodes[
                        replacement_index
                    ] = barcode
    fragment_sizes = np.asarray(
        fragment_sizes,
        dtype=np.int64
    )
    barcodes = np.asarray(
        barcodes,
        dtype=str
    )
    print(
        f"Total fragments: {total_fragments:,}"
    )
    print(
        f"Fragments sampled: {len(fragment_sizes):,}"
    )
    return (
        fragment_sizes,
        barcodes,
        total_fragments
    )
# ============================================================
# 1. Fragment size distribution
# ============================================================
def plot_fragment_size(
    fragment_sizes,
    outdir,
    max_fragment_size=700
):
    print("\nPlotting fragment size distribution...")
    fragment_sizes = fragment_sizes[
        fragment_sizes <= max_fragment_size
    ]
    bins = np.arange(
        0,
        max_fragment_size + 1
    )
    hist, edges = np.histogram(
        fragment_sizes,
        bins=bins
    )
    plt.figure(figsize=(5, 4))
    plt.plot(
        edges[:-1],
        hist,
        linewidth=1.2
    )
    plt.xlabel(
        "Fragment size (bp)"
    )
    plt.ylabel(
        "Count"
    )
    plt.xlim(
        0,
        max_fragment_size
    )
    plt.title(
        "Fragment Size Distribution"
    )
    save_figure(
        os.path.join(
            outdir,
            "fragment_size_distribution.png"
        )
    )
# ============================================================
# 2. TSS enrichment
# ============================================================
def calculate_tsse(adata, genome):
    print("\nCalculating TSS enrichment...")
    snap.metrics.tsse(
        adata,
        genome
    )
    return adata
def plot_tsse_density(
    adata,
    outdir
):
    print("\nPlotting TSS enrichment...")
    x = np.log10(
        adata.obs["n_fragment"].astype(float) + 1
    )
    y = adata.obs["tsse"].astype(float)
    valid = (
        np.isfinite(x)
        & np.isfinite(y)
    )
    x = x[valid]
    y = y[valid]
    if len(x) < 3:
        print(
            "Warning: insufficient data for KDE. "
            "Using ordinary scatter plot."
        )
        plt.figure(figsize=(5, 4))
        plt.scatter(
            x,
            y,
            s=6,
            alpha=0.5
        )
    else:
        xy = np.vstack([
            x,
            y
        ])
        try:
            z = gaussian_kde(xy)(xy)
            idx = z.argsort()
            x = x.iloc[idx]
            y = y.iloc[idx]
            z = z[idx]
            plt.figure(figsize=(5, 4))
            scat = plt.scatter(
                x,
                y,
                c=z,
                s=6,
                cmap="turbo",
                edgecolor="none"
            )
            plt.colorbar(
                scat,
                label="Density"
            )
        except Exception as e:
            print(
                f"Warning: KDE failed: {e}"
            )
            plt.figure(figsize=(5, 4))
            plt.scatter(
                x,
                y,
                s=6,
                alpha=0.5
            )
    plt.xlabel(
        "log10-transformed unique fragments"
    )
    plt.ylabel(
        "TSS enrichment score"
    )
    plt.title(
        "TSS Enrichment"
    )
    save_figure(
        os.path.join(
            outdir,
            "tss_enrichment.png"
        )
    )
# ============================================================
# 3. Preseq library complexity
# ============================================================
def plot_preseq_complexity(
    fragment_file,
    total_fragments,
    outdir,
    preseq_executable="preseq",
    extrap_factor=5.0,
    step_factor=0.01
):
    """Run Preseq lc_extrap on a BED fragment file and plot library complexity.
    Preseq lc_extrap accepts BED input. The input BED must be sorted by
    chromosome/start/end. Gzip input is temporarily decompressed because
    compatibility with gzipped BED varies across Preseq builds.
    """
    print("\nRunning Preseq library-complexity extrapolation...")
    if shutil.which(preseq_executable) is None:
        raise RuntimeError(
            f"Preseq executable not found: {preseq_executable}. "
            "Please install Preseq and ensure it is in PATH."
        )
    if total_fragments <= 0:
        print("Warning: no fragments available for Preseq.")
        return
    extrap = max(int(total_fragments * extrap_factor), total_fragments + 1)
    step = max(int(total_fragments * step_factor), 1)
    output_file = os.path.join(outdir, "preseq_library_complexity.tsv")
    with tempfile.TemporaryDirectory(prefix="preseq_atac_") as tmpdir:
        bed_for_preseq = os.path.join(tmpdir, "fragments.bed")
        # Preseq requires coordinate-sorted BED. The pipeline's sorted BED
        # should already satisfy this; decompress if needed.
        with open_fragment_file(fragment_file) as fin, open(bed_for_preseq, "w") as fout:
            shutil.copyfileobj(fin, fout)
        cmd = [
            preseq_executable,
            "lc_extrap",
            "-o", output_file,
            "-e", str(extrap),
            "-s", str(step),
            bed_for_preseq
        ]
        print("Preseq command:", " ".join(cmd))
        subprocess.run(cmd, check=True)
    # Preseq output: TOTAL_READS, EXPECTED_DISTINCT, lower CI, upper CI
    df = pd.read_csv(output_file, sep=r"\s+", comment="#")
    if df.shape[1] < 2:
        raise RuntimeError("Unexpected Preseq output format: " + output_file)
    x = pd.to_numeric(df.iloc[:, 0], errors="coerce")
    y = pd.to_numeric(df.iloc[:, 1], errors="coerce")
    plt.figure(figsize=(5, 4))
    plt.plot(x, y, linewidth=1.5, label="Expected distinct fragments")
    if df.shape[1] >= 4:
        lower = pd.to_numeric(df.iloc[:, 2], errors="coerce")
        upper = pd.to_numeric(df.iloc[:, 3], errors="coerce")
        valid = lower.notna() & upper.notna()
        plt.fill_between(
            x[valid], lower[valid], upper[valid],
            alpha=0.2, label="95% CI"
        )
    plt.xlabel("Sequencing depth (fragments)")
    plt.ylabel("Expected distinct fragments")
    plt.title("Preseq Library Complexity")
    plt.legend(frameon=False)
    save_figure(os.path.join(outdir, "preseq_library_complexity.png"))
    print(f"Preseq output: {output_file}")
# ============================================================
# 4. Estimated FRiP
# ============================================================
def plot_estimated_frip(
    summary,
    outdir
):
    print(
        "\nPlotting Estimated FRiP..."
    )
    estfrip_col = find_column(
        summary,
        [
            "estfrip",
            "estimatedfrip",
            "estimated_frip"
        ]
    )
    if estfrip_col is None:
        print(
            "Warning: estfrip column not found."
        )
        print(
            "Available columns:"
        )
        print(
            summary.columns.tolist()
        )
        return
    values = pd.to_numeric(
        summary[estfrip_col],
        errors="coerce"
    )
    values = values.dropna()
    plt.figure(figsize=(5, 4))
    plt.hist(
        values,
        bins=50,
        edgecolor="black",
        linewidth=0.3
    )
    plt.xlabel(
        "Estimated FRiP"
    )
    plt.ylabel(
        "Number of barcodes"
    )
    plt.title(
        "Estimated FRiP Distribution"
    )
    save_figure(
        os.path.join(
            outdir,
            "estimated_frip.png"
        )
    )
# ============================================================
# 5. Import fragments using SnapATAC2
# ============================================================
def import_fragments(
    fragment_file,
    genome,
    threads
):
    print(
        "\nImporting fragments into SnapATAC2..."
    )
    adata = snap.pp.import_fragments(
        fragment_file=fragment_file,
        chrom_sizes=genome,
        sorted_by_barcode=False,
        n_jobs=threads
    )
    print(
        "\nImported AnnData:"
    )
    print(
        adata
    )
    return adata
# ============================================================
# 6. Add tile matrix
# ============================================================
def create_tile_matrix(
    adata,
    bin_size,
    threads
):
    print(
        "\nCreating tile matrix..."
    )
    snap.pp.add_tile_matrix(
        adata,
        bin_size=bin_size,
        exclude_chroms=["chrM"],
        counting_strategy="fragment",
        n_jobs=threads
    )
    print(
        "\nTile matrix created:"
    )
    print(
        adata
    )
    return adata
# ============================================================
# 7. Add Chromap metadata
# ============================================================
def add_chromap_metadata(
    adata,
    summary
):
    print(
        "\nAdding Chromap metadata..."
    )
    barcode_col = find_column(
        summary,
        [
            "barcode",
            "barcodes",
            "cellbarcode"
        ]
    )
    if barcode_col is None:
        print(
            "Warning: barcode column not found in summary."
        )
        print(
            summary.columns.tolist()
        )
        return adata
    summary = summary.copy()
    summary[barcode_col] = (
        summary[barcode_col]
        .astype(str)
    )
    summary = summary.set_index(
        barcode_col
    )
    common_barcodes = (
        adata.obs_names
        .intersection(summary.index)
    )
    print(
        f"Matched barcodes: "
        f"{len(common_barcodes):,}"
    )
    column_mapping = {
        "total": [
            "total"
        ],
        "duplicate": [
            "duplicate",
            "duplicates"
        ],
        "unmapped": [
            "unmapped"
        ],
        "lowmapq": [
            "lowmapq",
            "lowmapquality"
        ],
        "fric": [
            "fric"
        ],
        "estfrip": [
            "estfrip",
            "estimatedfrip"
        ],
        "cachehit": [
            "cachehit"
        ],
        "numcacheslots": [
            "numcacheslots"
        ]
    }
    output_names = {
        "total": "chromap_total",
        "duplicate": "chromap_duplicate",
        "unmapped": "chromap_unmapped",
        "lowmapq": "chromap_lowmapq",
        "fric": "chromap_fric",
        "estfrip": "estimated_frip",
        "cachehit": "chromap_cachehit",
        "numcacheslots": "chromap_numcacheslots"
    }
    for metric, candidates in column_mapping.items():
        column = find_column(
            summary,
            candidates
        )
        if column is None:
            continue
        output_name = output_names[metric]
        adata.obs.loc[
            common_barcodes,
            output_name
        ] = pd.to_numeric(
            summary.loc[
                common_barcodes,
                column
            ],
            errors="coerce"
        ).values
    # Calculate duplicate rate if possible
    if (
        "chromap_total" in adata.obs.columns
        and "chromap_duplicate" in adata.obs.columns
    ):
        adata.obs[
            "duplicate_rate"
        ] = (
            adata.obs["chromap_duplicate"]
            /
            adata.obs["chromap_total"]
        )
    return adata
# ============================================================
# Main function
# ============================================================
def atac_qc_h5ad(
    fragment_file,
    summary_file,
    h5ad_file,
    outdir,
    bin_size=5000,
    threads=8,
    genome_name="mm10",
    max_fragment_size=700,
    max_sample_fragments=2000000,
    preseq_extrap_factor=5.0,
    preseq_step_factor=0.01,
    preseq_executable="preseq"
):
    # --------------------------------------------------------
    # Paths
    # --------------------------------------------------------
    os.makedirs(
        outdir,
        exist_ok=True
    )
    print("=" * 60)
    print(
        "ATAC QC and H5AD conversion"
    )
    print(
        "SnapATAC2 version:",
        snap.__version__
    )
    print("=" * 60)
    # --------------------------------------------------------
    # Read Chromap summary
    # --------------------------------------------------------
    summary = read_chromap_summary(
        summary_file
    )
    # --------------------------------------------------------
    # Sample fragments
    # --------------------------------------------------------
    (
        fragment_sizes,
        barcodes,
        total_fragments
    ) = sample_fragments(
        fragment_file,
        max_sample=max_sample_fragments,
        max_fragment_size=max_fragment_size
    )
    # --------------------------------------------------------
    # QC 1: Fragment size
    # --------------------------------------------------------
    plot_fragment_size(
        fragment_sizes,
        outdir,
        max_fragment_size=max_fragment_size
    )
    # --------------------------------------------------------
    # QC 2: Estimated FRiP
    # --------------------------------------------------------
    plot_estimated_frip(
        summary,
        outdir
    )
    # --------------------------------------------------------
    # QC 3: Preseq library complexity
    # --------------------------------------------------------
    plot_preseq_complexity(
        fragment_file,
        total_fragments,
        outdir,
        preseq_executable=preseq_executable,
        extrap_factor=preseq_extrap_factor,
        step_factor=preseq_step_factor
    )
    # --------------------------------------------------------
    # Import fragments
    # --------------------------------------------------------
    if genome_name == "mm10":
        genome = snap.genome.mm10
    else:
        raise ValueError(
            f"Unsupported genome: {genome_name}"
        )
    adata = import_fragments(
        fragment_file,
        genome,
        threads
    )
    # --------------------------------------------------------
    # QC 4: TSS enrichment
    # --------------------------------------------------------
    adata = calculate_tsse(
        adata,
        genome
    )
    plot_tsse_density(
        adata,
        outdir
    )
    # --------------------------------------------------------
    # Create tile matrix
    # --------------------------------------------------------
    adata = create_tile_matrix(
        adata,
        bin_size,
        threads
    )
    # --------------------------------------------------------
    # Add Chromap metadata
    # --------------------------------------------------------
    adata = add_chromap_metadata(
        adata,
        summary
    )
    # --------------------------------------------------------
    # Save H5AD
    # --------------------------------------------------------
    print(
        "\nSaving H5AD..."
    )
    adata.write(
        h5ad_file
    )
    # --------------------------------------------------------
    # Finish
    # --------------------------------------------------------
    print("\n" + "=" * 60)
    print(
        "Pipeline completed successfully!"
    )
    print(
        f"Fragments: {total_fragments:,}"
    )
    print(
        f"Barcodes: {adata.n_obs:,}"
    )
    print(
        f"Features: {adata.n_vars:,}"
    )
    print(
        f"QC output: {outdir}"
    )
    print(
        f"H5AD output: {h5ad_file}"
    )
    print("=" * 60)
