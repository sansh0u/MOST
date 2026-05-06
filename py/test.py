from toolkit.py.QC.scan_bc import scan
import typer
import os
import logging
from yaml_load import load_yaml, get_config, config_cal,setup_logger
from DBiT_RNA.run_zUMIs import zUMIs
from preprocess.qc import filter
from ATAC.chromap import chromap, sort_bed
from DBiT_RNA.stpipeline import stpipeline
from preprocess.bc_process import dbit_bc, atac_bc 
from pathlib import Path
import yaml
from preprocess.visual import detect_tissue_pixels
from preprocess.scan_bc import scan

def plot_fragment_distribution(config, out, bins=500, max_len=1000, log_y=False):
    df = pd.read_csv(
    config,
    sep="\t",
    header=None,
    names=["chr", "start", "end", "name", "score"],
    compression="infer"
    )

# 计算长度并过滤
    df["length"] = df["end"] - df["start"]
    df = df[df["length"] > 0]
    df = df[df["length"] <= max_len]

    plt.figure(figsize=(6, 4))

    counts, bin_edges = np.histogram(
        df["length"],
        bins=bins,
        range=(0, max_len)
    )

    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    
    plt.plot(bin_centers, counts)

    if log_y:
        plt.yscale("log")

    plt.xlim(0, max_len)

    plt.xlabel("Fragment size")
    plt.ylabel("Count")
    plt.title("Fragment Size Distribution")

    plt.tight_layout()

    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)

    plt.savefig(out, dpi=300)
    plt.close()


plot_fragment_distribution(df, "fragment.png")

config = [
    {"file2": "/mnt/d/prog/RNA2_1209_Illu_2.fq.gz"},
    {"Barcode": "/home/sanshou/project/toolkit/toolkit/barcode/20240614_2500barcode_AB_update.txt"}

]
scan(config)
app = typer.Typer(help = """
 pipeline toolkit

Commands:

  run     Run main pipeline
  zumis   Run zUMIs pipeline

Examples:

  toolkit run --config config.yaml

  toolkit zumis -dbit -in1 R1.fq -in2 R2.fq -out outdir
"""
, no_args_is_help = True)

@app.command(no_args_is_help = True)
def zumis(
    zpath: str = typer.Option(None, "--l", help="Path to zUMIs.sh"),
    config: str = typer.Option(None, "--config", help="Custom YAML")
):