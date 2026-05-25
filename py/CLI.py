
from numpy._core.defchararray import zfill
import typer
import os
import logging
from yaml_load import load_yaml,  config_cal,setup_logger
from DBiT_RNA.run_zUMIs import zUMIs, filled_yaml
from preprocess.qc import filter
from ATAC.chromap import chromap, sort_bed
from DBiT_RNA.stpipeline import stpipeline
from preprocess.bc_process import dbit_bc, atac_bc 
from pathlib import Path
import yaml
from preprocess.visual import detect_tissue_pixels
from DBiT_RNA.qc_adapt import qc_adapt
from preprocess.scan_bc import check_adapter

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

setup_logger()
logger = logging.getLogger("toolkit")

# =========================
# 主 pipeline
# =========================
@app.command(no_args_is_help = True)
def run(
    cfg_path: str = typer.Option(..., "--config", help="Pipeline config YAML"),
):
    """Run main pipeline"""

    print("Pipeline started")

    cfg = load_yaml(cfg_path)
    os.makedirs(cfg.out_dir, exist_ok=True)
    method = cfg.method
    cfg = config_cal(cfg)
    

    if method == "ATAC":
        #print(cfg)
        filter(cfg,method)
        atac_bc(cfg)
        chromap(cfg)
        sort_bed(cfg)

    elif method == "RNA":
        #print(cfg)
        qc_adapt(cfg)
        filter(cfg,method)
        dbit_bc(cfg)
        stpipeline(cfg)


@app.command(no_args_is_help=True)
def zumis(
    zpath: str = typer.Option(None, "--l", help="Path to zUMIs.sh"),
    cfg_path: str = typer.Option(..., "--config", help="YAML config file")
):
    """Run zUMIs pipeline"""
    print("zUMIs Pipeline started")
    
    BASE_DIR = Path(__file__).resolve().parent
    cfg_file = BASE_DIR / "config" / ".config.yaml"
    zcfg_path = BASE_DIR / "config"  / "RNA.yaml"
    
    # ========= zUMIs 路径处理 =========
    if zpath:
        zpath = Path(zpath)

        if zpath.is_dir():
            zpath = zpath / "zUMIs.sh"

        if not zpath.exists():
            raise typer.BadParameter(f"zUMIs not found: {zpath}")

        if zpath.name != "zUMIs.sh":
            raise typer.BadParameter(
                "Please provide zUMIs.sh or its directory"
            )

        cfg_file.parent.mkdir(parents=True, exist_ok=True)

        with open(cfg_file, "w") as f:
            yaml.safe_dump({"zumis_path": str(zpath)}, f)

        print("zUMIs path saved.")

    else:
        if cfg_file.exists():
            with open(cfg_file) as f:
                zcfg = yaml.safe_load(f) or {}

            if "zumis_path" in zcfg:
                zpath = Path(zcfg["zumis_path"])
            else:
                raise typer.BadParameter(
                    "No zUMIs path found, please provide --l once"
                )
        else:
            raise typer.BadParameter(
                "Please provide --l (zUMIs path) at least once"
            )

    if not zpath.exists():
        raise typer.BadParameter(f"zUMIs not found: {zpath}")

    print(f"Using zUMIs: {zpath}")

    # ========= config 检查 =========
    cfg_path = Path(cfg_path)

    if not cfg_path.exists():
        raise typer.BadParameter(f"Config file not found: {cfg_path}")

    print(f"Using config: {cfg_path}")

    # ========= 运行 =========
    cfg = load_yaml(cfg_path)
    filledcfg= cfg.out_dir / "filled.yaml"
    os.makedirs(cfg.out_dir, exist_ok=True)
    method = cfg.method
    cfg = config_cal(cfg)
    print(cfg)
    qc_adapt(cfg)
    check_adapter(cfg)
    filled_yaml(cfg, zcfg_path)
    zUMIs(zpath,filledcfg)

@app.command(no_args_is_help = True)
def astro(
    cfg_path: str = typer.Option(None, "--config", help="Custom YAML")
):
    """Run astro pipeline
    """
    print("Astro Pipeline started")
    """Run astro pipeline
        ASTRO --R1 R1.fq --R2 R2.fq \
    --barcode_file spatial_barcodes.txt \
    --gtffile hsa.no_piRNA.gtf --starref StarIndex/ \
    --PrimerStructure AAGCAGTGGTATCAACGCAGAGTGAATGGG_b_A{10}N{150} \
    --StructureUMI CAAGCGTTGGCTTCTCGCATCT_10 \
    --StructureBarcode 20_ATCCACGTGCTTGAGAGGCCAGAGCATTCG:ATCCACGTGCTTGAGAGGCCAGAGCATTCG...GTGGCCGATGTTTCGCATCGGCGTACGACT \
    --threadnum 16 \
    --steps 7 \
    --outputfolder output/"""
    cfg = load_yaml(cfg_path)
    os.makedirs(cfg.out_dir, exist_ok=True)
    method = cfg.method
    cfg = config_cal(cfg)
    

@app.command(no_args_is_help = True)
def matlab(
    input: str = typer.Option(None, "--in", help="Path to input file"),
    output: str = typer.Option(None, "--out", help="Path to output file")
):
    """Run matlab pipeline"""

@app.command(no_args_is_help = True)
def visual(
    image_path: str = typer.Option(..., "--in", help="Path to input image"),
    output_file: str = typer.Option("position.txt", "--out", help="Path to output file"),
    pixel: int = typer.Option(50, "--p", help="Pixel size"),
    threshold: int = typer.Option(0, "--t", help="Threshold")
):
    """Run detect_tissue_pixels pipeline"""
    detect_tissue_pixels(image_path, output_file, pixel, threshold)

# =========================
#  入口
# =========================
if __name__ == "__main__":
    app()

