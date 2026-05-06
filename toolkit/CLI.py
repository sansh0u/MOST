#!/user/bin/env python
import typer
import os
import logging
from yaml_load import load_yaml, get_config, config_cal
from logger import setup_logger
from DBiT_RNA.run_zUMIs import zUMIs
from qc import filter
from ATAC.bc_process import bc_pr
from ATAC.chromap import chromap, sort_bed
from DBiT_RNA.stpipeline import stpipeline
from DBiT_RNA.bc_pro import bc_pro 
from pathlib import Path
import yaml
from Pixel_identification import detect_tissue_pixels
from scan_bc import scan



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
    config_path: str = typer.Option(..., "--config", help="Pipeline config YAML"),
):
    """Run main pipeline"""

    print("Pipeline started")

    config = load_yaml(config_path)
    os.makedirs(get_config(config, "dir"), exist_ok=True)
    method = get_config(config, "Method")
    
    if method == "ATAC":
        bc2_loc, bc1_loc, read_len = scan(config)
        config = config_cal(config, bc2_loc, bc1_loc, read_len)
        print(config)
        filter(config)
        bc_pr(config)
        chromap(config)
        sort_bed(config)

    elif method == "RNA":
        
        bc2_loc, bc1_loc, read_len = scan(config)
        config = config_cal(config, bc2_loc, bc1_loc, read_len)
        print(config)
        #filter(config)
        #bc_pro(config)
        stpipeline(config)

@app.command(no_args_is_help = True)
def zumis(
    zpath: str = typer.Option(None, "--l", help="Path to zUMIs.sh"),
    dbit: bool = typer.Option(False, "--dbit", help="DBiT mode"),
    patho: bool = typer.Option(False, "--patho", help="Patho mode"),
    rna: bool = typer.Option(False, "--rna", help="RNA mode"),
    in1: str = typer.Option(None, "--in1"),
    in2: str = typer.Option(None, "--in2"),
    out: str = typer.Option(None, "--out"),
    config: str = typer.Option(None, "--config", help="Custom YAML"),
    illumina: bool = typer.Option(False, "--illumina"),
    pcr: bool = typer.Option(False, "--pcr"),   
):
    """Run zUMIs pipeline"""

    BASE_DIR = Path(__file__).resolve().parent
    CONFIG_FILE = BASE_DIR / "config" / ".config.yaml"
    # ========= zUMIs 路径处理 =========
    if zpath:
        zpath = Path(zpath)

        if zpath.is_dir():
            zpath = zpath / "zUMIs.sh"

        if not zpath.exists():
            raise typer.BadParameter(f"zUMIs not found: {zpath}")

        if zpath.name != "zUMIs.sh":
            raise typer.BadParameter("Please provide zUMIs.sh or its directory")

        CONFIG_FILE.parent.mkdir(parents=True, exist_ok=True)

        with open(CONFIG_FILE, "w") as f:
            yaml.safe_dump({"zumis_path": str(zpath)}, f)
        
        #  如果只是设置路径（没有任何运行参数），直接退出
        if not any([dbit, patho, rna, in1, in2, out, config]):
            print("zUMIs path saved. Ready to use.")
            raise typer.Exit()

    else:
        if CONFIG_FILE.exists():
            with open(CONFIG_FILE) as f:
                cfg = yaml.safe_load(f) or {}

            if "zumis_path" in cfg:
                zpath = Path(cfg["zumis_path"])
            else:
                raise typer.BadParameter("No zUMIs path found, please provide --l once")
        else:
            raise typer.BadParameter("Please provide --l (zUMIs path) at least once")

    if not zpath.exists():
        raise typer.BadParameter(f"zUMIs not found: {zpath}")

    print(f"Using zUMIs: {zpath}")

    # ========= mode 检查 =========
    if config:
        if any([dbit, patho, rna]):
            raise typer.BadParameter("--config cannot be used with --dbit/--patho/--rna")

    mode_count = sum([dbit, patho, rna])

    if mode_count == 0:
        raise typer.BadParameter("Select one mode: --dbit / --patho / --rna")

    if mode_count > 1:
        raise typer.BadParameter("Only one mode allowed")

    if not all([in1, in2, out]):
        raise typer.BadParameter("Require --in1 --in2 --out")

    # ========= 模式 =========
    if dbit:
        zUMIsconfig = BASE_DIR / "config" / "DBit.yaml"
        mode = "DBiT"
    elif patho:
        zUMIsconfig = BASE_DIR / "config" / "Patho.yaml"
        mode = "Patho"
    else:
        zUMIsconfig = BASE_DIR / "config" / "RNA.yaml"
        mode = "RNA"

    print(f"Running zUMIs in {mode} mode")


    zUMIs(zpath, zUMIsconfig, in1, in2, out)

@app.command(no_args_is_help = True)
def astro(
    apath: str = typer.Option(None, "--l", help="Path to astro"),
    config: str = typer.Option(None, "--config", help="Custom YAML")
):
    """Run astro pipeline"""

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

