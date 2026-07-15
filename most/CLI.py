import typer
import os
import subprocess
from most.yaml_load import load_yaml, config_cal

app = typer.Typer(help = """
 pipeline toolkit

Commands:

  run     Run main pipeline
  zumis   Run zUMIs pipeline

Examples:

  most run --config config.yaml

  most zumis -dbit -in1 R1.fq -in2 R2.fq -out outdir
"""
, no_args_is_help = True)

@app.command(no_args_is_help = True)
def run(config_file):
    cfg = load_yaml(config_file)
    method = cfg.method
    path,zpath = config_cal(cfg)
    success = False 
    try:
        subprocess.run(
    [
        "snakemake",
        "--snakefile",
        f"workflow/{method.lower()}.smk",
        "--configfile",
        str(path)
    ]
)
        success = True
    finally:
        if success:
            if os.path.exists(path):
                os.remove(path)
            if zpath and os.path.exists(zpath):
                os.remove(zpath)
            