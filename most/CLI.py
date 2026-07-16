import typer
import os
import subprocess
from most.yaml_load import load_yaml, config_cal
from importlib.metadata import version, PackageNotFoundError

import typer
from importlib.metadata import version

app = typer.Typer()


def get_version():
    return version("most")


@app.callback(invoke_without_command=True)
def main(
    ctx: typer.Context,
    config_file: str = typer.Option(None,"--config","-c"),
    version_flag: bool = typer.Option(False,"--version","-v",help="Show version and exit")
):

    if version_flag:
        typer.echo(f"most version {get_version()}")
        raise typer.Exit()

    if config_file is None:
        typer.echo(ctx.get_help())
        raise typer.Exit()

    cfg = load_yaml(config_file)
    method = cfg.method
    path, zpath = config_cal(cfg)

    try:
        subprocess.run(
            [
                "snakemake",
                "--snakefile",
                f"workflow/{method.lower()}.smk",
                "--configfile",
                str(path)
            ],
            check=True
        )

    finally:
        if os.path.exists(path):
            os.remove(path)
        if zpath and os.path.exists(zpath):
            os.remove(zpath)

            