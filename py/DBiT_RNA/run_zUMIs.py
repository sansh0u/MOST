import subprocess
import logging
import os
import typer
import yaml
from yaml_load import load_yaml
logger = logging.getLogger("toolkit")

def zUMIs(zumis_path, zumisconfig, in1, in2, out):
    '''
    调用zUMIs,要把in1,in2,out写进去
    '''
    if not zumisconfig.exists():
        raise typer.BadParameter("zUMIs config not found, please check the solutions on github")
    config = load_yaml(zumisconfig)
    

    base = os.path.basename(in1)
    if base.endswith("_1"):
        project_name = base[:-2]
    else:
        project_name = base
    
    config.setdefault("sequence_files", {}).setdefault("file1", {})["name"] = in1
    config.setdefault("sequence_files", {}).setdefault("file2", {})["name"] = in2
    config["out_dir"] = out
    config["project"] = project_name

    # 写回
    with open(config, "w") as f:
        yaml.safe_dump(config, f, sort_keys=False)

    print(f"project = {project_name}")
    
    cmd = [
        zumis_path,
        "-c", #-C是运行zumis自己的环境要重新下载，如果可以设置好conda
        "-y",
        config
    ]
    try: ####
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Error during zUMIs: {e}")
        raise

"""
def run_pipeline(config):
    stage = get_config(config, "which_Stage")

    if stage in ["Filtering"]:
        run_filtering(config)

    if stage in ["Filtering", "Mapping"]:
        run_mapping(config)

    if stage in ["Filtering", "Mapping", "Counting"]:
        run_counting(config)

    if stage in ["Filtering", "Mapping", "Counting", "Summarising"]:
        run_stats(config)
        """