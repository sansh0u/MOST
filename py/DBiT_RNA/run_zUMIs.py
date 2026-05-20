import subprocess
import logging
import os
import typer
import yaml
from yaml_load import load_yaml
from config_utils import get_config
logger = logging.getLogger("toolkit")
from pathlib import Path

def zUMIs(zpath, final_yaml):
    '''
    调用zUMIs,要把in1,in2,out写进去
    '''
    
    cmd = [
        zpath,
        "-c", 
        "-y",
        final_yaml
    ]
    try: ####
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Error during zUMIs: {e}")
        raise

def filled_yaml(cfg, cfg_path):
# -----------------------
# 填参数
# -----------------------
    bc2_start = get_config(cfg, "bc2_start")
    bc2_end = get_config(cfg, "bc2_end")
    bc1_start = get_config(cfg, "bc1_start")
    bc1_end = get_config(cfg, "bc1_end")
    umi_start = get_config(cfg, "umi_start")
    umi_len = get_config(cfg, "umi_len")
    out_dir = get_config(cfg, "Out_dir")
    out_dir = Path(out_dir)
    bc_str = (
    f"BC("
    f"{bc2_start+1}-{bc2_end},"
    f"{bc1_start+1}-{bc1_end}"
    f")"
)
    umi_str = (
    f"UMI("
    f"{umi_start+1}-"
    f"{umi_start+umi_len}"
    f")"
)
    
    with open(cfg_path) as f:
        zcfg = yaml.safe_load(f)
    zcfg["project"] = get_config(cfg, "Project")

    zcfg["sequence_files"]["file1"]["name"] = get_config(cfg, "file1")
    zcfg["sequence_files"]["file1"]["base_definition"] = [
        "cDNA(1,100)"
    ]

    zcfg["sequence_files"]["file2"]["name"] = get_config(cfg, "file2")
    zcfg["sequence_files"]["file2"]["base_definition"] = [
        bc_str,
        umi_str
    ]
    zcfg["reference"]["STAR_index"] = get_config(cfg, "star_index")
    zcfg["reference"]["GTF_file"] = get_config(cfg, "gtf_file")

    zcfg["out_dir"] = get_config(cfg, "Out_dir")

    zcfg["barcodes"]["barcode_file"] =  get_config(cfg, "barcode_file")

    # -----------------------
    # 输出新yaml
    # -----------------------

    with open(out_dir/"filled.yaml", "w") as f:
        yaml.dump(
            zcfg,
            f,
            default_flow_style=False,
            sort_keys=False
        )
