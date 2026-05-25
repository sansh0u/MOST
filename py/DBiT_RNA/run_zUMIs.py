import subprocess
import logging
import os
import typer
import yaml
from yaml_load import load_yaml
logger = logging.getLogger("toolkit")
from pathlib import Path

def zUMIs(zpath, filledcfg):
    '''
    调用zUMIs,要把in1,in2,out写进去
    '''
    
    cmd = [
        zpath,
        "-c", 
        "-y",
        filledcfg
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
    bc2_start = cfg.runtime.bc2_start
    bc2_end = cfg.runtime.bc2_end
    bc1_start = cfg.runtime.bc1_start
    bc1_end = cfg.runtime.bc1_end
    umi_start = cfg.runtime.umi_start
    umi_len = cfg.runtime.umi_len
    out_dir = cfg.out_dir
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
    zcfg["project"] = cfg.project

    zcfg["sequence_files"]["file1"]["name"] = cfg.sequence_file.file1
    zcfg["sequence_files"]["file1"]["base_definition"] = [
        "cDNA(1,100)"
    ]

    zcfg["sequence_files"]["file2"]["name"] = cfg.sequence_file.file2
    zcfg["sequence_files"]["file2"]["base_definition"] = [
        bc_str,
        umi_str
    ]
    zcfg["reference"]["STAR_index"] = cfg.reference.star_index
    zcfg["reference"]["GTF_file"] = cfg.reference.gtf_file

    zcfg["out_dir"] = cfg.out_dir

    zcfg["barcodes"]["barcode_file"] =  cfg.reference.barcode_file

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
