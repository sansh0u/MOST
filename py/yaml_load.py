"""
YAML configuration file loading and validation module
This module provides functionality for loading YAML configuration files,
validating configuration fields, and printing configuration summaries.
"""

import yaml
import logging
import os
from pathlib import Path
from preprocess.scan_bc import scan

logger = logging.getLogger("toolkit")



def setup_logger():
    """
    Setup logging configuration.
    """

    # 创建 logs 目录（如果不存在）
    os.makedirs("logs", exist_ok=True)

    # 创建 logger
    logger = logging.getLogger("toolkit")

    # 允许记录 INFO 以上日志
    logger.setLevel(logging.INFO)

    pipeline_handler = logging.FileHandler("logs/pipeline.log")
    pipeline_handler.setLevel(logging.INFO)

    error_handler = logging.FileHandler("logs/error.log")
    error_handler.setLevel(logging.ERROR)

    formatter = logging.Formatter(
        "%(asctime)s - %(levelname)s - %(message)s"
    )

    pipeline_handler.setFormatter(formatter)
    error_handler.setFormatter(formatter)

    # 添加 handler
    logger.addHandler(pipeline_handler)
    logger.addHandler(error_handler)

    return logger

def get_config(cfg, key, default=None):
    if isinstance(cfg, dict):
        if key in cfg and cfg[key] is not None:
            return cfg[key]
        for v in cfg.values():
            result = get_config(v, key, default)
            if result is not None:
                return result
    elif isinstance(cfg, list):
        for item in cfg:
            result = get_config(item, key, default)
            if result is not None:
                return result
    return default

def load_yaml(cfg_path):
    """
    Load YAML configuration file. 加载并判断YAML配置文件是否存在并可解析(已完成)
    """
    try:
        with open(cfg_path, "r", encoding="utf-8") as f:
            cfg = yaml.safe_load(f)
        
    except FileNotFoundError:
        logging.error(f"Configuration file not found: {cfg_path}")
        return None
    except yaml.YAMLError as e:
        logging.error(f"Error parsing YAML file {cfg_path}: {e}")
        return None
    except Exception as e:
        #logger.error(f"Unexpected error reading config: {e}")
        return None
    if cfg is None:
        #logger.error("YAML file is empty.")
        return None
    #直接输出config，要什么调用的时候自己取
    default_barcode = ( Path(__file__).resolve().parent.parent / "barcode" / "20240614_2500barcode_AB_update.txt")
    barcode_file = get_config(cfg, "barcode_file",default_barcode)
    cfg["barcode_file"] = barcode_file
    return cfg

def convert_range(r):
    # r like "23-32"
    start, end = map(int, r.split("-"))

    start0 = start - 1
    length = end - start + 1

    return start0, length    

def parse_pair(s):
    s = s.strip("()")
    r1, r2 = s.split(",")

    bc1_start, bc_len = convert_range(r1)
    bc2_start, bc_len = convert_range(r2)

    return bc1_start, bc_len, bc2_start, bc_len

def config_cal(cfg, method):
    """
    计算配置文件中的参数,如果advance里没有则默认
    """ 
    valid_methods = {"ATAC", "RNA"}
    if method not in valid_methods:
        raise ValueError(f"Unknown method: {method}")

    umi = get_config(cfg, "UMI")
    bc = get_config(cfg, "BC")

    if bc and umi is not None:
        umi_start, umi_len = convert_range(cfg["UMI"])
        bc1_start, bc_len, bc2_start, bc_len = parse_pair(cfg["BC"])
    else:
        print("No UMI or BC found in config file. Enter automatic mode.")
        result = scan(cfg, method)
        bc2_start = result["bc2"]
        bc1_start = result["bc1"]
        bc_len = result["bc_len"]
        read_len = result["read_len"]
        umi_start = result["umi_start"]
        umi_len = result["umi_len"]
    
    bc2_end = bc2_start + bc_len
    bc1_end = bc1_start + bc_len
    primer = get_config(cfg, "primer", "CAAGCGTTGGCTTCTCGCATCT")
    linker1 = get_config(cfg, "linker1", "GTGGCCGATGTTTCGCATCGGCGTACGACT")
    linker2 = get_config(cfg, "linker2", "ATCCACGTGCTTGAGAGGCCAGAGCATTCG")

    if read_len == 100:
        linker2 = ""
    
    k1 = len(linker1)
    k2 = len(linker2)

    restrictleft1 = bc1_end + k1 + umi_len
    restrictleft2 = bc2_end + k2 + umi_len
    seq_start = bc1_end + k1 + umi_len + 19
    cfg['Advanced'] = {
        'k1': k1,
        'k2': k2,
        'bc2_start': bc2_start,
        'bc2_end': bc2_end,
        'bc1_start': bc1_start,
        'bc1_end': bc1_end,
        'restrictleft1': restrictleft1,
        'restrictleft2': restrictleft2,
        'seq_start': seq_start,
        'umi_start': umi_start,
        'umi_len': umi_len,
        'linker1': linker1,
        'linker2': linker2,
        'primer': primer
    }
    
    return cfg


