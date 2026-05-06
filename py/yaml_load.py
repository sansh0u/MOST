"""
YAML configuration file loading and validation module
This module provides functionality for loading YAML configuration files,
validating configuration fields, and printing configuration summaries.
"""

import yaml
import logging
import os
from pathlib import Path

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

def get_config(config, key, default=None):
    if isinstance(config, dict):
        if key in config and config[key] is not None:
            return config[key]
        for v in config.values():
            result = get_config(v, key, default)
            if result is not None:
                return result
    elif isinstance(config, list):
        for item in config:
            result = get_config(item, key, default)
            if result is not None:
                return result
    return default

def load_yaml(config_path):
    """
    Load YAML configuration file. 加载并判断YAML配置文件是否存在并可解析(已完成)
    """
    try:
        with open(config_path, "r", encoding="utf-8") as f:
            config = yaml.safe_load(f)
        
    except FileNotFoundError:
        logging.error(f"Configuration file not found: {config_path}")
        return None
    except yaml.YAMLError as e:
        logging.error(f"Error parsing YAML file {config_path}: {e}")
        return None
    except Exception as e:
        #logger.error(f"Unexpected error reading config: {e}")
        return None
    if config is None:
        #logger.error("YAML file is empty.")
        return None
    #直接输出config，要什么调用的时候自己取
    BARCODE_FILE = get_config(config, "Barcode")
    if BARCODE_FILE is None:
        barcode_file = ( Path(__file__).resolve().parent.parent / "barcode" / "20240614_2500barcode_AB_update.txt")
        config["Barcode"] = barcode_file
    return config


def config_cal(config, result):
    """
    计算配置文件中的参数,如果advance里没有则默认
    """ 

    bc2_start = result["bc2"]
    bc1_start = result["bc1"]
    read_len = result["read_len"]
    bc2_end = bc2_start + 8
    bc1_end = bc1_start + 8 #8bp barcode
    restrictleft1 = bc1_end + 40
    restrictleft2 = bc2_end + 40
    seq_start = bc1_end + 40
    primer = len(get_config(config, "primer", "CAAGCGTTGGCTTCTCGCATCT"))
    linker1 = get_config(config, "linker1", "GTGGCCGATGTTTCGCATCGGCGTACGACT")
    linker2 = get_config(config, "linker2", "ATCCACGTGCTTGAGAGGCCAGAGCATTCG")
    if read_len == 100:
        linker1 = ""
    if bc2_start == 1:
        linker1 = "AGATGTGTATAAGAGACAGCATCGGCGTACGACT"
        linker2 = "CGAATGCTCTGGCCTCTCAAGCACGTGGAT"
    
    k1 = len(linker1)
    k2 = len(linker2)
    if bc2_start == primer :
        umi_start = bc1_end + k1 
    elif bc2_start == primer + 10:
        umi_start = primer 
    else:
        umi_start = get_config(config, "UMI") 
    restrictleft1 = bc1_end + k1 + 10
    restrictleft2 = bc2_end + k2 + 10
    seq_start = bc1_end + k1 + 19
    config['Advanced'] = {
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
        'linker1': linker1,
        'linker2': linker2
    }
    
    return config



