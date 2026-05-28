"""
YAML configuration file loading and validation module
This module provides functionality for loading YAML configuration files,
validating configuration fields, and printing configuration summaries.
"""

import yaml
import logging
import os
from pathlib import Path
from most.preprocess.scan_bc import scan,scan_len
from most.config_utils import Config
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


def load_yaml(cfg_path):
    """
    Load YAML configuration file. 加载并判断YAML配置文件是否存在并可解析(已完成)
    """
    try:

        with open(
            cfg_path,
            "r",
            encoding="utf-8"
        ) as f:

            data = yaml.safe_load(f)

        if data is None:

            raise ValueError(
                "YAML file is empty"
            )

        cfg = Config(**data)

        return cfg

    except FileNotFoundError:

        logging.error(
            f"Configuration file not found: "
            f"{cfg_path}"
        )

        raise

    except yaml.YAMLError as e:

        logging.error(
            f"YAML parse error: {e}"
        )

        raise

    except Exception as e:

        logging.error(
            f"Failed to load config: {e}"
        )

        raise

def convert_range(r):
    # r like "23-32"
    r = r.strip("()")
    start, end = map(int, r.split("-"))

    start0 = start - 1
    length = end - start + 1

    return start0, length    

def parse_pair(s):
    s = s.strip("()")
    r1, r2 = s.split(",")

    bc2_start, bc_len = convert_range(r1)
    bc1_start, bc_len = convert_range(r2)

    return bc2_start, bc1_start, bc_len

def config_cal(cfg):
    """
    计算配置文件中的参数,如果advance里没有则默认
    """ 
    method = cfg.method

    umi = cfg.advanced.UMI

    bc = cfg.advanced.BC

    if umi is not None and bc is not None:

        umi_start, umi_len = convert_range(umi)

        bc2_start, bc1_start, bc_len = parse_pair(bc)

        read_len = scan_len(cfg)

    else:

        print(
            "No UMI or BC found in config. Enter auto mode."
        )

        result = scan(cfg, method)

        bc2_start = result["bc2"]

        bc1_start = result["bc1"]

        bc_len = result["bc_len"]

        read_len = result["read_len"]

        umi_start = result["umi_start"]

        umi_len = result["umi_len"]

    bc2_end = bc2_start + bc_len

    bc1_end = bc1_start + bc_len

    primer = cfg.advanced.primer

    linker1 = cfg.advanced.linker1

    linker2 = cfg.advanced.linker2

    if read_len == 100:

        linker2 = ""

    k1 = len(linker1)

    k2 = len(linker2)

    restrictleft1 = bc1_end + k1 + umi_len

    restrictleft2 = bc2_end + k2 + umi_len

    if umi_start < bc2_start:

        seq_start = bc1_end + k1  + 19
    else:
        seq_start = bc1_end + k1 + umi_len + 19

    cfg.runtime.k1 = k1

    cfg.runtime.k2 = k2

    cfg.runtime.bc2_start = bc2_start

    cfg.runtime.bc2_end = bc2_end

    cfg.runtime.bc1_start = bc1_start

    cfg.runtime.bc1_end = bc1_end

    cfg.runtime.restrictleft1 = \
        restrictleft1

    cfg.runtime.restrictleft2 = \
        restrictleft2

    cfg.runtime.seq_start = seq_start

    cfg.runtime.umi_start = umi_start

    cfg.runtime.umi_len = umi_len

    return cfg
    