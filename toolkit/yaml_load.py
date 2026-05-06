"""
YAML configuration file loading and validation module
This module provides functionality for loading YAML configuration files,
validating configuration fields, and printing configuration summaries.
"""

import yaml
import logging

logger = logging.getLogger("toolkit")


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
    return config
"""
def method_check(config):
    """
    #检查配置文件中的method
"""
    #logger.info("Start method check.")
    method = get_config(config, "Method")
    method_dict = {
    "dbit": DBit_seq,
    "atac": ATAC_seq,
    "co_atac": co_ATAC,
    "co_rna": co_RNA,
    "patho_dbit": Patho_DBit,
    "patho_atac": Patho_ATAC
}
    if method not in method_dict:
        raise ValueError(f"Unknown Method: {method}")
    #logger.info("Start merge config.")
    config = merge_config(config, method_dict[method])

    if method in ["atac", "co_atac", "patho_atac","dbit", "patho_dbit"]:
        config_cal(config)
    elif method in ["co_rna"]:
        co_rna_cal(config)
    return config

#config, method = load_yaml("/home/sanshou/projects/tool/dbit/zUMIs.yaml")
#print(method)
    
def merge_config(config, default_config):
    for key, value in default_config.items():
        if key not in config or config[key] is None:
            config[key] = value
        elif isinstance(value, dict):
            config[key] = merge_config(config.get(key, {}), value)
    return config
"""


def config_cal(config, bc2_loc, bc1_loc, read_len):
    """
    计算配置文件中的参数,如果advance里没有则默认
    """ 
    bc2_start = bc2_loc 
    bc1_start = bc1_loc 
    bc2_end = bc2_loc + 8
    bc1_end = bc1_loc + 8 #8bp barcode
    restrictleft1 = bc1_end + 40
    restrictleft2 = bc2_end + 40
    seq_start = bc1_end + 40
    primer = len(get_config(config, "primer", "CAAGCGTTGGCTTCTCGCATCT"))
    linker1 = get_config(config, "linker1", "GTGGCCGATGTTTCGCATCGGCGTACGACT")
    linker2 = get_config(config, "linker2", "ATCCACGTGCTTGAGAGGCCAGAGCATTCG")
    if read_len == 100:
        linker1 = ""
    if bc2_loc == 1:
        linker1 = "AGATGTGTATAAGAGACAGCATCGGCGTACGACT"
        linker2 = "CGAATGCTCTGGCCTCTCAAGCACGTGGAT"
    
    k1 = len(linker1)
    k2 = len(linker2)
    if bc2_loc == primer :
        umi_start = bc1_end + k1 
    elif bc2_loc == primer + 10:
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

    """
    k1 = len(get_config(config, "linker1"))
    k2 = len(get_config(config, "linker2"))
    bc2_start = get_config(config, "primer") + get_config(config, "UMI")
    bc2_end = bc2_start + 8
    bc1_start = bc2_end + k2
    bc1_end = bc1_start + 8 #8bp barcode
    restrictleft1 = bc1_end + k2 + 10
    restrictleft2 = bc2_end + k1 + 10
    seq_start = bc1_end + k1 + 19
    umi_start = get_config(config, "primer") + 1
    config['preprocess'] = {
        'k1': k1,
        'k2': k2,
        'bc2_start': bc2_start,
        'bc2_end': bc2_end,
        'bc1_start': bc1_start,
        'bc1_end': bc1_end,
        'restrictleft1': restrictleft1,
        'restrictleft2': restrictleft2,
        'seq_start': seq_start,
        'umi_start': umi_start
    }
    """
    return config


