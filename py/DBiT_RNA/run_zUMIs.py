import subprocess
import logging
import os
import typer
import yaml
from yaml_load import load_yaml, get_config
logger = logging.getLogger("toolkit")

def zUMIs(zpath, zcfg_path):
    '''
    调用zUMIs,要把in1,in2,out写进去
    '''
    
    cmd = [
        zpath,
        "-c", #-C是运行zumis自己的环境要重新下载，如果可以设置好conda
        "-y",
        zcfg_path
    ]
    try: ####
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Error during zUMIs: {e}")
        raise
