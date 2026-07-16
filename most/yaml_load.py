import yaml
import os
import subprocess
from pathlib import Path
from most.preprocess.scan_bc import scan,scan_len
from most.model import Config
from importlib.resources import files,as_file

def unzip(file_path,thread):
    if not file_path.endswith(".gz"):
        return file_path
    
    out = file_path[:-3]
    with open(out, "wb") as f:
        subprocess.run(
            ["pigz","-p",thread,"-dc", file_path], stdout=f, check=True
        )

    return out

def load_yaml(cfg_path):
    """
    Load YAML configuration file. 加载并判断YAML配置文件是否存在并可解析(已完成)
    """
    try:

        with open(cfg_path, "r", encoding="utf-8") as f:
            data = yaml.safe_load(f)
        if data is None:
            raise ValueError("YAML file is empty")
        cfg = Config(**data)
        return cfg

    except FileNotFoundError:
       print(f"Configuration file not found: {cfg_path}")
       raise

    except yaml.YAMLError as e:
        print(f"YAML parse error: {e}")
        raise

    except Exception as e:
        print(f"Failed to load config: {e}")
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
    threads = cfg.threads
    out_dir = cfg.out_dir
    fa_file = cfg.reference.fa_file
    fa_file = unzip(fa_file,threads)
    cfg.reference.fa_file = fa_file
    gtf_file = cfg.reference.gtf_file
    gtf_file = unzip(gtf_file,threads)
    cfg.reference.gtf_file = gtf_file
    if umi is not None and bc is not None:
        umi_start, umi_len = convert_range(umi)
        bc2_start, bc1_start, bc_len = parse_pair(bc)
        read_len = scan_len(cfg)
    else:
        print("No UMI or BC found in config. Enter auto mode.")
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
    cfg.runtime.restrictleft1 = restrictleft1
    cfg.runtime.restrictleft2 = restrictleft2
    cfg.runtime.seq_start = seq_start
    cfg.runtime.umi_start = umi_start
    cfg.runtime.umi_len = umi_len
    
    path = out_dir / "runtime.yaml"
    zpath = None 
    if method == "ZUMIS":
        zpath = out_dir / "filled.yaml"
        bc_str = (f"BC("f"{bc2_start+1}-{bc2_end},"f"{bc1_start+1}-{bc1_end}"f")")
        umi_str = (f"UMI("f"{umi_start+1}-"f"{umi_start+umi_len}"f")")
        with as_file(files("most.config") / "RNA.yaml") as zcfg_path:
            with open(zcfg_path) as f:
                zcfg = yaml.safe_load(f)
        zcfg["project"] = cfg.project
        zcfg["sequence_files"]["file1"]["name"] = str(cfg.out_dir / "filtered_R1.fastq.gz")
        zcfg["sequence_files"]["file1"]["base_definition"] = [
            "cDNA(1-100)"
        ]
        zcfg["sequence_files"]["file2"]["name"] = str(cfg.out_dir / "filtered_R2.fastq.gz")
        zcfg["sequence_files"]["file2"]["base_definition"] = [
            bc_str,
            umi_str
        ]
        zcfg["reference"]["STAR_index"] = cfg.reference.star_index
        zcfg["reference"]["GTF_file"] = gtf_file
        zcfg["out_dir"] = cfg.out_dir
        zcfg['num_threads'] = threads
        zcfg["barcodes"]["barcode_file"] =  cfg.reference.barcode_file
        zcfg["reference"]["additional_STAR_params"]="--outFilterScoreMinOverLread 0.1 --outFilterMatchNminOverLread 0.1 --alignIntronMin 20 --alignIntronMax 1000000"
        with open(zpath, "w") as f:
            yaml.safe_dump(
                zcfg,
                f,
                default_flow_style=False,
                sort_keys=False
            )
    with open(path, "w", encoding="utf-8") as f:
        yaml.safe_dump(
            cfg.model_dump(),
            f,
            sort_keys=False
        )
    return path,zpath