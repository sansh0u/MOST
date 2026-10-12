import yaml
import subprocess
from pathlib import Path
from importlib.resources import files, as_file
from most.preprocess.scan_bc import scan, scan_len
from most.model import Config


def unzip(file_path, thread):
    file_path = str(file_path)
    if not file_path.endswith(".gz"):
        return file_path
    out = file_path[:-3]
    with open(out, "wb") as f:
        subprocess.run(
            ["pigz", "-p", str(thread), "-dc", file_path],
            stdout=f,
            check=True
        )
    return out


def load_yaml(cfg_path):
    with open(cfg_path, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f)
    if data is None:
        raise ValueError("Empty YAML file")
    return Config(**data)


def convert_range(r):
    r = r.strip("()")
    start, end = map(int, r.split("-"))
    if start < 1 or end < start:
        raise ValueError(f"Invalid coordinate range: {r}")
    return start - 1, end - start + 1


def parse_pair(s):
    s = s.strip("()")
    r1, r2 = s.split(",")
    bc2_start, bc_len1 = convert_range(r1)
    bc1_start, bc_len2 = convert_range(r2)
    if bc_len1 != bc_len2:
        raise ValueError("BC1 and BC2 length are different")
    return bc2_start, bc1_start, bc_len1


def dump_yaml(data, path):
    with open(path, "w", encoding="utf-8") as f:
        yaml.safe_dump(data, f, sort_keys=False)


def build_zumis_yaml(
    cfg, gtf_file, bc2_start, bc1_start, bc_len, umi_start, umi_len
):
    out = Path(cfg.out_dir)
    zpath = out / "filled.yaml"

    bc_str = (
        f"BC({bc2_start + 1}-{bc2_start + bc_len},"
        f"{bc1_start + 1}-{bc1_start + bc_len})"
    )

    base_definition = [bc_str]
    if umi_len > 0:
        umi_str = f"UMI({umi_start + 1}-{umi_start + umi_len})"
        base_definition.append(umi_str)

    with as_file(files("most.config") / "zUMIs.yaml") as template:
        with open(template, encoding="utf-8") as f:
            zcfg = yaml.safe_load(f)

    zcfg["project"] = cfg.project
    zcfg["sequence_files"]["file1"]["name"] = str(
        out / "filtered_R1.fastq.gz"
    )
    zcfg["sequence_files"]["file1"]["base_definition"] = ["cDNA(1-100)"]
    zcfg["sequence_files"]["file2"]["name"] = str(
        out / "filtered_R2.fastq.gz"
    )
    zcfg["sequence_files"]["file2"]["base_definition"] = base_definition
    zcfg["reference"]["STAR_index"] = cfg.reference.star_index
    zcfg["reference"]["GTF_file"] = gtf_file
    zcfg["out_dir"] = str(out)
    zcfg["num_threads"] = int(cfg.threads)
    zcfg["barcodes"]["barcode_file"] = cfg.reference.barcode_file
    zcfg["reference"]["additional_STAR_params"] = (
        "--outFilterScoreMinOverLread 0.1 "
        "--outFilterMatchNminOverLread 0.1 "
        "--alignIntronMin 20 "
        "--alignIntronMax 1000000"
    )

    dump_yaml(zcfg, zpath)
    return zpath


def config_cal(cfg):
    method = cfg.method
    out_dir = Path(cfg.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    threads = int(cfg.threads)

    cfg.reference.fa_file = unzip(cfg.reference.fa_file, threads)
    cfg.reference.gtf_file = unzip(cfg.reference.gtf_file, threads)

    umi = cfg.advanced.umi
    bc = cfg.advanced.bc
    umi_len_cfg = cfg.advanced.umi_len

    if umi_len_cfg < 0:
        raise ValueError("umi_len must be >= 0")

    # umi_len=0 means no UMI; an explicit UMI coordinate conflicts with it.
    if umi_len_cfg == 0 and umi:
        raise ValueError(
            "umi_len=0 disables UMI; do not specify advanced.umi."
        )

    if umi and not bc:
        raise ValueError(
            "advanced.umi was provided, but advanced.bc is missing."
        )

    if umi and bc:
        umi_start, umi_len = convert_range(umi)
        if umi_len_cfg == 0:
            raise ValueError("UMI coordinates conflict with umi_len=0.")
        if umi_len != umi_len_cfg:
            raise ValueError(
                f"UMI coordinate length ({umi_len}) does not match "
                f"advanced.umi_len ({umi_len_cfg})."
            )
        bc2_start, bc1_start, bc_len = parse_pair(bc)
        read_len = int(scan_len(cfg))

    elif bc and umi_len_cfg == 0:
        bc2_start, bc1_start, bc_len = parse_pair(bc)
        read_len = int(scan_len(cfg))
        umi_start, umi_len = 0, 0

    else:
        print("Auto detecting barcode and UMI positions...")
        result = scan(cfg, method)
        bc2_start = int(result["bc2"])
        bc1_start = int(result["bc1"])
        bc_len = int(result["bc_len"])
        read_len = int(result["read_len"])
        umi_start = int(result["umi_start"])
        umi_len = int(result["umi_len"])

    if umi_len == 0:
        umi_start = 0

    bc2_end = bc2_start + bc_len
    bc1_end = bc1_start + bc_len
    linker1 = cfg.advanced.linker1
    linker2 = cfg.advanced.linker2

    if read_len == 100:
        linker2 = ""

    k1 = len(linker1)
    k2 = len(linker2)

    restrictleft1 = bc1_end + k1 + umi_len
    restrictleft2 = bc2_end + k2 + umi_len

    if umi_len == 0 or umi_start < bc2_start:
        seq_start = bc1_end + k1 + 19
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

    runtime_path = out_dir / "runtime.yaml"
    dump_yaml(cfg.model_dump(mode="json"), runtime_path)

    zpath = None
    if method == "ZUMIS":
        zpath = build_zumis_yaml(
            cfg,
            cfg.reference.gtf_file,
            bc2_start,
            bc1_start,
            bc_len,
            umi_start,
            umi_len
        )

    return runtime_path, zpath