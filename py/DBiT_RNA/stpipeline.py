import logging
import subprocess
from yaml_load import get_config
import os

logger = logging.getLogger("toolkit")


def unzip(file_path):
    if not file_path.endswith(".gz"):
        return file_path

    out = file_path[:-3]

    subprocess.run(["pigz", "-dc", file_path], stdout=open(out, "wb"), check=True)

    return out

def stpipeline(cfg):
    """
    根据提供的配置对DBiT-seq数据进行chromap分析
    """
    stpipeline_id = get_config(cfg, 'Project', 'stpipeline')
    output_folder = f"{get_config(cfg, 'dir')}"
    temp_folder = f"{get_config(cfg, 'dir')}/temp"
    os.makedirs(temp_folder, exist_ok=True)
    output_file_R2 = f"{get_config(cfg, 'dir')}/linker2_R1.fastq.gz"
    output_file_R1 = f"{get_config(cfg, 'dir')}/output_R2.fastq.gz"
    b_file = f"{get_config(cfg, 'dir')}/linker2_R2.fastq.gz"
    star_index = get_config(cfg, 'star_index')
    gtf_file = get_config(cfg, 'gtf_file')
    bc_file = get_config(cfg, 'barcode_file')
    thread = str(get_config(cfg, 'Threads'))
    
    out_gtf = unzip(gtf_file)

    cmd = [ "st_pipeline_run", 
        "--output-folder", output_folder,
        "--temp-folder", temp_folder, 
        "--ids", bc_file,
        "--threads", thread,
        "--ref-map", star_index,
        "--ref-annotation", out_gtf,
        "--expName", stpipeline_id,
        "--log-file", f"{output_folder}/{stpipeline_id}_log.txt",
        "--htseq-no-ambiguous", "--demultiplexing-kmer", "5",
        "--umi-start-position", "16",
        "--umi-end-position", "26",
        "--demultiplexing-overhang", "0",
        "--min-length-qual-trimming", "18",
        "--no-clean-up", "--verbose",
        output_file_R2,
        output_file_R1
]
    print(cmd)
    try: ####
        subprocess.run(cmd)#, check=True)
        #logger.info("ATAC-seq filtering completed successfully.")
        #返回点东西让我知道成功了
    except subprocess.CalledProcessError as e:
        #logger.error(f"Error during chromap analysis: {e}")
        raise
    subprocess.run(["rm", "-r",output_file_R1, output_file_R2, b_file, out_gtf], check=True)
