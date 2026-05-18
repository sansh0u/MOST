import subprocess
import logging
from yaml_load import get_config

logger = logging.getLogger("toolkit")


def filter(cfg,method):
    """
    根据提供的配置对ATAC-seq数据进行过滤
    """
    # Placeholder for actual filtering logic需要校对logo信息
    #logger.info("Starting ATAC-seq quality control filtering...")
    if method == 'RNA':
        in1 = get_config(cfg, "Out_dir")+ "/filtered_R1.fastq"
        in2 = get_config(cfg, "Out_dir")+ "/filtered_R2.fastq"
    else:
        in1 = get_config(cfg, "file1")
        in2 = get_config(cfg, "file2")
    out_dir = get_config(cfg, "Out_dir")
   
    k1 = get_config(cfg, "k1")
    k2 = get_config(cfg, "k2")
    threads = get_config(cfg, "Threads")
    linker1 = get_config(cfg, "linker1")
    linker2 = get_config(cfg, "linker2")
    restrictleft1 = get_config(cfg, "restrictleft1")
    restrictleft2 = get_config(cfg, "restrictleft2")
    
    cmd1 =  [
        "bbduk",
        f"in={in1}", 
        f"in2={in2}",
        f"outm={out_dir}/linker1_R1.fastq.gz", 
        f"outm2={out_dir}/linker1_R2.fastq.gz", 
        f"hdist=3",
        f"k={k1}",
        f"literal={linker1}",
        f"threads={threads}",
        "mm=f", "rcomp=f", f"skipr1=t",
        f"restrictleft={restrictleft1}",
        f"stats={out_dir}/bbduk_stats_L1.txt"
    ]

    cmd2 =  [
        "bbduk",
        f"in={out_dir}/linker1_R1.fastq.gz",
        f"in2={out_dir}/linker1_R2.fastq.gz",
        f"outm={out_dir}/linker2_R1.fastq.gz", 
        f"outm2={out_dir}/linker2_R2.fastq.gz", 
        f"hdist=3",
        f"k={k2}",
        f"literal={linker2}",
        f"threads={threads}",
        "mm=f", "rcomp=f", f"skipr1=t",
        f"restrictleft={restrictleft2}",
        f"stats={out_dir}/bbduk_stats_L2.txt"
    ]
    
    cmd3 = [
        "bbduk",
        f"in={in1}", 
        f"in2={in2}",
        f"outm={out_dir}/linker2_R1.fastq.gz",
        f"outm2={out_dir}/linker2_R2.fastq.gz", 
        f"hdist=3",
        f"k={k1}",
        f"literal={linker2}",
        f"threads={threads}",
        "mm=f", "rcomp=f", f"skipr1=t",
        f"restrictleft={restrictleft2}",
        f"stats={out_dir}/bbduk_stats_L1.txt"
    ]
    try: ####
        if k2 != 0:
            subprocess.run(cmd1, check=True)
            subprocess.run(cmd2, check=True)
            subprocess.run(["rm", "-r",f"{out_dir}/linker1_R1.fastq.gz", f"{out_dir}/linker1_R2.fastq.gz"], check=True)
        else:
            subprocess.run(cmd3, check=True)
            
        #logger.info("ATAC-seq filtering completed successfully.")
        #返回点东西让我知道成功了
    except subprocess.CalledProcessError as e:
        logger.error(f"Error during ATAC-seq filtering: {e}")
        raise

    #logger.info("ATAC-seq quality control filtering completed.")
    
