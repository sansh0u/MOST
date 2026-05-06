import subprocess
import logging
from yaml_load import get_config

logger = logging.getLogger("toolkit")


def filter(config):
    """
    根据提供的配置对ATAC-seq数据进行过滤
    """
     # bbduk写入bin了，直接使用
    

    # Placeholder for actual filtering logic需要校对logo信息
    #logger.info("Starting ATAC-seq quality control filtering...")
    out_dir = get_config(config, "dir")
    in1 = get_config(config, "file1")
    in2 = get_config(config, "file2")
    k1 = get_config(config, "k1")
    k2 = get_config(config, "k2")
    threads = get_config(config, "Threads")
    linker1 = get_config(config, "linker1")
    linker2 = get_config(config, "linker2")
    restrictleft1 = get_config(config, "restrictleft1")
    restrictleft2 = get_config(config, "restrictleft2")

    cmd1 =  [
        "bbduk",
        f"in={in1}", 
        f"in2={in2}",
        f"outm={out_dir}/linker1_R1.fastq.gz", ####改名字
        f"outm2={out_dir}/linker1_R2.fastq.gz", ####
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
        f"outm={out_dir}/linker2_R1.fastq.gz", ####
        f"outm2={out_dir}/linker2_R2.fastq.gz", ####
        f"hdist=3",
        f"k={k2}",
        f"literal={linker2}",
        f"threads={config['Threads']}",
        "mm=f", "rcomp=f", f"skipr1=t",
        f"restrictleft={restrictleft2}",
        f"stats={out_dir}/bbduk_stats_L2.txt"
    ]
    
    cmd3 = [
        "bbduk",
        f"in={in1}", 
        f"in2={in2}",
        f"outm={out_dir}/linker2_R1.fastq.gz", ####改名字
        f"outm2={out_dir}/linker2_R2.fastq.gz", ####
        f"hdist=3",
        f"k={k2}",
        f"literal={linker2}",
        f"threads={threads}",
        "mm=f", "rcomp=f", f"skipr1=t",
        f"restrictleft={restrictleft2}",
        f"stats={out_dir}/bbduk_stats_L1.txt"
    ]
    try: ####
        if k1 != 0:
            subprocess.run(cmd1, check=True)
            subprocess.run(cmd2, check=True)
            subprocess.run(["rm", "-r",f"{out_dir}/linker1_R1.fastq.gz", f"{out_dir}/linker1_R2.fastq.gz"], check=True)
        else:
            subprocess.run(cmd3, check=True)
            print(cmd3)
        #logger.info("ATAC-seq filtering completed successfully.")
        #返回点东西让我知道成功了
    except subprocess.CalledProcessError as e:
        logger.error(f"Error during ATAC-seq filtering: {e}")
        raise

    #logger.info("ATAC-seq quality control filtering completed.")
    
