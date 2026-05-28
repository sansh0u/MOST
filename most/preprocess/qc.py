import subprocess
import logging

logger = logging.getLogger("toolkit")


def filter(cfg,method):
    """
    根据提供的配置对ATAC-seq数据进行过滤
    """
    print("Filtering started")
    # Placeholder for actual filtering logic需要校对logo信息
    #logger.info("Starting ATAC-seq quality control filtering...")
    if method == 'RNA':
        in1 = cfg.out_dir+ "/filtered_R1.fastq.gz"
        in2 = cfg.out_dir+ "/filtered_R2.fastq.gz"
    else:
        in1 = cfg.sequence_file.file1
        in2 = cfg.sequence_file.file2
    out_dir = cfg.out_dir
   
    k1 = cfg.runtime.k1
    k2 = cfg.runtime.k2
    threads = cfg.threads
    linker1 = cfg.advanced.linker1
    linker2 = cfg.advanced.linker2
    restrictleft1 = cfg.runtime.restrictleft1
    restrictleft2 = cfg.runtime.restrictleft2
    
    cmd1 =  [
        "bbduk.sh",
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
        "bbduk.sh",
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
        "bbduk.sh",
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
            #subprocess.run(["rm", "-r",f"{out_dir}/linker1_R1.fastq.gz", f"{out_dir}/linker1_R2.fastq.gz"], check=True)
        else:
            subprocess.run(cmd3, check=True)
            
        #logger.info("ATAC-seq filtering completed successfully.")
        #返回点东西让我知道成功了
    except subprocess.CalledProcessError as e:
        logger.error(f"Error during ATAC-seq filtering: {e}")
        raise

    #logger.info("ATAC-seq quality control filtering completed.")
    
