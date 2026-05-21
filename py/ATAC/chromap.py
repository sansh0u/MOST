import logging
import subprocess

logger = logging.getLogger("toolkit")

def chromap(cfg):
    """
    根据提供的配置对ATAC-seq数据进行chromap分析
    """
    output_file_R1 = f"{cfg.out_dir}/output_R1.fastq.gz"
    output_file_R2 = f"{cfg.out_dir}/linker2_R1.fastq.gz"
    subprocess.run(["rm", f"{cfg.out_dir}/linker2_R2.fastq.gz"], check=True)
    b_file = f"{cfg.out_dir}/output_R2.fastq.gz"
    index_file = cfg.reference.index_file
    fa_file = cfg.reference.fa_file
    output_file = f"{cfg.out_dir}/{cfg.project}.bed"
    bc_file = cfg.reference.barcode_file
    thread = str(cfg.threads)
    cmd = [ "chromap", 
        "--preset", "atac", "-x", index_file, 
        "-r", fa_file, 
        "-1", output_file_R1, 
        "-2", output_file_R2,
        "-b", b_file,
        "--barcode-whitelist", bc_file,
        "-t", thread, "-o", output_file
]
    try: ####
        subprocess.run(cmd, check=True)
        #logger.info("ATAC-seq filtering completed successfully.")
        #返回点东西让我知道成功了
    except subprocess.CalledProcessError as e:
        #logger.error(f"Error during chromap analysis: {e}")
        raise
    subprocess.run(["rm", "-r",output_file_R1, output_file_R2, b_file], check=True)


def sort_bed(cfg):
    """
    对chromap输出的bed文件进行排序
    """
    
    #subprocess.run(["sort", "-k1,1", "-k2,2n", "-k3,3n", "-k4,4", f"--parallel={cfg.threads}", "-S 36G", output_file +".bed", ">", output_file + "_sorted.bed"], check=True)
    threads = str(cfg.threads)
    output_file = f"{cfg.out_dir}/{cfg.project}"
    with open(output_file + "_sorted.bed", "w") as f:
        subprocess.run([
            "sort",
            "-k1,1",
            "-k2,2n",
            "-k3,3n",
            "-k4,4",
            f"--parallel={threads}",
            "-S", "36G",
            output_file + ".bed"
        ], stdout=f, check=True)

    subprocess.run([
        "bgzip",
        "-@",
        threads,
        output_file + "_sorted.bed"
    ], check=True)

    subprocess.run([
        "tabix",
        "-p", "bed",
        output_file + "_sorted.bed.gz"
    ], check=True)
    #subprocess.run(["bgzip", output_file + "_sorted.bed" ,"-@", cfg.threads], check=True)
    #subprocess.run(["tabix", "-p", "bed", output_file + "_sorted.bed.gz"], check=True)

#sort -k1,1 -k2,2n -k3,3n -k4,4 --parallel=12 -S 36G output_file.bed" > output_file_sorted.bed
#module load tabix
#bgzip output_file_sorted.bed
#tabix -p bed output_file_sorted.bed.gz
