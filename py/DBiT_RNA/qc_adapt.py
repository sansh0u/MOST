import subprocess
import logging
from preprocess.scan_bc import check_adapter

logger = logging.getLogger("toolkit")

def qc_adapt(cfg):
    Adapter = cfg.advanced.adapter
    CleanFq1 = cfg.out_dir + "/filtered_R1.fastq.gz"
    CleanFq2 = cfg.out_dir + "/filtered_R2.fastq.gz"
    fastq_intput_1 = cfg.sequence_file.file1
    fastq_intput_2 = cfg.sequence_file.file2
    score = check_adapter(cfg)
    """
    QC and adapt the primer to the fastq files.
    Args:
        config (dict): The configuration file.
    """

    cmd1 = [
    "cutadapt", "-m", "18", "-a", "A{10}N{150}",
    "--times", "4",
    "-g", Adapter,
    "-j", "12",
    "-o", CleanFq1,
    "-p", CleanFq2,
    fastq_intput_1, fastq_intput_2
]

    cmd2 = [
    "cutadapt", "-m", "18", "-a", "A{10}N{150}",
    "--times", "4",
    "-j", "12",
    "-o", CleanFq1,
    "-p", CleanFq2,
    fastq_intput_1, fastq_intput_2
]
    try: ####
        if score >= 0.5:
            print("Trimming Adapter")
            subprocess.run(cmd1, check=True)
            
        else:
            print("Trimming ployA")
            subprocess.run(cmd2, check=True)
            
    except subprocess.CalledProcessError as e:
            logger.error(f"Error during DBiT-seq filtering: {e}")
            raise

