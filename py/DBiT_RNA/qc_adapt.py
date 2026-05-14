import subprocess
import logging
from config_utils import get_config
from preprocess.scan_bc import check_adapter

logger = logging.getLogger("toolkit")

def qc_adapt(cfg):
    Adapter = get_config(cfg, "adapter","AAGCAGTGGTATCAACGCAGAGTGAATGGG")
    CleanFq1 = get_config(cfg, "file1")
    CleanFq2 = get_config(cfg, "file2")
    fastq_intput_1 = get_config(cfg, "dir")+ "/output_R1.fastq"
    fastq_intput_2 = get_config(cfg, "dir")+ "/output_R2.fastq"
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
            subprocess.run(cmd1, check=True)
            print("Trimming Adapter")
        else:
            subprocess.run(cmd2, check=True)
            print("Trimming ployA")
    except subprocess.CalledProcessError as e:
            logger.error(f"Error during DBiT-seq filtering: {e}")
            raise

