import subprocess
import logging
from yaml_load import get_config

logger = logging.getLogger("toolkit")

def qc_adapt(config,result):
    Adapter = get_config(config, "adapter","AAGCAGTGGTATCAACGCAGAGTGAATGGG")
    CleanFq1 = get_config(config, "file1")
    CleanFq2 = get_config(config, "file2")
    fastq_intput_1 = get_config(config, "dir")+ "/output_R1.fastq"
    fastq_intput_2 = get_config(config, "dir")+ "/output_R2.fastq"
    
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
        if rna == "illumina":
            subprocess.run(cmd1, check=True)
        else:
            subprocess.run(cmd2, check=True)
    except subprocess.CalledProcessError as e:
            logger.error(f"Error during DBiT-seq filtering: {e}")
            raise