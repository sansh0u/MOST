import subprocess
import logging
from yaml_load import get_config

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
            print("Trimming adapter")
        else:
            subprocess.run(cmd2, check=True)
            print("Trimming ployA")
    except subprocess.CalledProcessError as e:
            logger.error(f"Error during DBiT-seq filtering: {e}")
            raise

def check_adapter(cfg):################################################
    if method.upper() == "RNA":
            MAX_READS = 100000
            FASTQ = get_config(cfg, "file1")
            seqs1 = read_fastq_head(FASTQ, MAX_READS)
            print(f"reads loaded: {len(seqs1)}")
            Adapter = get_config(cfg, "adapter","AAGCAGTGGTATCAACGCAGAGTGAATGGG")
            mismatch = get_config(cfg, "adapter_mismatch", 2)

            if Adapter:
                hits = scan_adapter_positions(seqs1, Adapter, mismatch, 50)
                ratio = hits / len(seqs1)

                sub = ratio[:50]
                pos = sub.argmax()
                score = sub.max()

                print("\n=== Adapter result ===")
                print(f"position\t{pos}")
                print(f"ratio\t{score:.4f}")

                result = {
                    "score": score
                }
            else:
                print("\n=== Adapter skipped (no Adapter in config) ===")

    else:
            print("\n=== Adapter skipped (ATAC mode) ===")

    return result