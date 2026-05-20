from Bio.SeqIO.QualityIO import FastqGeneralIterator
from gzip import open as gzopen
import logging
import subprocess
from yaml_load import get_config
#import argparse

logger = logging.getLogger("toolkit")


def atac_bc(cfg):
    '''
    处理BCB_UMI格式的fastq文件,提取UMI和barcode,输出到新的fastq文件
    '''
    
    
    input_file = get_config(cfg, "Out_dir") + "/linker2_R2.fastq.gz"
    output_file_R1 = get_config(cfg, "Out_dir") + "/output_R1.fastq"
    output_file_R2 = get_config(cfg, "Out_dir") + "/output_R2.fastq"

    
    seq_start = get_config(cfg, "seq_start")
    bc2_start = get_config(cfg, "bc2_start")    
    bc2_end = get_config(cfg, "bc2_end")
    bc1_start = get_config(cfg, "bc1_start")
    bc1_end = get_config(cfg, "bc1_end")
    thread = get_config(cfg, "Thread")
    """
    logger.info(f"input_file: {input_file}")
    logger.info(f"output_file_R1: {output_file_R1}")
    logger.info(f"output_file_R2: {output_file_R2}")
    logger.info(f"seq_start: {seq_start}")
    logger.info(f"bc2_start: {bc2_start}")
    logger.info(f"bc2_end: {bc2_end}")
    logger.info(f"bc1_start: {bc1_start}")
    logger.info(f"bc1_end: {bc1_end}")
    """
    with gzopen(input_file, "rt") as in_handle_R1, open(output_file_R1, "w") as out_handle_R1, open(output_file_R2, "w") as out_handle_R2:
        #logger.info("Start BC processing")
        for title, seq, qual in FastqGeneralIterator(in_handle_R1):
            new_seq_R1 = seq[seq_start:]
            new_qual_R1 = qual[seq_start:]
            barcode = seq[bc2_start:bc2_end] + seq[bc1_start:bc1_end] # !!! BC2 + BC1
            new_qual_R2 = qual[bc2_start:bc2_end] + qual[bc1_start:bc1_end]        
            out_handle_R1.write("@%s\n%s\n+\n%s\n" % (title, new_seq_R1, new_qual_R1))
            out_handle_R2.write("@%s\n%s\n+\n%s\n" % (title, barcode, new_qual_R2))
    
    subprocess.run(["pigz", "-p", thread, "-f", output_file_R1], check=True)
    subprocess.run(["pigz", "-p", thread, "-f", output_file_R2], check=True)

def dbit_bc(cfg):
    """
    BC2,BC1,UMI
    """

    input_file = get_config(cfg, "Out_dir") + "/linker2_R2.fastq.gz"
    output_file = get_config(cfg, "Out_dir") + "/output_R2.fastq"
    
    umi_start = get_config(cfg, "umi_start")
    umi_len = get_config(cfg, "umi_len")
    bc2_start = get_config(cfg, "bc2_start")
    bc2_end = get_config(cfg, "bc2_end")
    bc1_start = get_config(cfg, "bc1_start")
    bc1_end = get_config(cfg, "bc1_end")
    thread = get_config(cfg, "Thread")


    with gzopen(input_file, "rt") as in_handle:
        with open(output_file, "w") as out_handle:
            for title, seq, qual in FastqGeneralIterator(in_handle):
                new_seq = seq[bc2_start:bc2_end] + seq[bc1_start:bc1_end] + seq[umi_start:umi_start+umi_len]  # BC2 + BC1 + UMI
                new_qual = qual[bc2_start:bc2_end] + qual[bc1_start:bc1_end] + qual[umi_start:umi_start+umi_len]
                out_handle.write("@%s\n%s\n+\n%s\n" % (title, new_seq, new_qual))
    
    subprocess.run(["pigz", "-p", thread, "-f", output_file], check=True)