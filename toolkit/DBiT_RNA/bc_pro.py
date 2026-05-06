from Bio.SeqIO.QualityIO import FastqGeneralIterator
from gzip import open as gzopen
import logging
import subprocess
from yaml_load import get_config
#import argparse

logger = logging.getLogger("toolkit")

def bc_pro(config):
    """
    BC2,BC1,UMI
    """
#ap = argparse.ArgumentParser()
#ap.add_argument("-i", "--input", required=True, help="Path to the inputfile")
#ap.add_argument("-o", "--output", required=True, help="Path to the outputfile")
#args = vars(ap.parse_args())

    input_file = get_config(config, "dir") + "/linker2_R2.fastq.gz"
    output_file = get_config(config, "dir") + "/output_R2.fastq"
    umi_start = get_config(config, "umi_start")
    bc2_start = get_config(config, "bc2_start")
    bc2_end = get_config(config, "bc2_end")
    bc1_start = get_config(config, "bc1_start")
    bc1_end = get_config(config, "bc1_end")


    with gzopen(input_file, "rt") as in_handle:
        with open(output_file, "w") as out_handle:
            for title, seq, qual in FastqGeneralIterator(in_handle):
                new_seq = seq[bc2_start:bc2_end] + seq[bc1_start:bc1_end] + seq[umi_start:umi_start+10]  # BC2 + BC1 + UMI
                new_qual = qual[bc2_start:bc2_end] + qual[bc1_start:bc1_end] + qual[umi_start:umi_start+10]
                out_handle.write("@%s\n%s\n+\n%s\n" % (title, new_seq, new_qual))
    
    subprocess.run(["pigz", "-p", "12", "-f", output_file], check=True)
