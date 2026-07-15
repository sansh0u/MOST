#!/usr/bin/env python3

import gzip
from Bio import SeqIO
from Bio.Seq import Seq
import subprocess
import logging


logger = logging.getLogger("toolkit")

def open_fastq(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    else:
        return open(path, "r")

def write_fastq_record(h, rid, seq, qual_ints):
    h.write("@{}\n".format(rid))
    h.write(str(seq) + "\n")
    h.write("+\n")
    h.write("".join(chr(q + 33) for q in qual_ints) + "\n")

def search_barcode85(read1_fastq, read2_fastq, CleanFq1,lambda_r1_out, lambda_r2_out, min_len=18, MOTIF="TTATTTTT"
    ):

    adapter_r2 = "AGATGTGTATAAGAGATAG"

    barcode85_reads = set()

    with open_fastq(read1_fastq) as r1, open_fastq(read2_fastq) as r2:

        for rec1, rec2 in zip(
                SeqIO.parse(r1, "fastq"),
                SeqIO.parse(r2, "fastq")):

            seq1 = str(rec1.seq)
            seq2 = str(rec2.seq)

            if MOTIF in seq1 or MOTIF in seq2:
                barcode85_reads.add(rec1.id)

    
    with open_fastq(CleanFq1) as r1, open_fastq(read2_fastq) as r2, \
         gzip.open(lambda_r1_out, "wt") as out1, \
         gzip.open(lambda_r2_out, "wt") as out2:

        for rec1, rec2 in zip(
                SeqIO.parse(r1, "fastq"),
                SeqIO.parse(r2, "fastq")):

            if rec1.id not in barcode85_reads:
                continue

            lambda_seq1 = rec1.seq
            lambda_qual1 = rec1.letter_annotations["phred_quality"]


            seq2 = str(rec2.seq)

            if adapter_r2 in seq2:

                trim_i = seq2.find(adapter_r2) + len(adapter_r2)

                lambda_seq2 = rec2.seq[trim_i:]
                lambda_qual2 = (
                    rec2.letter_annotations["phred_quality"][trim_i:]
                )

            else:

                lambda_seq2 = rec2.seq
                lambda_qual2 = (
                    rec2.letter_annotations["phred_quality"]
                )

            if len(lambda_seq1) < min_len:
                continue

            if len(lambda_seq2) < min_len:
                continue
            write_fastq_record(out1, rec1.id, lambda_seq1, lambda_qual1)
            write_fastq_record(out2, rec2.id, lambda_seq2, lambda_qual2)
            

def dmt_filter(cfg):
    
    CleanFq1 = cfg.out_dir + "/filtered_R1.fastq.gz"
    CleanFq2 = cfg.out_dir + "/filtered_R2.fastq.gz"
    read1_fastq = cfg.sequence_file.file1
    read2_fastq = cfg.sequence_file.file2
    threads = cfg.threads
    lambda_r1_out = cfg.out_dir + "/Lambda_R1.fastq.gz"
    lambda_r2_out = cfg.out_dir + "/Lambda_R2.fastq.gz"
    min_len=18
    linker1 = "GTGGTTGATGTTTTGTATTGGTGTATGATT"
    linker2 = "ATTTATGTGTTTGAGAGGTTAGAGTATTTG"
    #linker1 = cfg.advanced.linker1
    #linker2 = cfg.advanced.linker2
    out_dir = cfg.out_dir
   
    ######
    k1 = cfg.runtime.k1
    k2 = cfg.runtime.k2
    restrictleft1 = cfg.runtime.restrictleft1
    restrictleft2 = cfg.runtime.restrictleft2 
    #####

    
    """
    QC and adapt the primer to the fastq files.
    Args:
        config (dict): The configuration file.
    """
    
    cmd1 = ["cutadapt",  "-m", "18", "-a", "CTATCTCTTATA",
            "-a", "AGATGCGAGAAGCCAACGCTTG",
            "-j", threads,
            "-o", CleanFq1 ,
            read1_fastq
            ]
    
    cmd2 =  [
        "bbduk.sh",
        f"in={read1_fastq}", 
        f"in2={read2_fastq}",
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

    cmd3 =  [
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
    try: 
        subprocess.run(cmd1, check=True)
        subprocess.run(cmd2, check=True)
        subprocess.run(cmd3, check=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Error during DBiT-seq filtering: {e}")
        raise

    search_barcode85(read1_fastq, read2_fastq, CleanFq1, lambda_r1_out, 
                     lambda_r2_out, min_len, MOTIF="TTATTTTT")






