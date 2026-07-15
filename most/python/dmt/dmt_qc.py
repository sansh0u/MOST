#!/usr/bin/env python3
import gzip
import os
from fuzzysearch import find_near_matches
from Bio import SeqIO


def open_fastq_file(fp):
    return gzip.open(fp, "rt") if fp.endswith(".gz") else open(fp, "r")

def write_fastq_record(h, rid, seq, qual_ints):
    h.write(f"@{rid}\n")
    h.write(str(seq) + "\n")
    h.write("+\n")
    h.write("".join(chr(q + 33) for q in qual_ints) + "\n")


def load_whitelist(path, col=4):
    if not path:
        return None
    wl = set()
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= col:
                wl.add(parts[col-1])
    return wl


def assign_barcode_hamming1(bc, wl):
    if bc in wl:
        return bc
    bases = ("A","C","G","T")
    for i, x in enumerate(bc):
        for b in bases:
            if b == x:
                continue
            cand = bc[:i]+b+bc[i+1:]
            if cand in wl:
                return cand
    return None



def trim_adapters(cfg):

    r1_fastq = cfg.sequence_file.file1
    r2_fastq = cfg.sequence_file.file2
    out_prefix = cfg.out_dir

    adapters = ["CTATCTCTTATA","AGATGCGAGAAGCCAACGCTTG"]
    linker1 = "GTGGTTGATGTTTTGTATTGGTGTATGATT"
    linker2 = "ATTTATGTGTTTGAGAGGTTAGAGTATTTG"

    trim_end2 = "AGATGTGTATAAGAGATAG"
    min_len=10

    os.makedirs(os.path.dirname(out_prefix), exist_ok=True)

    r1_out = gzip.open(out_prefix+"_R1.fq.gz", "wt")
    r2_out = gzip.open(out_prefix+"_R2.fq.gz", "wt")

    BARCODE85_MOTIF="TTATTTTT"

    n_reads=0
    n_pass=0
    n_trim=0

    with open_fastq_file(r1_fastq) as h1,  open_fastq_file(r2_fastq) as h2:

        for rec1,rec2 in zip( SeqIO.parse(h1,"fastq"), SeqIO.parse(h2,"fastq")):

            n_reads+=1

            r1seq=rec1.seq
            r2seq=rec2.seq

            qual1=rec1.letter_annotations["phred_quality"]
            qual2=rec2.letter_annotations["phred_quality"]

            # =====================
            # trim R1 adapter
            # =====================

            hits=[]

            for a in adapters:

                hits.extend(find_near_matches(a, str(r1seq), max_l_dist=1))

            if hits:

                trim_i=min(x.start for x in hits)
                r1seq=r1seq[:trim_i]
                qual1=qual1[:trim_i]
                n_trim+=1

            # =====================
            # trim R2 linker
            # =====================

            te2=find_near_matches(trim_end2, str(r2seq), max_l_dist=1)

            if te2:
                idx=te2[0].end
                r2seq=r2seq[idx:]
                qual2=qual2[idx:]

            # =====================
            # linker barcode
            # =====================

            lk1=find_near_matches(linker1, str(rec2.seq), max_l_dist=2)

            lk2=find_near_matches(linker2, str(rec2.seq), max_l_dist=2)

            if not (len(lk1)==1 and len(lk2)==1 and len(te2)==1):
                continue

            b1=lk1[0].start-8
            b2=lk2[0].start-8

            barcode=str(
                rec2.seq[b1:lk1[0].start]
                +
                rec2.seq[b2:lk2[0].start]
            )

            if len(barcode)!=16:
                continue

            if len(r1seq)<min_len:
                continue

            if len(r2seq)<min_len:
                continue

            rid1=f"{barcode}_{rec1.id}"
            rid2=f"{barcode}_{rec2.id}"

            write_fastq_record(r1_out, rid1, r1seq, qual1)
            write_fastq_record(r2_out,rid2,r2seq,qual2)
            n_pass+=1

    r1_out.close()
    r2_out.close()

    with open(out_prefix+"_stats.txt","w") as f:
        f.write(f"Reads_total\t{n_reads}\n")
        f.write(f"Reads_passed\t{n_pass}\n")
        f.write(f"Reads_trimmed\t{n_trim}\n")
    print(f"trim adapters finished: {n_pass}/{n_reads}")