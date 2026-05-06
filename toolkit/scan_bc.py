import gzip
import numpy as np
from yaml_load import get_config
# -------------------------------
# config
# -------------------------------

MAX_READS = 100000
base2int = {'A':0, 'C':1, 'G':2, 'T':3}
# -------------------------------
# utils
# -------------------------------

def encode_kmer(seq):
    code = 0
    for c in seq:
        code = code * 4 + base2int.get(c, 0)
    return code


def load_barcodes(path):
    bc1_set = set()

    with open(path) as f:
        for line in f:

            parts = line.strip().split()
            seq = parts[0].upper()

            if len(seq) != 16:
                continue

            bc1_set.add(encode_kmer(seq[8:]))

    return bc1_set


def read_fastq_head(path, n):
    """only read first n reads"""
    seqs = []
    with gzip.open(path, "rt") as f:
        for i, line in enumerate(f):
            if i % 4 == 1:
                seqs.append(line.strip())
                if len(seqs) >= n:
                    break
    return seqs


# -------------------------------
# core (hash + rolling)
# -------------------------------
def scan_positions_hash(seqs, bc1_set):
    L = len(seqs[0])
    bc1_hits = np.zeros(L)
    mask = (4**7) - 1

    for seq in seqs:
        if len(seq) < 8:
            continue

        code = encode_kmer(seq[:8])

        for i in range(L - 7):

            if code in bc1_set:
                bc1_hits[i] += 1

            if i < L - 8:
                next_base = base2int.get(seq[i+8], 0)
                code = ((code & mask) * 4) + next_base

    return bc1_hits

def scan(config):
    """
    扫描 barcode
    """
    
    FASTQ = get_config(config, "file2")
    BARCODE_FILE = get_config(config, "Barcode")
    bc_set= load_barcodes(BARCODE_FILE)
    print(f"{len(bc_set)} bc loaded")
    seqs = read_fastq_head(FASTQ, MAX_READS)
    print(f"reads loaded: {len(seqs)}")
    read_len = len(seqs[0])
    bc_pos = scan_positions_hash(seqs, bc_set)
    total = len(seqs)
    ratio = bc_pos / total   

# 找最大的两个点
    idx = np.argsort(ratio)[-2:]

# 按位置排序（小→大）
    idx.sort()

    bc2_loc = idx[0] 
    bc1_loc = idx[1] 

    print("\n=== Peak result ===")
    print(f"bc2_location\t{bc2_loc}\t{ratio[bc2_loc]:.4f}")
    print(f"bc1_location\t{bc1_loc}\t{ratio[bc1_loc]:.4f}")
    print(f"read_length\t{read_len}")
    return bc2_loc, bc1_loc, read_len
