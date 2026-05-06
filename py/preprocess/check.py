import gzip
import numpy as np

MAX_READS = 100000

def read_fastq_head(path, n):
    seqs = []
    with gzip.open(path, "rt") as f:
        for i, line in enumerate(f):
            if i % 4 == 1:
                seqs.append(line.strip().upper())
                if len(seqs) >= n:
                    break
    return seqs


def scan_sequence_positions(seqs, query):
    """
    扫描 query 在每个位置的出现频率
    """
    query = query.upper()
    k = len(query)

    L = len(seqs[0])
    hits = np.zeros(L)

    for seq in seqs:
        if len(seq) < k:
            continue

        for i in range(L - k + 1):
            if seq[i:i+k] == query:
                hits[i] += 1

    return hits


def scan_query(fastq, query):
    seqs = read_fastq_head(fastq, MAX_READS)
    print(f"reads loaded: {len(seqs)}")

    hits = scan_sequence_positions(seqs, query)
    ratio = hits / len(seqs)

    # 找峰
    idx = np.argsort(ratio)[-5:]   # 取前5个位置
    idx.sort()

    print("\n=== Peak positions ===")
    for i in idx:
        print(f"pos {i+1}\t{ratio[i]:.4f}")

    return idx, ratio


# 用法
scan_query("/mnt/d/prog/RNA2_1209_Illu_1.fq.gz", "AAGCAGTGGTATCAACGCAGAGTGAATGGG")