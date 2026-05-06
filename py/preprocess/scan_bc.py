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
    seqs = []
    with gzip.open(path, "rt") as f:
        for i, line in enumerate(f):
            if i % 4 == 1:
                seqs.append(line.strip().upper())
                if len(seqs) >= n:
                    break
    return seqs


# -------------------------------
# barcode scan (你的原始逻辑)
# -------------------------------
def scan_positions_hash(seqs, fastq):
    L = len(seqs[0])
    bc1_hits = np.zeros(L)
    mask = (4**7) - 1

    for seq in seqs:
        if len(seq) < 8:
            continue

        code = encode_kmer(seq[:8])

        for i in range(L - 7):

            if code in fastq:
                bc1_hits[i] += 1

            if i < L - 8:
                next_base = base2int.get(seq[i+8], 0)
                code = ((code & mask) * 4) + next_base

    return bc1_hits


# -------------------------------
# TSO detection
# -------------------------------
def hamming(s1, s2):
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def scan_adapter_positions(seqs, query, max_mismatch, window):
    query = query.upper()
    k = len(query)

    L = len(seqs[0])
    hits = np.zeros(L)

    for seq in seqs:
        if len(seq) < k:
            continue

        # 只扫前 window
        max_i = min(window, len(seq) - k + 1)

        for i in range(max_i):
            if hamming(seq[i:i+k], query) <= max_mismatch:
                hits[i] += 1

    return hits





def scan(config, method):


    FASTQ = get_config(config, "file2")
    BARCODE_FILE = get_config(config, "Barcode")
    seqs = read_fastq_head(FASTQ, MAX_READS)
    print(f"reads loaded: {len(seqs)}")

    # ---------------- barcode ----------------
    bc_set = load_barcodes(BARCODE_FILE)
    print(f"{len(bc_set)} bc loaded")
    read_len = len(seqs[0])
    bc_hits = scan_positions_hash(seqs, bc_set)
    ratio = bc_hits / len(seqs)

    idx = np.argsort(ratio)[-2:]
    idx.sort()

    bc2_loc = idx[0]
    bc1_loc = idx[1]

    print("\n=== Barcode result ===")
    print(f"bc2_location\t{bc2_loc}\t{ratio[bc2_loc]:.4f}")
    print(f"bc1_location\t{bc1_loc}\t{ratio[bc1_loc]:.4f}")
    print(f"read_length\t{read_len}")
    result = {
        "bc2": bc2_loc,
        "bc1": bc1_loc,
        "read_len": read_len
    }

    # ---------------- adapter（可选） ----------------
    if method == "RNA":

        Adapter = get_config(config, "adapter","AAGCAGTGGTATCAACGCAGAGTGAATGGG")
        mismatch = get_config(config, "adapter_mismatch", 2)

        if Adapter:
            hits = scan_adapter_positions(seqs, Adapter, mismatch, 50)
            ratio = hits / len(seqs)

            sub = ratio[:50]
            pos = sub.argmax()
            score = sub.max()

            print("\n=== Adapter (TSO) result ===")
            print(f"position\t{pos}")
            print(f"ratio\t{score:.4f}")

            result["Adapter"] = {
                "pos": pos,
                "score": score
            }
        else:
            print("\n=== Adapter skipped (no Adapter in config) ===")

    else:
        print("\n=== Adapter skipped (ATAC mode) ===")

    return result
