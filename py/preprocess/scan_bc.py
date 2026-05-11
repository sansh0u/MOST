import gzip
import numpy as np
from yaml_load import get_config
from collections import Counter
import math


base2int = {'A':0, 'C':1, 'G':2, 'T':3}

# -------------------------------
# utils
# -------------------------------
def encode_kmer(seq):
    code = 0
    for c in seq:
        if c not in base2int:
            return None
        code = code * 4 + base2int[c]
    return code

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
# barcode
# -------------------------------
def load_barcodes(path):
    bc_set = set()
    bc_len = None

    with open(path) as f:
        for line in f:
            parts = line.strip().split()
            if not parts:
                continue

            seq = parts[0].upper()

            if bc_len is None:
                bc_len = len(seq) // 2

            if len(seq) != bc_len * 2:
                continue

            bc = seq[bc_len:]
            code = encode_kmer(bc)

            if code is not None:
                bc_set.add(code)

    return bc_set, bc_len

def scan_positions_hash(seqs, bc_set, bc_len):
    L = len(seqs[0])
    bc_hits = np.zeros(L)
    mask = (4 ** (bc_len - 1)) - 1

    for seq in seqs:
        if len(seq) < bc_len:
            continue

        code = encode_kmer(seq[:bc_len])
        if code is None:
            continue

        for i in range(L - bc_len + 1):
            if code in bc_set:
                bc_hits[i] += 1

            if i < L - bc_len:
                next_base = seq[i + bc_len]
                if next_base not in base2int:
                    break
                code = ((code & mask) * 4 + base2int[next_base])

    return bc_hits

def split_by_barcode(consensus, bc2_loc, bc1_loc, bc_len):

    left = (0, bc2_loc)

    mid = (bc2_loc + bc_len, bc1_loc)

    right = (bc1_loc + bc_len, len(consensus))

    return {
        "left": left,
        "mid": mid,
        "right": right
    }

# -------------------------------
# entropy
# -------------------------------
from collections import Counter
import math

def calc_entropy(column):
    """
    计算一个位置(column)的香农熵

    column:
        所有reads在同一位置的碱基列表
        例如:
            ['A','A','A','T','G']
    """

    counts = Counter(column)

    total = len(column)

    entropy = 0

    for count in counts.values():

        p = count / total

        entropy -= p * math.log2(p)

    return entropy


def global_entropy_profile(seqs):
    """
    计算:
        1. 每个位点的entropy
        2. consensus sequence

    参数
    ----
    seqs : list[str]

    返回
    ----
    entropy_list : list[float]
        每个位点entropy

    consensus : str
        consensus sequence
    """

    # 使用最短read长度
    L = min(len(s) for s in seqs)

    entropy_list = []

    consensus = []

    for i in range(L):

        # --------------------------------
        # 所有reads在第i个位点的碱基
        # --------------------------------
        column = []

        for s in seqs:

            base = s[i]

            if base in "ACGT":
                column.append(base)

        # 没有有效碱基
        if len(column) == 0:

            entropy_list.append(0)

            consensus.append("N")

            continue

        # --------------------------------
        # entropy
        # --------------------------------
        ent = calc_entropy(column)

        entropy_list.append(ent)

        # --------------------------------
        # consensus
        # --------------------------------
        counts = Counter(column)

        best_base = counts.most_common(1)[0][0]

        consensus.append(best_base)

    return entropy_list, "".join(consensus)

def find_high_entropy_regions(
        entropy,
        threshold=1.7,
        min_len=6):

    regions = []

    in_region = False

    start = 0

    for i, e in enumerate(entropy):

        if e >= threshold:

            if not in_region:

                start = i
                in_region = True

        else:

            if in_region:

                end = i

                if end - start >= min_len:

                    regions.append((start, end))

                in_region = False

    # tail
    if in_region:

        end = len(entropy)

        if end - start >= min_len:

            regions.append((start, end))

    return regions


# -------------------------------
# main
# -------------------------------
def scan(config,method):

    MAX_READS = 100000
    FASTQ = get_config(config, "file2")
    BARCODE_FILE = get_config(config, "Barcode")

    seqs = read_fastq_head(FASTQ, MAX_READS)
    print(f"reads loaded: {len(seqs)}")

    # barcode
    bc_set, bc_len = load_barcodes(BARCODE_FILE)
    bc_hits = scan_positions_hash(seqs, bc_set, bc_len)

    ratio = bc_hits / len(seqs)
    idx = np.argsort(ratio)[-2:]
    idx.sort()

    bc2_loc, bc1_loc = idx

    print("\n=== Barcode ===")
    print(f"bc2\t{bc2_loc}")
    print(f"bc1\t{bc1_loc}")

    entropy, consensus = global_entropy_profile(seqs)
    #plot_entropy(entropy,bc1_loc,bc2_loc)
    

    entropy = np.array(entropy)

    regions = find_high_entropy_regions(entropy,min_len=7)

    print(regions)
    

    # -----------------------------------------------------
    # output
    # -----------------------------------------------------

    print("\n=== Entropy Profile ===")

    #for i, (e, b) in enumerate(zip(entropy, consensus)):

        #print(f"{i}\t{e:.3f}\t{b}")
    return {
        "bc2": bc2_loc,
        "bc1": bc1_loc,
        "bc_len": bc_len
    }


import matplotlib.pyplot as plt

def plot_entropy(entropy, bc1_loc , bc2_loc, title="Global Entropy Profile" ):

    plt.figure(figsize=(14, 5))

    plt.plot(entropy)

    plt.xlabel("Position")

    plt.ylabel("Shannon entropy")
    plt.axvline(bc1_loc, linestyle="--")
    plt.axvline(bc1_loc + 8, linestyle="--")
    plt.axvline(bc2_loc, linestyle="--")
    plt.axvline(bc2_loc + 8, linestyle="--")
    plt.title(title)

    plt.ylim(0, 2.1)

    plt.grid(alpha=0.3)

    plt.tight_layout()

    plt.show()

