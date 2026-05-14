import gzip
import numpy as np
from yaml_load import get_config
from collections import Counter
import math
import matplotlib.pyplot as plt

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


# -------------------------------
# entropy
# -------------------------------


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
# -----------------------------------------------------
    # plot entropy profile
# -----------------------------------------------------

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

def filtered(regions,bc1_loc,bc2_loc):
    filtered_regions = []
    for start, end in regions:

        if end >= 117:
            continue

        overlap_bc1 = not (end <= bc1_loc or start >= bc1_loc+8)

        overlap_bc2 = not (end <= bc2_loc or start >= bc2_loc+8)

        if overlap_bc1 or overlap_bc2:
            continue

        filtered_regions.append((start, end))
    return  filtered_regions

def hamming(s1, s2):

    if len(s1) != len(s2):
        raise ValueError(
            "Sequences must have same length"
        )

    return sum(
        c1 != c2
        for c1, c2 in zip(s1, s2)
    )

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
# -------------------------------
# main
# -------------------------------
def scan(cfg,method):

    MAX_READS = 100000
    fastq = get_config(cfg, "file2")
    barcode_file = get_config(cfg, "barcode_file")

    seqs = read_fastq_head(fastq, MAX_READS)
    print(f"reads loaded: {len(seqs)}")
    read_len = len(seqs[0])
    # barcode
    bc_set, bc_len = load_barcodes(barcode_file)
    bc_hits = scan_positions_hash(seqs, bc_set, bc_len)

    ratio = bc_hits / len(seqs)
    idx = np.argsort(ratio)[-2:]
    idx.sort()

    bc2_loc, bc1_loc = idx

    print("\n=== Barcode ===")
    print(f"bc2 starts at bp {bc2_loc+1}")
    print(f"bc1 starts at bp {bc1_loc+1}")

    entropy, _ = global_entropy_profile(seqs)
    #plot_entropy(entropy,bc1_loc,bc2_loc)
    #print("\n=== Entropy Profile ===")

    #for i, (e, b) in enumerate(zip(entropy, consensus)):

        #print(f"{i}\t{e:.3f}\t{b}")
    entropy = np.array(entropy)

    print("\n=== UMI ===")
    regions = find_high_entropy_regions(entropy,min_len=6)
    
    for start, end in regions:
        print(f"UMI candidate region: bp {start+1}-{end}")
        
    print("\n=== Filtered UMI regions ===")
    filtered_regions = filtered(regions,bc1_loc,bc2_loc)
    umi_start = umi_end = umi_len = 0
    
    if method == "ATAC":
        if not filtered_regions:
            print("No UMI region found.")
            
        else:
            
            for start, end in filtered_regions:
                print(f"UMI may be located between bp {start+1}-{end}")
            umi_start, umi_end = filtered_regions[0]
            
    elif method == "RNA":
        if len(filtered_regions) != 1:
            raise ValueError(
                f"Expected exactly 1 UMI region, found {len(filtered_regions)}: {filtered_regions}"
            )
        
        for start, end in filtered_regions:
            print(f"UMI may be located between bp {start+1}-{end}")
        umi_start, umi_end = filtered_regions[0]
   
    umi_len = umi_end - umi_start 
    

    # -----------------------------------------------------
    # output
    # -----------------------------------------------------

   
    return {
        "bc2": bc2_loc,
        "bc1": bc1_loc,
        "bc_len": bc_len,
        "read_len": read_len,
        'umi_start': umi_start,
        'umi_len': umi_len
    }

def check_adapter(cfg):
    
    MAX_READS = 100000
    fastq = get_config(cfg, "file1")
    seqs1 = read_fastq_head(fastq, MAX_READS)
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
                
    else:
        print("\n=== Adapter skipped (no Adapter in config) ===")

    return score




    
