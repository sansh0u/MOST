import gzip
import numpy as np
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


def filtered(regions, bc1_loc, bc2_loc, umi_len):
    filtered_regions = []

    for start, end in regions:

        # 只保留指定长度的 UMI
        # 不符合长度的候选区域直接跳过
        if end - start != umi_len:
            print(
                f"Skip UMI candidate region: "
                f"({start}, {end}), "
                f"length={end-start}, "
                f"expected={umi_len}"
            )
            continue

        # 限制 UMI 区域不能超过指定位置
        if end >= 117:
            continue

        # barcode overlap
        overlap_bc1 = not (
            end <= bc1_loc or
            start >= bc1_loc + 8
        )

        overlap_bc2 = not (
            end <= bc2_loc or
            start >= bc2_loc + 8
        )

        if overlap_bc1 or overlap_bc2:
            continue

        filtered_regions.append((start, end))

    return filtered_regions

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
def scan(cfg, method):

    MAX_READS = 100000

    fastq = cfg.sequence_file.file2
    barcode_file = cfg.reference.barcode_file

    # -----------------------------------------------------
    # UMI length comes from YAML
    # -----------------------------------------------------
    umi_len = cfg.advanced.umi_len

    if umi_len is None:
        raise ValueError(
            "umi_len is not provided in YAML/config"
        )

    umi_len = int(umi_len)

    print(f"Configured UMI length: {umi_len} bp")

    # -----------------------------------------------------
    # Read FASTQ
    # -----------------------------------------------------
    seqs = read_fastq_head(fastq, MAX_READS)

    print(f"reads loaded: {len(seqs)}")

    if not seqs:
        raise ValueError("No reads found in FASTQ")

    read_len = len(seqs[0])

    # -----------------------------------------------------
    # Barcode
    # -----------------------------------------------------
    bc_set, bc_len = load_barcodes(barcode_file)

    bc_hits = scan_positions_hash(
        seqs,
        bc_set,
        bc_len
    )

    ratio = bc_hits / len(seqs)

    idx = np.argsort(ratio)[-2:]
    idx.sort()

    bc2_loc, bc1_loc = idx

    print("\n=== Barcode ===")
    print(f"bc2 starts at bp {bc2_loc + 1}")
    print(f"bc1 starts at bp {bc1_loc + 1}")

    # -----------------------------------------------------
    # Entropy
    # -----------------------------------------------------
    entropy, _ = global_entropy_profile(seqs)

    entropy = np.array(entropy)

    print("\n=== UMI ===")

    # -----------------------------------------------------
    # Find high entropy regions
    # -----------------------------------------------------
    regions = find_high_entropy_regions(
        entropy,
        min_len=6
    )

    for start, end in regions:
        print(
            f"UMI candidate region: "
            f"bp {start + 1}-{end} "
            f"(length={end-start})"
        )

    # -----------------------------------------------------
    # Filter UMI regions
    # -----------------------------------------------------
    print("\n=== Filtered UMI regions ===")

    filtered_regions = filtered(
        regions,
        bc1_loc,
        bc2_loc,
        umi_len
    )

    # 不要再把 umi_len 改成 0
    umi_start = 0
    umi_end = 0

    # -----------------------------------------------------
    # ATAC / DMT
    # -----------------------------------------------------
    if method == "ATAC" or method == "DMT":

        if not filtered_regions:

            print(
                f"No valid {umi_len} bp UMI region found."
            )

        else:

            for start, end in filtered_regions:

                print(
                    f"UMI may be located between "
                    f"bp {start + 1}-{end}"
                )

            umi_start, umi_end = filtered_regions[0]

    # -----------------------------------------------------
    # RNA
    # -----------------------------------------------------
    elif method == "RNA":

        if len(filtered_regions) != 1:

            raise ValueError(
                f"Expected exactly 1 UMI region, "
                f"found {len(filtered_regions)}: "
                f"{filtered_regions}"
            )

        start, end = filtered_regions[0]

        print(
            f"UMI may be located between "
            f"bp {start + 1}-{end}"
        )

        umi_start = start
        umi_end = end

    # -----------------------------------------------------
    # Output
    # -----------------------------------------------------
    return {
        "bc2": int(bc2_loc),
        "bc1": int(bc1_loc),
        "bc_len": int(bc_len),
        "read_len": int(read_len),
        "umi_start": int(umi_start),
        "umi_len": int(umi_len)
    }

def check_adapter(fastq,Adapter,mismatch):
    
    MAX_READS = 100000
    seqs1 = read_fastq_head(fastq, MAX_READS)

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

def scan_len(cfg):
    fastq = cfg.sequence_file.file2
    seqs = read_fastq_head(fastq, 1)
    print(f"reads loaded: {len(seqs)}")
    read_len = len(seqs[0])
    return read_len



    
