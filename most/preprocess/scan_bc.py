import gzip
import math
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter

base2int = {"A": 0, "C": 1, "G": 2, "T": 3}


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
    if bc_len is None:
        raise ValueError(f"No valid barcodes found in {path}")
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
                code = (code & mask) * 4 + base2int[next_base]
    return bc_hits


def calc_entropy(column):
    counts = Counter(column)
    total = len(column)
    return -sum((count / total) * math.log2(count / total)
                for count in counts.values())


def global_entropy_profile(seqs):
    L = min(len(s) for s in seqs)
    entropy_list = []
    consensus = []
    for i in range(L):
        column = [s[i] for s in seqs if s[i] in "ACGT"]
        if not column:
            entropy_list.append(0)
            consensus.append("N")
            continue
        entropy_list.append(calc_entropy(column))
        consensus.append(Counter(column).most_common(1)[0][0])
    return entropy_list, "".join(consensus)


def find_high_entropy_regions(entropy, threshold=1.7, min_len=6):
    regions = []
    in_region = False
    start = 0
    for i, e in enumerate(entropy):
        if e >= threshold:
            if not in_region:
                start = i
                in_region = True
        elif in_region:
            if i - start >= min_len:
                regions.append((start, i))
            in_region = False
    if in_region and len(entropy) - start >= min_len:
        regions.append((start, len(entropy)))
    return regions


def filtered(regions, bc1_loc, bc2_loc, umi_len=None):
    filtered_regions = []
    for start, end in regions:
        if umi_len is not None and end - start != umi_len:
            print(f"Skip UMI candidate region: ({start}, {end}), "
                  f"length={end-start}, expected={umi_len}")
            continue
        if end >= 117:
            continue
        overlap_bc1 = not (end <= bc1_loc or start >= bc1_loc + 8)
        overlap_bc2 = not (end <= bc2_loc or start >= bc2_loc + 8)
        if overlap_bc1 or overlap_bc2:
            continue
        filtered_regions.append((start, end))
    return filtered_regions


def hamming(s1, s2):
    if len(s1) != len(s2):
        raise ValueError("Sequences must have same length")
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def scan_adapter_positions(seqs, query, max_mismatch, window):
    query = query.upper()
    k = len(query)
    L = len(seqs[0])
    hits = np.zeros(L)
    for seq in seqs:
        if len(seq) < k:
            continue
        max_i = min(window, len(seq) - k + 1)
        for i in range(max_i):
            if hamming(seq[i:i + k], query) <= max_mismatch:
                hits[i] += 1
    return hits


def scan(cfg, method):
    MAX_READS = 100000
    fastq = cfg.sequence_file.file2
    barcode_file = cfg.reference.barcode_file
    umi_len = cfg.advanced.umi_len
    umi_len = 10 if umi_len is None else int(umi_len)
    if umi_len < 0:
        raise ValueError("umi_len must be >= 0")
    print(f"Configured UMI length: {umi_len} bp")

    seqs = read_fastq_head(fastq, MAX_READS)
    print(f"reads loaded: {len(seqs)}")
    if not seqs:
        raise ValueError("No reads found in FASTQ")
    read_len = len(seqs[0])

    bc_set, bc_len = load_barcodes(barcode_file)
    bc_hits = scan_positions_hash(seqs, bc_set, bc_len)
    ratio = bc_hits / len(seqs)
    if len(ratio) < 2:
        raise ValueError("Read length is too short to detect two barcode positions")
    idx = np.argsort(ratio)[-2:]
    idx.sort()
    bc2_loc, bc1_loc = idx

    print("\n=== Barcode ===")
    print(f"bc2 starts at bp {bc2_loc + 1}")
    print(f"bc1 starts at bp {bc1_loc + 1}")

    # umi_len=0: retain barcode positions and skip UMI detection
    if umi_len == 0:
        print("\n=== UMI ===")
        print("umi_len = 0; UMI detection skipped.")
        return {
            "bc2": int(bc2_loc),
            "bc1": int(bc1_loc),
            "bc_len": int(bc_len),
            "read_len": int(read_len),
            "umi_start": 0,
            "umi_len": 0
        }

    entropy, _ = global_entropy_profile(seqs)
    entropy = np.array(entropy)
    print("\n=== UMI ===")
    regions = find_high_entropy_regions(entropy, min_len=6)

    for start, end in regions:
        print(f"UMI candidate region: bp {start + 1}-{end} "
              f"(length={end-start})")

    print("\n=== Filtered UMI regions ===")
    filtered_regions = filtered(regions, bc1_loc, bc2_loc, umi_len)
    umi_start = 0
    umi_end = 0

    if method in ("ATAC", "DMT"):
        if not filtered_regions:
            print(f"No valid {umi_len} bp UMI region found.")
        else:
            for start, end in filtered_regions:
                print(f"UMI may be located between bp {start + 1}-{end}")
            umi_start, umi_end = filtered_regions[0]

    elif method == "RNA":
        if len(filtered_regions) != 1:
            raise ValueError(
                f"Expected exactly 1 UMI region, found "
                f"{len(filtered_regions)}: {filtered_regions}"
            )
        start, end = filtered_regions[0]
        print(f"UMI may be located between bp {start + 1}-{end}")
        umi_start, umi_end = start, end

    return {
        "bc2": int(bc2_loc),
        "bc1": int(bc1_loc),
        "bc_len": int(bc_len),
        "read_len": int(read_len),
        "umi_start": int(umi_start),
        "umi_len": int(umi_len)
    }


def check_adapter(fastq, Adapter, mismatch):
    MAX_READS = 100000
    seqs1 = read_fastq_head(fastq, MAX_READS)
    if not seqs1:
        raise ValueError("No reads found in FASTQ")
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
        score = 0.0
        print("\n=== Adapter skipped (no Adapter in config) ===")
    return score


def scan_len(cfg):
    fastq = cfg.sequence_file.file2
    seqs = read_fastq_head(fastq, 1)
    if not seqs:
        raise ValueError("No reads found in FASTQ")
    print(f"reads loaded: {len(seqs)}")
    return len(seqs[0])