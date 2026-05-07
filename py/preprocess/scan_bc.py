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


def calc_entropy(column):
    counts = Counter(column)
    total = len(column)
    ent = 0
    for v in counts.values():
        p = v / total
        ent -= p * math.log2(p)
    return ent


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
# entropy + consensus
# -------------------------------
def global_entropy_profile(seqs):
    L = len(seqs[0])

    entropy_list = []
    consensus = []

    for i in range(L):
        column = [s[i] for s in seqs]
        counts = Counter(column)

        ent = calc_entropy(column)
        entropy_list.append(ent)

        base = counts.most_common(1)[0][0]
        consensus.append(base)

    return entropy_list, "".join(consensus)


# -------------------------------
# 🔥 区块切分（核心新增）
# -------------------------------
def segment_by_entropy(entropy, threshold=1.0,jump=0.5, min_len=5):
    """
    根据entropy切分
    低熵 = 结构区barcode 
    高熵 = 随机区UMI / RNA
    """

    segments = []
    start = 0
    state = entropy[0] < threshold  # True = 低熵
    
    for i in range(1, len(entropy)):
        cur_state = entropy[i] < threshold
        if (cur_state != state) or (entropy[i] - entropy[i-1] > jump):

            if i - start >= min_len:
                segments.append((start, i, state))

            start = i
            state = cur_state

    if len(entropy) - start >= min_len:
        segments.append((start, len(entropy), state))

    return segments
        


# -------------------------------
# main pipeline
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

    # entropy
    entropy, consensus = global_entropy_profile(seqs)

    # segmentation
    raw_segments = segment_by_entropy(entropy)

    segments = []
    for s, e, low_entropy in raw_segments:
        refined = refine_segment(seqs, s, e)
        for rs, re in refined:
            segments.append((rs, re, low_entropy))

    print("\n=== Segments ===")
    for s, e, low_entropy in segments:
        label = "STRUCTURE" if low_entropy else "RANDOM"
        print(f"{s}-{e}\t{label}\t{consensus[s:e]}")

    return {
        "bc2": bc2_loc,
        "bc1": bc1_loc,
        "entropy": entropy,
        "consensus": consensus,
        "segments": segments
    }


def get_freq(bases):
    total = len(bases)
    freq = {b:0 for b in 'ACGT'}
    for b in bases:
        if b in freq:
            freq[b] += 1
    for b in freq:
        freq[b] /= total
    return freq

def entropy_from_freq(freq):
    H = 0
    for p in freq.values():
        if p > 0:
            H -= p * math.log2(p)
    return H

def jsd(freq1, freq2):
    M = {b:(freq1[b]+freq2[b])/2 for b in freq1}
    return entropy_from_freq(M) - 0.5*(entropy_from_freq(freq1)+entropy_from_freq(freq2))

def find_best_split(seqs, start, end):
    best_k = None
    best_score = -1

    for k in range(start+2, end-2):  # 避免太短
        left = []
        right = []

        for s in seqs:
            left.extend(s[start:k])
            right.extend(s[k:end])

        f1 = get_freq(left)
        f2 = get_freq(right)

        score = jsd(f1, f2)

        if score > best_score:
            best_score = score
            best_k = k

    return best_k, best_score

def refine_segment(seqs, start, end, jsd_thres=0.08, min_len=6):
    k, score = find_best_split(seqs, start, end)

    if k is None or score < jsd_thres or (end - start) < min_len:
        return [(start, end)]

    left = refine_segment(seqs, start, k, jsd_thres, min_len)
    right = refine_segment(seqs, k, end, jsd_thres, min_len)

    return left + right