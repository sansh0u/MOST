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



def segment_by_entropy(entropy,
                       threshold=1.0,
                       gradient_thres=0.15,
                       min_len=5,
                       smooth_window=3):

    """
    使用:
    1. smoothed entropy
    2. entropy gradient
    进行区块切分
    """

    entropy = smooth_entropy(entropy, smooth_window)

    # entropy梯度
    grad = np.gradient(entropy)

    segments = []

    start = 0
    state = entropy[0] < threshold

    for i in range(1, len(entropy)-1):

        cur_state = entropy[i] < threshold

        # -------------------------
        # 条件1:
        # 低熵/高熵状态改变
        # -------------------------
        state_change = (cur_state != state)

        # -------------------------
        # 条件2:
        # entropy快速变化
        # -------------------------
        sharp_drop = abs(grad[i]) > gradient_thres

        if state_change or sharp_drop:

            # 避免过短segment
            if i - start >= min_len:
                segments.append((start, i, state))

                start = i
                state = cur_state

    # last segment
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
    raw_entropy, consensus = global_entropy_profile(seqs)

# -------------------------------
# smoothed entropy
# -------------------------------
    smooth_ent = smooth_entropy(raw_entropy, window=5)
    """
    print("\n=== Entropy Profile ===")

    for i, (e, b) in enumerate(zip(raw_entropy, consensus)):

        print(f"{i}\t{e:.3f}\t{b}")
    """
    # -------------------------------
    # split by BC anchors
    # -------------------------------
    regions = split_by_barcode(
        consensus,
        bc2_loc,
        bc1_loc,
        bc_len
    )

    print("\n=== Regions ===")

    for name, (s, e) in regions.items():

        seq = consensus[s:e]

        print(f"{name}\t{s}-{e}\t{seq}")

    # -------------------------------
    # MID = linker2
    # -------------------------------
    mid_s, mid_e = regions["mid"]

    is_linker, score = detect_linker(
        smooth_ent,
        mid_s,
        mid_e
    )

    print("\n=== Linker2 ===")
    print(f"{mid_s}-{mid_e}")
    print(f"mean entropy = {score:.3f}")

    # -------------------------------
    # RIGHT: search UMI
    # -------------------------------
    right_s, right_e = regions["right"]

    umi = find_umi_region(
    raw_entropy,
    right_s,
    right_e,
    umi_len=10,
    jump_thres=0.6
)

    print("\n=== RIGHT ===")

    if umi:

        us, ue, score = umi

        print(f"UMI\t{us}-{ue}\tentropy={score:.3f}")

        linker1 = (right_s, us)

        print(f"LINKER1\t{linker1[0]}-{linker1[1]}")
        print(consensus[linker1[0]:linker1[1]])

        print(f"UMI_SEQ\t{consensus[us:ue]}")

        print(f"RNA\t{ue}-{right_e}")

    else:

        print("No UMI detected")

    # -------------------------------
    # LEFT: optional UMI
    # -------------------------------
    left_s, left_e = regions["left"]

    left_umi = find_umi_region(
        raw_entropy,
        left_s,
        left_e,
        umi_len=10
    )

    print("\n=== LEFT ===")

    if left_umi:

        us, ue, score = left_umi

        print(f"LEFT_UMI\t{us}-{ue}")

        print(f"PRIMER\t{left_s}-{us}")

    else:

        print(f"PRIMER\t{left_s}-{left_e}")
    # segmentation

    return {
        "bc2": bc2_loc,
        "bc1": bc1_loc,
        "raw_entropy": raw_entropy,
        "smoothed_entropy": smooth_ent,
        "consensus": consensus
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
# entropy smoothing
# -------------------------------
def smooth_entropy(entropy, window=5):

    smoothed = []

    for i in range(len(entropy)):

        s = max(0, i - window)
        e = min(len(entropy), i + window + 1)

        smoothed.append(np.mean(entropy[s:e]))

    return np.array(smoothed)


# -------------------------------
# detect umi
# -------------------------------
def find_umi_region(entropy,
                    start,
                    end,
                    umi_len=10,
                    jump_thres=0.25):

    """
    找:
    entropy突然升高的位置
    """

    best_pos = None
    best_jump = -1

    for i in range(start + 1, end - umi_len):

        jump = entropy[i] - entropy[i - 1]

        if jump > best_jump:

            best_jump = jump
            best_pos = i

    if best_jump >= jump_thres:

        return (
            best_pos,
            best_pos + umi_len,
            best_jump
        )

    return None

# -------------------------------
# linker detection
# -------------------------------
def detect_linker(entropy,
                  start,
                  end,
                  low_thres=1.0):

    region = entropy[start:end]

    mean_entropy = np.mean(region)

    if mean_entropy < low_thres:
        return True, mean_entropy

    return False, mean_entropy