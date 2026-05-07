import gzip
import numpy as np
from yaml_load import get_config
from collections import Counter

# -------------------------------
# config
# -------------------------------

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

    bc_set = set()
    bc_len = None

    with open(path) as f:

        for line in f:

            parts = line.strip().split()

            if not parts:
                continue

            seq = parts[0].upper()

            # BC2+BC1 总长度
            if bc_len is None:
                bc_len = len(seq) // 2

            # 长度不一致跳过
            if len(seq) != bc_len * 2:
                continue

            # 取后半段（BC1）
            bc = seq[bc_len:]

            code = encode_kmer(bc)

            if code is not None:
                bc_set.add(code)

    return bc_set, bc_len



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

                code = (((code & mask) * 4)+ base2int[next_base])

    return bc_hits


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

    MAX_READS = 100000
    FASTQ = get_config(config, "file2")
    BARCODE_FILE = get_config(config, "Barcode")
    seqs = read_fastq_head(FASTQ, MAX_READS)
    print(f"reads loaded: {len(seqs)}")

    # ---------------- barcode ----------------
    bc_set, bc_len = load_barcodes(BARCODE_FILE)

    print(f"{len(bc_set)} barcode loaded")
    print(f"barcode_length\t{bc_len}")

    read_len = len(seqs[0])

    bc_hits = scan_positions_hash(seqs, bc_set, bc_len)

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
        "bc_len": bc_len,
        "read_len": read_len
    }
    linker2 = find_linker_region(
    seqs,
    start=bc2_loc + bc_len,
    end=bc1_loc
    )
    if linker2 is not None:

        result["linker2"] = linker2
# 在 BC1 后面找 linker1
    linker1 = find_linker_region(
    seqs,
    start=bc1_loc + bc_len,
    end=bc1_loc + bc_len + 50
    )
    if linker1 is not None:

        result["linker1"] = linker1


        
    # ---------------- adapter（可选） ----------------
    if method.upper() == "RNA":
        MAX_READS = 100000
        FASTQ = get_config(config, "file1")
        seqs1 = read_fastq_head(FASTQ, MAX_READS)
        print(f"reads loaded: {len(seqs1)}")
        Adapter = get_config(config, "adapter","AAGCAGTGGTATCAACGCAGAGTGAATGGG")
        mismatch = get_config(config, "adapter_mismatch", 2)

        if Adapter:
            hits = scan_adapter_positions(seqs1, Adapter, mismatch, 50)
            ratio = hits / len(seqs1)

            sub = ratio[:50]
            pos = sub.argmax()
            score = sub.max()

            print("\n=== Adapter result ===")
            print(f"position\t{pos}")
            print(f"ratio\t{score:.4f}")

            result = {
                "score": score
            }
        else:
            print("\n=== Adapter skipped (no Adapter in config) ===")

    else:
        print("\n=== Adapter skipped (ATAC mode) ===")

    return result

from collections import Counter
import numpy as np


def calc_entropy(column):

    total = len(column)

    counts = Counter(column)

    ent = 0.0

    for v in counts.values():

        p = v / total

        ent -= p * np.log2(p)

    return ent


def find_linker_region(
    seqs,
    start,
    end,
    entropy_threshold=1.0,
    jump_threshold=0.5,
    min_conserve=10
):
    """
    根据entropy寻找linker
    使用entropy jump自动截断
    """

    region_seqs = []

    for seq in seqs:

        if len(seq) < end:
            continue

        region_seqs.append(seq[start:end])

    if len(region_seqs) == 0:
        return None

    length = end - start

    entropy_list = []

    consensus = []

    # --------------------------------
    # 计算entropy
    # --------------------------------
    for i in range(length):

        column = [s[i] for s in region_seqs]

        counts = Counter(column)

        ent = calc_entropy(column)

        entropy_list.append(ent)

        base = counts.most_common(1)[0][0]

        consensus.append(base)

    consensus = "".join(consensus)

    # --------------------------------
    # 找低entropy起点
    # --------------------------------
    best_start = None

    for i, ent in enumerate(entropy_list):

        if ent < entropy_threshold:

            best_start = i

            break

    if best_start is None:

        print("\n=== No linker detected ===")

        return None

    # --------------------------------
    # entropy jump截断
    # --------------------------------
    best_end = length

    for i in range(best_start, length - 1):

        jump = entropy_list[i + 1] - entropy_list[i]

        # entropy突然升高
        if jump > jump_threshold:

            best_end = i + 1

            break

    # --------------------------------
    # 长度过滤
    # --------------------------------
    if best_end - best_start < min_conserve:

        print("\n=== Linker too short ===")

        return None

    linker_seq = consensus[best_start:best_end]

    linker_len = len(linker_seq)

    linker_start = start + best_start

    linker_end = linker_start + linker_len

    # --------------------------------
    # 输出
    # --------------------------------
    print("\n=== Linker detected ===")

    print(f"start\t{linker_start}")

    print(f"end\t{linker_end}")

    print(f"length\t{linker_len}")

    print(f"sequence\t{linker_seq}")

    # 调试用：显示entropy jump
    for i in range(best_start, best_end - 1):

        jump = entropy_list[i + 1] - entropy_list[i]

        if jump > jump_threshold:

            print(
                f"entropy_jump\t{i}->{i+1}\t"
                f"{entropy_list[i]:.3f} -> "
                f"{entropy_list[i+1]:.3f}"
            )

    return {
        "start": int(linker_start),
        "end": int(linker_end),
        "length": int(linker_len),
        "sequence": linker_seq
    }

def global_entropy_profile(seqs):
    MAX_READS = 100000
    FASTQ = get_config(config, "file2")
    
    seqs = read_fastq_head(FASTQ, MAX_READS)
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

    consensus = "".join(consensus)

    return entropy_list, consensus