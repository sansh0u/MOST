#!/usr/bin/env python3
"""
dmt_preprocess.py

Refactored entry point for DMT preprocessing.
This module exposes run() for importing from a Typer/Click CLI,
while keeping argparse compatibility when executed directly.

The complete processing logic should be placed inside run().
"""

import argparse, gzip, os, sys
from fuzzysearch import find_near_matches
from Bio import SeqIO
from Bio.Seq import Seq



dmux_handles_R1 = {}
dmux_handles_R2 = {}

OPEN_HANDLE_CAP = 7000
def open_fastq_file(fp):
    return gzip.open(fp, "rt") if fp.endswith(".gz") else open(fp, "r")


def write_fastq_record(h, rid, seq, qual_ints):
    h.write("@{}\n".format(rid))
    h.write(str(seq) + "\n")
    h.write("+\n")
    h.write("".join(chr(q + 33) for q in qual_ints) + "\n")


def load_whitelist(path, col=4):
    if not path:
        return None
    wl = set()
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= col:
                wl.add(parts[col - 1])
    return wl


def assign_barcode_hamming1(bc, wl, ambig_policy="drop"):
    if bc in wl:
        return bc
    hits = []
    bases = ("A","C","G","T")
    bc_list = list(bc)
    for i, orig in enumerate(bc_list):
        for b in bases:
            if b == orig: continue
            cand = bc[:i] + b + bc[i+1:]
            if cand in wl:
                hits.append(cand)
                if ambig_policy == "first":
                    return cand
    if not hits:
        return None
    if ambig_policy == "first":
        return sorted(hits)[0]
    return None

def get_dmux_handles(barcode):
    if barcode not in dmux_handles_R1:
        # If you want LRU safety, call _maybe_evict_lru() here
        r1p = os.path.join(args.dmux_dir, f"{barcode}_R1.{args.chunk_id}.fq")
        r2p = os.path.join(args.dmux_dir, f"{barcode}_R2.{args.chunk_id}.fq")
        # Use large buffers to cut syscalls (1 MiB)
        dmux_handles_R1[barcode] = open(r1p, "a", buffering=1024*1024)
        dmux_handles_R2[barcode] = open(r2p, "a", buffering=1024*1024)
    return dmux_handles_R1[barcode], dmux_handles_R2[barcode]

def run(read1_fastq, read2_fastq, output_prefix="out",
        adapters=None,
        linker1="GTGGTTGATGTTTTGTATTGGTGTATGATT",
        linker2="ATTTATGTGTTTGAGAGGTTAGAGTATTTG",
        trim_end2="AGATGTGTATAAGAGATAG",
        min_len=10,
        whitelist=None,
        whitelist_col=4,
        dmux_dir=None,
        chunk_id="",
        barcode_len=16,
        bc_mismatch=0,
        bc_ambig="drop"):

    if adapters is None:
        adapters=["CTATCTCTTATA","AGATGCGAGAAGCCAACGCTTG"]

    os.makedirs(os.path.dirname(output_prefix) or ".", exist_ok=True)

    wl = load_whitelist(whitelist, whitelist_col) if whitelist else None
    do_dmux = wl is not None and dmux_dir is not None

    if do_dmux:
        os.makedirs(dmux_dir, exist_ok=True)

    r1_out = gzip.open(output_prefix + "_R1.fq.gz", "wt")
    r2_out = gzip.open(output_prefix + "_R2.fq.gz", "wt")
    lambda_r1_out = gzip.open(output_prefix + "_Lambda_R1.fq.gz", "wt")
    lambda_r2_out = gzip.open(output_prefix + "_Lambda_R2.fq.gz", "wt")

    n_reads = n_pass = n_trim = n_trim_bases = n_reads_barcode85 = n_written_barcode85 = 0
    BARCODE85_MOTIF = "TTATTTTT"

    with open_fastq_file(read1_fastq) as r1h, open_fastq_file(read2_fastq) as r2h:
        for rec1, rec2 in zip(SeqIO.parse(r1h, "fastq"), SeqIO.parse(r2h, "fastq")):
            n_reads += 1

            r1seq = rec1.seq
            r2seq = rec2.seq
            is_barcode85 = BARCODE85_MOTIF in str(r1seq) or BARCODE85_MOTIF in str(r2seq)
            if is_barcode85:
                n_reads_barcode85 += 1
            lambda_seq1 = r1seq
            lambda_qual1 = rec1.letter_annotations["phred_quality"]

            hits = []
            for a in adapters:
                hits.extend(find_near_matches(a, str(r1seq), max_l_dist=1))

            if hits:
                trim_i = min(m.start for m in hits)
                if trim_i < len(lambda_seq1):
                    n_trim += 1
                    n_trim_bases += (len(lambda_seq1) - trim_i)
                    lambda_seq1 = lambda_seq1[:trim_i]
                    lambda_qual1 = lambda_qual1[:trim_i]

            te2 = find_near_matches(trim_end2, str(r2seq), max_l_dist=1)

            if te2:
                trim2_i = te2[0].end
                lambda_seq2 = r2seq[trim2_i:]
                lambda_qual2 = rec2.letter_annotations["phred_quality"][trim2_i:]
            else:
                lambda_seq2 = r2seq
                lambda_qual2 = rec2.letter_annotations["phred_quality"]

            if is_barcode85 and len(lambda_seq1) >= args.min_len and len(lambda_seq2) >= args.min_len:
                write_fastq_record(lambda_r1_out, rec1.id, lambda_seq1, lambda_qual1)
                write_fastq_record(lambda_r2_out, rec2.id, lambda_seq2, lambda_qual2)
                n_written_barcode85 += 1
            lk1 = find_near_matches(linker1, str(r2seq), max_l_dist=2)
            lk2 = find_near_matches(linker2, str(r2seq), max_l_dist=2)
            if not (len(lk1)==1 and len(lk2)==1 and len(te2)==1):
                continue

            l1s = lk1[0].start; l2s = lk2[0].start
            b1s = max(0, l1s-8); b2s = max(0, l2s-8)
            barcode = str(r2seq[b1s:l1s] + r2seq[b2s:l2s])
            
            if barcode in barcode_cnt:
                barcode_cnt[barcode] += 1
            else:
                barcode_cnt[barcode] = 1
                    
            if len(barcode) != args.barcode_len:
                continue

            new_seq1 = lambda_seq1
            new_qual1 = lambda_qual1
            new_seq2 = lambda_seq2
            new_qual2 = lambda_qual2

            if len(new_seq1) < args.min_len or len(new_seq2) < args.min_len:
                continue

            # Write trimmed chunk-level outputs (keep sample-wide merge compatible)
            rid1 = f"{barcode}_{rec1.id}"
            rid2 = f"{barcode}_{rec2.id}"
            write_fastq_record(r1_out, rid1, new_seq1, new_qual1)
            write_fastq_record(r2_out, rid2, new_seq2, new_qual2)
            n_pass += 1
            
            assigned = None
            if do_dmux and wl is not None:
                if args.bc_mismatch <= 0:
                    if barcode in wl:
                        assigned = barcode
                else:
                # allow 1 mismatch (you can allow >1 but not recommended)
                    if barcode in wl:
                        assigned = barcode
                    elif args.bc_mismatch >= 1:
                        assigned = assign_barcode_hamming1(barcode, wl, args.bc_ambig)

            if assigned is not None:
                h1, h2 = get_dmux_handles(assigned)
                h1.write(f"@{rid1}\n{str(new_seq1)}\n+\n{''.join(chr(q+33) for q in new_qual1)}\n")
                h2.write(f"@{rid2}\n{str(new_seq2)}\n+\n{''.join(chr(q+33) for q in new_qual2)}\n")
                dmux_counts_R1[assigned] = dmux_counts_R1.get(assigned, 0) + 1
                dmux_counts_R2[assigned] = dmux_counts_R2.get(assigned, 0) + 1
            elif do_dmux:
                h1, h2 = get_dmux_handles("UNMATCHED")
                h1.write(f"@{rid1}\n{str(new_seq1)}\n+\n{''.join(chr(q+33) for q in new_qual1)}\n")
                h2.write(f"@{rid2}\n{str(new_seq2)}\n+\n{''.join(chr(q+33) for q in new_qual2)}\n")
                dmux_counts_R1["UNMATCHED"] = dmux_counts_R1.get("UNMATCHED", 0) + 1
                dmux_counts_R2["UNMATCHED"] = dmux_counts_R2.get("UNMATCHED", 0) + 1



    r1_out.close()
    r2_out.close()
    lambda_r1_out.close()
    lambda_r2_out.close()