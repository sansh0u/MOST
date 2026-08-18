import subprocess
from most.preprocess.scan_bc import check_adapter


def qc_adapt(fastq_intput_1,fastq_intput_2,CleanFq1,CleanFq2,threads,Adapter,mismatch):
    print("START check_adapter", flush=True)
    score = check_adapter(fastq_intput_1,Adapter,mismatch)
    print(f"adapter score = {score}", flush=True)
    cmd1 = [
    "cutadapt", "-m", "18", "-a", "A{10}N{150}",
    "--times", "4",
    "-g", Adapter,
    "--pair-filter=any",
    "-j", str(threads),
    "-o", CleanFq1,
    "-p", CleanFq2,
    fastq_intput_1, fastq_intput_2
]

    cmd2 = [
    "cutadapt", "-m", "18", "-a", "A{10}N{150}",
    "--times", "4",
    "--pair-filter=any",
    "-j", str(threads),
    "-o", CleanFq1,
    "-p", CleanFq2,
    fastq_intput_1, fastq_intput_2
]
    try: 
        if score >= 0.5:
            print("Trimming Adapter", flush=True)
            print("CMD:", " ".join(cmd1), flush=True)
            subprocess.run(cmd1, check=True)
            
        else:
            print("Trimming ployA", flush=True)
            print("CMD:", " ".join(cmd2), flush=True)
            subprocess.run(cmd2, check=True)
            
    except subprocess.CalledProcessError as e:
        raise
    print("QC finished", flush=True)