import subprocess
from most.preprocess.scan_bc import check_adapter


def qc_adapt(fastq_intput_1,fastq_intput_2,CleanFq1,CleanFq2,threads,Adapter,mismatch):
    score = check_adapter(fastq_intput_1,Adapter,mismatch)
    cmd1 = [
    "cutadapt", "-m", "18", "-a", "A{10}N{150}",
    "--times", "4",
    "-g", Adapter,
    "-j", str(threads),
    "-o", CleanFq1,
    "-p", CleanFq2,
    fastq_intput_1, fastq_intput_2
]

    cmd2 = [
    "cutadapt", "-m", "18", "-a", "A{10}N{150}",
    "--times", "4",
    "-j", str(threads),
    "-o", CleanFq1,
    "-p", CleanFq2,
    fastq_intput_1, fastq_intput_2
]
    try: 
        if score >= 0.5:
            print(cmd1)
            print("Trimming Adapter")
            subprocess.run(cmd1, check=True)
            
        else:
            print(cmd2)
            print("Trimming ployA")
            subprocess.run(cmd2, check=True)
            
    except subprocess.CalledProcessError as e:
        raise

