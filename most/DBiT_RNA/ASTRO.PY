from pathlib import Path
import subprocess
import logging

def build_astro_command(cfg):

    """
    Convert internal config to ASTRO command
    """

    r1 = cfg.sequence_file.file1
    r2 = cfg.sequence_file.file2

    barcode_file = \
        cfg.reference.barcode_file

    gtf_file = \
        cfg.reference.gtf_file

    star_index = \
        cfg.reference.star_index

    threads = str(cfg.threads)

    out_dir = cfg.out_dir

    primer = cfg.advanced.primer

    linker1 = cfg.advanced.linker1

    linker2 = cfg.advanced.linker2

    umi_len = cfg.runtime.umi_len

    # -------------------------
    # UMI structure
    # -------------------------

    structure_umi = \
        f"{primer}_{umi_len}"

    # -------------------------
    # Barcode structure
    # -------------------------

    structure_barcode = \
        f"20_{linker2}:{linker2}...{linker1}"

    # -------------------------
    # PrimerStructure
    # -------------------------

    primer_structure = \
        f"{cfg.advanced.adapter}_b_A{{10}}N{{150}}"

    cmd = [

        "ASTRO",

        "--R1",
        r1,

        "--R2",
        r2,

        "--barcode_file",
        barcode_file,

        "--gtffile",
        gtf_file,

        "--starref",
        star_index,

        "--PrimerStructure",
        primer_structure,

        "--StructureUMI",
        structure_umi,

        "--StructureBarcode",
        structure_barcode,

        "--threadnum",
        threads,

        "--steps",
        "7",

        "--outputfolder",
        out_dir
    ]
    try: ####
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Error during ASTRO : {e}")
        raise
    