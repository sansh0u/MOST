rule all:
    input:
        f"{config['out_dir']}/astro_output"

rule astro:
    input:
        r1=config["sequence_file"]["file1"],
        r2=config["sequence_file"]["file2"],
        barcode_file=config["reference"]["barcode_file"],
        gtf_file=config["reference"]["gtf_file"],
        star_index=config["reference"]["star_index"]
    output:
        r1=directory(f"{config['out_dir']}/astro_output")
    params:
        primer_structure=f"{config["advanced"]["adapter"]}_b_A{{10}}N{{150}}"
        structure_umi=f"{config["advanced"]["primer"]}_{config["runtime"]["umi_len"]}"
        structure_barcode=f"20_{config["advanced"]["linker2"]}:{config["advanced"]["linker2"]}...{config["advanced"]["linker1"]}"
    threads:
        config["threads"]
    shell:
        """
        ASTRO --R1 {input.r1} --R2 {input.r2} \
            --barcode_file {input.barcode_file} \
            --gtf_file {input.gtf_file} \
            --starref {input.star_index} \
            --PrimerStructure {params.primer_structure} \
            --StructureUMI {params.structure_umi} \
            --StructureBarcode {params.structure_barcode} \
            --threadnum {threads} \
            --steps 7 \
            --outputfolder {output.r1}
        """