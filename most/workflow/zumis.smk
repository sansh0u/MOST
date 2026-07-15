rule all:
    input:
        f"{config['out_dir']}/zUMIS_output"


rule qc:
    input:
        r1=config["sequence_file"]["file1"],
        r2=config["sequence_file"]["file2"]
    output:
        r1=temp(f"{config['out_dir']}/filtered_R1.fastq.gz"),
        r2=temp(f"{config['out_dir']}/filtered_R2.fastq.gz")
    params:
        Adapter=config["advanced"]["adapter"],
        mismatch=config["advanced"]["adapter_mismatch"]
    threads:
        config["threads"]
    run:
        from most.python.rna.qc_adapt import qc_adapt
        qc_adapt(
            input.r1,
            input.r2,
            output.r1,
            output.r2,
            threads,
            params.Adapter,
            params.mismatch
        )


rule run_zumis:
    input:
        zumis=config["tools"]["zumis"],
        zcfg=f"{config['out_dir']}/filled.yaml",
        r1=f"{config['out_dir']}/filtered_R1.fastq.gz",
        r2=f"{config['out_dir']}/filtered_R2.fastq.gz"
    output:
        directory(f"{config['out_dir']}/zUMIS_output")
    shell:
        """
        {input.zumis} \
        -c\
        -y {input.zcfg}
        """