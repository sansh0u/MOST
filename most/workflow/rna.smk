rule all:
    input:
        (f"{config['out_dir']}/temp")


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
        from most.preprocess.qc_adapt import qc_adapt
        qc_adapt(
            input.r1,
            input.r2,
            output.r1,
            output.r2,
            threads,
            params.Adapter,
            params.mismatch
        )


rule bbduk_linker1:
    input:
        r1=f"{config['out_dir']}/filtered_R1.fastq.gz",
        r2=f"{config['out_dir']}/filtered_R2.fastq.gz"
    output:
        r1=temp(f"{config['out_dir']}/linker1_R1.fastq.gz"),
        r2=temp(f"{config['out_dir']}/linker1_R2.fastq.gz"),
        stats=(f"{config['out_dir']}/bbduk_stats_L1.txt")
    params:
        k=config["runtime"]["k1"],
        linker=config["advanced"]["linker1"],
        restrictleft=config["runtime"]["restrictleft1"]
    threads:
        config["threads"]
    shell:
        """
        bbduk.sh \
        in={input.r1} \
        in2={input.r2} \
        outm={output.r1} \
        outm2={output.r2} \
        hdist=3 \
        k={params.k} \
        literal={params.linker} \
        threads={threads} \
        mm=f \
        rcomp=f \
        skipr1=t \
        restrictleft={params.restrictleft} \
        stats={output.stats}
        """


rule bbduk_linker2:
    input:
        r1=f"{config['out_dir']}/linker1_R1.fastq.gz",
        r2=f"{config['out_dir']}/linker1_R2.fastq.gz"
    output:
        r1=temp(f"{config['out_dir']}/linker2_R1.fastq.gz"),
        r2=temp(f"{config['out_dir']}/linker2_R2.fastq.gz"),
        stats=(f"{config['out_dir']}/bbduk_stats_L2.txt")
    params:
        k=config["runtime"]["k2"],
        linker=config["advanced"]["linker2"],
        restrictleft=config["runtime"]["restrictleft2"]
    threads:
        config["threads"]
    shell:
        """
        bbduk.sh \
        in={input.r1} \
        in2={input.r2} \
        outm={output.r1} \
        outm2={output.r2} \
        hdist=3 \
        k={params.k} \
        literal={params.linker} \
        threads={threads} \
        mm=f \
        rcomp=f \
        skipr1=t \
        restrictleft={params.restrictleft} \
        stats={output.stats}
        """


rule dbit_qc:
    input:
        r1=f"{config['out_dir']}/linker2_R1.fastq.gz"
    output:
        r1=temp(f"{config['out_dir']}/output_R2.fastq.gz")
    params:
        umi_start=config["runtime"]["umi_start"],
        umi_len=config["runtime"]["umi_len"],
        bc2_start=config["runtime"]["bc2_start"],
        bc2_end=config["runtime"]["bc2_end"],
        bc1_start=config["runtime"]["bc1_start"],
        bc1_end=config["runtime"]["bc1_end"]
    run:
        from most.preprocess.bc_process import dbit_bc
        dbit_bc(
            input.r1,
            output.r1,
            params.umi_start,
            params.umi_len,
            params.bc2_start,
            params.bc2_end,
            params.bc1_start,
            params.bc1_end
        )

    
rule gzip_file:
    input:
        r1=f"{config['out_dir']}/output_R2.fastq"
    output:
        r1=temp(f"{config['out_dir']}/output_R2.fastq.gz")
    threads:
        config["threads"]
    shell:
        """
        pigz -p {threads} -f {input.r1}
        """


rule stpipeline:
    input:
        r1=f"{config['out_dir']}/output_R2.fastq.gz",
        r2=f"{config['out_dir']}/linker2_R1.fastq.gz",
        b_file=f"{config['out_dir']}/linker2_R2.fastq.gz",
        star_index=config["reference"]["star_index"],
        gtf_file=config["reference"]["gtf_file"],
        bc_file=config["reference"]["barcode_file"]
    output:
        r1={config['out_dir']},
        r2=(f"{config['out_dir']}/temp")
    params:
        stpipeline_id=config["project"],
        log_file=f"{config['out_dir']}/{config['project']}_log.txt"
    threads:
        config["threads"]
    shell:
        """
        st_pipeline_run \
        --output-folder {output.r1} \
        --temp-folder {output.r2} \
        --ids {input.bc_file} \
        --threads {threads} \
        --ref-map {input.star_index} \
        --ref-annotation {input.gtf_file} \
        --expName {params.stpipeline_id} \
        --log-file {params.log_file} \
        --htseq-no-ambiguous --demultiplexing-kmer 5 \
        --umi-start-position 16 \
        --umi-end-position 26 \
        --demultiplexing-overhang 0 \
        --min-length-qual-trimming 18 \
        --verbose \
        {input.r1} \
        {input.r2}
        """
    