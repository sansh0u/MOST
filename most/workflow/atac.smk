rule all:
    input:
        f"{config['out_dir']}/{config['project']}_sorted.bed.gz",
        f"{config['out_dir']}/{config['project']}_sorted.bed.gz.tbi"


rule bbduk_linker1:
    input:
        r1=config["sequence_file"]["file1"],
        r2=config["sequence_file"]["file2"]
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


rule atac_bc:
    input:
        r1=f"{config['out_dir']}/linker2_R2.fastq.gz"
    output:
        r1=temp(f"{config['out_dir']}/output_R1.fastq"),
        r2=temp(f"{config['out_dir']}/output_R2.fastq")
    params:
        seq_start=config['runtime']['seq_start']
        bc2_start=config['runtime']['bc2_start']   
        bc2_end=config['runtime']['bc2_end'] 
        bc1_start=config['runtime']['bc1_start']
        bc1_end=config['runtime']['bc1_end']
    run:
        from most.preprocess.bc_process import atac_bc
        atac_bc(
            input.r1,
            input.r2,
            output.r1,
            output.r2,
            params.seq_start,
            params.bc2_start,
            params.bc2_end,
            params.bc1_start,
            params.bc1_end
        )


rule gzip_file:
    input:
        r1=f"{config['out_dir']}/output_R1.fastq",
        r2=f"{config['out_dir']}/output_R2.fastq"
    output:
        r1=temp(f"{config['out_dir']}/output_R1.fastq.gz"),
        r2=temp(f"{config['out_dir']}/output_R2.fastq.gz")
    threads:
        config["threads"]
    shell:
        """
        pigz -p {threads} {input.r1}
        pigz -p {threads} {input.r2}
        """


rule chromap:
    input:
        r1=f"{config['out_dir']}/output_R1.fastq.gz",
        r2=f"{config['out_dir']}/linker2_R1.fastq.gz",
        b_file=f"{config['out_dir']}/output_R2.fastq.gz",
        index_file=config['reference']['chromap_index'],
        fa_file=config['reference']['fa_file'],
        bc_file=config['reference']['barcode_file']
    output:
        bed=temp(f"{config['out_dir']}/{config['project']}.bed")
    threads:
        config["threads"]
    shell:
        """
        chromap \
        --preset atac \
        -x {input.index_file} \
        -r {input.fa_file} \
        -1 {input.r1} \
        -2 {input.r2} \
        -b {input.b_file} \
        --barcode-whitelist {input.bc_file} \
        -t {threads} \
        -o {output.bed}
        """

rule sort_bed:
    input:
        bed=f"{config['out_dir']}/{config['project']}.bed"
    output:
        bed=f"{config['out_dir']}/{config['project']}_sorted.bed"
    threads:
        config["threads"]
    shell:
        """
        sort \
        -k1,1 \
        -k2,2n \
        -k3,3n \
        -k4,4 \
        --parallel {threads} \
        -S 36G \
        {input.bed} \
        > {output.bed}
        """


rule bgzip_bed:
    input:
        bed=f"{config['out_dir']}/{config['project']}_sorted.bed"
    output:
        gz=f"{config['out_dir']}/{config['project']}_sorted.bed.gz"
    threads:
        config["threads"]
    shell:
        """
        bgzip \
        -@ {threads} \
        {input.bed}
        """


rule tabix_bed:
    input:
        gz=f"{config['out_dir']}/{config['project']}_sorted.bed.gz"
    output:
        tbi=f"{config['out_dir']}/{config['project']}_sorted.bed.gz.tbi"
    shell:
        """
        tabix \
        -p bed \
        {input.gz}
        """
