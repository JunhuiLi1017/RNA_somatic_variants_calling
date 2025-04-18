rule map_raw_sort:
    message:
        "use sambamba instead of samtools sort due to the large mem required by samtools"
    input:
        unpack(get_unsort_bam)
    output:
        sort_bam=temp("{outpath}/02_map/bwa/{sample}/{sample}.sort.bam"),
        sort_stat="{outpath}/02_map/bwa/{sample}/{sample}.sort.stat",
        sort_index="{outpath}/02_map/bwa/{sample}/{sample}.sort.bam.bai"
    log:
        "{outpath}/02_map/logs/{sample}.raw.sort.log"
    params:
        tmpdir="{outpath}/02_map/bwa/{sample}/tmpdir_{sample}",
        command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 1000) // threads
    threads:
        resource['resource']['medium']['threads']
    resources:
        mem_mb=resource['resource']['medium']['mem_mb']
    conda:
        "../envs/sambamba1.0.1.yaml"
    shell:
        """
        mkdir -p {params.tmpdir} && \
        sambamba sort \
            -t {threads} \
            -m {params.command_mem}M \
            -o {output.sort_bam} \
            --tmpdir {params.tmpdir} \
            {input.unsort_bam} > {log} 2>&1
        samtools stats {output.sort_bam} > {output.sort_stat}
        samtools index {output.sort_bam}
        """

rule mark_dup:
    input:
        bam="{outpath}/02_map/bwa/{sample}/{sample}.sort.bam",
        bai="{outpath}/02_map/bwa/{sample}/{sample}.sort.bam.bai"
    output:
        bam=temp("{outpath}/02_map/dup/{sample}/{sample}.sort.markdup.bam"),
        bai="{outpath}/02_map/dup/{sample}/{sample}.sort.markdup.bai",
        mtx="{outpath}/02_map/dup/{sample}/{sample}.sort.markdup.matrix"
    log:
        "{outpath}/02_map/logs/{sample}.sort.markdup.log"
    threads:
        resource['resource']['high']['threads']
    resources:
        mem_mb=resource['resource']['high']['mem_mb']
    conda:
        "../envs/gatk4.6.1.0.yaml"
    params:
        tmpdir="{outpath}/02_map/dup/{sample}/tmpdir_{sample}",
        gatk=config['gatk_current_using'],
        command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 4000)
    shell:
        """
        source ~/anaconda3/etc/profile.d/conda.sh; conda activate gatk4.6.1.0
        mkdir -p {params.tmpdir}
        java -Xms{params.command_mem}m -XX:ParallelGCThreads={threads} -jar {params.gatk} \
        MarkDuplicates \
        --TMP_DIR {params.tmpdir} \
        --INPUT {input.bam} \
        --METRICS_FILE {output.mtx} \
        --OUTPUT {output.bam} \
        --CREATE_INDEX true > {log} 2>&1
        conda deactivate
        """

rule split_NCigarReads:
    input:
        bam="{outpath}/02_map/dup/{sample}/{sample}.sort.markdup.bam",
        bai="{outpath}/02_map/dup/{sample}/{sample}.sort.markdup.bai"
    output:
        bam=temp("{outpath}/02_map/splitNCigar/{sample}/{sample}.sort.markdup.splitNCigar.bam"),
        bai="{outpath}/02_map/splitNCigar/{sample}/{sample}.sort.markdup.splitNCigar.bai"
    log:
        "{outpath}/02_map/logs/{sample}.sort.markdup.splitNCigar.log"
    params:
        gatk=config['gatk_current_using'],
        ref=config['reference'],
        command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 4000)
    threads:
        resource['resource']['very_high']['threads']
    resources:
        mem_mb=resource['resource']['very_high']['mem_mb']
    conda:
        "../envs/gatk4.6.1.0.yaml"
    shell:
        """
        source ~/anaconda3/etc/profile.d/conda.sh; conda activate gatk4.6.1.0
        java -Xms{params.command_mem}m -XX:ParallelGCThreads={threads} -jar {params.gatk} \
        SplitNCigarReads \
        -R {params.ref} \
        -I {input.bam} \
        -O {output.bam} > {log} 2>&1
        conda deactivate
        """

rule split_NCigarReads_sort:
    input:
        bam="{outpath}/02_map/splitNCigar/{sample}/{sample}.sort.markdup.splitNCigar.bam",
        bai="{outpath}/02_map/splitNCigar/{sample}/{sample}.sort.markdup.splitNCigar.bai"
    output:
        sort_bam="{outpath}/02_map/splitNCigar_sort/{sample}/{sample}.markdup.splitNCigar.sort.bam",
        sort_bai="{outpath}/02_map/splitNCigar_sort/{sample}/{sample}.markdup.splitNCigar.sort.bam.bai",
        sort_stat="{outpath}/02_map/splitNCigar_sort/{sample}/{sample}.markdup.splitNCigar.sort.bam.stats.txt"
    log:
        "{outpath}/02_map/logs/{sample}.markdup.splitNCigar.sort.log"
    params:
        tmpdir="{outpath}/02_map/splitNCigar_sort/{sample}/tmpdir_{sample}",
        command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 1000) // threads
    threads:
        resource['resource']['medium']['threads']
    resources:
        mem_mb=resource['resource']['medium']['mem_mb']
    conda:
        "../envs/sambamba1.0.1.yaml"
    shell:
        """
        mkdir -p {params.tmpdir} && \
        sambamba sort \
            -t {threads} \
            -m {params.command_mem}M \
            -o {output.sort_bam} \
            --tmpdir {params.tmpdir} \
            {input.bam} > {log} 2>&1
        samtools stats {output.sort_bam} > {output.sort_stat}
        samtools index {output.sort_bam}
        """

rule base_reca_librator:
    input:
        bam="{outpath}/02_map/splitNCigar/{sample}/{sample}.sort.markdup.splitNCigar.bam",
        bai="{outpath}/02_map/splitNCigar/{sample}/{sample}.sort.markdup.splitNCigar.bai"
    output:
        table="{outpath}/02_map/bqsr/{sample}/{sample}.recal_data.table"
    log:
        "{outpath}/02_map/logs/{sample}.recal_data.log"
    params:
        gatk=config['gatk_current_using'],
        ref=config['reference'],
        dpsnp138=config['dpsnp138'],
        known_indels=config['known_indels'],
        Mills_and_1000G=config['Mills_and_1000G'],
        command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000)
    threads:
        resource['resource']['very_high']['threads']
    resources:
        mem_mb=resource['resource']['very_high']['mem_mb']
    conda:
        "../envs/gatk4.6.1.0.yaml"
    shell:
        """
        source ~/anaconda3/etc/profile.d/conda.sh; conda activate gatk4.6.1.0
        java -Xms{params.command_mem}m -XX:ParallelGCThreads={threads} \
        -jar {params.gatk} BaseRecalibrator \
         -I {input.bam} \
         -O {output.table} \
         -R {params.ref} \
         --known-sites {params.dpsnp138} \
         --known-sites {params.known_indels} \
         --known-sites {params.Mills_and_1000G} > {log} 2>&1
         conda deactivate
        """

rule apply_bqsr:
    input:
        bam="{outpath}/02_map/splitNCigar/{sample}/{sample}.sort.markdup.splitNCigar.bam",
        table="{outpath}/02_map/bqsr/{sample}/{sample}.recal_data.table"
    output:
        bam="{outpath}/02_map/bqsr/{sample}/{sample}.sort.markdup.splitNCigar.bqsr.bam",
        bai="{outpath}/02_map/bqsr/{sample}/{sample}.sort.markdup.splitNCigar.bqsr.bai",
        flag="{outpath}/02_map/bqsr/{sample}/{sample}.sort.markdup.splitNCigar.bqsr.stat.txt"
    log:
        "{outpath}/02_map/logs/{sample}.applyBQSR.log"
    params:
        ref=config['reference'],
        gatk=config['gatk_current_using'],
        prefix="{outpath}/02_map/bqsr/{sample}/{sample}",
        command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000)
    threads:
        resource['resource']['very_high']['threads']
    resources:
        mem_mb=resource['resource']['very_high']['mem_mb']
    conda:
        "../envs/gatk4.6.1.0.yaml"
    shell:
        """
        source ~/anaconda3/etc/profile.d/conda.sh; conda activate gatk4.6.1.0
        java -Xms{params.command_mem}m -XX:ParallelGCThreads={threads} \
        -jar {params.gatk} ApplyBQSR \
        -I {input.bam}\
        -O {output.bam} \
        -R {params.ref} \
        --use-original-qualities \
        --bqsr-recal-file {input.table} > {log} 2>&1
        samtools flagstat {output.bam} > {output.flag}
        conda deactivate
        """
