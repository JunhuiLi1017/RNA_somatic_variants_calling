rule map_raw:
    input:
        get_fastqc_fastp
    output:
        raw_bam="{outpath}/02_map/01_raw/{sample}.{library}.{flowlane}.Aligned.sortedByCoord.out.bam"
    log:
        "{outpath}/02_map/logs/{sample}.raw.log"
    params:
        ref=config['reference_dir'],
        prefix="{sample}.{library}.{flowlane}",
        temp_dir="{outpath}/02_map/01_raw/",
        command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 1000) // threads
    threads:
        resource['resource']['high']['threads']
    resources:
        mem_mb=resource['resource']['high']['mem_mb']
    container:
        config["STAR"]
    shell:
        """
        mkdir -p {params.temp_dir}
        cd {params.temp_dir}
        STAR --runThreadN {threads} \
        --runMode alignReads \
        --genomeDir {params.ref} \
        --readFilesIn {input} \
        --outFileNamePrefix {params.prefix} \
        --outSAMtype BAM SortedByCoordinate
        """

rule map_raw_index:
    input:
        "{outpath}/02_map/01_raw/{sample}.{library}.{flowlane}.Aligned.sortedByCoord.out.bam"
    output:
        "{outpath}/02_map/01_raw/{sample}.{library}.{flowlane}.Aligned.sortedByCoord.out.bam.bai"
    log:
        "{outpath}/02_map/logs/{sample}.raw.index.log"
    threads:
        resource['resource']['medium']['threads']
    resources:
        mem_mb=resource['resource']['medium']['mem_mb']
    container:
        config["samtools"]
    shell:
        """
        samtools index {input}
        """

rule mark_dup:
    message:
        "Note: we assumed the RNA sample is pcr-based sample"
    input:
        lambda wildcards: [
            f"{wildcards.outpath}/02_map/01_raw/{wildcards.sample}.{u.library}.{u.flowlane}.Aligned.sortedByCoord.out.bam"
            for u in units[units["sample"] == wildcards.sample].itertuples()
        ]
    output:
        bam=temp("{outpath}/02_map/02_dup/{sample}.markdup.bam"),
        bai="{outpath}/02_map/02_dup/{sample}.markdup.bai",
        mtx="{outpath}/02_map/02_dup/{sample}.markdup.matrix"
    log:
        "{outpath}/02_map/logs/{sample}.markdup.log"
    threads:
        resource['resource']['high']['threads']
    resources:
        mem_mb=resource['resource']['high']['mem_mb']
    params:
        tmpdir="{outpath}/02_map/02_dup/{sample}/tmpdir_{sample}",
        input_args=lambda wildcards, input: " ".join(f"--INPUT {bam}" for bam in input),
        command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 4000)
	container:
        config["gatk_4.1.8.1"]
    shell:
        """
        mkdir -p {params.tmpdir}
        gatk --java-options "-Xms{params.command_mem}m -XX:ParallelGCThreads={threads}" \
        MarkDuplicates \
        {params.input_args} \
        --TMP_DIR {params.tmpdir} \
        --METRICS_FILE {output.mtx} \
        --OUTPUT {output.bam} \
        --CREATE_INDEX true > {log} 2>&1
        """

rule split_NCigarReads:
    input:
        bam="{outpath}/02_map/02_dup/{sample}.markdup.bam",
        bai="{outpath}/02_map/02_dup/{sample}.markdup.bai"
    output:
        bam="{outpath}/02_map/03_splitNCigar/{sample}.splitNCigar.bam",
        bai="{outpath}/02_map/03_splitNCigar/{sample}.splitNCigar.bai"
    log:
        "{outpath}/02_map/logs/{sample}.splitNCigar.log"
    params:
        gatk=config['gatk_current_using'],
        ref=config['reference'],
        command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 4000)
    threads:
        resource['resource']['very_high']['threads']
    resources:
        mem_mb=resource['resource']['very_high']['mem_mb']
    container:
        config["gatk_4.1.8.1"]
    shell:
        """
        gatk --java-options "-Xms{params.command_mem}m -XX:ParallelGCThreads={threads}" \
        SplitNCigarReads \
        -R {params.ref} \
        -I {input.bam} \
        --CREATE_INDEX true \
        -O {output.bam} > {log} 2>&1
        """

rule NCigarReads_sort:
    input:
        bam="{outpath}/02_map/03_splitNCigar/{sample}.splitNCigar.bam",
        bai="{outpath}/02_map/03_splitNCigar/{sample}.splitNCigar.bai"
    output:
        sort_bam="{outpath}/02_map/04_splitNCigar_sort/{sample}.markdup.splitNCigar.sort.bam",
        sort_bai="{outpath}/02_map/04_splitNCigar_sort/{sample}.markdup.splitNCigar.sort.bai"
    log:
        "{outpath}/02_map/logs/{sample}.splitNCigar.sort.log"
    params:
        tmpdir="{outpath}/02_map/04_splitNCigar_sort/{sample}/tmpdir_{sample}",
        command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 1000) // threads
    threads:
        resource['resource']['medium']['threads']
    resources:
        mem_mb=resource['resource']['medium']['mem_mb']
    container:
        config["sambamba_1.0.1"]
    shell:
        """
        mkdir -p {params.tmpdir} && \
        sambamba sort \
            -t {threads} \
            -m {params.command_mem}M \
            -o {output.sort_bam} \
            --tmpdir {params.tmpdir} \
            {input.bam} > {log} 2>&1
        sambamba index {output.sort_bam} {output.sort_bai}
        """

rule base_reca_librator:
    input:
        bam="{outpath}/02_map/03_splitNCigar/{sample}.splitNCigar.bam",
        bai="{outpath}/02_map/03_splitNCigar/{sample}.splitNCigar.bai"
    output:
        table="{outpath}/02_map/04_bqsr/{sample}.recal_data.table"
    log:
        "{outpath}/02_map/logs/{sample}.recal_data.log"
    params:
        gatk=config['gatk_current_using'],
        ref=config['reference'],
        dbsnp138=config['dbsnp138'],
        g1000_known_indels=config['g1000_known_indels'],
        mills_and_1000g=config['mills_and_1000g'],
        command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000)
    threads:
        resource['resource']['very_high']['threads']
    resources:
        mem_mb=resource['resource']['very_high']['mem_mb']
    container:
        config["gatk_4.1.8.1"]
    shell:
        """
        gatk --java-options "-Xms{params.command_mem}m -XX:ParallelGCThreads={threads}" \
        BaseRecalibrator \
        -I {input.bam} \
        -O {output} \
        -R {params.ref} \
        --known-sites {params.dbsnp138} \
        --known-sites {params.g1000_known_indels} \
        --known-sites {params.mills_and_1000g} > {log} 2>&1
        """

rule apply_bqsr:
    input:
        bam="{outpath}/02_map/03_splitNCigar/{sample}.splitNCigar.bam",
        table="{outpath}/02_map/04_bqsr/{sample}.recal_data.table"
    output:
        bam="{outpath}/02_map/05_apply_bqsr/{sample}.bqsr.bam",
        bai="{outpath}/02_map/05_apply_bqsr/{sample}.bqsr.bai"
    log:
        "{outpath}/02_map/logs/{sample}.applyBQSR.log"
    params:
        ref=config['reference'],
        gatk=config['gatk_current_using'],
        prefix="{outpath}/02_map/05_apply_bqsr/{sample}",
        command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000)
    threads:
        resource['resource']['very_high']['threads']
    resources:
        mem_mb=resource['resource']['very_high']['mem_mb']
    container:
        config["gatk_4.1.8.1"]
    shell:
        """
        gatk --java-options "-Xms{params.command_mem}m -XX:ParallelGCThreads={threads}" \
        ApplyBQSR \
        -I {input.bam}\
        -O {output.bam} \
        -R {params.ref} \
        --use-original-qualities \
        --create-output-bam-index true
        --bqsr-recal-file {input.table} > {log} 2>&1
        """

rule apply_bqsr_stat:
    input:
        bam="{outpath}/02_map/05_apply_bqsr/{sample}.bqsr.bam",
        bai="{outpath}/02_map/05_apply_bqsr/{sample}.bqsr.bai",
    output:
        flag="{outpath}/02_map/06_stats/{sample}.bqsr.flatstat.txt",
        idxstat="{outpath}/02_map/06_stats/{sample}.bqsr.idxstat",
        stat="{outpath}/02_map/06_stats/{sample}.bqsr.stat.txt"
    log:
        "{outpath}/02_map/logs/{sample}.flagstat.log"
    params:
        ref=config['reference'],
        gatk=config['gatk_current_using'],
        prefix="{outpath}/02_map/05_apply_bqsr/{sample}",
        command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000)
    threads:
        resource['resource']['high']['threads']
    resources:
        mem_mb=resource['resource']['high']['mem_mb']
    container:
        config["samtools_1.20"]
    shell:
        """
        samtools flagstat {output.bam} > {output.flag}
        samtools stats {output.bam} > {output.stat}
        samtools idxstat {output.bam} > {output.idxstat}
        """
        
rule multiqc_mapping:
	input:
		# Collect all statistics files for MultiQC
		samtools_stats=lambda wildcards: expand(
			"{outpath}/02_map/06_stats/{sample}.bqsr.stat.txt",
			outpath=wildcards.outpath,
			sample=Sample
		),
		samtools_idxstats=lambda wildcards: [
			f"{wildcards.outpath}/02_map/06_stats/{sample}.samtools_idxstats.txt"
			for sample in Sample
		],
		samtools_flagstats=lambda wildcards: expand(
			"{outpath}/02_map/06_stats/{sample}.bqsr.flatstat.txt",
			outpath=wildcards.outpath,
			sample=Sample
		)
	output:
		report_map="{outpath}/02_map/multiqc_report_map.html",
		data_dir="{outpath}/02_map/06_stats/multiqc_data"
	params:
		outdir="{outpath}/02_map",
		title="Mapping Quality Report"
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	container:
		config['multiqc_1.22.3']
	shell:
		"""
		multiqc \
			--outdir {params.outdir} \
			--title "{params.title}" \
			--filename multiqc_report_map.html \
			--data-dir {output.data_dir} \
			--force \
			{params.outdir}
		"""