rule Mutect2:
	input:
		bam="{outpath}/02_map/bqsr/{sample}/{sample}.sort.markdup.splitNCigar.bqsr.bam"
	output:
		raw_vcf="{outpath}/03_variants/02_mutect2/sub/{sample}/{sample}.{chrno}.mt2pon.vcf.gz",
		raw_stat="{outpath}/03_variants/02_mutect2/sub/{sample}/{sample}.{chrno}.mt2pon.vcf.gz.stats"
	log:
		"{outpath}/03_variants/logs/{sample}.{chrno}.mutect2.log"
	params:
		pon=config['pon'],
		af_only_gnomad=config['af_only_gnomad'],
		ref=config['reference'],
		sample="{sample}",
		interval_list=intervals_dir + "/{chrno}.intervals.list",
		dbsnp138=config['dpsnp138'],
		gatk=config['gatk_current_using'],
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000)
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	conda:
		"../envs/gatk4.6.1.0.yaml"
	shell:
		'''
		source ~/anaconda3/etc/profile.d/conda.sh; conda activate gatk4.6.1.0
		java -Xms{params.command_mem}m -XX:ParallelGCThreads={threads} \
		-jar {params.gatk} Mutect2 \
		-R {params.ref} \
		-I {input.bam} \
		--pon {params.pon} \
		-tumor {params.sample} \
		--germline-resource {params.af_only_gnomad} \
		-L {params.interval_list} \
		--interval-padding 100 \
		-O {output.raw_vcf} > {log} 2>&1
		conda deactivate
		'''

rule vcf_merge:
	input:
		lambda wildcards: [f"{outpath}/03_variants/02_mutect2/sub/{wildcards.sample}/{wildcards.sample}.{chr}.mt2pon.vcf.gz" for chr in chromosomes]
	output:
		args="{outpath}/03_variants/02_mutect2/{sample}/{sample}.args",
		vcf="{outpath}/03_variants/02_mutect2/{sample}/{sample}.mt2pon.vcf.gz",
		tbi="{outpath}/03_variants/02_mutect2/{sample}/{sample}.mt2pon.vcf.gz.tbi"
	log:
		"{outpath}/03_variants/logs/{sample}.args.logs"
	shell:
		'''
		ls {input} > {output.args}
		vcf-concat -f {output.args} | bgzip > {output.vcf}
		tabix -p vcf {output.vcf}
		'''

rule merge_mutectstats:
	input:
		lambda wildcards: [f"{wildcards.outpath}/03_variants/02_mutect2/sub/{wildcards.sample}/{wildcards.sample}.{chr}.mt2pon.vcf.gz.stats" for chr in chromosomes]
	output:
		"{outpath}/03_variants/02_mutect2/{sample}/{sample}.mt2pon.vcf.gz.stats"
	params:
		list_para = lambda wildcards: ' '.join([f"--stats {wildcards.outpath}/03_variants/02_mutect2/sub/{wildcards.sample}/{wildcards.sample}.{chr}.mt2pon.vcf.gz.stats" for chr in chromosomes]),
		gatk=config['gatk_current_using'],
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 1000)
	log:
		"{outpath}/03_variants/logs/{sample}.MergeMutectStats.log"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	conda:
		"../envs/gatk4.6.1.0.yaml"
	shell:
		'''
		source ~/anaconda3/etc/profile.d/conda.sh; conda activate gatk4.6.1.0
		java -Xms{params.command_mem}m -XX:ParallelGCThreads={threads} \
		-jar {params.gatk} MergeMutectStats \
		{params.list_para} \
		-O {output} > {log} 2>&1
		conda deactivate
		'''

rule FilterMutectCall:
	input:
		vcf="{outpath}/03_variants/02_mutect2/{sample}/{sample}.mt2pon.vcf.gz",
		stats="{outpath}/03_variants/02_mutect2/{sample}/{sample}.mt2pon.vcf.gz.stats",
		index="{outpath}/03_variants/02_mutect2/{sample}/{sample}.mt2pon.vcf.gz.tbi"
	output:
		"{outpath}/03_variants/02_mutect2/{sample}/{sample}.mt2pon.filter.vcf.gz"
	log:
		"{outpath}/03_variants/logs/{sample}.filter.log"
	params:
		ref=config['reference'],
		gatk=config['gatk_current_using'],
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000)
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	conda:
		"../envs/gatk4.6.1.0.yaml"
	shell:
		'''
		source ~/anaconda3/etc/profile.d/conda.sh; conda activate gatk4.6.1.0
		java -Xms{params.command_mem}m -XX:ParallelGCThreads={threads} \
		-jar {params.gatk} FilterMutectCalls \
		-R {params.ref} \
		-V {input.vcf} \
		-O {output} > {log} 2>&1
		conda deactivate
		'''

rule anno_mutect2:
	input:
		vcf_gz="{outpath}/03_variants/02_mutect2/{sample}/{sample}.mt2pon.filter.vcf.gz"
	output:
		input_anno="{outpath}/03_variants/02_mutect2/{sample}/{sample}.input.anno.{ref_version}.txt",
		vcf_anno="{outpath}/03_variants/02_mutect2/{sample}/{sample}.anno.{ref_version}_multianno.txt"
	params:
		ref_version=config["ref_version"],
		data_anno="/pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/03Database/database_fromHPCC/dataset/ANNOVAR/annovar",
		annovar="/pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/05softwares/ANNOVAR/annovar/table_annovar.pl",
		outdir="{outpath}/03_variants/02_mutect2/{sample}/",
		outputanno="{outpath}/03_variants/02_mutect2/{sample}/{sample}.anno"
	log:
		"{outpath}/logs/anno_mutect2/{sample}.{ref_version}.log"
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	shell:
		'''
		zcat {input.vcf_gz} | grep -v "#" | awk '{{if($7=="PASS"){{print $0}}}}' | awk '{{ print $1"\\t"$2"\\t"$2+length($4)-1"\\t"$4"\\t"$5 }}' > {output.input_anno}
		perl {params.annovar} \
		{output.input_anno} \
		{params.data_anno}/humandb_{ref_version} \
		-buildver {ref_version} \
		-out {params.outputanno} \
		-remove \
		-protocol refGene,dbnsfp42a \
		-operation g,f \
		-nastring . 
		'''

