rule haplotypecaller:
	input:
		bam="{outpath}/02_map/05_apply_bqsr/{sample}.bqsr.bam"
	output:
		vcf="{outpath}/03_variants/haplotypecaller/01_raw_sub/{sample}.{chr}.vcf.gz",
		tbi="{outpath}/03_variants/haplotypecaller/01_raw_sub/{sample}.{chr}.vcf.gz.tbi"
	log:
		"{outpath}/03_variants/logs/{sample}.{chr}.haplotypecaller.log"
	params:
		ref=config['reference'],
		dbsnp138=config['dbsnp138'],
		interval=intervals_dir + "/{chr}.intervals.list"
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
		HaplotypeCaller \
		-O {output.vcf} \
		-R {params.ref} \
		-I {input.bam}\
		-L {params.interval} \
		--dont-use-soft-clipped-bases \
		--dbsnp {params.dbsnp138} > {log} 2>&1
		"""

rule mergevcfs:
	input:
		vcf=lambda wildcards: expand(
			"{outpath}/03_variants/haplotypecaller/01_raw_sub/{sample}.{chr}.vcf.gz", 
			outpath=wildcards.outpath, 
			sample=wildcards.sample, 
			chr=CHROMOSOMES
		),
		tbi=lambda wildcards: expand(
			"{outpath}/03_variants/haplotypecaller/01_raw_sub/{sample}.{chr}.vcf.gz.tbi", 
			outpath=wildcards.outpath, 
			sample=wildcards.sample, 
			chr=CHROMOSOMES
		)
	output:
		vcf_gz="{outpath}/03_variants/haplotypecaller/02_merge/{sample}.vcf.gz",
		vcf_tbi="{outpath}/03_variants/haplotypecaller/02_merge/{sample}.vcf.gz.tbi"
	params:
		gatk=config['gatk_current_using'],
		input_args=lambda wildcards, input: " ".join(f"--INPUT {vcf}" for vcf in input.vcf),
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 1000)
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	container:
		config["gatk_4.1.8.1"]
	shell:
		'''
		gatk --java-options "-Xms{params.command_mem}m -XX:ParallelGCThreads={threads}" \
		MergeVcfs \
		{params.input_args} \
		--OUTPUT {output.vcf_gz}
		'''

rule haplotypecaller_anno:
	input:
		vcf_gz="{outpath}/03_variants/haplotypecaller/02_merge/{sample}.vcf.gz",
		vcf_tbi="{outpath}/03_variants/haplotypecaller/02_merge/{sample}.vcf.gz.tbi"
	output:
		vcf_anno="{outpath}/03_variants/haplotypecaller/03_annotation/{sample}.anno.{ref_version}_multianno.txt"
	params:
		ref_version=config['ref_version'],
		annovar_dir=config['annovar_dir'],
		outputanno="{outpath}/03_variants/haplotypecaller/03_annotation/{sample}.anno"
	log:
		"{outpath}/03_variants/logs/{sample}.{ref_version}.log"
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	container:
		config["terra_perl_anno"]
	shell:
		"""
		perl {params.annovar_dir}/table_annovar.pl \
		{input.vcf_gz} \
		{params.annovar_dir}/humandb_{params.ref_version} \
		-buildver {params.ref_version} \
		-out {params.outputanno} \
		-remove \
		-protocol refGene,dbnsfp42a,clinvar_20240917,gnomad41_genome,gnomad41_exome \
		-operation g,f,f,f,f \
		-nastring . \
		-vcfinput
		"""

rule haplotypecaller_summary:
	input:
		expand("{outpath}/03_variants/haplotypecaller/03_annotation/{sample}.anno.{ref_version}_multianno.txt",
		outpath=Outpath,
		sample=Sample,
		ref_version=Ref_version)
	output:
		"{outpath}/03_variants/haplotypecaller/summary.txt"
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	shell:
		"""
		echo "haplotypecaller step is finished" > {output}
		"""