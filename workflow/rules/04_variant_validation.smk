rule samtools_mpileup:
	input:
		bam="{outpath}/02_map/05_apply_bqsr/{sample}.bqsr.bam",
		bai="{outpath}/02_map/05_apply_bqsr/{sample}.bqsr.bai",
		site_bed=get_site_bed
	output:
		mpileup="{outpath}/04_validation/01_samtools_mpileup/{sample}.txt"
	log:
		log="{outpath}/04_validation/logs/{sample}.samtools_mpileup.log"
	params:
		ref=config['reference'],
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000)
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	container:
		config["samtools_1.20"]
	shell:
		"""
		samtools mpileup -l {input.site_bed} -f {params.ref} {input.bam} > {output.mpileup} 2> {log}
		"""

rule reads_count:
	input:
		"{outpath}/04_validation/01_samtools_mpileup/{sample}.txt"
	output:
		"{outpath}/04_validation/02_reads_count/{sample}.reads_count.txt"
	log:
		"{outpath}/04_validation/logs/{sample}.reads_count.log"
	params:
		site_bed=get_site_bed,
		reads_count="/pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/00script/00_pipeline/target_sequence_analysis/workflow/bin/reads_count_v1.0.py"
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	container:
		config["python_alpine3.21"]
	shell:
		"""
		python {params.reads_count} --input {input} --bed {params.site_bed} --output {output}
		"""

rule summary_validation:
	input:
		expand("{outpath}/04_validation/02_reads_count/{sample}.reads_count.txt",
		outpath=Outpath,
		sample=Sample
		)
	output:
		"{outpath}/04_validation/summary.txt"
	log:
		"{outpath}/04_validation/logs/summary.log"
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	shell:
		"""
		echo "validation is finished" > {output}
		"""