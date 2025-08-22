rule splice_ai:
	input:
		vcf_gz="{outpath}/03_variants/haplotypecaller/02_merge/{sample}.vcf.gz"
	output:
		vcf="{outpath}/05_splice_variant/spliceai/01_raw/{sample}.vcf.gz"
	log:
		"{outpath}/05_splice_variant/logs/{sample}.spliceai.log"
	params:
		ref=config['reference'],
		anno="grch38",
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000)
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	container:
		container_image['spliceai_1.3.1']
	shell:
		"""
		spliceai -I {input.vcf_gz} -O {output.vcf} -R {params.ref} -A {params.anno} > {log} 2>&1
		"""

rule spliceai_anno:
	input:
		vcf_gz="{outpath}/05_splice_variant/spliceai/01_raw/{sample}.vcf.gz"
	output:
		vcf_anno="{outpath}/05_splice_variant/spliceai/02_annotation/{sample}.anno.{ref_version}_multianno.txt"
	params:
		ref_version=config['ref_version'],
		annovar_dir=config['annovar_dir'],
		outputanno="{outpath}/05_splice_variant/spliceai/02_annotation/{sample}.anno"
	log:
		"{outpath}/05_splice_variant/logs/{sample}.{ref_version}.log"
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	container:
		container_image["terra_perl_anno"]
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

rule summary_spliceai:
	input:
		expand("{outpath}/05_splice_variant/spliceai/02_annotation/{sample}.anno.{ref_version}_multianno.txt",
		outpath=Outpath,
		sample=Sample,
		ref_version=Ref_version
		)
	output:
		"{outpath}/05_splice_variant/summary.txt"
	log:
		"{outpath}/05_splice_variant/logs/summary.log"
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	shell:
		"""
		echo "splice_variant is finished" > {output}
		"""