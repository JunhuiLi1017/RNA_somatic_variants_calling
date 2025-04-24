# Rule 7: Filter common SNPs using dbSNP, gnomAD, and 1000 Genomes
rule filter_common_snps:
	input:
		vcf="{outpath}/03_variants/01_haplotypecaller/{sample}/{sample}.vcf.gz",
		tbi="{outpath}/03_variants/01_haplotypecaller/{sample}/{sample}.vcf.gz.tbi"
	output:
		vcf="{outpath}/03_variants/01_haplotypecaller/04_sub_vcf_variant_commonSNP_filter/{sample}.common_filter.vcf.gz",
		tbi="{outpath}/03_variants/01_haplotypecaller/04_sub_vcf_variant_commonSNP_filter/{sample}.common_filter.vcf.gz.tbi"
	log:
		"{outpath}/03_variants/logs/{sample}.filter_indels.log"
	params:
		gnomad=config["af_only_gnomad"],
		pon=config["pon"]
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	conda:
		"../envs/gatk4.6.1.0.yaml"
	shell:
		"""
		source ~/anaconda3/etc/profile.d/conda.sh; conda activate bcftools	
		# Filter out common variants (dbSNP rsID, or present in 1000G PoN)
		bcftools view -e 'ID != "."' {input.vcf} | \
		bcftools view -T ^{params.pon} | bgzip > {output.vcf}
		tabix -p vcf {output.vcf}
		conda deactivate
		"""

rule annotate_clinvar:
	input:
		vcf="{outpath}/03_variants/01_haplotypecaller/04_sub_vcf_variant_commonSNP_filter/{sample}.common_filter.vcf.gz"
	output:
		txt="{outpath}/03_variants/01_haplotypecaller/06_annovar/{sample}.{ref_version}_multianno.txt"
	log:
		"{outpath}/03_variants/logs/{sample}.{ref_version}.06_annovar.log"
	params:
		ref_version=config['ref_version'],
		annovar_dir=config['annovar_dir'],
		outputanno="{outpath}/03_variants/01_haplotypecaller/06_annovar/{sample}",
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000)
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	shell:
		"""
		perl {params.annovar_dir}/table_annovar.pl \
		{input.vcf} \
		{params.annovar_dir}/humandb_{ref_version} \
		-buildver {ref_version} \
		-out {params.outputanno} \
		-remove \
		-protocol refGene,dbnsfp42a,clinvar_20240917,gnomad30_genome,collins_dosage \
		-operation g,f,f,f,f \
		-nastring . \
		-vcfinput
		"""

