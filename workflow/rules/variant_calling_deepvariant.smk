rule run_deepvariant:
	input:
		bam="{outpath}/02_map/splitNCigar_sort/{sample}/{sample}.markdup.splitNCigar.sort.bam"
	output:
		vcf="{outpath}/03_variants/03_deepvariant/{sample}/{sample}.{chr}.{ref_version}.output.vcf.gz"
	log:
		"{outpath}/03_variants/logs/{sample}.{chr}.{ref_version}.deepvariant.log"
	params:
		model_type="WES",
		extra_args="split_skip_reads=true,channels=''",
		ref=config['reference'],
		interval_list=intervals_dir + "/{chr}.intervals.list",
		intermediate_dir="{outpath}/03_variants/03_deepvariant/{sample}/{sample}.{chr}.{ref_version}.intermediate_results_dir",
		log_dir="{outpath}/03_variants/logs/deepvariant/{chr}",
		chr="{chr}",
		model="/pi/michelle.kelliher-umw/Junhui/RNAseq_Bella_mutant_calling_03_26_2025/03_database/deepvariant/model/model.ckpt"
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	container:
		"library://junhuili/deepvariant/deepvariant:1.8.0"
	shell:
		"""
		run_deepvariant \
		  --model_type={params.model_type} \
		  --customized_model={params.model} \
		  --ref={params.ref} \
		  --reads={input.bam} \
		  --output_vcf={output.vcf} \
		  --num_shards={threads} \
		  --logging_dir={params.log_dir} \
		  --make_examples_extra_args={params.extra_args} \
		  --regions {params.chr} \
		  --intermediate_results_dir={params.intermediate_dir} > {log} 2>&1
		"""