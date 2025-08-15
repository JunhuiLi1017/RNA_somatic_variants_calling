# RNA Somatic Variants Calling Pipeline

A comprehensive Snakemake-based pipeline for calling somatic variants from RNA sequencing data. This pipeline supports multiple variant callers and provides quality control, mapping, variant calling, and validation workflows.

## Overview

This pipeline performs end-to-end analysis of RNA sequencing data for somatic variant detection, including:

- **Quality Control**: FastQC, FastP, and MultiQC for comprehensive QC reports
- **Read Mapping**: Best practice of GATK
- **Variant Calling**: Multiple callers including HaplotypeCaller, DeepVariant, DeepSomatic, and Mutect2
- **Variant Validation**: Comprehensive validation and filtering of called variants
- **Annotation**: Functional annotation using ANNOVAR

## Pipeline Structure

```
workflow/
├── Snakefile                    # Main workflow file
├── Snakefile.common_snp         # Workflow for common SNP filtering
├── rules/                       # Individual workflow rules
│   ├── 00_common.smk           # Common functions and configurations
│   ├── 01_qc.smk               # Quality control rules
│   ├── 02_map_reads.smk        # Read mapping rules
│   ├── 03_variant_calling_haplotypecaller.smk  # HaplotypeCaller rules
│   ├── 04_variant_validation.smk # Variant validation rules
│   └── ...                     # Additional variant calling rules
├── envs/                        # Conda environment files
├── schemas/                     # Configuration schemas
└── submit.sh                    # Submission script
```

## Workflow Steps

### 1. Quality Control (`01_qc.smk`)
- **FastQC**: Raw read quality assessment
- **FastP**: Adapter trimming and quality filtering
- **MultiQC**: Aggregated QC reports

### 2. Read Mapping (`02_map_reads.smk`)
- **STAR**: Read alignment to reference genome
- **SAMtools**: BAM file processing and sorting
- **GATK BQSR**: Base quality score recalibration

### 3. Variant Calling (`03_variant_calling_*.smk`)
Supported variant callers:
- **HaplotypeCaller**: GATK's haplotype-based caller
- **DeepVariant**: Google's deep learning-based caller
- **DeepSomatic**: Somatic variant detection with DeepVariant
- **Mutect2**: GATK's somatic variant caller

### 4. Variant Validation (`04_variant_validation.smk`)
- reads count for the candidate sites
- Summary statistics generation

## Configuration

### Main Configuration Files

- `config/config_hg38.yaml`: Configuration for hg38 reference genome
- `config/config_hg19.yaml`: Configuration for hg19 reference genome
- `config/units.tsv`: Sample information and file paths
- `config/resources.yaml`: Computational resources allocation

### Required Input Files

The `units.tsv` file should contain the following columns:
- `sample`: Sample identifier
- `library`: Library identifier
- `flowlane`: Flow cell and lane information
- `platform`: Sequencing platform (e.g., ILLUMINA)
- `fq1`: Path to forward reads (R1)
- `fq2`: Path to reverse reads (R2) [optional for single-end]
- `trim_front1/2`: Number of bases to trim from front of reads
- `trim_tail1/2`: Number of bases to trim from end of reads
- `site_bed`: BED file defining target regions

### Key Configuration Parameters

- `ref_version`: Reference genome version (hg38/hg19)
- `reference`: Path to reference genome FASTA
- `interval`: Directory containing interval files for targeted sequencing
- `callers`: List of variant callers to use
- `outpath`: Output directory for all results

## Usage

### Prerequisites

1. **Snakemake**: Install Snakemake with conda
   ```bash
   conda install -c bioconda snakemake
   ```

2. **Docker/Singularity**: For containerized execution
   - Docker containers are specified in the config file
   - Singularity containers are also supported

3. **Reference Data**: Ensure all reference files are available:
   - Reference genome FASTA
   - Known variants (dbSNP, Mills, 1000G)
   - Interval files for targeted regions

### Running the Pipeline

1. **Configure the pipeline**:
   ```bash
   # Edit configuration files
   vim config/config_hg38.yaml
   vim config/units.tsv
   ```

2. **Run the pipeline**:
   ```bash
   # Using the provided submission script
   bash workflow/submit.sh
   
   # Or run directly with Snakemake
   snakemake -s workflow/Snakefile --configfile config/config_hg38.yaml -p -j 99
   ```

3. **For common SNP filtering workflow**:
   ```bash
   #bsub -W 72:00 -q long 'source ~/anaconda3/etc/profile.d/conda.sh; conda activate snakemake; bash submit.sh &> submit.log'
    source ~/anaconda3/etc/profile.d/conda.sh; conda activate snakemake
    snakemake -s Snakefile --configfile config/config.yaml -p -j 99 --latency-wait 500 --default-resources mem_mb=10000 disk_mb=10000 --cluster 'bsub -q long -o main.log -R "rusage[mem={resources.mem_mb}]" -n {threads} -R span[hosts=1] -W 72:00'
   ```

### Output Structure

```
{outpath}/
├── 01_multiqc/                 # Quality control reports
│   ├── fastqc/                 # FastQC results
│   │   ├── {sample}_{library}_{flowlane}.R1_fastqc.html
│   │   ├── {sample}_{library}_{flowlane}.R1_fastqc.zip
│   │   ├── {sample}_{library}_{flowlane}.R2_fastqc.html
│   │   └── {sample}_{library}_{flowlane}.R2_fastqc.zip
│   ├── fastp/                  # FastP results
│   │   ├── {sample}_{library}_{flowlane}.json
│   │   └── {sample}_{library}_{flowlane}.html
│   └── multiqc_report_fq.html  # Aggregated QC report
├── 02_map/                     # Mapping results
│   ├── 01_raw/                 # STAR alignment files
│   │   ├── {sample}.{library}.{flowlane}.Aligned.sortedByCoord.out.bam
│   │   └── {sample}.{library}.{flowlane}.Aligned.sortedByCoord.out.bam.bai
│   ├── 02_dup/                 # MarkDuplicates results
│   │   ├── {sample}.markdup.bam
│   │   ├── {sample}.markdup.bai
│   │   └── {sample}.markdup.matrix
│   ├── 03_splitNCigar/         # SplitNCigarReads results
│   │   ├── {sample}.splitNCigar.bam
│   │   └── {sample}.splitNCigar.bai
│   ├── 04_splitNCigar_sort/    # Sorted SplitNCigar results
│   │   ├── {sample}.markdup.splitNCigar.sort.bam
│   │   └── {sample}.markdup.splitNCigar.sort.bai
│   ├── 04_bqsr/                # Base recalibration tables
│   │   └── {sample}.recal_data.table
│   ├── 05_apply_bqsr/          # Final BQSR processed BAMs
│   │   ├── {sample}.bqsr.bam
│   │   └── {sample}.bqsr.bai
│   ├── 06_stats/               # Mapping statistics
│   │   ├── {sample}.bqsr.flatstat.txt
│   │   ├── {sample}.bqsr.idxstat
│   │   ├── {sample}.bqsr.stat.txt
│   │   └── multiqc_data/       # MultiQC data directory
│   └── multiqc_report_map.html # Mapping QC report
├── 03_variants/                # Variant calling results
│   └── haplotypecaller/        # HaplotypeCaller results
│       ├── 01_raw_sub/         # Per-chromosome VCFs
│       │   ├── {sample}.{chr}.vcf.gz
│       │   └── {sample}.{chr}.vcf.gz.tbi
│       ├── 02_merge/           # Merged VCFs
│       │   ├── {sample}.vcf.gz
│       │   └── {sample}.vcf.gz.tbi
│       ├── 03_annotation/      # ANNOVAR annotation
│       │   └── {sample}.anno.{ref_version}_multianno.txt
│       └── summary.txt         # HaplotypeCaller summary
└── 04_validation/              # Validation results
    ├── 01_samtools_mpileup/    # SAMtools mpileup results
    │   └── {sample}.txt
    ├── 02_reads_count/         # Read count analysis
    │   └── {sample}.reads_count.txt
    └── summary.txt             # Validation summary
```

## Supported Reference Genomes

- **hg38**: Human genome assembly GRCh38
- **hg19**: Human genome assembly GRCh37

## Variant Callers

### HaplotypeCaller
- GATK's haplotype-based germline variant caller
- Suitable for both germline and somatic variant detection
- Provides high sensitivity for complex variants

### DeepVariant
- Google's deep learning-based variant caller
- Excellent accuracy for SNVs and indels
- Requires GPU for optimal performance

### DeepSomatic
- Specialized for somatic variant detection
- Built on DeepVariant architecture
- Optimized for tumor-normal comparisons

### Mutect2
- GATK's somatic variant caller
- Designed specifically for somatic variant detection
- Includes extensive filtering capabilities

## Quality Control

The pipeline includes comprehensive quality control measures:

- **Read Quality**: FastQC assessment of raw reads
- **Adapter Contamination**: FastP adapter trimming
- **Mapping Quality**: BWA-MEM alignment statistics
- **Coverage Analysis**: Depth and uniformity assessment
- **Variant Quality**: Multiple quality metrics and filters

## Computational Resources

The pipeline is designed for high-performance computing environments:

- **Memory**: Configurable per rule (default: 10GB)
- **CPU**: Multi-threaded execution
- **Storage**: Efficient intermediate file management
- **Cluster Support**: LSF job submission included

## Citation
## Support

For issues and questions:
1. Check the log files in `workflow/logs/`
2. Review the Snakemake documentation
3. Ensure all dependencies and reference files are properly configured

## License

This pipeline is provided as-is for research purposes. Please ensure compliance with the licenses of all included tools and reference data.
