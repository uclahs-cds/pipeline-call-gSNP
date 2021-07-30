# Regenotype gSNP

- [regenotype-gSNP](#regenotype-gsnp)
  - [Overview](#overview)
  - [How To Run](#how-to-run)
  - [Flow Diagram](#flow-diagram)
  - [Pipeline Steps](#pipeline-steps)
    - [1. GATK HaplotypeCaller](#1-gatk-haplotypecaller)
    - [2. GATK Split Intervals](#2-gatk-split-intervals)
    - [3. GATK GenomicsDBImport](#3-gatk-genomicsdbimport)
    - [4. GATK GenotypeGVCFs](#4-gatk-genotypegvcfs)
    - [5. GATK SortVcf](#5-gatk-sortvcf)
    - [6. GATK/Picard MergeVcfs](#6-gatkpicard-mergevcfs)
    - [7. GATK VariantRecalibrator](#7-gatk-variantrecalibrator)
    - [8. GATK ApplyVQSR](#8-gatk-applyvqsr)
    - [9. RTG vcfstats](#9-rtg-vcfstats)

  - [Inputs](#inputs)
  - [Outputs](#outputs)
  - [Testing and Validation](#testing-and-validation)
    - [Test Data Set](#test-data-set)
    - [Validation](#validation)
    - [Validation Tool](#validation-tool)
  - [References](#references)

## Overview

The regenotype gSNP pipeline performs joint genotyping on a cohort.
This allows population level information to be used to have greater sensitivity for low-frequency variants and filter out false positives.
This pipeline takes in single-sample BAMs or GVCFs (output from the [call-gSNP pipeline](https://github.com/uclahs-cds/pipeline-call-gSNP)) and then performs joint-genotyping on the entire cohort.

---

## How To Run

1. Update the params section of the .config file

2. Update the input csv

3. See the submission script, [here](https://github.com/uclahs-cds/tool-submit-nf), to submit your pipeline

---

## Flow Diagram

[Source: GATK Documentation](https://gatk.broadinstitute.org/hc/en-us/articles/360035890431-The-logic-of-joint-calling-for-germline-short-variants)

![Joint Calling Flow](https://drive.google.com/uc?id=15MxnVvtt8yR8nn7_wChzsZz_vqFHFmON)

[Source GATK Documentation](https://gatk.broadinstitute.org/hc/en-us/articles/360035890411)

![Joint Calling Flow](https://drive.google.com/uc?id=1EOTFKQPBLApJkAHx1GU6VImq_xeGac1w)

---

## Pipeline Steps

### 1. GATK HaplotypeCaller

Call variants individually per-sample. Use the HaplotypeCaller to call variants individually per-sample in `-ERC GVCF` mode.

### 2. GATK Split Intervals
Scatter the input intervals for parallelization.

### 3. GATK GenomicsDBImport
Combine all of the samples into one database, per interval.

### 4. GATK GenotypeGVCFs
Perform a joint genotyping analysis of the gVCFs produced for all samples in a cohort from the genomics database.

### 5. GATK SortVcf
Sort each of the VCF files before merging.

### 6. GATK/Picard MergeVcfs
Combine all of the scattered interval vcfs into one vcf file.

### 7. GATK VariantRecalibrator
Build the variant recalibration model for SNPs and INDELs separately.

### 8. GATK ApplyVQSR
Apply the variant recalibration model from step 5 to the GVCF file.

### 9. RTG vcfstats
Generate statistics about the content of the output VCF file(s) from step 7.

---

## Inputs

### Input CSV

 Field | Type | Description |
| ------------ | ------------ | ------------------------ |
| sample |  String | Sample ID |
| input_file | Path | Path to the BAM or GVCF file. Both BAMs and GVCFs can be passed as inputs. BAMs with be processed with HaplotypeCaller. |
| chr | String | (Optional) Omit if BAM/GVCF is over entire genome. Specify the single chromosome for the given BAM/GVCF file. |

### Config
| Parameter | Description |
| --- | --- |
| dataset_id | ID of the data set |
| blcds_registered_dataset | Boolean indicating if you want the output to be registered |
| sge_scheduler | Boolean indicating if we are on SGE. Use `false` if on slurm. |
| input_csv | Path to input.csv |
| intervals | A list of interval to run in parallel. Defaults to canonical chromosomes in [hg38_decoy_chromosomes_canonical.txt](pipeline/config/single-node/hg38_decoy_chromosomes_canonical.txt) |
| save_intermediate_files | true to save all intermediate files |
| reference_fasta |   path to genome genome.fa |
| reference_fasta_fai |  path to genome genome.fa.fai |
| reference_dict |  path to genome genome.dict |
| bundle_v0_dbsnp138_vcf_gz |  path to resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf.gz |
| bundle_mills_and_1000g_gold_standard_indels_vcf_gz |  path to Mills_and_1000G_gold_standard.indels.hg38.vcf.gz |
| bundle_known_indels_vcf_gz |  path to Homo_sapiens_assembly38.known_indels.vcf.gz |
| bundle_v0_dbsnp138_vcf_gz |  path to resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf.gz |
| bundle_hapmap_3p3_vcf_gz |  path to hapmap_3.3.hg38.vcf.gz |
| bundle_omni_1000g_2p5_vcf_gz |  path to 1000G_omni2.5.hg38.vcf.gz |
| bundle_phase1_1000g_snps_high_conf_vcf_gz |  path to 1000G_phase1.snps.high_confidence.hg38.vcf.gz |
| bundle_contest_hapmap_3p3_vcf_gz |  path to hapmap_3.3.hg38.BIALLELIC.PASS.vcf.gz |
| output_dir |  path to outputs  |
| output_log_directory |  path to logs |
| temp_dir |  path to temp, e.g /scratch |
| scatter_count | 50 |
| split_intervals_extra_args | Extra arguments to pass into SplitIntervals |
| max_num_intervals_to_import_in_parallel | Max number of intervals to import in parallel; higher values may improve performance, but require more memory and a higher number of file descriptors open at the same time. |
| reader_threads | How many simultaneous threads to use when opening VCFs in batches; higher values may improve performance when network latency is an issue |
| batch_size | Batch size controls the number of samples for which readers are open at once and therefore provides a way to minimize memory consumption. If batch_size > 100, use --consolidate True in extra args. Start with `batch_size = 50` for large datasets |
| genomics_db_import_extra_args | Extra arguments for GenomicsDBImport |

---

## Outputs

| Output | Description |
| ------------ | ------------------------ |
| regenotyped_merged_recalibrated_SNP.vcf.gz | Recalibrated germline SNPs |
| regenotyped_merged_recalibrated_SNP_AND_INDEL.vcf.gz | Recalibrate germline SNPs and INDELs |
| regenotyped_merged_recalibrated_SNP_vcfstats.txt | vcfstats for germline SNP vcf output |
| regenotyped_merged_recalibrated_SNP_AND_INDEL_vcfstats.txt | vcfstats for germline SNP and INDEL vcf output |

---

## Testing and Validation

### Test Data Sets

#### Soragni-Eng-cNF-MNFO (cNF-organoid) WGS
A small dataset used for testing (5 samples) is the cNF-organoid data from Alice Soragni's lab. The input.csv used for testing was
```
sample,input_file
SECNMNFO000001,/hot/users/stefaneng/project-cNF-organoid/pipeline/pipeline-call-gSNP/output/SECNMNFO000001/call-gSNP-20210220-190250/SECNMNFO000001/Picard-2.23.3/merged_raw_gvcf/SECNMNFO000001_merged_raw_variants.g.vcf.gz
SECNMNFO000002,/hot/users/stefaneng/project-cNF-organoid/pipeline/pipeline-call-gSNP/output/SECNMNFO000002/call-gSNP-20210220-200223/SECNMNFO000002/Picard-2.23.3/merged_raw_gvcf/SECNMNFO000002_merged_raw_variants.g.vcf.gz
SECNMNFO000003,/hot/users/stefaneng/project-cNF-organoid/pipeline/pipeline-call-gSNP/output/SECNMNFO000003/call-gSNP-20210220-190250/SECNMNFO000003/Picard-2.23.3/merged_raw_gvcf/SECNMNFO000003_merged_raw_variants.g.vcf.gz
SECNMNFO000004,/hot/users/stefaneng/project-cNF-organoid/pipeline/pipeline-call-gSNP/output/SECNMNFO000004/call-gSNP-20210220-190219/SECNMNFO000004/Picard-2.23.3/merged_raw_gvcf/SECNMNFO000004_merged_raw_variants.g.vcf.gz
SECNMNFO000005,/hot/users/stefaneng/project-cNF-organoid/pipeline/pipeline-call-gSNP/output/SECNMNFO000005/call-gSNP-20210220-190222/SECNMNFO000005/Picard-2.23.3/merged_raw_gvcf/SECNMNFO000005_merged_raw_variants.g.vcf.gz
```

This run is quick
```
Duration    : 48m 42s
CPU hours   : 6.0
```

|process_name                                          |realtime              |cpu    |peak_rss |peak_vmem |
|:-----------------------------------------------------|:---------------------|:------|:--------|:---------|
|GATK_merge_vcfs                                       |116s (~1.93 minutes)  |96.9%  |394.6 MB |16.1 GB   |
|GATK_run_GenomicsDBImport_import_gvcfs                |170s (~2.83 minutes)  |190.7% |2.5 GB   |17.1 GB   |
|GATK_run_GenotypeGVCFs                                |344s (~5.73 minutes)  |98.5%  |789.4 MB |6.6 GB    |
|GATK_sort_vcf                                         |18s                   |95.1%  |1.2 GB   |16.1 GB   |
|GATK_split_intervals                                  |5.7s                  |94.2%  |608 MB   |32.1 GB   |
|recalibrate_indels:GATK_run_ApplyVQSR                 |127s (~2.12 minutes)  |97.2%  |917.2 MB |32.1 GB   |
|recalibrate_indels:GATK_run_VariantRecalibrator_indel |202s (~3.37 minutes)  |97.4%  |1.8 GB   |35.2 GB   |
|recalibrate_snps:GATK_run_ApplyVQSR                   |146s (~2.43 minutes)  |97.2%  |948.7 MB |32.1 GB   |
|recalibrate_snps:GATK_run_VariantRecalibrator_snp     |1032s (~17.2 minutes) |97.9%  |3.6 GB   |35.3 GB   |

#### COVID19 WGS
The [COVID dataset](https://github.com/uclahs-cds/dataset-register-file/tree/master/COVID19) is the largest test so far for regenotype-gSNP with 702 `.g.vcf` files. The current settings in the config file provided were sufficient for this size of dataset but not guaranteed for other datasets. The run time was `6d 23h 23m 13s` (CPU hours `7'864.4`) on `M64` node with scratch space > 4T. The config, input csv, timeline.html, report.html and trace.txt are available on SLURM (10.18.82.18) `/hot/pipelines/development/slurm/regenotype-gSNP/benchmarks/COVID`

|process_name                                          |realtime              |cpu   |peak_rss |peak_vmem |
|:-----------------------------------------------------|:---------------------|:-----|:--------|:---------|
|GATK_merge_vcfs                                       |15159s (~4.21 hours)  |97.9% |978.3 MB |35.1 GB   |
|GATK_run_GenomicsDBImport_import_gvcfs                |117840s (~1.36 days)  |99.3% |14.3 GB  |33.4 GB   |
|GATK_run_GenotypeGVCFs                                |294203s (~3.41 days)  |98.3% |3.4 GB   |9.2 GB    |
|GATK_sort_vcf                                         |1182s (~19.7 minutes) |96.9% |33.7 GB  |35.3 GB   |
|GATK_split_intervals                                  |7.6s                  |91.1% |743.9 MB |32.1 GB   |
|recalibrate_indels:GATK_run_ApplyVQSR                 |10651s (~2.96 hours)  |97.2% |969.2 MB |32.1 GB   |
|recalibrate_indels:GATK_run_VariantRecalibrator_indel |6817s (~1.89 hours)   |98%   |3.8 GB   |34.9 GB   |
|recalibrate_snps:GATK_run_ApplyVQSR                   |11049s (~3.07 hours)  |97.7% |1 GB     |32.1 GB   |
|recalibrate_snps:GATK_run_VariantRecalibrator_snp     |10084s (~2.8 hours)   |92.9% |23.5 GB  |34.9 GB   |

#### Zlotta-Gebo-PRAD-RGK6 WGS

The output and logs are available on SLURM here: `/hot/projects/diseases/prostate-cancer/tag/zlotta/output_regenotype-gSNP/`
The run time was `1d 23h 12m 48s` with `1'421.9` CPU hours.

|process_name                                          |realtime               |cpu   |peak_rss |peak_vmem |
|:-----------------------------------------------------|:----------------------|:-----|:--------|:---------|
|GATK_merge_vcfs                                       |2931s (~48.85 minutes) |97%   |599.6 MB |35.1 GB   |
|GATK_run_GenomicsDBImport_import_gvcfs                |67845s (~18.85 hours)  |99.3% |16.1 GB  |31.4 GB   |
|GATK_run_GenotypeGVCFs                                |71702s (~19.92 hours)  |99.3% |2.8 GB   |8.8 GB    |
|GATK_sort_vcf                                         |321s (~5.35 minutes)   |96.7% |15.4 GB  |35.3 GB   |
|GATK_split_intervals                                  |7s                     |93.5% |619 MB   |32.1 GB   |
|recalibrate_indels:GATK_run_ApplyVQSR                 |2169s (~36.15 minutes) |96.4% |964.4 MB |32.1 GB   |
|recalibrate_indels:GATK_run_VariantRecalibrator_indel |1638s (~27.3 minutes)  |95.6% |2.7 GB   |34.9 GB   |
|recalibrate_snps:GATK_run_ApplyVQSR                   |2295s (~38.25 minutes) |96.6% |979.5 MB |32.1 GB   |
|recalibrate_snps:GATK_run_VariantRecalibrator_snp     |3582s (~59.7 minutes)  |98.3% |12.6 GB  |34.9 GB   |

### Validation

### Validation Tool

Included is a template for validating your input files. For more information on the tool check out: https://github.com/uclahs-cds/tool-validate-nf

---

## References
1. [Old Perl Pipeline](https://github.com/uclahs-cds/BoutrosLab/blob/master/CPC-GENE/GermlineSomatic/genome-wide/branches/GWAS/gwas_pipelines/regenotype_sig_snps.pl)
2. [The logic of joint calling for germline short variants](https://gatk.broadinstitute.org/hc/en-us/articles/360035890431-The-logic-of-joint-calling-for-germline-short-variants)
3. [Calling variants on cohorts of samples using the HaplotypeCaller in GVCF mode](https://gatk.broadinstitute.org/hc/en-us/articles/360035890411)
4. [GVCF - Genomic Variant Call Format](https://gatk.broadinstitute.org/hc/en-us/articles/360035531812)
