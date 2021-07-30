# Changelog
All notable changes to the pipeline-regenotype-gSNP pipeline.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

## [1.3.0] - 2021-07-30
### Added
- Add process for RTG vcfstats on VQSRed vcf files
### Changed
- Update validate docker image to 2.1.5
- Standardize process names
- Config directory structure updated

## [1.2.0] - 2021-07-13
### Changed
- Copy log `.command.*` files into `params.log_output_dir`.
- `params.save_intermediate_files` parameter enables temporary files to be copied to the final output directory
- Adds scattered parallelization via GATK's SplitIntervals (taken from Caden's call-sSNV) #18
  - Required addition of Picard/GATK's SortVcf
- Add task memory for all GATK processes #18
  - Set java heap memory to `task.memory - params.gatk_command_mem_diff` as in call-sSNV #16
- Upgrade GATK to 4.2.0.0 #14
- Add Genomic DB arguments #7 #18
  - `--batch-size` (to decrease memory usage)
  - `--reader-threads` (for faster network IO)
  - `--tmp-dir` to help resolve memory issues
  - `--max-num-intervals-to-import-in-parallel` to improve import speed
- Replace bcftools concatenate + tabix generate index with Picard/GATK's MergeVcfs
- Replace java tmp directory with `/scratch` in all remaining places #2

---

## [1.1.0] - 2021-04-29
  - Auto detect bam and gvcf files.
  - Fix intermediate files log
  - Apply INDEL VQSR recalibration to SNP recalibrated vcf

## [1.0.0] - 2021-02-08
Initial release of the regenotype-gSNP pipeline

### Added
  -  Regenotype / Joint genotyping pipeline
