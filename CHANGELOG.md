# Changelog
All notable changes to the pipeline-call-gSNP pipeline.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]
### Changed
- Reorganize repo with pipeline entrypoint at root of repo

---

## [8.0.0] - 2022-03-14
### Added
- BETA: Support for paired inputs with a single normal sample and multiple tumour samples

### Changed
- Switch to SAMtools for indexing BAMs
- Use sample ID and intervals as identifiers for log output directories
- Standardize config structure
- Partially revert BQSR parallelization and group ApplyBQSR by interval
- Parallelize BaseRecalibrator per sample
- Save VQSR output for QC
- Save SNP+INDEL VQSRed VCF to output

---

## [7.2.1] - 2022-01-25
### Changed
- Parallelize BQSR
- Update .gitignore

### Security
- Update GATK to 4.2.4.1 to address Log4j vulnerabilities (https://github.com/advisories/GHSA-8489-44mv-ggj8, https://github.com/advisories/GHSA-p6xc-xr62-6r2g)
- Update Picard version to 2.26.10 to address Log4j vulnerabilities (https://github.com/advisories/GHSA-8489-44mv-ggj8)

---

## [7.2.0] - 2021-12-17
### Changed
- Enable threading for MergeSamFiles
- Parallelize reheadering and indexing processes
- Update reheadering to use -c option
- Modularize workflows for different modes (single vs. paired, WGS vs targeted)
- Update GATK to 4.2.4.0 to address Log4j critical vulnerability (https://github.com/advisories/GHSA-jfh8-c2jp-5v3q)
- Update Picard to 2.26.8 to address Log4j critical vulnerability (https://github.com/advisories/GHSA-jfh8-c2jp-5v3q) 

---

## [7.1.0] - 2021-11-09
### Added
- Parallelize IR and BQSR in WXS/WES mode

### Fixed
- Fix targeted, single sample mode bugs

---

## [7.0.0] - 2021-10-27
### Added
- Update call-gSNP to DSL2
- Add GPL2 license
- Parallelize MergeVcfs
- Parallelize MergeSamFiles
- Standardize output and log directories
- Add process to remove intermediate files when save_intermediate_files is disabled
- Parallelize GetPileupSummaries, CalculateContamination, and DepthOfCoverage processes
- Split HaplotypeCaller process into process for VCF and GVCF modes
- Parallelize GVCF HC process
- Extract genome intervals from reference dictionary

### Changed
- Adjust static resource allocation to be more efficient
- Auto-detect reference fasta dictionary
- Rename ".bai" output files to ".bam.bai"
- Auto-detect when in targeted mode and when in WGS mode
