# Changelog
All notable changes to the pipeline-call-gSNP pipeline.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]
### Changed
- Enable threading for MergeSamFiles
- Parallelize reheadering and indexing processes

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
