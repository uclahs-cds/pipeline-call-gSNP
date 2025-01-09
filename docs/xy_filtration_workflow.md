# Filter XY calls from a germline VCF file

## Steps:
1. Extract autosomes and chrX/Y variants from input VCF
2. Filter chrX/Y variants
3. Merge autosomal and filtered chrX/Y variants

## chrX/Y Filter Criteria:
- Extract chrX/Y calls
- Extract chrX/Y calls overlapping with Pseudo-Autosomal Regions (PARs)
- For non-PAR chrX/Y calls
    - if `sample_sex` is `XY`:
        - Filter out heterozygous `GT` calls in chrX and chrY
        - Transform homozygous `GT=1/1` to hemizygous `GT=1`
    - if `sample_sex` is `XX`:
        - Filter out `chrY` calls