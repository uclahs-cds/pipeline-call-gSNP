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

## Pseudo-Autosomal Regions (PARs)
### GRCh38
| CHROM | START | END | PAR | REGION | REFERENCE |
|---|---|---|---|---|---|
| chrX | 10001 | 2781479 | PAR1 | Xp22 | EMSEMBL |
| chrX | 91434839 | 91438584 | PAR3/XTR | Xq21.3 | PMID:23708688 |
| chrX | 155701383 | 156030895 | PAR2 | Xq28 | ENSEMBL |
| chrY | 10001 | 10300000 | PAR1+PAR3/XTR | Yp11 | ENSEMBL +PMID:23708688 |
| chrY | 56887903 | 57217415 | PAR2 | Yq12 | ENSEMBL |
