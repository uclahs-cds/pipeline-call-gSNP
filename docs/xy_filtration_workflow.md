# Filter XY calls from a germline VCF file
> **Note**: Currently, the XY Filter would only work on species with XY/XX sex chromosome model.
## Steps:
1. Extract autosomes and chrX/Y variants from input VCF
2. Filter chrX/Y variants
3. Merge autosomal and filtered chrX/Y variants

## chrX/Y Filter Criteria:
- Extract chrX/Y calls
- Extract chrX/Y calls overlapping with Pseudo-Autosomal Regions (PARs)
- if `genetic_sex` is `XY`:
    - For non-PAR chrX/Y calls
        - Filter out heterozygous `GT` calls in chrX and chrY
        - Transform homozygous `GT=1/1` to hemizygous `GT=1`
- If `genetic_sex` is `XX`:
    - Filter out any `chrY` calls regardless of PAR or non-PAR

## Pseudo-Autosomal Regions (PARs)
### GRCh38 - GIAB XY Stratitification PAR [v3.5](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.5/GRCh38@all/XY/)
| CHROM | START | END | PAR |
|---|---|---|---|---|
| chrX | 10000 | 44821 | PAR1 |
| chrX | 94821 | 133871 | PAR1 |
| chrX | 222346 | 226276 | PAR1 |
| chrX | 226351 | 1949345 | PAR1 |
| chrX | 2132994 | 2137388 | PAR1 |
| chrX | 2137488 | 2781479 | PAR1 |
| chrX | 89140830 | 93428068 | PAR3/XTR |
| chrX | 155701383 | 156030895 | PAR2 |
| chrY | 10001 | 44821 | PAR1+PAR3/XTR |
| chrY | 94821 | 133871 | PAR1+PAR3/XTR |
| chrY | 222346 | 226276 | PAR1+PAR3/XTR |
| chrY | 226351 | 1949345 | PAR1+PAR3/XTR |
| chrY | 2132994 | 2137388 | PAR1+PAR3/XTR |
| chrY | 2137488 | 2781479 | PAR1+PAR3/XTR |
| chrY | 3050044 | 6235111 | PAR1+PAR3/XTR |
| chrY | 6532906 | 6748713 | PAR1+PAR3/XTR |
| chrY | 56887902 | 57217415 | PAR2 |

### GRCm39
| CHROM | START | END | PAR | REFERENCE |
|---|---|---|---|---|
| X | 168752755 | 169376592 | PAR | GRCm39 |
| Y | 4072168 | 4161965 | CEN | GRCm39 |
| Y | 90757114 | 91355967 | PAR | GRCm39 |
