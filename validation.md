Validation and benchmarking with A-full and CPCG0196-B1 - single WGS mode
1. [A-full](#a-full)
2. [CPCG0196-B1](#cpcg0196-b1)

## A-full

Summary tables comparing call-gSNP v7.0.0 to call-gSNP-DSL1 v5.4.5

Input: `/hot/pipeline/development/slurm/call-gSNP/7.0.0-rc1/pre_standardization/single_A-full/A-full.csv`

- **call-gSNP v7.0.0**: `/hot/pipeline/development/slurm/call-gSNP/7.0.0-rc1/pre_standardization/single_A-full`
    - Runtime: 2d 4h 14m 2s
- **call-gSNP-DSL1 v5.4.5**: `/hot/pipeline/development/slurm/call-gSNP/DSL1_outputs/v5.4.5_for_comparison`
    - Runtime: 2d 15h 9m 53s

Raw VCF comparison:
|Version|VCF (Common)|VCF (Unique)|
|:---:|:---:|:---:|
|7.0.0|4813793|74 (0.0015%)|
|DSL1-5.4.5|4813793|62 (0.0013%)|

Raw GVCF comparison:
|Version|VCF (Common)|VCF (Unique)|
|:---:|:---:|:---:|
|7.0.0|51962094|3234 (0.0062%)|
|DSL1-5.4.5|51962094|3454 (0.0066%)|

Filtered SNV GVCF:
|Version|VCF (Common)|VCF (Unique)|
|:---:|:---:|:---:|
|7.0.0|3265039|24550 (0.74%)|
|DSL1-5.4.5|3265039|17728 (0.54%)|

Filtered INDEL GVCF:
|Version|VCF (Common)|VCF (Unique)|
|:---:|:---:|:---:|
|7.0.0|811214|2875 (0.35%)|
|DSL1-5.4.5|811214|6501 (0.79%)|

Realigned and recalibrated BAM:
|SAMPLE|LEFT_FILE|RIGHT_FILE|MAPPINGS_MATCH|MAPPINGS_DIFFER|UNMAPPED_BOTH|UNMAPPED_LEFT|UNMAPPED_RIGHT|MISSING_LEFT|MISSING_RIGHT|DUPLICATE_MARKINGS_DIFFER|ARE_EQUAL|
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
|Normal|DSL1-5.4.5|7.0.0|1723823909|11|4585406|0|0|3753|0|0|N|

## CPCG0196-B1

Summary tables comparing call-gSNP v7.0.0 to call-gSNP-DSL1 v5.4.3

Input: `/hot/pipeline/development/slurm/call-gSNP/7.0.0-rc1/pre_standardization/single_CPCG0196-B1/CPCG0196-B1.csv`

- **call-gSNP v7.0.0**: `/hot/pipeline/development/slurm/call-gSNP/7.0.0-rc1/pre_standardization/single_CPCG0196-B1`
    - Runtime: 2d 13h 30m 29s
- **call-gSNP-DSL1 v5.4.3**: `/hot/pipeline/development/slurm/call-gSNP/DSL1_outputs/call-gSNP_v5.4.3/outputs_CPCG0196-B1/call-gSNP-20210324-210359/`
    - Runtime: 3d 0h 9m 53s

Raw VCF comparison:
|Version|VCF (Common)|VCF (Unique)|
|:---:|:---:|:---:|
|7.0.0|4863231|133768 (2.68%)|
|DSL1-5.4.3|4863231|577 (0.012%)|

Raw GVCF comparison:
|Version|VCF (Common)|VCF (Unique)|
|:---:|:---:|:---:|
|7.0.0|51588596|3653326 (6.61%)|
|DSL1-5.4.3|51588596|7506 (0.014%)|

Filtered SNV GVCF:
|Version|VCF (Common)|VCF (Unique)|
|:---:|:---:|:---:|
|7.0.0|3035126|130500 (4.02%)|
|DSL1-5.4.3|3035126|77095 (2.38%)|

Filtered INDEL GVCF:
|Version|VCF (Common)|VCF (Unique)|
|:---:|:---:|:---:|
|7.0.0|844167|6042 (0.71%)|
|DSL1-5.4.3|844167|361 (0.042%)|

Realigned and recalibrated BAM:
|SAMPLE|LEFT_FILE|RIGHT_FILE|MAPPINGS_MATCH|MAPPINGS_DIFFER|UNMAPPED_BOTH|UNMAPPED_LEFT|UNMAPPED_RIGHT|MISSING_LEFT|MISSING_RIGHT|DUPLICATE_MARKINGS_DIFFER|ARE_EQUAL|
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
|Normal|DSL1-5.4.3|7.0.0|1935641125|182|16050267|0|0|21372|0|0|N|