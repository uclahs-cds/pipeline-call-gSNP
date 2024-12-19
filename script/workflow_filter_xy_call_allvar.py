#!/usr/bin/env python3
"""
Filter XY calls from call-gSNP single sample VCF file

Steps:
- Extract autosomes and chrX/Y variants from input VCF
- Filter chrX/Y variants
- Merge autosomal and filtered chrX/Y variants

Filter criteria:
- Extract XY calls
- Extract XY calls overlapping with Pseudo-Autosomal Regions (PARs)
- For non-PAR
    - Male sample: Filter out heterozygous GT calls in chrX
    - Female sample: Filter out chrY calls

Dependencies:
- Python 3
- HAIL python library (pip install hail)

Note:
- Do not export VCF to a path that is being read from in the same pipeline,\
based on HAIL recommendation
"""

import os
import argparse
import hail as hl

script_dir = os.getcwd()

parser = argparse.ArgumentParser()
parser.add_argument(
    '--sample_name',
    dest='sample_name',
    help = 'Sample name',
    required=True
    )
parser.add_argument(
    '--input_vcf',
    dest='input_vcf',
    help = 'Input single sample VCF file path',
    required=True
    )
parser.add_argument(
    '--sample_sex',
    dest='sample_sex',
    help = 'Sample sex, XY or XX',
    required=True
    )
parser.add_argument(
    '--par_bed',
    dest='par_bed',
    help = 'Input BED file path for Pseudo-Autosomal Regions (PAR)',
    required=True
    )
parser.add_argument(
    '--output_dir',
    dest='output_dir',
    help = 'Output path where filtered XY variant VCF will be written',
    required=True
    )

args = parser.parse_args()

sample_name = args.sample_name
sample_sex = args.sample_sex
vcf_file = args.input_vcf
par_bed_file = args.par_bed
output_dir = args.output_dir

#Import PAR BED file
par = hl.import_bed(
    path = args.par_bed,
    reference_genome = 'GRCh38',
    skip_invalid_intervals = True
    )

#Extract VCF file header
vcf_header = hl.get_vcf_metadata(vcf_file)
vcf_source = output_dir + '/call-gSNP_caller_source_VCF_header.txt'

#Import VCF file into a hail MatrixTable
vcf_matrix = hl.import_vcf(
    path = vcf_file,
    reference_genome = 'GRCh38',
    force_bgz = True
    )

#Filter XY calls
##Extract XY calls
X_contig = vcf_matrix.locus.contig.startswith('chrX') | vcf_matrix.locus.contig.startswith('X')
Y_contig = vcf_matrix.locus.contig.startswith('chrY') | vcf_matrix.locus.contig.startswith('Y')
extract_condition = (X_contig) | (Y_contig)
vcf_XY = vcf_matrix.filter_rows(extract_condition)
print('chrX/Y variants before XY filtration:', vcf_XY.count())

##Extract autosomes
vcf_autosomes = vcf_matrix.filter_rows(~extract_condition)

##Extract PAR and non-PAR regions
par_variants = vcf_XY.filter_rows(hl.is_defined(par[vcf_XY.locus]))
non_par_variants = vcf_XY.filter_rows(hl.is_missing(par[vcf_XY.locus]))

if sample_sex == 'XY':
    #If MALE (XY), remove heterozygous non-PAR chrX calls
    non_par_filtered_variants = non_par_variants.filter_rows(
        hl.agg.all(
            non_par_variants.GT.is_diploid() & non_par_variants.GT.is_hom_var()
            )
        )
    non_par_filtered_variants = non_par_filtered_variants.annotate_entries(
        GT = hl.call(non_par_filtered_variants.GT[0])
        )

elif sample_sex == 'XX':
    #If Female (XX), remove non-PAR chrY calls
    non_par_filtered_variants = non_par_variants.filter_rows(
        non_par_variants.locus.contig.startswith('chrX') | non_par_variants.locus.contig.startswith('X')
        )

#Combine PAR and filtered non-PAR regions
par_non_par = [par_variants, non_par_filtered_variants]
filterXY = hl.MatrixTable.union_rows(*par_non_par)
print('chrX/Y variant counts after XY filtration:', filterXY.count())

#Combine filtered X/Y + autosomal variants
autosomes_XYfiltered = [vcf_autosomes, filterXY]
output_vcf = hl.MatrixTable.union_rows(*autosomes_XYfiltered)

#Export MatrixTable to VCF
output_file = output_dir + '/' + sample_name + '_filterXY.vcf.bgz'

hl.export_vcf(
    dataset = filterXY,
    output = output_file,
    tabix = True,
    metadata = vcf_header,
    append_to_header = vcf_source
    )
