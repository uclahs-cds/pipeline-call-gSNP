#!/usr/bin/env python3
Task
    - add sample name input
    - make script work if VCF has multiple samples
Note
    - script works on both single sample and multi sample vcf
    - filtration can be done using hom GT, not necessarily AF
"""
Filter XY calls from call-gSNP single sample VCF file

Filter criteria:
- Extract XY calls
- Extract XY calls overlapping with Pseudo-Autosomal Regions (PARs)
- For non-PAR
    - Male sample: Filter out chrX calls where VAF < 80%
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

sample = args.sample_name
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
vcf_source = script_dir + '/call-gSNP_caller_source_VCF_header.txt'

#Import VCF file into a hail MatrixTable
vcf_matrix = hl.import_vcf(
    path = vcf_file,
    reference_genome = 'GRCh38',
    force_bgz = True
    )

###Get Sample ID
sample_name = vcf_matrix.s.collect()[0] #list has one sample ID

#Filter XY calls
##Extract XY calls
X_contig = vcf_matrix.locus.contig.startswith('chrX')
Y_contig = vcf_matrix.locus.contig.startswith('chrY')
extract_condition = (X_contig) | (Y_contig)
vcf_XY = vcf_matrix.filter_rows(extract_condition)
print('variants in chrX/Y:', vcf_XY.count())

##Remove calls with DP=0
depth_field = sample_name + '.DP'
depth = vcf_XY.make_table()
zero_depth = depth.filter(depth[depth_field] == 0)
zero_depth_count = zero_depth.count()
print('variants with zero depth:', zero_depth_count)
zero_depth_contig = zero_depth.locus.contig.collect()
zero_depth_pos = zero_depth.locus.position.collect()

for allele in range(zero_depth_count):
    zero_depth_match = (vcf_XY.locus.contig == zero_depth_contig[allele]) & (vcf_XY.locus.position == zero_depth_pos[allele])
    vcf_XY = vcf_XY.filter_rows(
        ~(zero_depth_match)
    )

print('variants after filtering zero depth:', vcf_XY.count())

##Extract PAR and non-PAR regions
par_filtered = vcf_XY.filter_rows(hl.is_defined(par[vcf_XY.locus]))
non_par_filtered = vcf_XY.filter_rows(hl.is_missing(par[vcf_XY.locus]))

###For non-PAR regions, extract VQSR PASS calls. (note: Hail parses PASS as an empty set {})
#non_par_filtered = non_par_filtered.filter_rows(
#    hl.len(non_par_filtered.filters) == 0
#    )

##Predict SEX of the sample
#imputed_sex = hl.impute_sex(vcf_file.GT)
#temp place holder for SEX
SEX = 'XY'

if SEX == 'XY':
    #If MALE (XY), remove non-PAR chrX calls with AF=0.5
    filter_non_par_call = non_par_filtered.filter_rows(
        non_par_filtered.info.AF[0] != 0.5
        )
elif SEX == 'XX':
    #If Female (XX), remove non-PAR chrY calls
    filter_non_par_call = non_par_filtered.filter_rows(
        non_par_filtered.locus.contig.startswith('chrX')
        )

#Combine PAR and filtered non-PAR regions
par_non_par = [par_filtered, filter_non_par_call]
filterXY = hl.MatrixTable.union_rows(*par_non_par)

#Export MatrixTable to VCF
output_file = output_dir + '/' + sample_name + '_filterXY.vcf.bgz'

hl.export_vcf(
    dataset = filterXY,
    output = output_file,
    tabix = True,
    metadata = vcf_header,
    append_to_header = vcf_source
    )
