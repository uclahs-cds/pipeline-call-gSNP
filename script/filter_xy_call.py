#!/usr/bin/env python3
"""
Filter XY calls from a germline VCF file

Steps:
- Extract autosomes and chrX/Y variants from input VCF
- Filter chrX/Y variants
- Merge autosomal and filtered chrX/Y variants

Filter criteria:
- Extract XY calls
- Extract XY calls overlapping with Pseudo-Autosomal Regions (PARs)
- Male sample (XY):
    - For non-PAR
        - Filter out heterozygous GT calls in chrX and chrY
        - Transform homozygous GT=1/1 to hemizygous GT=1
- Female sample (XX): Filter out any chrY calls

Dependencies:
- Python 3
- HAIL python library (pip install hail)

Note:
- Do not export VCF to a path that is being read from in the same pipeline,\
based on HAIL recommendation
"""

import argparse
import os
import sys
import tempfile
import hail as hl

script_dir = os.getcwd()

parser = argparse.ArgumentParser()
parser.add_argument(
    '--input_vcf',
    dest='input_vcf',
    help = 'Input single sample VCF file path',
    required=True
    )
parser.add_argument(
    '--vcf_source_file',
    dest='vcf_source_file',
    help = 'A TXT file containing variant caller source details (eg. ##source=HaplotypeCaller)',
    required=True
    )
parser.add_argument(
    '--genetic_sex',
    dest='genetic_sex',
    help = 'Sample genetic sex, XY or XX',
    required=True
    )
parser.add_argument(
    '--par_bed',
    dest='par_bed',
    help = 'Input BED file path for Pseudo-Autosomal Regions (PAR)',
    required=True
    )
parser.add_argument(
    '--ref-genome',
    dest='ref_genome',
    help = 'Reference genome used to produce input VCF',
    required=True
    )
parser.add_argument(
    '--ref-index',
    dest='ref_index',
    help = 'Reference genome index file',
    required=True
    )
parser.add_argument(
    '--cpu',
    dest='cpu',
    help = 'Number of CPUs to initialize hail',
    required=False
    )
parser.add_argument(
    '--output_name',
    dest='output_name',
    help = 'Output file name',
    required=True
    )
parser.add_argument(
    '--output_dir',
    dest='output_dir',
    help = 'Output path where filtered XY variant VCF will be written',
    required=True
    )

args = parser.parse_args()

genetic_sex = args.genetic_sex
vcf_file = args.input_vcf
vcf_source_file = args.vcf_source_file
par_bed = args.par_bed
reference_genome = args.ref_genome
reference_index = args.ref_index
output_name = args.output_name
output_dir = args.output_dir

if args.cpu is not None:
    #Initialize CPUs
    cpu = args.cpu
    hl.init(master=f'local[{cpu}]', backend='spark')

#args: matrix table, 'XY' or 'XX', 'before' or 'after'
def get_xy_counts(matrix_table, input_sex, state):
    """Function to get variant counts in chrX and chrY"""
    chr_x = matrix_table.filter_rows(
        matrix_table.locus.contig.startswith('chrX') | matrix_table.locus.contig.startswith('X')
        )
    chr_y = matrix_table.filter_rows(
        matrix_table.locus.contig.startswith('chrY') | matrix_table.locus.contig.startswith('Y')
        )
    result = print(
        f'chrX variant counts {state} {input_sex} filtration:', chr_x.count_rows(),
        f'\nchrY variant counts {state} {input_sex} filtration:', chr_y.count_rows()
        )
    return result

#Get reference for input VCF
input_reference = hl.ReferenceGenome.from_fasta_file(
    "input_reference",
    reference_genome,
    reference_index
    )

#Extract VCF file header
vcf_header = hl.get_vcf_metadata(vcf_file)

#Add script system command to VCF source
SCRIPT_COMMAND = ' '.join(sys.argv)

with open(vcf_source_file, 'r', encoding='utf-8') as vcf_source:
    vcf_source_content = vcf_source.read()

script_command_entry = f'##XYFiltration=<CommandLine={SCRIPT_COMMAND}>'
vcf_source = vcf_source_content + script_command_entry
temp_file_path = os.path.join(tempfile.gettempdir(), 'temp_file.txt')
with open(temp_file_path, 'w', encoding='utf-8') as temp_file:
    temp_file.write(vcf_source)

#Import PAR BED file
par = hl.import_bed(
    path = par_bed,
    reference_genome = input_reference,
    skip_invalid_intervals = True
    )

#Import VCF file into a hail MatrixTable
vcf_matrix = hl.import_vcf(
    path = vcf_file,
    reference_genome = input_reference,
    force_bgz = True
    )

#Filter XY calls
##Extract XY calls
X_contig = vcf_matrix.locus.contig.startswith('chrX') | vcf_matrix.locus.contig.startswith('X')
Y_contig = vcf_matrix.locus.contig.startswith('chrY') | vcf_matrix.locus.contig.startswith('Y')
extract_condition = (X_contig) | (Y_contig)
vcf_XY = vcf_matrix.filter_rows(extract_condition)
##Show chrX/Y counts before XY filter
get_xy_counts(vcf_XY, genetic_sex, 'before')

##Extract autosomes
vcf_autosomes = vcf_matrix.filter_rows(~extract_condition)

##Extract PAR and non-PAR regions
par_variants = vcf_XY.filter_rows(hl.is_defined(par[vcf_XY.locus]))
non_par_variants = vcf_XY.filter_rows(hl.is_missing(par[vcf_XY.locus]))

if genetic_sex == 'XY':
    #If MALE (XY), remove heterozygous non-PAR chrX calls
    non_par_filtered_variants = non_par_variants.filter_rows(
        hl.agg.all(
            non_par_variants.GT.is_diploid() & non_par_variants.GT.is_hom_var()
            )
        )
    non_par_filtered_variants = non_par_filtered_variants.annotate_entries(
        GT = hl.call(non_par_filtered_variants.GT[0])
        )
    #Combine PAR and filtered non-PAR regions
    par_non_par = [par_variants, non_par_filtered_variants]
    filterXY = hl.MatrixTable.union_rows(*par_non_par)
    #Show chrX/Y counts after XY filter
    get_xy_counts(filterXY, genetic_sex, 'after')
    #Combine filtered XY + autosomal variants
    autosomes_sex_filtered = [vcf_autosomes, filterXY]

elif genetic_sex == 'XX':
    #If Female (XX), keep ONLY chrX calls
    filterXX = vcf_XY.filter_rows(
        vcf_XY.locus.contig.startswith('chrX') | \
            vcf_XY.locus.contig.startswith('X')
        )
    #Show chrX/Y counts after XX filter
    get_xy_counts(filterXX, genetic_sex, 'after')
    #Combine filtered XX + autosomal variants
    autosomes_sex_filtered = [vcf_autosomes, filterXX]

#Export MatrixTable to VCF
output_vcf = hl.MatrixTable.union_rows(*autosomes_sex_filtered)
OUTPUT_FILE = f'{output_dir}/{output_name}.vcf.bgz'

hl.export_vcf(
    dataset = output_vcf,
    output = OUTPUT_FILE,
    tabix = True,
    metadata = vcf_header,
    append_to_header = temp_file_path
    )
