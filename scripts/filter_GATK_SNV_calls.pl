#!/usr/bin/env perl
#
### filter_GATK_SNV_calls.pl
##################################################################################
# Filter only confident germline calls from the VCFs created by Richard's somatic clasifier script
# 
### HISTORY
#######################################################################################
# Version               Date                Coder                     Comments
# 0.01                  2014-06-15          Takafumi Yamaguchi        Initial development
# 0.02                  2015-02-24          Takafumi Yamaguchi        Added --ref option
# 0.03                  2015-03-16          Takafumi Yamaguchi        Added INDEL filter
# 0.04                  2015-06-11          sespiritu                 Added option to skip ambiguous SNPs/indels, and splitting SNP/indel calls
###INCLUDES
######################################################################################

use strict;
use warnings;
use Getopt::Long;
use feature qw(say);
use Carp;
use Pod::Usage;
use Path::Class;
use YAML qw(LoadFile);
use File::Path qw(make_path);
use File::Spec;
use File::Basename;
use List::MoreUtils qw(first_index);
use Data::Dumper;
### COMMAND LINE DEFAULT ARGUMENTS ################################################################
# list of arguments and default values go here as hash key/value pairs
our %opts = (
	input => undef,
	output_dir => undef,
	sample => undef,
	normal => undef,
	tumour => undef,
	ref => undef,
	filter_somatic => 'Y',
	filter_ambiguous => 'Y',
	split_calls => 'Y'
	);

### MAIN CALLER ###################################################################################
my $result = main();
exit($result);

### FUNCTIONS #####################################################################################

### main ##########################################################################################
# Description:
#   Main subroutine for program
# Input Variables:
#   %opts = command line arguments
# Output Variables:
#   N/A

sub main {
	# get the command line arguments
	GetOptions(
	\%opts,
	"help|?",
	"man",
	"input=s" => \$opts{'input'},
	"output_dir:s" => \$opts{'output_dir'},
	"ref=s" => \$opts{'ref'},
	"sample:s" => \$opts{'sample'},
	"normal:s" => \$opts{'normal'},
	"tumour:s" => \$opts{'tumour'},
	"filter_somatic:s" => \$opts{'filter_somatic'},
	"filter_ambiguous:s" => \$opts{'filter_ambiguous'},
	"split_calls:s" => \$opts{'split_calls'},
	"dir|d:s"
	) or pod2usage(64);

	if ($opts{'help'}) { pod2usage(1) };
	if ($opts{'man'}) { pod2usage(-exitstatus => 0, -verbose => 2) };

	while(my ($arg, $value) = each(%opts)) {
		if (!($arg=~/\:/) and !defined $value) {
			print "ERROR: Missing argument $arg\n";
			pod2usage(2);
			}
		}		

	my ($data_snv, $data_indel, $data_combined);

	#read vcf
	my($vcf_filename, $vcf_directory, $suffix) = fileparse($opts{'input'}, ".vcf");
	print "Input VCF file is ".$opts{'input'}."\n";

	if (defined $opts{output_dir}) {
		$vcf_directory = $opts{output_dir};
		}
	
	# define the output VCF file name here 
	my $output_files = {};
	$output_files->{snp} = file($vcf_directory, join("_", 'filtered_germline_snv', $opts{'sample'}.'_nosomatic.vcf'));
	$output_files->{indel} = file($vcf_directory, join("_", 'filtered_germline_indel', $opts{'sample'}.'_nosomatic.vcf'));
	$output_files->{combined} = file($vcf_directory, join("_", 'filtered_germline_all', $opts{'sample'}.'.vcf'));

	my $vcf_fh;
	if ($opts{'input'} =~ /.gz$/) {
		open($vcf_fh, "gunzip -c $opts{'input'} |") || die "Can't open pipe to $opts{'input'}!\n";
		}
	else {
		open ($vcf_fh, $opts{'input'}) || die "Can't open file $opts{'input'}!\n";
		}

	my @metadata;
	my ($normal_index, $tumour_index, $genotype, $class);
	while (<$vcf_fh>) {
		my $line = $_;
		chomp($line);
		if ($line=~/^##/) {
			push @metadata, $line;
			next;
			} #Skip meta data

		my @info = split /\t/, $line;
		my $filter = $info[6];
				
		if ($line=~/^#CHROM/ and $opts{filter_somatic} and $opts{normal} and $opts{tumour}) {
			# pull out GT for normal and tumour (FORMAT should be 9th)
			$normal_index = first_index { $_ eq $opts{normal} } @info;
			$tumour_index = first_index { $_ eq $opts{tumour} } @info;

			push @metadata, $line;
			next;
			} # skip somatic 

		if ($filter ne 'PASS') { next; }	# skip non-PASSED calls
		
		my $chr = $info[0];
		my $pos = $info[1];
		my $id = $info[2];
		my $ref = $info[3];
		my $alt = $info[4];

		# -----------------------------------------------------------------------------------------
		# variant classfication
		if ($opts{filter_somatic} and $opts{normal} and $opts{tumour}) {
			# get indices for GT, normal and tumour columns
			my @format_elements = split /\:/, $info[8];
			my $gt_index = first_index { $_ eq 'GT' } @format_elements;
			
			my @format_normal = split /\:/, $info[($normal_index)];
			my @format_tumour = split /\:/, $info[($tumour_index)];

			if (!$format_tumour[$gt_index] or !$format_normal[$gt_index]) {
				croak "NO GT field found although GT filed is indicated in the FORMAT column\n";
				}
			
			my @ref_entries = split /,/, $ref;
			my @alt_entries = split /,/, $alt;
			my $type;

			# determine variant type (SNP/INDEL)
			foreach my $elm (@ref_entries, @alt_entries) {
				#INDEL - at least one entry should be length > 1.
				if(length($elm) > 1){
					$type = 'INDEL';
					last;
					}
					$type = 'SNV';
				}

			# no info for either normal or tumour to classify variant
			if ('./.' eq $format_tumour[$gt_index] or './.' eq $format_normal[$gt_index]) {
				$class->{$type}->{unknown}++;
				}

			if ($format_normal[$gt_index] eq $format_tumour[$gt_index]) {
				$class->{$type}->{germline}++;
				}
			# SOMATIC or LOH
			else {
				my @normal_gt = split /\//, $format_normal[$gt_index];
				my @tumour_gt = split /\//, $format_tumour[$gt_index];

				# TO DO: handle more than 3 allels? It seems there's no 3 allels in GATK VCF due to the ploidy setting?
				if (scalar @normal_gt > 2 or scalar @tumour_gt > 2) {
					print "CNA?\n";
					}
				
				if(	($normal_gt[0] eq $normal_gt[1] and $tumour_gt[0] ne $tumour_gt[1]) or
					($normal_gt[0] eq $normal_gt[1] and $tumour_gt[0] eq $tumour_gt[1])
					) {
					# First case e.g.
					#	0/0 -> 0/1
					#	1/1 -> 0/1
					# Second case e.g.
					#	0/0 -> 1/1
					$class->{$type}->{SOMATIC}++;
					}
				elsif ($normal_gt[0] ne $normal_gt[1] and $tumour_gt[0] eq $tumour_gt[1]) {
					# e.g.
					#	0/1 -> 0/0
					#	0/1 -> 1/1
					$class->{$type}->{LOH}++;
					}
				elsif ($normal_gt[0] ne $normal_gt[1] and $tumour_gt[0] ne $tumour_gt[1]) {
					# e.g.
					#	0/1 -> 0/2
					#	1/2 -> 0/2
					#	2/3 -> 1/3
					$class->{$type}->{SOMATIC}++;
					}
				}

			$genotype->{$type}->{$format_normal[$gt_index].':'.$format_tumour[$gt_index]}++;
			
			# removing non-germline SNV/INDELs and unknown variants in normal
			# these cases are not germline or unknown -> Skip!
			if ($format_normal[$gt_index] eq '0/0' or $format_normal[$gt_index] eq './.') { next; }
			}	# variant classification
	
		# -----------------------------------------------------------------------------------------
		# filtering SNP and indel
		my $ref_length = length($ref);
		my $alt_length = length($alt);

		if((1 ne $ref_length) or (1 ne $alt_length)) {
			if(('Y' eq $opts{'filter_ambiguous'}) and ($alt =~/\,/)) {
				next;	# skip ambiguous SNPs/indel
				}

			# multiple ALT alleles -> could be an SNP, indel, or both
			# e.g.
			#	REF = A; ALT = C,G			(SNP)
			#	REF = A; ALT = C,GT			(SNP/indel)
			#	REF = A; ALT = CC,GG		(indel)
			elsif((1 eq $ref_length) and ($alt =~/\,/)) {
				# split by ALT allele and calculate their lengths, sorted in increasing order
				my @alleles = split(',', $alt);
				my @alleles_length = sort(map{ length($_) } @alleles);
				
				# if there exists an ALT allele of length 1 by checking min length
				if($alleles_length[0] == 1) {
					$data_snv->{$chr}->{$pos} = $line;

					# if there exists an ALT allele of length > 1 by checking max length
					if($alleles_length[-1] > 1) {
						$data_indel->{$chr}->{$pos} = $line;
						}
					$data_combined->{$chr}->{$pos} = $line;
					}
				}
			# e.g.
			#	REF = A; ALT = CG
			#	REF = AC; ALT = G
			#	REF = AC; ALT = GTT
			#	REF = AC; ALT = G,TTT
			else {
				$data_indel->{$chr}->{$pos} = $line;
				$data_combined->{$chr}->{$pos} = $line;
				}
			}
		# e.g.	REF = A; ALT = C
		else {
			$data_snv->{$chr}->{$pos} = $line;
			$data_combined->{$chr}->{$pos} = $line;
			}
		}	# while reading VCF

	close $vcf_fh;

	# -----------------------------------------------------------------------------------------
	# output genotype and variant class counts
	if ($opts{filter_somatic} and $opts{normal} and $opts{tumour}) {
		# define the output VCF file name here 
		my $output_genotype_count = file(
			$vcf_directory,
			join("_",
				'filtered_germline_genotype_count',
				$opts{'sample'}.'.tsv'
				) );
		my $output_variant_class_count = file(
			$vcf_directory,
			join("_",
				'filtered_germline_variant_class_count',
				$opts{'sample'}.'.tsv'
				) );
		
		# printing genotype counts
		open (my $fh_genotype, '>', $output_genotype_count);
		foreach my $type (sort {$a cmp $b} keys %{$genotype}){
			foreach my $genotype_cmb (sort {$a cmp $b} keys %{$genotype->{$type}}){
				print $fh_genotype join("\t", $type, $genotype_cmb, $genotype->{$type}->{$genotype_cmb})."\n";
				}
			}
		close($fh_genotype);

		# print variant counts
		open (my $fh_class, '>', $output_variant_class_count);
		foreach my $type (sort {$a cmp $b} keys %{$class}){
			foreach my $variant_class (sort {$a cmp $b} keys %{$class->{$type}}){
				print $fh_class join("\t", $type, $variant_class, $class->{$type}->{$variant_class})."\n";
				}
			}
		close($fh_class);
		}

	my $index_file = $opts{'ref'}.'.fai';
	open (my $fh_index, $index_file);
	my $metadata = join("\n", @metadata);

	if('Y' eq $opts{'split_calls'}) {
		print "Output VCF file for SNV is ".$output_files->{snp}."\n";
		print "Output VCF file for INDEL is ".$output_files->{indel}."\n";

		open (my $fh_snv_output, ">", $output_files->{snp});
		open (my $fh_indel_output, ">", $output_files->{indel});
		print $fh_snv_output $metadata."\n";
		print $fh_indel_output $metadata."\n";

		while(<$fh_index>){
			my $line = $_;
			my @info_chr = split /\t/, $line;
			my $chr = $info_chr[0];
			
			foreach my $pos (sort {$a <=> $b} keys (%{$data_snv->{$chr}})){ #sorting positions
				print $fh_snv_output $data_snv->{$chr}->{$pos}."\n";
				}#end of foreach (pos)
	
			foreach my $pos (sort {$a <=> $b} keys (%{$data_indel->{$chr}})){ #sorting positions
				print $fh_indel_output $data_indel->{$chr}->{$pos}."\n";
				}
			}#end of while
		close $fh_snv_output;
		close $fh_indel_output;
		}
	# else combine them into one file
	else {
		print "Output VCF file for combined SNV/INDEL is ".$output_files->{combined}."\n";

		open (my $fh_combined_output, ">", $output_files->{combined});
		print $fh_combined_output $metadata."\n";

		while(<$fh_index>){
			my $line = $_;
			my @info_chr = split /\t/, $line;
			my $chr = $info_chr[0];

			foreach my $pos (sort {$a <=> $b} keys (%{$data_combined->{$chr}})){ #sorting positions
				print $fh_combined_output $data_combined->{$chr}->{$pos}."\n";
				}#end of foreach (pos)
			}#end of while
		close $fh_combined_output;
		}

	close $fh_index;

	# -----------------------------------------------------------------------------------------
	# bgzipping and indexing VCF files
	foreach my $type (keys %{$output_files}) {
		if (-e $output_files->{$type}) {
			system(join(' ; ',
				'bgzip -f ' . $output_files->{$type},
				'tabix -p vcf ' . $output_files->{$type} . '.gz'
				) );
			unlink $output_files->{$type};
			}
		}
	return 0;
	}

__END__


=head1 NAME

filter_GATK_SNV_calls.pl

=head1 SYNOPSIS

B<filter_GATK_SNV_calls.pl> [options] [file ...]

	Options:
	--help               brief help message
	--man                full documentation
	--input              input VCF (required)
	--sample             sample name (required)
	--normal             normal id (required)
	--tumour             tumour id (required)
	--ref                reference genome (required)
	--filter_somatic     filter out any misclassified somatic mutations (default: Y)
	--filter_ambiguous   filter out "ambiguous" SNPs/indels (default: Y)
	--split_calls        split germline calls into SNPs and indels (default: Y)

=head1 OPTIONS

=over 8

=item B<--help>

Print a brief help message and exit.

=item B<--man>

Print the manual page.

=item B<--input>

Input VCF file.

=item B<--sample>

Sample name.

=item B<--normal>

Normal id.

=item B<--sample>

Tumour id.

=item B<--ref>

Reference genome.

=item B<--filter_somatic>

Whether to filter out misclassified somatic mutations.

=item B<--filter_ambiguous>

Whether to filter out germline calls with multiple ALT alleles.

=item B<--split_calls>

Whether to split germline calls into SNPs and indels.

=back

=head1 DESCRIPTION

B<filter_GATK_SNV_calls.pl> Filter only confident germline calls from the VCFs created by Richard's somatic clasifier script

=head1 EXAMPLE

filter_GATK_SNV_calls.pl --input input.vcf --sample sample.name --normal normal.id --tumour tumour.id

=head1 AUTHOR

Takafumi Yamaguchi  -- Boutros Lab

With modifications by Shadrielle Melijah G. Espiritu -- Boutros Lab

The Ontario Institute for Cancer Research

=head1 ACKNOWLEDGEMENTS

Paul Boutros, PhD, PI - Boutros Lab
