#!/usr/bin/perl
use strict;
use warnings;
use v5.14;

## Sandesh Shrestha ##
## sshrest1@ksu.edu ##
## Oct 2019 ##

$| = 1;

my $fq1 = $ARGV[0];
my $fq2 = $ARGV[1];
my $i1 = $ARGV[2];
my $i2 = $ARGV[3];
my $barcode_file = $ARGV[4];

unless (defined $ARGV[0] and defined $ARGV[1] and defined $ARGV[2] and defined $ARGV[3] and defined $ARGV[4]){
  say "$0 file.R1.fastq file.R2.fastq i1.fastq i2.fastq barcode.txt";
  say "\nBarcode file:\nSample1 i1_i2\nSample2 AACCACTC_AAGACTGG";
  exit;
}


if ($fq1 =~ /\.gz$/i){
open INFQ1, "gunzip -c $fq1 |" or die "File not found in the folder";
}
else{
open INFQ1, $fq1 or die "File not found in the folder";
}

if ($fq2 =~ /\.gz$/i){
open INFQ2, "gunzip -c $fq2 |" or die "File not found in the folder";
}
else{
open INFQ2, $fq2 or die "File not found in the folder";
}

if ($i1 =~ /\.gz$/i){
open INI1, "gunzip -c $i1 |" or die "File not found in the folder";
}
else{
open INI1, $i1 or die "File not found in the folder";
}

if ($i2 =~ /\.gz$/i){
open INI2, "gunzip -c $i2 |" or die "File not found in the folder";
}
else{
open INI2, $i2 or die "File not found in the folder";
}

open( BC_FILE, $barcode_file ) or die "No barcode file";

my ( @sample_file_names, @barcode );
foreach (<BC_FILE>) {
	chomp $_;

	#say $_;
	my @elements = split( "\t", $_ );
	
	push @sample_file_names, $elements[0];
	#say $elements[0];
	
	push @barcode, $elements[1];
    	#say $elements[1];
}


my @fh_array_R1;
my @fh_array_R2;
foreach (@sample_file_names) {
	#say $_;
	local *OUT_R1;
	local *OUT_R2;
	my $r1_fn = "$_\_R1.fq.gz";
	my $r2_fn = "$_\_R2.fq.gz";
	open( OUT_R1, "| gzip -1 > $r1_fn") or die "cannot write file";
	open( OUT_R2, "| gzip -1 > $r2_fn") or die "cannot write file";
	
	#open( OUT_R1, ">", "$_\_R1.fq" ) or die "cannot write file";
	#open( OUT_R2, ">", "$_\_R2.fq" ) or die "cannot write file";
	
	push @fh_array_R1, *OUT_R1;
	push @fh_array_R2, *OUT_R2;        
}

# make hash to count reads matched to sample
my %count_sample_reads;
foreach (@sample_file_names) {

	# initialiaing all sample reads as zero
	$count_sample_reads{$_} = 0;
}

# for unknown-reads count
my $unknown_reads = 0;

# total count reads in a file
my $total_reads_file = 0;

# unknown barcode file
open( UNKNOWN_R1, "| gzip -1 > unknown-barcode-R1_$fq1-$barcode_file.fq.gz" ) or die "cannot create unknown-barcode-R1_$fq1.fq";
open( UNKNOWN_R2, "| gzip -1 > unknown-barcode-R2_$fq2-$barcode_file.fq.gz" ) or die "cannot create unknown-barcode-R2_$fq2.fq";

#open( UNKNOWN_R1, ">unknown-barcode-R1_$fq1-$barcode_file.fq" ) or die "cannot create unknown-barcode-R1_$fq1.fq";
#open( UNKNOWN_R2, ">unknown-barcode-R2_$fq2-$barcode_file.fq" ) or die "cannot create unknown-barcode-R2_$fq2.fq";

while ( my $p1_first_line  = <INFQ1>, my $p1_second_line = <INFQ1>, my $p1_third_line  = <INFQ1>, my $p1_fourth_line = <INFQ1>,
		my $p2_first_line  = <INFQ2>, my $p2_second_line = <INFQ2>, my $p2_third_line  = <INFQ2>, my $p2_fourth_line = <INFQ2>,
		my $i1_first_line  = <INI1>, my $i1_second_line = <INI1>, my $i1_third_line  = <INI1>, my $i1_fourth_line = <INI1>,
		my $i2_first_line  = <INI2>, my $i2_second_line = <INI2>, my $i2_third_line  = <INI2>, my $i2_fourth_line = <INI2>
) {

	chomp( $p1_first_line, $p1_second_line, $p1_third_line, $p1_fourth_line, $p2_first_line, $p2_second_line, $p2_third_line, $p2_fourth_line, $i1_second_line, $i2_second_line);
	
	# concat indices
	my $i1i2 = $i1_second_line ."_" . $i2_second_line;
	#say $i1i2;
	
	my $matched_R1 = "$p1_first_line\n$p1_second_line\n$p1_third_line\n$p1_fourth_line\n";
	my $matched_R2 = "$p2_first_line\n$p2_second_line\n$p2_third_line\n$p2_fourth_line\n";

	#print $no_matched_R1;

# not-match count to avoid duplication inside loop: to compare with total barcodes, if equal, then no match
	my $nomatch_count = 0;
	for ( my $j = 0 ; $j < scalar @barcode ; $j++ ) {

		if ( $i1i2 eq $barcode[$j] ) {
			
			print { $fh_array_R1[$j] } $matched_R1;
			print { $fh_array_R2[$j] } $matched_R2;
			
			$count_sample_reads{ $sample_file_names[$j] }++;
			last;
		}
		else {
			$nomatch_count++;

			#print single unknown
			if ( $nomatch_count == scalar @barcode) {

				#print to unknown;
				print UNKNOWN_R1 $matched_R1;
				print UNKNOWN_R2 $matched_R2;

				# count unknown reads
				$unknown_reads++;

			}
		}
	}
	$total_reads_file++;
}

open( REPORT, ">summary_report-$fq1-$barcode_file.txt" );

my $sample_total_read = 0;
foreach my $key ( keys %count_sample_reads ) {

	#print $key;
	$sample_total_read += $count_sample_reads{$key};

	#print $key, "\t", $count_sample_reads{$key}, "\n";
	print REPORT "$key\t$count_sample_reads{$key}\n";
}

print REPORT "Total_reads_in_samples:\t$sample_total_read\n";
print REPORT "Total_unknown_reads:\t$unknown_reads\n";
print REPORT "Total_reads_in_file:\t$total_reads_file\n";

exit;
