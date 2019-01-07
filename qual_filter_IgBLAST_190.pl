#!/usr/bin/perl

use warnings;
use strict;
use IO::File;
use IO::Uncompress::Bunzip2 qw(bunzip2 $Bunzip2Error);
use IO::Compress::Bzip2 qw(bzip2 $Bzip2Error);

if (@ARGV != 2) {
	print "usage: $0 <input> <output>\n";
	print "where:\n\t<input> igblast, bzip2 compressed, outfm 19 (AIRR tab-delim)\n";
	print "\t<output> same format as input, with filtering applied to remove poor alignments\n";
	exit;
}

chomp (my $inputFile = $ARGV[0]);
chomp (my $outputFile = $ARGV[1]);

#open and uncompress the input file
my $input = new IO::File $inputFile, or die "There was a problem opening the input file: $inputFile\n$!\n";
my $in = new IO::Uncompress::Bunzip2 $input, or die "There was a problem uncompressing the input file: $inputFile\n$Bunzip2Error\n";

#prepare compressed output file
my $output = new IO::File ">$outputFile", or die "There was a problem opening the output file: $outputFile\n$!\n";
my $out = new IO::Compress::Bzip2 $output, or die "There was a problem compressing the output file: $outputFile\n$Bzip2Error\n";

#work through the input file
my $hdr = <$in>;
chomp $hdr;

#use the header to find the fields that are to be used for the filtering
my %fields = ();
my @colNames = split(/\t/, $hdr, -1);
for (my $i = 0 ; $i < scalar @colNames; $i++) {
	my $colName = $colNames[$i];

	$fields{$colName} = $i;
}

#use the header from the input file for the output file, nothing to change as not adding/removing any cols
print $out $hdr . "\n";

while (my $l = <$in>) {
	chomp $l;
	my @data = split(/\t/, $l, -1);

	#flag to track whether or not to include the current sequence in the output file
	my $keep = 1;

	#start with checking if the sequence contains stop codons, is in-frame and is productive
	#stop_codon	vj_in_frame	productive
	my $stopCodons = $data[($fields{stop_codon})];
	my $inFrame = $data[($fields{vj_in_frame})];
	my $productive = $data[($fields{productive})];
	if ($stopCodons =~ m/T/) {
		$keep = 0;
	} 
	if ($inFrame =~ m/F/) {
		$keep = 0;
	}
	if ($productive =~ m/F/) {
		$keep = 0;
	}

	#confirm there's at least a v and j call
	#v_call	d_call	j_call
	my $v = $data[($fields{v_call})];
	my $j = $data[($fields{j_call})];
	if ((length $v) == 0) {
		$keep = 0;
	}
	if ((length $j) == 0) {
		$keep = 0;
	}

	#calculate the length of the IGHV alignment and check it passes a minimum threshold
	#v_sequence_alignment
	my $vStr = $data[($fields{v_sequence_alignment})];
	my $vLen = length $vStr;
	if ($vLen < 180) {
		$keep = 0;
	}

	#check the percent identity for the IGHV alignment and check it passes minimum threshold
	#v_identity
	my $vPID = $data[($fields{v_identity})];
	if ((length $vPID) == 0) {
		$vPID = 0;
	}
	if ($vPID < 70) {
		$keep = 0;
	}
	#check that the CDR3 found meets a minimum length
	#cdr3_aa
	my $cdr3AA = $data[($fields{cdr3_aa})];
	my $cdr3Len = length $cdr3AA;
	if ($cdr3Len < 4) {
		$keep = 0;
	}

	#check the flag and output if it it still true
	if ($keep == 1) {
		print $out $l . "\n";
	}

}
#close input and output file
close $in;
close $out;

#done
exit;
