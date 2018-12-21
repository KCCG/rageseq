#!/usr/bin/perl

use strict;
use warnings;
use IO::File;
use IO::Uncompress::Bunzip2 qw(bunzip2 $Bunzip2Error);

if (@ARGV != 2) {
	print "usage: $0 <input> <output>\n";
	print "where:\n\t<input> bzip2 compressed, igblast 1.9.0 output in AIRR tab-delim (outfmt 19)\n";
	print "\t<output> fasta formatted file with the CDR3s intended for clustering to call clonal relationships, description lines of FASTA tracking id, v, j usage\n";
	exit;
}

chomp (my $inputFile = $ARGV[0]);
chomp (my $outputFile = $ARGV[1]);

#open the input file, prepare the output file
my $input = new IO::File $inputFile, or die "There was a problem opening the input file: $inputFile\n$!\n";
my $in = new IO::Uncompress::Bunzip2 $input, or die "There was a problem uncompressing the input file: $inputFile\n$Bunzip2Error\n";

open (OUT, ">$outputFile"), or die "There was a problem opening the output file: $outputFile\n$!\n";

#work through the input file exctracting the details for the description line and the CDR3 sequence
#header
my $hdr = <$in>;
chomp $hdr;
my %fields = ();
my @colNames = split(/\t/, $hdr, -1);
for (my $i = 0; $i < scalar @colNames; $i++) {
	my $colName = $colNames[$i];

	$fields{$colName} = $i;
}

#output is FASTA format therefore no header needed for the output file
while (my $l = <$in>) {
	chomp $l;

	my @data = split(/\t/, $l, -1);

	my $id = $data[($fields{sequence_id})];
	my $v = $data[($fields{v_call})];
	#remove the allele from the IGHV
	if ($v =~ m/\*/) {
		($v) = split(/\*/, $v, -1);
	}
	my $j = $data[($fields{j_call})];
	#remove the allele from the IGHJ
	if ($j =~ m/\*/) {
		($j) = split(/\*/, $j, -1);
	}
	my $cdr3Nt = $data[($fields{cdr3})];
	$cdr3Nt = lc $cdr3Nt;

	#format the data fields into the FASTA output
	my $outStr = ">" . $id . " " . $v . " " . $j . "\n";
	$outStr .= $cdr3Nt . "\n";

	print OUT $outStr;
}
#close the input/output files
close $in;
close OUT;

#done
exit;