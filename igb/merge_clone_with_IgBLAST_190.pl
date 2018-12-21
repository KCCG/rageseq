#!/usr/bin/perl

use strict;
use warnings;
use IO::File;
use IO::Uncompress::Bunzip2 qw(bunzip2 $Bunzip2Error);
use IO::Compress::Bzip2 qw(bzip2 $Bzip2Error);


if (@ARGV != 3) {
	print "usage: $0 <igblast> <clustering> <output>\n";
	print "where:\n\t<igblast> bzip2 compressed igblast tab delim with AIRR tab-delim format from IgBLAST 1.9.0 (outfmt 19)\n";
	print "\t<clustering> clustering detatils, tab-delim, parsed from cd-hit output, compressed\n";
	print "\t<output> igblast input with extra col, joined on sequence_id col, compressed\n";
	exit;
}

chomp (my $inputFile = $ARGV[0]);
chomp (my $clusterFile = $ARGV[1]);
chomp (my $outputFile = $ARGV[2]);

#hash the cluster details
my $input = new IO::File $clusterFile, or die "There was a problem opening the cluster file: $clusterFile\n$!\n";
my $in = new IO::Uncompress::Bunzip2 $input, or die "There was a problem uncompressing the cluster file: $clusterFile\n$Bunzip2Error\n";

my %clusters = ();

#no header for this file
while (my $l = <$in>) {
	chomp $l;

	(my $label, my $id, my $status) = split(/\t/, $l, -1);

	#hash it
	$clusters{$id} = $label;
}
#close the cluster details file
close $in;

#open the IgBLAST file and prepare the output file
$input = new IO::File $inputFile, or die "There was a problem opening the input file: $inputFile\n$!\n";
$in = new IO::Uncompress::Bunzip2 $input, or die "There was a problem uncompressing the input file: $inputFile\n$Bunzip2Error\n";

my $output = new IO::File ">$outputFile", or die "There was a problem opening the output file: $outputFile\n$!\n";
my $out = new IO::Compress::Bzip2 $output, or die "There was a problem compressing the output file: $outputFile\n$Bzip2Error\n";

#work through the input file, match the clone label from the hash (if there is one), output
#get the header
my $hdr = <$in>;
chomp $hdr;

#use the header from input file as base for output file header, adding an extra column
print $out "$hdr\tcloneLabel\n";

while (my $l = <$in>) {
	chomp $l;

	my @data = split(/\t/, $l, -1);

	my $seqId = $data[0];

	my $cloneLabel = "";
	if (exists $clusters{$seqId}) {
		$cloneLabel = $clusters{$seqId};
	}

	print $out $l . "\t" . $cloneLabel . "\n";
}
#close the file
close $in;
close $out;

#done
exit;