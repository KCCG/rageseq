#!/usr/bin/perl

use warnings;
use strict;
use IO::File;
use IO::Uncompress::Bunzip2 qw(bunzip2 $Bunzip2Error);

if (@ARGV != 2) {
	print "usage: $0 <igblast> <output>\n";
	print "where:\n\t<igblast> igblast, bzip2 compressed, clones added\n";
	print "\t<output> summary of the V, J and VJ usage in the input file at the read and clone level\n";
	exit;
}

chomp (my $inputFile = $ARGV[0]);
chomp (my $outputFile = $ARGV[1]);

#uncompress the input file
my $input = new IO::File $inputFile, or die "There was a problem opening the input file: $inputFile\n$!\n";
my $in = new IO::Uncompress::Bunzip2 $input, or die "There was a problem uncompressing the input file: $inputFile\n$Bunzip2Error\n";

#work through the input file
my $hdr = <$in>;
chomp $hdr;
my @colNames = split(/\t/, $hdr, -1);
my %fields = ();
for (my $i = 0; $i < scalar @colNames; $i++) {
	my $colName = $colNames[$i];
	$fields{$colName} = $i;
}
my %readCount = ();
my %cloneCount = ();
while (my $l = <$in>) {
	chomp $l;

	my @data = split(/\t/, $l, -1);

	#track the v, j and v+j usage for total reads and total clones
	my $clone = $data[($fields{cloneLabel})];
	my $v = $data[($fields{v_call})];
	my $j = $data[($fields{j_call})];
	my $locus = $data[($fields{locus})];
	my $id = $data[($fields{sequence_id})];

	#strip the allele off of the gene name
	if ($v =~ m/\*/) {
		($v) = split(/\*/, $v, -1);
	}
	if ($j =~ m/\*/) {
		($j) = split(/\*/, $j, -1);
	}

	#track the read level usage
	$readCount{$locus}{$v}{$j}{$id} = "";

	#track the clone level usage
	$cloneCount{$locus}{$v}{$j}{$clone} = "";
}
close $in;

#get the total number of reads and clones to calculate frequences
my %totals = ();
#reads
foreach my $locus (keys %readCount) {
	foreach my $v (keys %{$readCount{$locus}}) {
		foreach my $j (keys %{$readCount{$locus}{$v}}) {
			my @idLst = keys %{$readCount{$locus}{$v}{$j}};

			$totals{$locus}{reads} += scalar @idLst; 
		}
	}
}
#clones
foreach my $locus (keys %cloneCount) {
	foreach my $v (keys %{$cloneCount{$locus}}) {
		foreach my $j (keys %{$cloneCount{$locus}{$v}}) {
			my @idLst = keys %{$cloneCount{$locus}{$v}{$j}};

			$totals{$locus}{clones} += scalar @idLst; 
		}
	}
}

#prepare the output file
open (OUT, ">$outputFile"), or die "There was a problem opening the output file: $outputFile\n$!\n";
#make a header for the output file
print OUT "locus\ttarget\ttype\ttotalCount\tfreq\n";
foreach my $locus (keys %readCount) {
	#keep track of the js
	my %js = ();

	foreach my $v (keys %{$readCount{$locus}}) {
		#keep track of the totals for the v
		my $vTotal = 0;

		foreach my $j (keys %{$readCount{$locus}{$v}}) {
			my @idLst = keys %{$readCount{$locus}{$v}{$j}};
			my $count = scalar @idLst;

			$vTotal += $count; 
			$js{$j} += $count;

			my $freq = $count / ($totals{$locus}{reads});

			#output the count for the v+j
			print OUT "$locus\t$v" . "_" . $j . "\treads\t" . $count . "\t" . $freq . "\n";
		}

		my $vFreq = $vTotal / ($totals{$locus}{reads});
		#output the count for the v
		print OUT "$locus\t$v\treads\t" . $vTotal . "\t" . $vFreq . "\n";

	}

	foreach my $j (keys %js) {
		my $count = $js{$j};
		my $freq = $count / ($totals{$locus}{reads});
		#output the count for the v
		print OUT "$locus\t$j\treads\t" . $count . "\t" . $freq . "\n";
	}
}
foreach my $locus (keys %cloneCount) {
	#keep track of the js
	my %js = ();

	foreach my $v (keys %{$cloneCount{$locus}}) {
		#keep track of the totals for the v
		my $vTotal = 0;

		foreach my $j (keys %{$cloneCount{$locus}{$v}}) {
			my @idLst = keys %{$cloneCount{$locus}{$v}{$j}};
			my $count = scalar @idLst;

			$vTotal += $count; 
			$js{$j} += $count;

			my $freq = $count / ($totals{$locus}{clones});

			#output the count for the v+j
			print OUT "$locus\t$v" . "_" . $j . "\tclones\t" . $count . "\t" . $freq . "\n";
		}

		my $vFreq = $vTotal / ($totals{$locus}{clones});
		#output the count for the v
		print OUT "$locus\t$v\tclones\t" . $vTotal . "\t" . $vFreq . "\n";

	}

	foreach my $j (keys %js) {
		my $count = $js{$j};
		my $freq = $count / ($totals{$locus}{clones});
		#output the count for the v
		print OUT "$locus\t$j\tclones\t" . $count . "\t" . $freq . "\n";
	}
}
close OUT;

#done
exit;

