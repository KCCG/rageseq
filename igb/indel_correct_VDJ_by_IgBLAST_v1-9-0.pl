#!usr/bin/perl

### usage example ###
# for processing IgBLAST v1.9.0 AIRR formatted output to revert indels to germline gene segment matches
# REQUIRED: 
# * an install of IgBLAST at /data/apps/ncbi-igblast-1.9.0/
# * BLAST databases (makeblastdb) for the germline reference dataasets at /data/apps/ncbi-igblast-1.9.0/database/
#
# if FASTA files contain mix of TCR and BCR sequences need to run on instance of IgBLAST for each locus 
# (example for processing all *.fa files in a directory):
#
## BCR
# export IGDATA=/data/apps/ncbi-igblast-1.9.0/
# for f in *.fa; do echo $f; /data/apps/ncbi-igblast-1.9.0/bin/igblastn \
# -germline_db_V /data/apps/ncbi-igblast-1.9.0/database/human_V_imgt.fa \
# -germline_db_D /data/apps/ncbi-igblast-1.9.0/database/human_D_imgt.fa \
# -germline_db_J /data/apps/ncbi-igblast-1.9.0/database/human_J_imgt.fa \
# -auxiliary_data /data/apps/ncbi-igblast-1.9.0/optional_file/human_gl.aux \
# -domain_system imgt -ig_seqtype Ig \
# -num_alignments_V 1 -num_alignments_D 1 -num_alignments_J 1 -organism human -outfmt 19 \
#-query $f -out ${f%.fa}.igblast.txt; done;
#
## TCR
# for f in *.fa; do echo $f; /data/apps/ncbi-igblast-1.9.0/bin/igblastn \
# -germline_db_V /data/apps/ncbi-igblast-1.9.0/database/human_tcr_v.fa \
# -germline_db_D /data/apps/ncbi-igblast-1.9.0/database/human_tcr_d.fa \
# -germline_db_J /data/apps/ncbi-igblast-1.9.0/database/human_tcr_j.fa \
# -auxiliary_data /data/apps/ncbi-igblast-1.9.0/optional_file/human_gl.aux \
# -domain_system imgt -ig_seqtype TCR \
# -num_alignments_V 1 -num_alignments_D 1 -num_alignments_J 1 -organism human -outfmt 19 \
# -query $f -out ${f%.fa}.tcr.igblast.txt; done;
#
## running the indel correction to generate new FASTA files, passing "sample" label is optional, if "" passed won't append the sample label to the IDs
# for f in *.igblast.txt; do echo $f; perl indel_correct_VDJ_by_IgBLAST_v1-9-0.pl $f ${f%.igblast.txt}l.indel_corr.fa "sample"
#
#####################


use warnings;
use strict;
use IO::File;
use IO::Uncompress::Bunzip2 qw(bunzip2 $Bunzip2Error);

if (@ARGV != 3) {
	print "usage: $0 <input> <output> <sample>\n";
	print "where:\n\t<input> bzip2 compressed IgBLAST output file, tab-delim\n";
	print "\t<output> fasta formatted version of VDJs that have been error corrected based on IgBLAST alignments\n";
	print "\t<sample> sample label to be added to ID\n";
	print "\tnote: this is for post-processing of IgBLAST v1.9.0 AIRR format output (-outfmt 19) with single IGHV, D and J alignments (-num_alignments 1)\n";
	exit;
}

chomp (my $inputFile = $ARGV[0]);
chomp (my $outputFile = $ARGV[1]);
chomp (my $sample = $ARGV[2]);

#open and uncompress the input file
my $input = new IO::File $inputFile, or die "There was a problem opening the input file: $inputFile\n$!\n";
my $in = new IO::Uncompress::Bunzip2 $input, or die "There was a problem uncompressing the input file: $inputFile\n$Bunzip2Error\n";

#open the output file
open (OUT, ">$outputFile"), or die "There was a problem opening the output file: $outputFile\n$!\n";

#get the header from the input file
my $hdr = <$in>;
chomp $hdr;

#use the header to get the position of the cols with the germline and query sequences
my $qSeqIndex = -1;
my $glSeqIndex = -1;
my $idIndex = -1;
my $vIndex = -1;
my $jIndex = -1;
my $prodIndex = -1;

my @colNames = split(/\t/, $hdr, -1);

for (my $i = 0; $i < scalar @colNames; $i++) {
	my $colName = $colNames[$i];

	if ($colName =~ m/^sequence_alignment$/) {
		$qSeqIndex = $i;
	} elsif ($colName =~ m/^germline_alignment$/) {
		$glSeqIndex = $i;
	} elsif ($colName =~ m/^sequence_id$/) {
		$idIndex = $i;
	} elsif ($colName =~ m/^v_call$/) {
		$vIndex = $i;
	} elsif ($colName =~ m/^j_call$/) {
		$jIndex = $i;
	} elsif ($colName =~ m/^productive$/) {
		$prodIndex = $i;
	} 
}
#check
if (($glSeqIndex == -1) || ($qSeqIndex == -1) || ($idIndex == -1) || ($vIndex == -1) || ($jIndex == -1) || ($prodIndex == -1)) {
	print "There was a problem finding one or more of the column indicies in the header.\n";
	exit;
}

#work through the input file
while (my $l = <$in>) {
	chomp $l;

	my @data = split(/\t/, $l, -1);

	#get the id
	my $id = $data[$idIndex];
	($id) = split(/\s/, $id, -1);

	my $v = $data[$vIndex];
	my $j = $data[$jIndex];
	my $prod = $data[$prodIndex];

	if ((length $v) == 0) {
		$v = "None";
	}
	if ((length $j) == 0) {
		$j = "None";
	}

	if (($v !~ m/None/) && ($j !~ m/None/)) { 
		#build the query and germline strings to correct the indels for
		my $gl = $data[$glSeqIndex];
		my $vdj = $data[$qSeqIndex];

		if ((length $gl) != (length $vdj)) {
			print $id . "\n";
			print "VDJ: $vdj\n";
			print "GL:  $gl\n";
			exit;
		}

		#check the the V-REGION as built form the fr1,2,3 and cdr1,2 is at least 270 nts
		if (((length $gl) > 250) && ((length $vdj) > 250)) {
			my @ntsVDJ = split(//, $vdj, -1);
			my @ntsGL = split(//, $gl, -1);

			my $correctedVDJ = "";
			for (my $i = 0; $i < scalar @ntsVDJ; $i++) {
				my $ntVDJ = $ntsVDJ[$i];
				my $ntGL = $ntsGL[$i];

				if ($ntGL =~ m/\-/) {
					#insertion relative the germline
					#remove the insertion, so don't add any nts to the corrected string at this position

				} elsif ($ntVDJ =~ m/\-/) {
					#deletion relative to the germline
					#re-instate the germline nt to correct the gap
					$correctedVDJ .= $ntGL;
				} else {
					#in all other cases (that is, not a gap, use the query sequence for the rearranged VDJ)
					$correctedVDJ .= $ntVDJ;
				}

			}

			#can now output the 'corrected' vdj string
			#make everything lowercase
			$correctedVDJ = lc($correctedVDJ);
			#remove the 'reversed' flag from the ID
			$id =~ s/^reversed\|//;
			if ((length $sample) > 0) {
				print OUT ">" . $id . "_" . $sample . "\n$correctedVDJ\n";
			} else {
				print OUT ">" . $id . "\n$correctedVDJ\n";
			}
		}
	}
}
#close files
close $in;
close OUT;

#done
exit;
