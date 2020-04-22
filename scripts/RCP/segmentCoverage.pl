#!/bin/env perl

###
### segmentCoverage.pl uses HMMSeg to segment a genome's normalized coverage
###
### Copyright (C) 2014 by Gustavo Glusman, Institute for Systems Biology
### Contact: Gustavo@SystemsBiology.org
### 
### This is free software; you can redistribute it and/or modify
### it under the same terms as Perl itself, either Perl version 5.8.8 or,
### at your option, any later version of Perl 5 you may have available.
###

use strict;
use File::Basename;

our $dirname = dirname (__FILE__);

$|=1;

### some configuration may be required here
my $hmmseg = "$dirname/../../bin/HMMseg"; # change this to point to the *directory* where you have HMMSeg installed
chomp(my $java = `which java`); # change this if needed

my $input  = shift @ARGV;
my $outfile = shift @ARGV;
die "
Usage:\t% segmentCoverage.pl input [output]
\t\tinput is a genome's normalized coverage
\t\toutfile is where you want the output to go, defaults to location of input
e.g.:\t% segmentCoverage.pl coverage/GS000012345

" unless $input;
my $inext = "norm";
my $outext = "seg";
$input =~ s/\.$inext(\.gz)?$//;
$outfile ||= $input;

die "
FATAL: no $input.$inext.gz file found
Try running the command: bin/normalizeCoverage.pl $input

" unless -s "$input.$inext.gz";

# find tabix
my $tabixDir = `which tabix`;
chomp $tabixDir;
die "
FATAL: segmentCoverage.pl requires tabix
http://sourceforge.net/projects/samtools/files/tabix/

" unless $tabixDir;
$tabixDir =~ s/\/tabix$//;

# make sure HMMSeg is there
die "
FATAL: segmentCoverage.pl requires HMMSeg.jar, but it doesn't seem to be present in $hmmseg
http://noble.gs.washington.edu/proj/hmmseg/

" unless -e "$hmmseg/HMMSeg.jar";

# determine which chromosomes are present
open CHROMLIST, "$tabixDir/tabix -l $input.$inext.gz |";
my @chromlist = <CHROMLIST>;
close CHROMLIST;
chomp(@chromlist);

# create the HMM model file
my $model = "$outfile.hmm";
createModel($model);

# final preparations
my $binsize = 1000;
my $tmp = "$outfile.tmp";
$outfile =~ /\/(.+)$/;
open FL, ">$outfile.list";
print FL basename($1) . ".tmp";
close FL;

open OUTF, ">$outfile.$outext";
print OUTF "#binsize\t$binsize\n";
print OUTF join("\t", "#chrom", "start", "end", "state", "length (kb)"), "\n";

# segment each chromosome
foreach my $chrom (@chromlist) {

	# Skipping chrX and chrY
	if ($chrom eq "chrX" || $chrom eq "chrY" || $chrom eq "chr23" || $chrom eq "chr24" ) {
		next;
	}

	open TMP, ">$tmp";
	my $norm = readNormalizedCoverage("$input.$inext.gz", $chrom);
	foreach my $bin (1..$#{$norm}) {
		print TMP join("\t", $bin, $norm->[$bin] || 0), "\n";
	}
	close TMP;
	
	print " INFO: Segmenting $chrom\n";

	`$java -jar $hmmseg/HMMSeg.jar --model $model --column 2 $outfile.list >/dev/null 2>&1`;

	my($prevState, $end);
	my $bin = 1;
	my $start = 1;
	open VIT, "$tmp.vit";
	while (<VIT>) {
		chomp;
		if ($_ != $prevState) {
			print OUTF join("\t", $chrom, ($start-1)*$binsize, $end*$binsize, $prevState, $end-$start+1), "\n" if defined $prevState;
			$start = $bin;
		}
		$prevState = $_;
		$end = $bin;
		$bin++;
	}
	close VIT;
	print OUTF join("\t", $chrom, ($start-1)*$binsize, $end*$binsize, $prevState, $end-$start+1), "\n";
}
close OUTF;

# cleanup
unlink $tmp;
unlink "$outfile.list";
unlink "$tmp.vit";

`$tabixDir/bgzip $outfile.$outext -f`;
`$tabixDir/tabix $outfile.$outext.gz -s 1 -b 2 -e 3 -f`;



sub readNormalizedCoverage {
	my($file, $chrom) = @_;
	
	my $fix = sqrt(1.05);
	my @norm;
	open NORM, "$tabixDir/tabix $file $chrom |";
	while (<NORM>) {
		chomp;
		my($c, $start, $end, $norm) = split /\t/;
		next unless $c eq $chrom;
		$norm[sprintf("%.0f", $end/$binsize)] = $norm/$fix;
	}
	close NORM;
	return \@norm;
}

sub createModel {
	my($file) = @_;
	
	open MODEL, ">$file";
	print MODEL "HMMSEG 3
State Start
Emissions none
Transitions 0.5 0.0 0.5 0.0 0.0
State 0
Emissions 0 50
Transitions 0.99998 0.00001 0.00001 0.0 0.0
State 1
Emissions  50 50
Transitions 0.00001 0.99998 0.00001 0.0 0.0
State 2
Emissions  100 300
Transitions 0.00000002 0.00000008 0.9999998 0.00000008 0.00000002
State 3
Emissions  150 50
Transitions 0.0 0.0 0.00001 0.99998 0.00001
State 4
Emissions  200 100
Transitions 0.0 0.0 0.00001 0.00001 0.99998
";
	close MODEL;	
}


