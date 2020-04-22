#!/bin/env perl

use strict;
use warnings;

my $infile = $ARGV[0];
my $outfile= $ARGV[1];
my $chromo = $ARGV[2];

if (@ARGV < 2) {
	print "Usage: $0 <coverages.bed> <outfile>\n"; exit;
}

my $ext = "cov";
my $overext = "over";
my $binsize = 1000;
my @encoded = (1..8, 12, 14..47, 49..255);
my @encodedChr = map chr($encoded[$_]), (0..249);
my($total, %total, $auto);
my $switch = 200;
my @recoded;
my @decoded = ("null", 0..7, "tab", "newline", "verttab", 8, "carriagereturn", 9..42, "zero", 43..249);

# find tabix
my $tabixDir = `which tabix`;
chomp $tabixDir;
die "
FATAL: condenseCoverage.pl requires tabix
http://sourceforge.net/projects/samtools/files/tabix/

" unless $tabixDir;
$tabixDir =~ s/\/tabix$//;

my $maxrecoded = (250-$switch)**2+$switch;
for (my $i=$switch;$i<=$maxrecoded;$i++) {
	$recoded[$i] = int(sqrt($i-$switch)+$switch+0.5+1e-6);
}

if (-e "$outfile.$ext" && "$outfile.$overext") {
	#
}

if ($infile =~/.gz/) {
	open (IN, "gunzip -c $infile |");
}
else {
	open (IN, "<", $infile);
}
open (OUT, ">", "$outfile.$ext");
open OVERFLOW, ">", "$outfile.$overext";

while (my $line=<IN>) {
    chomp $line;
    my @tmp = split (/\t/, $line);

	if ($chromo) {
		if ($tmp[0]!~/chr/) {
			$chromo =~s/chr//;
		}
		if ($tmp[0] ne $chromo) {
			next;
		} 

		# Skipping chrX and chrY
		if ($chromo eq "chrX" || $chromo eq "chrY" || $chromo eq "chr23" || $chromo eq "chr24" ) {
			next;
		}
	}

    my $counts = $tmp[3];
    my $total_bases = $counts*101;
    my $ref_len = $tmp[2]-$tmp[1];

    my $coverage = $total_bases/$ref_len;
	if ($coverage>$switch) {
		my $recoded = ($recoded[$coverage] // (sqrt($coverage-$switch)+$switch));
		if ($recoded>=250) {
            print OVERFLOW "$tmp[0]\t$tmp[1]\t$coverage\n";
			$coverage = 249;
		} else {
			$coverage = $recoded;
		}
	}

    my $enc = $encodedChr[$coverage];
	my $cov = $decoded[ord($enc)];
    print OUT "$tmp[0]\t$tmp[1]\t$enc\n";

}

`$tabixDir/bgzip $outfile.cov -f`;
`$tabixDir/tabix $outfile.cov.gz -s 1 -b 2 -e 2 -f`;
`$tabixDir/bgzip $outfile.$overext -f`;
`$tabixDir/tabix $outfile.$overext.gz -s 1 -b 2 -e 2 -f`;