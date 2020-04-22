#!/bin/env perl

###
### totalCoverageByGC.pl computes summary statistics from a genome's condensed coverage and the reference's GC profile
###
### Copyright (C) 2014 by Gustavo Glusman, Institute for Systems Biology
### Contact: Gustavo@SystemsBiology.org
### 
### This is free software; you can redistribute it and/or modify
### it under the same terms as Perl itself, either Perl version 5.8.8 or,
### at your option, any later version of Perl 5 you may have available.
###

use strict;
$|=1;
my $source  = shift @ARGV;
my $refGC   = shift @ARGV;
my $outDir  = shift @ARGV;
my $outName = shift @ARGV;
my $chromo  = shift @ARGV;

# read referenceGC buckets profile
my %refGC;
open GC, "gunzip -c $refGC |";
$_ = <GC>;
chomp;
my(undef, $binsize) = split /\t/;
die "Unexpected reference GC profile format: could not identify binsize\n" unless $binsize+1-1;

$_ = <GC>;
while (<GC>) {
	chomp;
	my($chrom, $bin, $bucket) = split /\t/;

	if ($chromo) {
		if ($chrom!~/chr/) {
			$chromo =~s/chr//;
		}
		if ($chrom ne $chromo) {
			next;
		} 
		# Skipping chrX and chrY
		if ($chromo eq "chrX" || $chromo eq "chrY" || $chromo eq "chr23" || $chromo eq "chr24" ) {
			next;
		}
	}

	$refGC{$chrom}[$bin] = $bucket;
}
close GC;
# read overflows
my %overflow;
open OVER, "gunzip -c $source.over.gz |";
while (<OVER>) {
	chomp;
	my($chrom, $pos, $cov) = split /\t/; #pos base 1
	if ($chromo) {
		if ($chrom!~/chr/) {
			$chromo =~s/chr//;
		}
		if ($chrom ne $chromo) {
			next;
		} 
	}	
	$overflow{$chrom}{$pos} = $cov;
}
close OVER;

# read and process condensed coverage
my @decoded = ("null", 0..7, "tab", "newline", "verttab", 8, "carriagereturn", 9..42, "zero", 43..249);
my $switch = 200;
my %gchist;
open COV, "gunzip -c $source.cov.gz |";
while (<COV>) {
	chomp;
	my($chrom, $bin, $map) = split /\t/;
	if ($chromo) {
		if ($chrom!~/chr/) {
			$chromo =~s/chr//;
		}
		if ($chrom ne $chromo) {
			next;
		} 
		# Skipping chrX and chrY
		if ($chromo eq "chrX" || $chromo eq "chrY" || $chromo eq "chr23" || $chromo eq "chr24" ) {
			next;
		}

	}

	my $pos    = $bin; #pos base 0
	my $posbin = int($pos/$binsize)+1; #posbin base 1
	my $bucket = $refGC{$chrom}[$posbin];

	next unless defined $bucket;
	my $cov = $decoded[ord($map)];
	#print "$chrom\t$pos\t$cov\t$posbin\t$bucket\n";

	my $type = 0;
	if ($cov>$switch) {
		my $overflow = $overflow{$chrom}{$bin*$binsize+1};
		if ($overflow) {
			$cov = $overflow;
			$type = 2;
		} else {
			$cov = ($cov-$switch)**2+$switch;
			$type = 1;
		}
	}
	$gchist{$chrom}[$bucket][$type] += $cov;
}
close COV;

# output
my(%line, %total);
my($total, $auto, @total, @auto);
foreach my $chrom (sort keys %gchist) {
	$line{$chrom} = $chrom;
	foreach my $i (0..$#{$gchist{$chrom}}) {
		$line{$chrom} .= "\t";
		if (defined $gchist{$chrom}[$i]) {
			foreach my $j (0..2) {
				my $v = $gchist{$chrom}[$i][$j]*2;

				if ($v) {
					$gchist{$chrom}[$i][$j]=$v
				}
				$total{'all'}[$i][$j] += $v;
				$total[$j] += $v;
				$total += $v;
				if ($chrom !~ /chr[XYM]/) {
					$auto += $v;
					$auto[$j] += $v;
					$total{'autosomes'}[$i][$j] += $v;
				}
			}
			$line{$chrom} .= join(",", @{$gchist{$chrom}[$i]});
		}
	}
}
foreach my $chrom ('all', 'autosomes') {
	$line{$chrom} = $chrom;
	foreach my $i (0..$#{$total{$chrom}}) {
		$line{$chrom} .= "\t";
		if (defined $total{$chrom}[$i]) {
			$line{$chrom} .= join(",", @{$total{$chrom}[$i]});
		}
	}
}

my $outfile .= "$outDir/$outName.totalByGC";
open OUTF, ">$outfile";
print OUTF join("\t", ("#total", $total, join(",", @total))), "\n";
print OUTF join("\t", ("#autosomes", $auto, join(",", @auto))), "\n";
	
foreach my $chrom ('all', 'autosomes', sort keys %gchist) {
	print OUTF $line{$chrom}, "\n";
}
close OUTF;

