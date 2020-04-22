#!/bin/env perl

###
### ReferenceCoverageProfile.pl computes, from the ScaledCoverageProfile based on a collection of individual genomes, the reference value at each genomic bin.
###
### Copyright (C) 2014 by Gustavo Glusman, Institute for Systems Biology
### Contact: Gustavo@SystemsBiology.org
###
### This is free software; you can redistribute it and/or modify
### it under the same terms as Perl itself, either Perl version 5.8.8 or,
### at your option, any later version of Perl 5 you may have available.###
###

use strict;
$|=1;

my($scpfile, $outfile) = @ARGV;
die "
Usage:\t% ReferenceCoverageProfile.pl scpfile outfile
\t\tscpfile points at the ScaledCoverageProfile file.
\t\toutfile is where you want the output to go
e.g.:\t% ReferenceCoverageProfile.pl mySet.scp.gz mySet.rcp

" unless $scpfile && $outfile;

my $binsize = 1000;
my @targets = (0, .5, 1, 1.5, 2, 2.5, 3, 3.5, 4);

# find tabix
my $tabixDir = `which tabix`;
chomp $tabixDir;
die "
FATAL: ReferenceCoverageProfile.pl requires tabix
http://sourceforge.net/projects/samtools/files/tabix/

" unless $tabixDir;
$tabixDir =~ s/\/tabix$//;

my $targetsFile = $scpfile;
$targetsFile =~ s/gz$/targetSizes/;
my @targets;
open TARGETS, $targetsFile;
while (<TARGETS>) {
	chomp;
	my($bucket, $total) = split/\t/;
	$targets[$bucket] = $total;
}
close TARGETS;

open OUTF, ">$outfile";
print OUTF "#binsize\t$binsize\n";
print OUTF join("\t", "#targets", @targets), "\n";
print OUTF join("\t", "#CHROM", "START", "END", "GC", "RCP", "MAD"), "\n";

open SCP, "gunzip -c $scpfile |";
while (<SCP>) {
	next if /^#/;
	chomp;
	my($chrom, $bin, $bucket, $hist) = split /\t/;
	my @v = split /,/, $hist;
	next unless $#v>0;
	my ($median, $total) = medianFromHist(\@v);
	next unless $median;
	my ($mad1) = madFromHist(\@v, $median);
	
	my $best;
	my $bestd = 10000;
	my @bestHit;
	my $bestInitial;
	my $jump = $mad1/10;
	$jump = 1 if $jump<1;
	for (my $tm=int($median/2);$tm<$median*2.5;$tm+=$jump) {
		next unless $tm;
		my $totald;
		my @hit = ('-','-','-','-','-');
		foreach my $i (0..$#v) {
			next unless my $individuals = $v[$i];
			my $d = 100000;
			my $bestTarget;
			my $bestDev = 10000000;
			foreach my $t (0..$#targets) {
				my $dev = abs($tm*$targets[$t]-$i)/$tm/2;
				if ($dev<$bestDev) {
					$bestDev = $dev;
					$bestTarget = $t;
				}
			}
			$totald += $bestDev*$individuals;
			$hit[$bestTarget]+=$individuals;
		}
		
		my $initialTotal = $totald;
		my @stresses = ($totald);
		if ($hit[0]+$hit[1]) {
			my $delpop = $hit[0]+$hit[1]+$hit[2];
			my $f = (2*$hit[0]+$hit[1])/$delpop/2;
			my $expnull = $delpop*$f**2;
			my $exphemi = 2*$delpop*$f;
			$totald += ($stresses[1] = abs($expnull-$hit[0])/4);
			$totald += ($stresses[2] = abs($exphemi-$hit[1])/4);
		}
		if ($hit[3]+$hit[4]) {
			my $exppop = $hit[4]+$hit[3]+$hit[2];
			my $f = (2*$hit[4]+$hit[3])/$exppop/2;
			my $expfour = $exppop*$f**2;
			my $expthree = 2*$exppop*$f;
			$totald += ($stresses[3] = abs($expfour-$hit[4])/4);
			$totald += ($stresses[4] = abs($expthree-$hit[3])/4);
		}
		
		if ($totald<$bestd) {
			$bestd = $totald;
			$bestInitial = $initialTotal;
			$best = $tm;
			@bestHit = @hit;
		} elsif ($bestd && $totald/$bestd > 10) {
			last;
		}
	}
	my $prevRoundBest = $best;
	for (my $tm=$prevRoundBest-.9*$jump;$tm<$prevRoundBest+$jump;$tm+=.1) {
		my $totald;
		my @hit = ('-','-','-','-','-');
		foreach my $i (0..$#v) {
			next unless my $individuals = $v[$i];
			my $d = 100000;
			my $bestTarget;
			my $bestDev = 10000000;
			foreach my $t (0..$#targets) {
				my $target = $targets[$t];
				my $dev = abs($tm*$targets[$t]-$i)/$tm/2;
				if ($dev<$bestDev) {
					$bestDev = $dev;
					$bestTarget = $t;
				}
			}
			$totald += $bestDev*$individuals;
			$hit[$bestTarget]+=$individuals;
		}
		
		my $initialTotal = $totald;
		my @stresses = ($totald);
		if ($hit[0]+$hit[1]) {
			my $delpop = $hit[0]+$hit[1]+$hit[2];
			my $f = (2*$hit[0]+$hit[1])/$delpop/2;
			my $expnull = $delpop*$f**2;
			my $exphemi = 2*$delpop*$f;
			$totald += ($stresses[1] = abs($expnull-$hit[0])/4);
			$totald += ($stresses[2] = abs($exphemi-$hit[1])/4);
		}
		if ($hit[3]+$hit[4]) {
			my $exppop = $hit[4]+$hit[3]+$hit[2];
			my $f = (2*$hit[4]+$hit[3])/$exppop/2;
			my $expfour = $exppop*$f**2;
			my $expthree = 2*$exppop*$f;
			$totald += ($stresses[3] = abs($expfour-$hit[4])/4);
			$totald += ($stresses[4] = abs($expthree-$hit[3])/4);
		}
		
		if ($totald<$bestd) {
			$bestd = $totald;
			$bestInitial = $initialTotal;
			$best = $tm;
			@bestHit = @hit;
		} elsif ($bestd && $totald/$bestd > 10) {
			last;
		}
	}
	my ($mad2) = madFromHist(\@v, $best);
	$best ||= 1;
	$median ||= 1;
	print OUTF join("\t", 
		$chrom, $bin*$binsize-$binsize, $bin*$binsize,
		$bucket,
		sprintf("%.1f", $best), sprintf("%.2f", 100*$mad2/$best),
		), "\n";
}
close SCP;
close OUTF;
`$tabixDir/bgzip $outfile -f`;
`$tabixDir/tabix $outfile.gz -s 1 -b 2 -e 3 -f`;




sub medianFromHist {
	my($h) = @_;
	
	my $total;
	foreach my $v (@{$h}) {
		$total += $v;
	}
	$total /= 2;
	my $cumul;
	foreach my $i (0..$#{$h}) {
		$cumul += $h->[$i];
		return $i, $total if $cumul>=$total;
	}
	die "failed to identify median from histogram ".join(",", @{$h})."\n";
}

sub madFromHist {
	my($h, $median) = @_;
	my @dev;
	foreach my $i (0..$#{$h}) {
		$dev[abs($i-$median)] += $h->[$i];
	}
	return medianFromHist(\@dev);
}

