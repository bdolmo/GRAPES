#!/bin/env perl

###
### ScaledCoverageProfile.pl computes a coverage profile from a collection of individual genomes.
### For each genome, it expects to find the condensedCoverage (extensions .cov and .over) and the totalCoverageByGC statistics (extension .totalByGC).
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

my($filelist, $outfile) = @ARGV;
die "
Usage:\t% ScaledCoverageProfile.pl filelist outfile
\t\tfilelist is a file enumerating, one per line, the genome assemblies to be collected - including paths (absolute or relative).
\t\toutfile is where you want the output to go
e.g.:\t% ScaledCoverageProfile.pl mySet.list mySet.scp

" unless $filelist && $outfile;

my $frz = 'hg19';
my $chromSizesFile = "chromInfo.$frz.gz";
my $parFile = "par.$frz";
my $GCbucketsFile = "GCbuckets.$frz.gz";
my $buckets = 25;
my $ws = 20;
my $binsize = 1000;
my $chunksize = 1e6;
my $covext = "cov";
my $overext = "over";
my $tbgext = "totalByGC";

# find tabix
my $tabixDir = `which tabix`;
chomp $tabixDir;
die "
FATAL: ScaledCoverageProfile.pl requires tabix
http://sourceforge.net/projects/samtools/files/tabix/

" unless $tabixDir;
$tabixDir =~ s/\/tabix$//;

my $chromSizes = readChromSizes($chromSizesFile);
my $par = readPar($parFile);
#get GC profile for reference genome, in 1 kb bins, in buckets
my $GCbuckets = readGCbuckets($GCbucketsFile);

my(%doubleStart, %doubleEnd);
foreach my $chrom (qw/chrX chrY/) {
	$doubleStart{$chrom} = int($par->{$chrom}{'PAR1'}{'end'}/$binsize);
	$doubleEnd{$chrom}   = int($par->{$chrom}{'PAR2'}{'start'}/$binsize);
}

my @assemblies;
open F, $filelist;
while (<F>) {
	next if /^#/;
	chomp;
	my($asm) = split /\t/;
	push @assemblies, $asm;
}
close F;

my $ngenomes = scalar @assemblies;
print "$ngenomes assemblies\n";

#prepare coverage for each assembly - and enumerate chromosomes observed in each
my(%asmObj, %logs, @totalLog, %obsChr, %yPresent);
foreach my $assembly (@assemblies) {
	my $tot = loadTotalByGC("$assembly.$tbgext");
	
	if (defined $tot->{'chrY'}) {
		my $total;
		foreach $_ (@{$tot->{'chrY'}}) {
			my($low) = split /,/;
			$total += $low;
		}
		$yPresent{$assembly} = 1 if $total>3000000; #hard-coded cutoff based on observed distributions...
	}
	while (my($key) = each %$tot) {
		$obsChr{$key}{$assembly}++ if $key =~ /^chr/;
	}
	
	my $parts = $tot->{'autosomes'};
	foreach my $bucket (0..$buckets-1) {
		@_ = split /,/, $parts->[$bucket];
		unless ($_[0]+$_[1]) {
			die join("\t", "processing:", $assembly, $bucket, @_), "\n";
		}
		
		my $l = log($_[0]+$_[1]);
		$logs{$assembly}[$bucket] = $l;
		$totalLog[$bucket] += $l;
	}
}

#compute bucket-wise target sizes for scaling assemblies
my @avgLog;
open TARGET, ">$outfile.targetSizes";
foreach my $bucket (0..$buckets-1) {
	$avgLog[$bucket] = $totalLog[$bucket] / $ngenomes;
	print TARGET $bucket, "\t", exp($avgLog[$bucket]), "\n";
}
close TARGET;
	
#compute bucket-wise scaling factor for each assembly
my %scaling;
foreach my $assembly (@assemblies) {
	foreach my $bucket (0..$buckets-1) {
		$scaling{$assembly}[$bucket] = exp($logs{$assembly}[$bucket]-$avgLog[$bucket]);
	}
}
	
#ready to start processing!
open OUTF, ">$outfile";
print OUTF join("\t", "#binsize", $binsize), "\n";
print OUTF join("\t", qw/#chrom bin bucket SCP/), "\n";
foreach my $chrom (sort keys %obsChr) {
	#print join("\t", $chrom, scalar keys %{$obsChr{$chrom}}), "\n";
	#next unless $chrom eq 'chrX';
	my $chromSize = $chromSizes->{$chrom};
	foreach my $meg (0..int($chromSize/$chunksize)) {
		my $start = $meg*$chunksize+1;
		my $end = ($meg+1)*$chunksize;
		$end = $chromSize if $chromSize<$end;
		print "$chrom:$start-$end ";
		
		my %obs;
		#read coverages from all assemblies, and scale
		foreach my $g (0..$#assemblies) {
			my $assembly = $assemblies[$g];
			unless (defined $obsChr{$chrom}{$assembly}) {
				print "missing $chrom for $assembly\n";
			}
			
			print "." unless ($g+1) % 10;
			print " " unless ($g+1) % 100;
			my($asCov, $lastBin) = coverage($assembly, $chrom, $start, $end, $binsize);
			
			while (my($bin, $cov) = each %{$asCov}) {
				$cov *= 2 if $yPresent{$assembly}
				&& ($chrom eq 'chrX' || $chrom eq 'chrY')
				&& $bin>=$doubleStart{$chrom} && $bin<=$doubleEnd{$chrom};
				$obs{$bin}[$g] = $cov;
			}
		}
		print "\n";
		
		foreach my $bin (sort {$a<=>$b} keys %obs) {
			next unless defined $obs{$bin};
			my $bucket = $GCbuckets->{$chrom}[$bin];
			my @hist;
			
			foreach my $g (0..$#assemblies) {
				next unless defined $obs{$bin}[$g];
				next if $chrom eq 'chrY' && !$yPresent{$assemblies[$g]};
				my $v = $obs{$bin}[$g] / $scaling{$assemblies[$g]}[$bucket];
				
				$hist[$v+.5]++; ## dithering may make more sense, but using integer counts is more convenient
			}
			print OUTF join("\t", $chrom, $bin, $bucket, join(",", @hist)), "\n";
		}
	}
}
close OUTF;
`$tabixDir/bgzip $outfile -f`;
`$tabixDir/tabix $outfile.gz -s 1 -b 2 -e 2 -f`;

#############
sub coverage {
	my($asm, $chrom, $start, $end) = @_;
	my $rescale = $binsize/$ws;
	
	my @decoded = ("null", 0..7, "tab", "newline", "verttab", 8, "carriagereturn", 9..42, "zero", 43..249);
	my $switch = 200;

	#default to whole chromosome
	$start ||= 1;
	$end ||= 3e8;
	
	my $startbin = sprintf("%.0f", ($start-1)/$binsize+1);
	my $endbin = sprintf("%.0f", ($end-1)/$binsize+1);
	
	#read overflows
	my %overflow;
	my $qstart = int(($start-1)/$ws+1); #base 1
	my $qend = int(($end-1)/$ws+1); #base 1
	open OVER, "$tabixDir/tabix $asm.$overext.gz $chrom:$qstart-$qend |";
	while (<OVER>) {
		chomp;
		my($chrom, $pos, $cov) = split /\t/; #pos base 1
		$overflow{$chrom}{$pos} = $cov;
	}
	close OVER;
	
	#read coverage and integrate
	my(%summed);
	my $lastBin;
	open COV, "$tabixDir/tabix $asm.$covext.gz $chrom:$qstart-$qend |";
	while (my $line = <COV>) {
		chomp $line;
		my(undef, $bin, $map) = split /\t/, $line;
		$bin--;
		foreach my $i (0..length($map)-1) {
			my $pos = $bin*$binsize + $i*$ws+1; #pos base 1
			next if $pos<$start;
			last if $pos>$end;
			my $cov = $decoded[ord(substr($map, $i, 1))];
			if ($cov>$switch) {
				my $overflow = $overflow{$chrom}{$bin*$chunksize+$i+1};
				if ($overflow) {
					$cov = $overflow;
				} else {
					$cov = ($cov-$switch)**2+$switch;
				}
			}
			my $outputBin = sprintf("%.0f", ($bin*$chunksize+$i)/$rescale+1); #base 1
			$summed{$outputBin} += $cov/$rescale;
			$lastBin = $outputBin if $outputBin>$lastBin; #base 1
		}
	}
	close COV;
	return \%summed, $lastBin;	
}

sub readChromSizes {
	my($file) = @_;
	my %chromSize;
	open F, "gunzip -c $file |";
	while (<F>) {
		chomp;
		my($chrom, $size) = split /\t/;
		$chromSize{$chrom} = $size;
	}
	close F;
	return \%chromSize;
}
sub readPar {
	my($file) = @_;
	my %par;
	open F, $file;
	while (<F>) {
		next if /^#/;
		chomp;
		my($chrom, $start, $end, $name) = split /\t/;
		$par{$chrom}{$name}{'start'} = $start;
		$par{$chrom}{$name}{'end'} = $end;
	}
	close F;
	return \%par;
}
sub readGCbuckets {
	my($file) = @_;
	my %buckets;
	open F, "gunzip -c $file |";
	while (<F>) {
		chomp;
		my($chrom, $bin, $bucket) = split /\t/;
		$buckets{$chrom}[$bin] = $bucket;
	}
	close F;
	return \%buckets;
}
sub loadTotalByGC {
	my($file) = @_;
	my %tot;
	open F, $file;
	while (<F>) {
		chomp;
		my($field, @values) = split /\t/;
		$tot{$field} = \@values;
	}
	close F;
	return \%tot;
}

