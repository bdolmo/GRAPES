#!/bin/env perl

###
### normalizeCoverage.pl computes, from a genome's condensed coverage and a Reference Coverage Profile, the percentage of diploid coverage at each locus.
###
### Copyright (C) 2014 by Gustavo Glusman, Institute for Systems Biology
### Contact: Gustavo@SystemsBiology.org
### 
### This is free software; you can redistribute it and/or modify
### it under the same terms as Perl itself, either Perl version 5.8.8 or,
### at your option, any later version of Perl 5 you may have available.
###

use strict;
use Sort::Key::Natural qw(natsort);

$|=1;

my $input   = $ARGV[0];
my $rcpfile = $ARGV[1];
my $outDir  = $ARGV[2];
my $outName = $ARGV[3];
my $chromo  = $ARGV[4];

if (@ARGV < 2) {
	print " Usage: $0 coverage/GS000012345 RCPs/CGI-10.rcp.gz\n"; exit;
}

my $inext   = 'cov';
my $overext = 'over';
my $outext  = 'norm';
$input =~ s/\.$inext(\.gz)$//;

# find tabix
my $tabixDir = `which tabix`;
chomp $tabixDir;
die "
FATAL: normalizeCoverage.pl requires tabix
http://sourceforge.net/projects/samtools/files/tabix/

" unless $tabixDir;
$tabixDir =~ s/\/tabix$//;

die "
FATAL: no $input.totalByGC file found
Try running the command: bin/totalCoverageByGC.pl $input reference/GCprofile.hg19.gz

" unless -s "$input.totalByGC";

# read RCP
my %rcp;
my %bucket;
open RCP, "gunzip -c $rcpfile |";
$_ = <RCP>;
chomp;
my(undef, $binsize) = split /\t/;
$_ = <RCP>;
chomp;
my(undef, @targetSizes) = split /\t/;
$_ = <RCP>;
while (<RCP>) {
	chomp;
	my($chrom, $start, $end, $bucket, $ref, $mad) = split /\t/;
	my $bin = sprintf("%.0f", $end/$binsize);

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

	$rcp{$chrom}[$bin] = $ref;
	$bucket{$chrom}[$bin] = $bucket;
}
close RCP;

# determine scaling per GC bucket only on Autosomes
my @totalAuto;
open(TOTAL, "<$input.totalByGC") || die "Can't open $input.totalByGC\n";
while (my $line =<TOTAL>) {
	chomp $line;
	my($field, @values) = split /\t/, $line;
	next unless $field eq 'autosomes';

	@totalAuto = @values;
	last;
}
close TOTAL;
my @scalings;
foreach my $bucket (0..24) {
	my @tmp = split (/,/, $totalAuto[$bucket]);
	push @scalings, ($tmp[0]+$tmp[1]+$tmp[2])/$targetSizes[$bucket]; 
}


# determine which chromosomes are present
open CHROMLIST, "$tabixDir/tabix -l $input.$inext.gz |";
my @chromlist = <CHROMLIST>;
close CHROMLIST;
chomp(@chromlist);

if ($chromo) {
	my @chromlist = ("$chromo");
}

# process and output
open OUTF, ">$outDir/$outName.$outext";
print OUTF "#binsize\t$binsize\n";
print OUTF "#profile\t$rcpfile\n";
print OUTF join("\t", "#scalings", @scalings), "\n";

foreach my $chrom (@chromlist) {


		# Skipping chrX and chrY
		if ($chrom eq "chrX" || $chrom eq "chrY" || $chrom eq "chr23" || $chrom eq "chr24" ) {
			next;
		}

	print " INFO: Normalizing $chrom\n";
	my $chromCov = condensedCoverage($chrom, 0, 0, $binsize);
	foreach my $bin (natsort keys %$chromCov) { #base 1
		if (!$rcp{$chrom}[$bin]) {
			next;
		}

		my $ch = $rcp{$chrom}[$bin];
		next unless $ch>0;
		my $bincov = int ($chromCov->{$bin}->{COV}/$chromCov->{$bin}->{NUM});

		my $start = $bin*$binsize;
		my $end   = ($bin*$binsize)+$binsize;

		my $scaledcov = $bincov / $scalings[$bucket{$chrom}[$bin]];
		#print "$chrom\t$start\t$end\t$bincov\t$ch\n";

		print OUTF join("\t", $chrom, $bin*$binsize, ($bin*$binsize)+$binsize, sprintf("%.0f", 100*$scaledcov/$ch)), "\n";
	}
}
close OUTF;

`$tabixDir/bgzip $outDir/$outName.$outext -f`;
`$tabixDir/tabix $outDir/$outName.$outext.gz -s 1 -b 2 -e 3 -f`;


sub condensedCoverage {
	my($chrom, $start, $end, $reqbinsize) = @_;
	
	my @decoded = ("null", 0..7, "tab", "newline", "verttab", 8, "carriagereturn", 9..42, "zero", 43..249);
	my $switch = 200;
	my $ws = 20;
	$reqbinsize ||= $ws;
	my $rescale = $reqbinsize/$ws;
	
	#default to whole chromosome
	$start ||= 1;
	$end ||= 3e8;

	my %overflow;
	my $qstart = int(($start-1)/$ws+1); #base 1
	my $qend = int(($end-1)/$ws+1); #base 1

	open OVER, "$tabixDir/tabix $input.$overext.gz $chrom:$start-$end |";
	while (<OVER>) {
		chomp;
		my($chrom, $pos, $cov) = split /\t/; #pos base 1
		$overflow{$pos} = $cov;
	}
	close OVER;
	
	my $chunksize = 1000;
	my $binsize = $chunksize*$ws;
	my $startbin = sprintf("%.0f", ($start-1)/$binsize+1);
	my $endbin = sprintf("%.0f", ($end-1)/$binsize+1);

	my %summed;
	open COV, "$tabixDir/tabix $input.$inext.gz $chrom:$start-$end |";
	while (<COV>) {
		chomp;
		my($chromoSome, $bin, $map) = split /\t/;
		$bin--;

		my $pos = $bin; #pos base 1
		next if $pos<$start;
		last if $pos>$end;
		my $cov = $decoded[ord($map)];

		if ($cov > $switch) {
			my $overflow = $overflow{$bin*$chunksize+1};
			if ($overflow) {
				$cov = $overflow;
			} else {
				$cov = ($cov-$switch)**2+$switch;
			}
		}
		#my $outputBin = int(($bin*$chunksize)/$rescale+1+1e-9); #base 1
		#print "$chromoSome\t$bin\t$outputBin\n";
		my $rounded  = int($bin/1000+1)*1000;
		$rounded = sprintf("%.0f", $rounded/1000);
		#print "$bin\t$rounded\n";
		if (!$summed{$rounded}{NUM}) {
			$summed{$rounded}{NUM} = 0;
		}
		$summed{$rounded}{COV}+= $cov;
		$summed{$rounded}{NUM}++;		
	}
	close COV;
	
	return \%summed;
}

sub roundup
{
   my $number = shift;
   my $round = shift;

   my @tmp = split (//, $number);

   if ($tmp[-3] >5) {

   }
	else {
		
	}


}