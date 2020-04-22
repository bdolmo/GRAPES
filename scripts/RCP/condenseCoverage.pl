#!/bin/env perl

###
### condenseCoverage.pl transforms the depth of coverage information in a CGI REF directory or a BAM file into a condensed format.
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
my $source = shift @ARGV;
my $outfile = shift @ARGV;
die "
Usage:\t% condenseCoverage.pl input outfile
\t\tinput can be a Complete Genomics REF directory or a BAM file
\t\toutfile is where you want the output to go
e.g.:\t% condensedCoverage.pl myCGIdata/GS000012345/GS12345-DNA_A01/ASM/REF coverage/GS000012345
e.g.:\t% condensedCoverage.pl myIlluminaData/SS1234567/genome/bam/SS1234567.bam coverage/SS1234567

" unless $source && $outfile;

my %todo;
foreach (fulldirlist($source)) {
	$todo{$1} = $_ if /^coverageRefScore-chr(.+?)-/;
}
#print join("\t", "to do:", sort {$a<=>$b || $a cmp $b} keys %todo), "\n";

# find tabix
my $tabixDir = `which tabix`;
chomp $tabixDir;
die "
FATAL: condenseCoverage.pl requires tabix
http://sourceforge.net/projects/samtools/files/tabix/

" unless $tabixDir;
$tabixDir =~ s/\/tabix$//;

my $ext = "cov";
my $overext = "over";
my $ws = 20;
my $binsize = 1000;
my @encoded = (1..8, 12, 14..47, 49..255);
my @encodedChr = map chr($encoded[$_]), (0..249);
my $nullmap = chr(1) x $binsize;
my($total, %total, $auto);
my $switch = 200;
my @recoded;
my $maxrecoded = (250-$switch)**2+$switch;
for (my $i=$switch;$i<=$maxrecoded;$i++) {
	$recoded[$i] = int(sqrt($i-$switch)+$switch+0.5+1e-6);
}

foreach my $shit (@encodedChr) {
	print "$shit\n";
}


# processing a BAM file
	# find samtools
	my $samtools = `which samtools`;
	chomp $samtools;
	die "
	FATAL: condenseCoverage.pl requires samtools for processing BAM files
	http://samtools.sourceforge.net/
	
	" unless $samtools;
	
	my $prev_chrom;
	open RD, "$samtools depth $source 2> /dev/null |" || die "condensedCoverage can't open samtools depth pipe\n";
	
	open OUTF, ">$outfile.$ext";
	open OVERFLOW, ">$outfile.$overext";
	print "BAM file: $source\n";
	my($prevbin, $count, @values);
	my $prevw = -1;
	my $map = $nullmap;
	while (<RD>) {
		chomp;
		my($chrom, $pos, $gccov) = split /\t/;
		$chrom = "chr$chrom" if $chrom !~ /^chr/o;
		$chrom = "chrM" if $chrom eq "chrMT";
		
		if ($prev_chrom && $chrom ne $prev_chrom) {
			# finish off previous chromosome
			saveValue($prev_chrom, \$map, $prevw, \@values, $binsize);
			print OUTF join("\t", $prev_chrom, $prevbin+1, $map), "\n" if $count;
			print " $prev_chrom";
			
			# init new chromosome
			$map = $nullmap;
			$prevw = -1;
			$prevbin = $count = 0;
			@values = ();
		}
		$prev_chrom = $chrom;
		my $w = int($pos/$ws);
		if ($w>$prevw) {
			if ($prevw>=0 && scalar @values) {
				saveValue($chrom, \$map, $prevw, \@values, $binsize);
			}
			$prevw = $w;
			@values = ();
		}
		
		my $bin = int($w/$binsize);
		if ($bin>$prevbin) {
			print OUTF join("\t", $chrom, $prevbin+1, $map), "\n" if $count;
			$count = 0;
			$prevbin = $bin;
			$map = $nullmap;
		}
		
		push @values, $gccov;
		$count++;
	}
	
	close RD;
	#finish off last chromosome
	if ($prevw>=0 && scalar @values) {
		saveValue($prev_chrom, \$map, $prevw, \@values, $binsize);
	}
	print OUTF join("\t", $prev_chrom, $prevbin+1, $map), "\n" if $count;
	print " $prev_chrom\n";

`$tabixDir/bgzip $outfile.$ext -f`;
`$tabixDir/tabix $outfile.$ext.gz -s 1 -b 2 -e 2 -f`;
`$tabixDir/bgzip $outfile.$overext -f`;
`$tabixDir/tabix $outfile.$overext.gz -s 1 -b 2 -e 2 -f`;


sub saveValue {
	my($chrom, $mapref, $prevw, $values, $binsize) = @_;
	
	my $avg;
	foreach my $i (@{$values}) {
		$avg += $i;
	}
	$avg /= scalar @{$values};
	
	if ($avg>$switch) {
		my $recoded = ($recoded[$avg] // (sqrt($avg-$switch)+$switch));
		if ($recoded>=250) {
			print OVERFLOW join("\t", $chrom, $prevw+1, $avg), "\n";
			$avg = 249;
		} else {
			$avg = $recoded;
		}
	}
	substr(${$mapref}, $prevw % $binsize, 1) = $encodedChr[$avg];
}


sub fulldirlist {
	my($dir) = @_;
	opendir (DIR, $dir);
	my @files = grep /^[^.]/, readdir DIR;
	closedir DIR;
	return @files;
}

