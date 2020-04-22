#!/bin/env perl

#This script computes the coverage level distribution from a set of HMM-segmented individual genome assemblies.
#
# It expects as parameters the fileglob of segmentCoverage output files.
#
# It produces to the standard output, to be saved in the file of your choice, the map of CNV frequencies.
# Out columns:
#  CHROM: chromosome
#  BINSTART: base-1 first bin in segment (binsize depends on input, typically 1 kb)
#  BINEND: base-1 last bin in segment
#  SKEW: -1 indicates absent in everybody, 0 means diploid in everybody (on average), 1 would mean tetraploid in everybody
#  COUNTS: five columns indicating how many genomes had 0, 1, 2, 3, 4+ copies of segment.


my(%level, %total, $ngenomes);

foreach my $file (@ARGV) {
	if ($file =~ /\.gz$/) {
		open SEGOUT, "gunzip -c $file |";
	} else {
		open SEGOUT, $file;
	}
	open(SEGOUT, "gunzip -c $segout.gz |");
	while (<SEGOUT>) {
		next if /^#/;
		chomp;
		my($chrom, $start, $end, $state) = split /\t/;
		$start--;
		$level{$chrom}{$start}[$state]++;
		$level{$chrom}{$end}[$state]--;
		$total{$chrom}{$start} += $state;
		$total{$chrom}{$end} -= $state;
	}
	close SEGOUT;
	$ngenomes++;
}

print join("\t", "#assemblies", $ngenomes), "\n";
print join("\t", qw/#CHROM BINSTART BINEND SKEW COUNTS/), "\n";
foreach my $chrom (sort keys %level) {
	my($c, @c);
	my @pos = sort {$a<=>$b} keys %{$level{$chrom}};
	my $chromEnd = $#pos;
	foreach my $i (0..$#pos-1) {
		my $pos = $pos[$i];
		$c += $total{$chrom}{$pos};
		foreach my $state (0..$#{$level{$chrom}{$pos}}) {
			$c[$state] += $level{$chrom}{$pos}[$state];
		}
		
		print join("\t", $chrom, $pos+1, $pos[$i+1], sprintf("%.4f", $c/$ngenomes/2-1), @c), "\n";
	}
}

