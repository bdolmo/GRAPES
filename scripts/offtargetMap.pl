#!/usr/bin/env perl


use strict;
use warnings;

my %hash = ();

my $A;
my $B;
my $flag = 0;
my $sum;

my %regional = ();

while(my $line=<STDIN>) {
	chomp $line;


	my @tmp = split (/\t/, $line);

	$regional{"$tmp[0]\t$tmp[1]\t$tmp[2]"}+=$tmp[6];

	if ($flag == 0) {
		$A = join ("\t", @tmp[0..2]);
		$flag = 1;
		$sum+=$tmp[6];
		next;
	}
	if ($flag == 1) {
		$B = join ("\t", @tmp[0..2]);
		if ($A eq $B) {
			$sum+=$tmp[6];
		}
		else {
			my @tmp2 = split (/\t/, $A);
			#print "$A\t$tmp2[0]:$tmp2[1]-$tmp2[2]\t$sum\n";
			$sum = $tmp[6];
			$flag = 0;
		}
	}
}

foreach my $regions (sort keys %regional) {
	my @tmp = split (/\t/, $regions);
	#chr1	10000	16345	chr1:10000-16345	28.120215837
	$regional{$regions} = sprintf "%.3f",$regional{$regions};
	print "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[0]:$tmp[1]-$tmp[2]\t$regional{$regions}\n";
}
