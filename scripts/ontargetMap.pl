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

	$regional{"$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]"}+=$tmp[6];

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
	print "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t$regional{$regions}\n";
}



#chr1	813014	839814	pdwindow_2;piece_5	0.333333	1	0.00124378
#chr1	813014	839814	pdwindow_2;piece_5	0.166667	4	0.00248757
#chr1	813014	839814	pdwindow_2;piece_5	0.5	2	0.00373134
#chr1	813014	839814	pdwindow_2;piece_5	1	3321	12.3918
#chr1	813014	839814	pdwindow_2;piece_5	0.166667	1	0.000621892
#chr1	813014	839814	pdwindow_2;piece_5	0.125	2	0.000932836
#chr1	813014	839814	pdwindow_2;piece_5	0.25	4	0.00373134
#chr1	813014	839814	pdwindow_2;piece_5	1	1917	7.15299

#chr1	851756	851845	pdwindow_2;piece_6	1	89	100

#chr1	856796	930946	pdwindow_2;piece_7	1	6696	9.03034
#chr1	856796	930946	pdwindow_2;piece_7	0.5	7	0.00472016
#chr1	856796	930946	pdwindow_2;piece_7	0.125	1	0.000168577
#chr1	856796	930946	pdwindow_2;piece_7	1	8	0.0107889
#chr1	856796	930946	pdwindow_2;piece_7	0.333333	3	0.00134862
#chr1	856796	930946	pdwindow_2;piece_7	0.047619	2	0.00012844
