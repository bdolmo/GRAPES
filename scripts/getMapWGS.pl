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

	$regional{"$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]"}{SUM}+=int($tmp[4]);
	$regional{"$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]"}{NUM}++;
}

foreach my $regions (sort keys %regional) {
	next if $regions eq "";
	my @tmp = split (/\t/, $regions);
	#chr1	10000	16345	chr1:10000-16345	28.120215837
	my @info = split (/;/, $tmp[3]);
	my ($CN) = grep ($_=~/CN=/, @info);
	$CN =~s/CN=//;
	my $ratio = 0;
	if ($CN > 0) {
		$ratio = $CN/2;
	}
	my $mapq = int ($regional{$regions}{SUM}/$regional{$regions}{NUM});

	print "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3];MAPQ=$mapq;KDIV=0;BREAKREADS=0;ASSEMBLED=0;PE=0;CSDISC=0;NINS=0;RDsupp=1;LOHsupp=.\n";

}

sub regla3 {
	my $num = shift;
	return int(($num*60)/100);
}