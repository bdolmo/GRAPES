#!/usr/bin/env perl


use strict;
use warnings;
use File::Basename;

my $breakPointCalls = $ARGV[0];
my $cnvCalls = $ARGV[1];

my $bedtools = `which bedtools`;
chomp $bedtools;

if (!$bedtools) {
    print " ERROR: bedtools was not found. Please add it on PATH\n";
    exit;
}

if ( (!-e $breakPointCalls || -z $breakPointCalls) && (-e $cnvCalls) ) {
    open (IN, "<", $cnvCalls) || die " ERROR: cannot open $cnvCalls\n";
    while (my $line =<IN>) {
        chomp $line;
        print "$line\n";
    }
    close IN;
}
elsif ( (!-e $cnvCalls || -z $cnvCalls) && (-e $breakPointCalls) ) {
    open (IN, "<", $breakPointCalls) || die " ERROR: cannot open $breakPointCalls\n";
    while (my $line =<IN>) {
        chomp $line;
        print "$line\n";
    }
    close IN;
}

if (-e $breakPointCalls && !-z $breakPointCalls && -e $cnvCalls && !-z $cnvCalls) {
 
    # Variants that do not intersect
    my $noStrA = `$bedtools intersect -a $cnvCalls -b $breakPointCalls -v`;
    chomp $noStrA;
    my @tmpNoStrA = split (/\n/, $noStrA);
    foreach my $var (@tmpNoStrA) {
	print "$var\n";
    }


    my $noStrB = `$bedtools intersect -b $cnvCalls -a $breakPointCalls -v`;
    chomp $noStrB;
    my @tmpNoStrB = split (/\n/, $noStrB);
    foreach my $var (@tmpNoStrB) {
	print "$var\n";
    }
    my $str = `$bedtools intersect -a $cnvCalls -b $breakPointCalls -f 0.2 -r -wao`;
    chomp $str;
    my @tmpStr = split (/\n/, $str);
   
    my $svtype_A;
    my $svtype_B;

    foreach my $variant ( @tmpStr ) {

        my @tmp = split (/\t/, $variant);
        if ($variant =~/\t-1\t/) {
            my $joinedLine = join ("\t", @tmp[0..3]);
            print "$joinedLine\n";
        }
        else {
            my $infoA = $tmp[3];
            my $infoB = $tmp[7];

            my @tmpA = split (/;/, $infoA);
            my @tmpB = split (/;/, $infoB);

            my @genesA   = grep ($_ =~/GENE=/, @tmpA); 
            my @genesB   = grep ($_ =~/GENE=/, @tmpB); 

            my @ratiosA  = grep ($_ =~/RRD=/, @tmpA);
            my @ratiosB  = grep ($_ =~/RRD=/, @tmpB);

            my @regionsA = grep ($_ =~/REGIONS=/, @tmpA);
            my @regionsB = grep ($_ =~/REGIONS=/, @tmpB);

            my ($ciposA) = grep ($_=~/CIPOS=/, @tmpA);
            $ciposA = "CIPOS=." if !$ciposA;

            my ($ciposB) = grep ($_=~/CIPOS=/, @tmpB);
            $ciposB = "CIPOS=." if !$ciposB;

            my ($ciendA) = grep ($_=~/CIEND=/, @tmpA);
            $ciendA = "CIEND=." if !$ciendA;

            my ($ciendB) = grep ($_=~/CIEND=/, @tmpB);
            $ciendB = "CIEND=." if !$ciendB;

            my @breakReadsA = grep ($_ =~/BREAKREADS=/, @tmpA);
            my @breakReadsB = grep ($_ =~/BREAKREADS=/, @tmpB);

            my @assembledA = grep ($_ =~/ASSEMBLED=/, @tmpA);
            my @assembledB = grep ($_ =~/ASSEMBLED=/, @tmpB);

            my @AF_A = grep ($_ =~/AF=/, @tmpA);
            my @AF_B = grep ($_ =~/AF=/, @tmpB);

            my @PE_A = grep ($_ =~/PE=/ && $_ !~ /SVTYPE=/, @tmpA);
            my @PE_B = grep ($_ =~/PE=/ && $_ !~ /SVTYPE=/, @tmpB);

            my @s2n_A = grep ($_ =~/SNR=/, @tmpA);
            my @s2n_B = grep ($_ =~/SNR=/, @tmpB);

            my @zscore_A = grep ($_ =~/ZSCORE=/, @tmpA);
            my @zscore_B = grep ($_ =~/ZSCORE=/, @tmpB);

            my @GC_A = grep ($_ =~/GC=/, @tmpA);
            my @GC_B = grep ($_ =~/GC=/, @tmpB);

            my @MAP_A = grep ($_ =~/MAP=/, @tmpA);
            my @MAP_B = grep ($_ =~/MAP=/, @tmpB);

            if ($infoA =~/=DEL/ && $infoB =~/=DEL/) {
                if ($infoA =~/PRECISE/ &&  $infoB =~/IMPRECISE/) {
                    print "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3];$ciposB;$ciendB;$genesB[0];$ratiosB[0];$regionsB[0];$s2n_B[0];$zscore_B[0];$GC_B[0];$MAP_B[0]\n";
                }
                if ($infoB =~/PRECISE/ &&  $infoA =~/IMPRECISE/) {
                    print "$tmp[4]\t$tmp[5]\t$tmp[6]\t$tmp[7];$ciposA;$ciendA;$genesA[0];$ratiosA[0];$regionsA[0];$s2n_A[0];$zscore_A[0];$GC_A[0];$MAP_A[0]\n";
                }
            }
            if ($infoA =~/=DUP/ && $infoB =~/=DUP/) {
                if ($infoA =~/PRECISE/ &&  $infoB =~/IMPRECISE/) {
                    print "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3];$ciposB;$ciendB;$genesB[0];$ratiosB[0];$regionsB[0];$s2n_B[0];$zscore_B[0];$GC_B[0];$MAP_B[0]\n";
                }
                if ($infoB =~/PRECISE/ &&  $infoA =~/IMPRECISE/) {
                    print "$tmp[4]\t$tmp[5]\t$tmp[6]\t$tmp[7];$ciposA;$ciendA;$genesA[0];$ratiosA[0];$regionsA[0];$s2n_A[0];$zscore_A[0];$GC_A[0];$MAP_A[0]\n";
                }
            }
            # Complex variant
            if ($infoA =~/=INV/ && ( $infoB =~/=DUP/ || $infoB =~/=DEL/ )) {
                if ($infoA =~/PRECISE/ &&  $infoB =~/IMPRECISE/) {
                    print "$tmp[0]\t$tmp[1]\t$tmp[2]\tPRECISE;SVTYPE=CX;$ciposA;$ciendA;$breakReadsA[0];$assembledA[0];$AF_A[0];$PE_A[0];$ratiosB[0];$regionsB[0];$s2n_B[0];$zscore_B[0]\n";
                }
            }
            if ($infoB =~/=INV/ && ( $infoA =~/=DUP/ || $infoA =~/=DEL/) ) {
                if ($infoB =~/PRECISE/ &&  $infoA =~/IMPRECISE/) {
                    print "$tmp[4]\t$tmp[5]\t$tmp[6]\tPRECISE;SVTYPE=CX;$ciposB;$ciendB;$breakReadsB[0];$assembledB[0];$AF_B[0];$PE_B[0];$ratiosA[0];$regionsA[0];$s2n_A[0];$zscore_A[0]\n";
                }
            }
        } 
    }
}

