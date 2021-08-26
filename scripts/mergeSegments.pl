#!/usr/bin/env perl

use strict;
use warnings;
use List::MoreUtils qw(first_index indexes );
use List::Util qw/first/;
use List::Util qw(sum);
use Statistics::Descriptive;
use List::MoreUtils qw(uniq);

my $input = $ARGV[0];
my $bed   = $ARGV[1];

# Aim: Concatenate adjacent CNVs that have not been segmented throught DNAcopy.
# Input: segmented CNVs and single-exon CNVs.
# Output: final CNV file

# Rules
# 1. Current exon has to be immediately adjacent to the previous affected exon
# 2. Current exon has to present the same kind of CNV (del or dup) at similar ratios
# 3. Current exon's signal to noise has to be equal or higher than 8.5
# 4. Current exon has to be located on the same chromosome that previous exon
# 5. Gap between current exon and previous exon has to be smaller than 1mb

sub printHelp {
    print "Usage: perl $0 <CNV_SAMPLE.BED> <BED_REGIONS>\n";
    exit;
}

if( !$input ) {
    print " ERROR: no input file was introduced\n";
    printHelp();
}
if (!$bed) {
    print " ERROR: no BED regions file was introduced\n";
    printHelp();
}

my $max_gap_size = 10e6;

my $tmpBed = `cat $bed | grep -v 'chr\tstart\t'`;
chomp $tmpBed;
my @arrBed = split (/\n/, $tmpBed);

my $tmpInput = `cat $input`;
chomp $tmpInput;
my @arrInput = split (/\n/, $tmpInput);

my %segments = ();
my %seen = ();
my $ratio_A  = '.';
my $s2n_A    = '.';
my $SNRC_A   = '.';
my $zscore_A = '.';
my $MAD_A    = '.';
my $GC_A     = '.';
my $MAP_A    = '.';

my $ratio_B  = '.';
my $s2n_B    = '.';
my $SNRC_B   = '.';
my $zscore_B = '.';
my $MAD_B    = '.';
my $GC_B     = '.';
my $MAP_B    = '.';

foreach my $call_A (@arrInput) {

    my @tmp = split (/\t/, $call_A);

    my ($chr_A, $start_A, $end_A, $regions_A, $cipos_A, $ciend_A, $svtype_A,
    $gene_A, $ratio_A, $SNR_A, $SNRC_A, $MAD_A, $zscore_A, $GC_A, $MAP_A) = parseCall($call_A);

    my $coord_A = "$chr_A\t$start_A\t$end_A";

    my @ratios = ();
    my @SNR    = ();
    my @SNRC   = ();
    my @zscores= ();
    my @genes  = ();
    my @MADS   = ();
    my $regions = $regions_A;
    my @GC     = ();
    my @MAP    = ();

    $seen{$coord_A}++;
    next if $seen{$coord_A} > 1;

    push @genes,  $gene_A;
    push @ratios, $ratio_A if $ratio_A ne '.';
    push @MADS, $MAD_A if $MAD_A ne '.';
    push @SNR, $SNR_A if $SNR_A ne '.';
    push @SNRC, $SNRC_A if $SNRC_A ne '.';
    push @zscores, $zscore_A if $zscore_A ne '.';
    push @GC, $GC_A if $GC_A ne '.';
    push @MAP, $MAP_A if $MAP_A ne '.';

    my $index_A = first_index { /$tmp[0]\t$tmp[1]\t/ } @arrBed;
    my $outEnd = $end_A;

    foreach my $call_B (@arrInput) {

        next if $call_A eq $call_B;
        my ($chr_B, $start_B, $end_B, $regions_B, $cipos_A, $ciend_A, $svtype_B,
        $gene_B, $ratio_B, $SNR_B, $SNRC_B, $MAD_B, $zscore_B, $GC_B, $MAP_B) = parseCall($call_B);
        my $gap_size = $start_B - $end_A;

        my $coord_B = "$chr_B\t$start_B\t$end_B";

        my $index_B = first_index { /$chr_B\t$start_B\t/ } @arrBed;
        my $idx_diff = abs($index_B-$index_A);

        if ($chr_A eq $chr_B && $gap_size < $max_gap_size && $gap_size > 0 && $svtype_A eq $svtype_B && $idx_diff < 4 ) {
            $outEnd = $end_B;
            $regions+= $regions_B;

            push @SNR, $SNR_B if $SNR_B ne '.';
            push @SNRC, $SNRC_B if $SNRC_B ne '.';
            push @genes, $gene_B;
            push @ratios, $ratio_B if $ratio_B ne '.';
            push @MADS, $MAD_B if $MAD_B ne '.';
            push @zscores, $zscore_B if $zscore_B ne '.';
            push @GC, $GC_B if $GC_B ne '.';
            push @MAP, $MAP_B if $MAP_B ne '.';

            $seen{$coord_B}++;
            $chr_A   = $chr_B;
            $end_A   = $end_B;
            $index_A = first_index { /$chr_B\t$start_B\t/ } @arrBed;
        }
        else {
            next;
        }
    }

    my $out_ratios = '.';
    if (@ratios > 1) {
        $out_ratios = meanArray (@ratios);
    }
    elsif (@ratios == 1) {
        $out_ratios = $ratio_B eq '.' ? $ratio_A : $ratio_B;
    }

    my $out_MAD = '.';
    if (@MADS > 1) {
        $out_MAD = meanArray (@MADS);
    }
    elsif (@MADS == 1) {
        $out_MAD = $MAD_B eq '.' ? $MAD_A : $MAD_B;
    }

    my $out_SNR = '.';
    if (@SNR > 1) {
        $out_SNR = meanArray (@SNR);
    }
    elsif (@SNR == 1) {
        $out_SNR = $SNR_A;
    }
    my $out_SNRC = '.';
    if (@SNRC > 1){
        $out_SNRC = meanArray(@SNRC);
    }
    elsif (@SNRC == 1) {
        $out_SNRC = $SNRC_A;
    }

    my $out_zscores = '.';
    if (@zscores > 1) {
        $out_zscores = meanArray (@zscores);
    }
    elsif (@zscores == 1) {
        $out_zscores = $zscore_A;
    }

    my $out_GC = '.';
    if (@GC > 1) {
        $out_GC = meanArray (@GC);
    }
    elsif (@GC == 1) {
        $out_GC = $GC_A;
    }

    my $out_MAP = '.';
    if (@MAP > 1) {
        $out_MAP = meanArray (@MAP);
    }
    elsif (@MAP == 1) {
        $out_MAP = $MAP_A;
    }

    my @uniq = uniq(@genes);
    my $gene_str = join(",", @uniq);
    my $size = $outEnd-$start_A;
	  my $copyNumber = int ($out_ratios*2+0.5);

    print "$chr_A\t$start_A\t$outEnd\tIMPRECISE;SVTYPE=$svtype_A;CIPOS=$cipos_A;CIEND=$ciend_A;SVLEN=$size;GC=$out_GC;MAP=$out_MAP;GENE=$gene_str;REGIONS=$regions;RRD=$out_ratios;MADRD=$out_MAD;CN=$copyNumber;SNR=$out_SNR;SNRC=$out_SNRC;ZSCORE=$out_zscores\n";
}

#############################
sub parseCall {

    my $line = shift;

    my @tmp   = split (/\t/, $line);
    my $chr   = $tmp[0];
    my $start = $tmp[1];
    my $end   = $tmp[2];

    my @tmpInfo   = split (";", $tmp[3]);
    my ($regions) = grep ($_=~/REGIONS=/, @tmpInfo);
    $regions=~s/REGIONS=//;

    my ($svtype)  = grep ($_=~/SVTYPE=/, @tmpInfo);
    $svtype =~s/SVTYPE=//;

    my ($cipos) = grep ($_=~/CIPOS=/, @tmpInfo);
    $cipos =~s/CIPOS=// if $cipos;
    $cipos = '.' if !$cipos;

    my ($ciend) = grep ($_=~/CIEND=/, @tmpInfo);
    $ciend =~s/CIEND=// if $ciend;
    $ciend = '.' if !$ciend;

    my ($gene) =  grep ($_=~/GENE=/, @tmpInfo);
    $gene =~s/GENE=//;

    my ($ratio) = grep ($_=~/RRD=/, @tmpInfo);
    $ratio =~s/RRD=//;
    $ratio = '.' if !$ratio;

    my ($MAD) = grep ($_=~/MADRD=/, @tmpInfo);
    $MAD =~s/MADRD=//;
    $MAD = '.' if !$MAD;

    my ($SNR) = grep ($_=~/SNR=/, @tmpInfo);
    $SNR =~s/SNR=//;
    $SNR = '.' if !$SNR;

    my ($SNRC) = grep ($_=~/SNRC=/, @tmpInfo);
    $SNRC =~s/SNRC=//;
    $SNRC = '.' if !$SNRC;

    my ($GC) = grep ($_=~/GC=/, @tmpInfo);
    $GC =~s/GC=//;
    $GC = '.' if !$GC;

    my ($MAP) = grep ($_=~/MAP=/, @tmpInfo);
    $MAP =~s/MAP=//;
    $MAP = '.' if !$MAP;

    my ($zscore) = grep ($_=~/ZSCORE=/, @tmpInfo);
    $zscore =~s/ZSCORE=//;
    $zscore = '.' if !$zscore;

    return ($chr, $start, $end, $regions, $cipos, $ciend, $svtype, $gene, $ratio, $SNR, $SNRC, $MAD, $zscore, $GC, $MAP);

}

sub meanArray {
 my @data = @_;
 my $stat = Statistics::Descriptive::Full->new();
 $stat -> add_data(@data);
 my $mean = $stat->mean();
 return sprintf "%.3f",$mean;
}
