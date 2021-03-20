#!/usr/bin/env perl

package VCF;

use strict;
use warnings;
use File::Basename;


my @arrFields = ( "SVTYPE", "CIPOS", "CIEND", "SVLEN", "EV", "REGIONS", 
"GENE", "MAPQ", "KDIV", "GC", "MAP", "BR", "ASBR", "PE", "PPE", "RRD", 
"MADRD", "CN", "SNR", "SNRC", "ZSCORE", "PRD", "NSNV", "BAF", "CONF_SCORE");

sub generate {

    my $inputBed = shift;
    my $VCF = $inputBed;
    $VCF =~s/.bed/.vcf/;

    open (IN, "<", $inputBed) || die " ERROR: Unable to open $inputBed\n";
    open (VCF, ">", $VCF)     || die " ERROR: Unable to open $VCF\n";
    my $header = printHeader();
    print VCF "$header";

    while (my $line=<IN>) {
        chomp $line;

        my @tmp = split (/\t/, $line);
        my %Info = parse($line);
        my @Arr = ();
        foreach my $field (@arrFields) {
            push @Arr, $Info{$field};
        }
        #$Info{FILTER} =~s/FILTER=//;
        $Info{FILTER} = 'PASS'; # 

        $Info{SVTYPE} =~s/SVTYPE=//;
        print VCF "$tmp[0]\t$tmp[1]\t.\tN\t<$Info{SVTYPE}>\t.\t$Info{FILTER}\t$Info{PRECISION};END=$tmp[2];". join(";", @Arr) . "\tGT:CN\t./.:.\n";
    }
    close IN;
    close VCF;
}

###########################
sub parse {
    my $line = shift;

    my @t = split (/\t/, $line);
    my @tmp  = split (/;/, $t[3]);
    my $size = $t[2]-$t[1];

    my %hash = (
        PRECISION => $tmp[0],
        SVTYPE => "SVTYPE=" . fetchPattern("SVTYPE", \@tmp),
        CIPOS  => "CIPOS=" . fetchPattern("CIPOS", \@tmp),
        CIEND  => "CIEND=" . fetchPattern("CIEND", \@tmp),
        SVLEN  => "SVLEN=$size",
        EV     => "EV=RD",
        REGIONS=> "REGIONS=". fetchPattern("REGIONS", \@tmp),
        GENE   => "GENE=" . fetchPattern("GENE", \@tmp),
        MAPQ   => "MAPQ=" . fetchPattern("MAP", \@tmp),
        KDIV   => "KDIV=" . fetchPattern("KDIV", \@tmp),
        GC     => "GC=" . fetchPattern("GC", \@tmp),
        MAP    => "MAP=" . fetchPattern("MAP", \@tmp),
        BR     => "BR=" . fetchPattern("BREAKREADS", \@tmp),
        ASBR   => "ASBR=" . fetchPattern("ASSEMBLED", \@tmp),
        PE     => "PE=" . fetchPattern("PE", \@tmp),
        PPE    => "PPE=" . fetchPattern("PPE", \@tmp),
        RRD    => "RRD=" . fetchPattern("RRD", \@tmp),
        CN     => "CN=" . fetchPattern("CN", \@tmp),
        MADRD  => "MADRD=" . fetchPattern("MADRD", \@tmp),
        SNR    => "SNR=" . fetchPattern("SNR", \@tmp),
        SNRC   => "SNRC=" . fetchPattern("SNRC",\@tmp),
        ZSCORE => "ZSCORE=" . fetchPattern("ZSCORE", \@tmp),
        PRD    => "PRD=" . fetchPattern("PRD", \@tmp),
        NSNV   => "NSNV=" .fetchPattern("NSNV", \@tmp),
        BAF    => "BAF=" . fetchPattern("BAF", \@tmp),
        CONF_SCORE => "CONF_SCORE=". fetchPattern("CONF_SCORE",\@tmp)
    );
    return %hash;
}


###########################
sub fetchPattern {
    my $pattern = shift;
    my $arr  = shift;
    my @array = @$arr;
    #print "@array\n";

    my ($field) = grep ($_ =~/^$pattern=/, @array);
    $field = "$pattern=." if !$field;

    $field =~s/$pattern=//;
    return $field;
}

###########################
sub printHeader {

    my $refGenome = basename ($::genome);
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
    my $yyyymmdd = sprintf "%.4d%.2d%.2d", $year+1900, $mon+1, $mday;

    my $header = "##fileformat=VCFv4.3
##fileDate=$yyyymmdd
##reference=$refGenome
##source=GRAPES
##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">
##INFO=<ID=PRECISE,Number=0,Type=Flag,Description=\"Precise structural variation\">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">
##INFO=<ID=CIPOS,Number=1,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">
##INFO=<ID=CIEND,Number=1,Type=Integer,Description=\"Confidence interval around END for imprecise variants\">
##INFO=<ID=REGIONS,Number=1,Type=Integer,Description=\"Total number of affected regions\">
##INFO=<ID=GENE,Number=1,Type=Integer,Description=\"Gene/s affected at 5prime and 3prime\">
##INFO=<ID=MAPQ,Number=1,Type=Integer,Description=\"Mean mapping quality of variant supporting reads\">
##INFO=<ID=KDIV,Number=1,Type=Float,Description=\"Kmer diversity of assembled break-reads\">
##INFO=<ID=GC,Number=1,Type=Integer,Description=\"GC content\">
##INFO=<ID=MAP,Number=1,Type=Integer,Description=\"Mappability percentage extracted from GEM's mappability track\">
##INFO=<ID=BR,Number=1,Type=Integer,Description=\"Total number of break-reads\">
##INFO=<ID=ASBR,Number=1,Type=Integer,Description=\"Total assembled break-reads into contig/s\">
##INFO=<ID=PE,Number=1,Type=Integer,Description=\"Total discordant read pairs\">
##INFO=<ID=PPE,Number=1,Type=Integer,Description=\"Phred-scale p-value using Poisson distribution on clustered discordant pairs\">
##INFO=<ID=RRD,Number=1,Type=Float,Description=\"Read-depth ratio at flanking coordinates\">
##INFO=<ID=MADRD,Number=1,Type=Float,Description=\"Median Absoute Deviation of the segmented read-depth\">
##INFO=<ID=SNR,Number=1,Type=Float,Description=\"Signal to noise of the read depth ratios for the case sample\">
##INFO=<ID=SNRC,Number=1,Type=Float,Description=\"Signal to noise of the read depth ratios for all control samples\">
##INFO=<ID=ZSCORE,Number=1,Type=Float,Description=\"Number of std.devs the target is from the reference mean\">
##INFO=<ID=PRD,Number=1,Type=Integer,Description=\"Phred-scale p-value using Fisher's exact test on flanking read-depth\">
##INFO=<ID=CN,Number=1,Type=Integer,Description=\"Copy number of segment with breakend support\">
##INFO=<ID=NSNV,Number=1,Type=Float,Description=\"Total SNVs found within the SV\">
##INFO=<ID=BAF,Number=1,Type=Float,Description=\"Mean B-allele frequency of the SNVs found within the SV\">
##INFO=<ID=CONF_SCORE,Number=1,Type=Float,Description=\"Confidence score of the variant call (0-1)\">
##FILTER=<ID=PASS,Description=\"Variant passes filtering criteria\">
##FILTER=<ID=ProcessedPseudogene,Description=\"Variant probably comes from a processed pseudogene\">
##ALT=<ID=DEL,Description=\"Deletion\">
##ALT=<ID=DUP,Description=\"Duplication\">
##ALT=<ID=INV,Description=\"Inversion\">
##ALT=<ID=INS,Description=\"Insertion\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Copy number genotype for imprecise events\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tSAMPLE\n";
    return $header;
}
