#!/usr/bin/env perl

package VCF;

use strict;
use warnings;
use File::Basename;

# Given a BAM file and a BED file with SV records, generate a VCF
my @arrFields = ( "SVTYPE", "CNVR", "SVLEN", "MQ", "EVIDENCE", "KDIV", "GC", "BR", "ASBR", "PE", "PPE", "RRD", "MADRD", "PRD", "NSNV", "HSNVR");

sub generate {
    my $inputBed = shift;
    my $VCF = $inputBed;
    $VCF =~s/.filtered.bed/.vcf/;
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
        $Info{FILTER} =~s/FILTER=//;
        $Info{SVTYPE} =~s/SVTYPE=//;

        print VCF "$tmp[0]\t$tmp[1]\t.\tN\t<$Info{SVTYPE}>\t.\t$Info{FILTER}\t$Info{PRECISION};END=$tmp[2];". join(";", @Arr) . "\tGT:CN\t./.:.\n";
    
    }
    close IN;
    close VCF;

    #bgzip and index VCF
   my $cmd = "$::bgzip -c $VCF > $VCF.gz";
   system $cmd;

   $cmd = "$::tabix -p vcf $VCF.gz";
   system $cmd;

}

sub parse {
    my $line = shift;
    my @t = split (/\t/, $line);

    my @tmp = split (/;/, $t[3]);
    #IMPRECISE;SVTYPE=DEL;SVLEN=435;MAPQ=27;KDIV=0;GC=55;BREAKREADS=0;ASSEMBLED=0;PE=9;phredPE=99;ratioRD=0.994113;madRD=0.1;phredRD=0;numSNV=3;homPctSNV=0.666667;
    #IMPRECISE;FILTER=LowQual;SVTYPE=DUP;CNVR=no;SVLEN=5000;MQ=34;KDIV=0;GC=0.52;BR=0;ASBR=0;PE=0;PPE=0;RRD=1.54;MADRD=0.09;PRD=0;NSNV=0;percHomSNV=0;5pRRD=0.00;3pRRD=0.00;
    my %hash = (
        PRECISION => $tmp[0],
        FILTER => $tmp[1],
        SVTYPE => $tmp[2],
        CNVR  => $tmp[3],
        SVLEN  => $tmp[4],
        MQ     => $tmp[5],
        EVIDENCE => $tmp[6],
        KDIV   => $tmp[7],
        GC     => $tmp[8],
        BR     => $tmp[9],
        ASBR   => $tmp[10],
        PE     => $tmp[11],
        PPE    => $tmp[12],
        RRD    => $tmp[13],
        MADRD  => $tmp[14],
        PRD    => $tmp[15],
        NSNV   => $tmp[16],
        HSNVR => $tmp[17],
    );
    return %hash;
}

sub getContigs  {

    my $str = `$::samtools view -HS $::bam | $::grep 'SN:' | $::awk '{ gsub(\"SN:\", \"\"); gsub(\"LN:\", \"\");print \"##contig=<ID=\"\$2\",length=\"\$3\">\" }'`;
    chomp $str;
    return $str;
}

sub printHeader {

    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
    my $fileDate = $year . $mon . $mday;
    my $refGenome = basename ($::genome);

    my $contigs = getContigs();

    my $header = "##fileformat=VCFv4.2
##fileDate=$fileDate
##FILTER=<ID=PASS,Description=\"All filters passed\">
##reference=$refGenome
##source=GRAPES
##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">
##INFO=<ID=PRECISE,Number=0,Type=Flag,Description=\"Precise structural variation\">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">
##INFO=<ID=MQ,Number=1,Type=Integer,Description=\"Mean mapping quality of all variant-supporting reads\">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">
##INFO=<ID=CNVR,Number=1,Type=String,Description=\"Copy number variable retion\">
##INFO=<ID=EVIDENCE,Number=1,Type=String,Description=\"SV signals present\">
##INFO=<ID=BR,Number=1,Type=Integer,Description=\"Total number of break-reads\">
##INFO=<ID=ASBR,Number=1,Type=Integer,Description=\"Total assembled break-reads into contig/s\">
##INFO=<ID=KDIV,Number=1,Type=Float,Description=\"Kmer diversity of assembled break-reads\">
##INFO=<ID=PE,Number=1,Type=Integer,Description=\"Total discordant read pairs\">
##INFO=<ID=PPE,Number=1,Type=Integer,Description=\"Phred-scale p-value using Poisson distribution on clustered discordant pairs\">
##INFO=<ID=RRD,Number=1,Type=Float,Description=\"Read-depth ratio at flanking coordinates\">
##INFO=<ID=MADRD,Number=1,Type=Float,Description=\"Median Absoute Deviation of the segmented read-depth\">
##INFO=<ID=PRD,Number=1,Type=Integer,Description=\"Phred-scale p-value using Fisher's exact test on flanking read-depth\">
##INFO=<ID=CN,Number=1,Type=Integer,Description=\"Copy number of segment containing breakend\">
##INFO=<ID=GC,Number=1,Type=Integer,Description=\"GC content\">
##INFO=<ID=NSNV,Number=1,Type=Float,Description=\"Total SNVs found within a deletion\">
##INFO=<ID=HSNVR,Number=1,Type=Float,Description=\"Ratio of homozygous SNVs within a deletion\">
##ALT=<ID=DEL,Description=\"Deletion\">
##ALT=<ID=DUP,Description=\"Duplication\">
##ALT=<ID=INV,Description=\"Inversion\">
##ALT=<ID=INS,Description=\"Insertion\">
##FILTER=<ID=LowQual,Description=\"Filtered by low-quality\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Copy number genotype for imprecise events\">
$contigs
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$::outName\n";
    return $header;
}