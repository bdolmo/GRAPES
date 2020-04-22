#!/usr/bin/env perl

package HardFilter;

use strict;
use warnings;
use Sort::Key::Natural qw(natsort);

# Hard Filter in default mode:
# a) Size lower than 1Mb
# b) Clustered Discordant pairs must be significant
# c) DELs flanking ratio < 0.7
# d) DUPs flanking ratio > 1.3
# e) CNVs, flanking read-depth must be significant
# f) GC is lower than 20% or greater than 80%
# g) Number of assembled breakreads must be greater than 50%
# h) MAPQ greater than 10
# NB: Filtered SVs are not removed, just tagged as <LowQual>

my @arrFields = ( "FILTER", "SVTYPE","CNVR", "SVLEN", "MQ", "EV", "KDIV", "GC", "BR", "ASBR", "PE", "PPE", "RRD", "MADRD", "PRD", "NSNV", "HSNVR", "5pRRD", "3pRRD");
#my @arrFields = ( "FILTER", "SVTYPE","CNVR", "SVLEN", "MQ", "KDIV", "GC", "BR", "ASBR", "PE", "PPE", "RRD", "MADRD", "PRD", "NSNV", "percHomSNV", "5pRRD", "3pRRD", "inSegDup");

###################
sub annotSegDup { # Deprecated
    my $infile = shift;

    my $cmd = "$::bedtools intersect -a $infile -b $::centromeres -v | $::bedtools intersect -a stdin -b $::segDups > $::outDir/$::outName.annot.segdup.bed -loj";
    system $cmd;

    my %Hash = ();

    open (IN ,"<", "$::outDir/$::outName.annot.segdup.bed") || die " ERROR: Unable to open $::outDir/$::outName.annot.segdup.bed\n";
    while (my $line=<IN>) {
        chomp $line;
        next if $line =~/^#Chr/;
        my @tmp = split (/\t/, $line);
        my $inSegDup = "inSegDup";
#chr1	53564935	53566982	PRECISE	DEL	no	2047	28	0.31	56	48	0	1.0176	0.1	0	1	4	1	53.8086	0.744767	0.785025	.	-1	-1
        if ($tmp[22] eq '-1') {
            $inSegDup = "noSegDup";
        }
        $line = join ("\t", @tmp[0..20]);
        $line.= "\t$inSegDup";
        $Hash{$line}++;
    }
    close IN;
    return %Hash;
}

###################
sub apply {

    my $inputBed = shift;
    my $outputBed = $inputBed;

    $outputBed=~s/annFeatures/filtered/;
    open (OUT, ">", $outputBed) || die " ERROR: Unable to open $outputBed\n";
    open (IN, "$::bedtools intersect -a $inputBed -b $::centromeres -v |")   || die " ERROR: Unable to open $inputBed\n";

    #foreach my $line (natsort keys %Hash) {
    while (my $line=<IN>) {
        chomp $line;
        next if $line =~/^#Chr/;
        #Chr	Start	End	Precision	SVTYPE	SVLEN	MAPQ	Kmer_diversity	Breakreads	Assembled	Discordants	Ratio_RD	MAD_RD	Phred_discordant	Phred_RD	Total_SNV	Homozygous_ratio	GC_content
        my @tmp = split (/\t/, $line);
        my %Info = getInfo($line);

        $Info{RRD} = 0 if $Info{RRD} eq "inf";

        next if $Info{SVLEN} > 1000000;

        $Info{FILTER} = "PASS";

        ###################
        #  Filters  CNVs  #
        ###################
        # Now checking deletion filters
        if ($Info{SVTYPE} eq 'DEL') {
            
            $Info{FILTER} = "LowQual" if $Info{RRD} > 10;

            # Fore RD-only calls 
            if ($Info{BR} == 0 && $Info{PE} == 0) {
                $Info{FILTER} = "LowQual" if $Info{MQ} < 35;
                if ( $Info{CNVR} eq 'no') {
                    $Info{FILTER} = "LowQual" if $Info{NSNV} > 2 && $Info{HSNVR} < 0.75;
                }
            }
            if ($Info{BR} > 0) {
                my $ratioAsmbled = $Info{ASBR}/$Info{BR};
                $Info{FILTER} = "LowQual" if $ratioAsmbled < 0.50;
            }
            if ($Info{SVLEN} >= 200) {

                $Info{FILTER} = "LowQual" if $Info{RRD} > 0.75 && $Info{CNVR} eq 'no';
                $Info{FILTER} = "LowQual" if $Info{RRD} > 1.5 && $Info{CNVR} eq 'yes';

                #$Info{FILTER} = "LowQual" if $Info{PRD} < 10 && $Info{CNVR} eq 'no';
                if ( $Info{"5pRRD"} <  1.5 && $Info{"3pRRD"} < 1.5 && $Info{CNVR} eq 'no') {
                    $Info{FILTER} = "LowQual" if $Info{NSNV} > 2 && $Info{HSNVR} < 0.75;
                }
            }
            my $meanFlankingDepth = ( $Info{"5pRRD"}+$Info{"3pRRD"} )/2;

            if ( $Info{"5pRRD"} >  3 ||  $Info{"3pRRD"} > 3) {
                $Info{FILTER} = "LowQual";
            }            
        }
        # Now checking duplication filters
        elsif ($Info{SVTYPE} eq 'DUP' && $Info{SVLEN} >= 1000) {
            $Info{FILTER} = "LowQual" if $Info{RRD} < 1.25;
            $Info{FILTER} = "LowQual" if $Info{PRD} < 20;
        }
        if ( $Info{MQ} < 15 ) {
            $Info{FILTER} = "LowQual";
        }

        # We make sure that all descriptors have default values
        foreach my $field ( @arrFields ) {
            $Info{$field} = defined $Info{$field} ? $Info{$field} : '.';
        }

        # We transform the hash to an array
        # with this format "key=$value" to facilitate print on "join"
        my @Data = map { $_ . '=' . $Info{$_} } @arrFields;

        # Printing to VCF
        # Note that re-use $outline
        print OUT join("\t", @tmp[0..2]) . "\t$tmp[3];", join( ";", @Data ), ";\n";

    }
    close IN;
    close OUT;
}

    ##INFO=<ID=BR,Number=1,Type=Integer,Description=\"Total number of break-reads\">
    ##INFO=<ID=ASBR,Number=1,Type=Integer,Description=\"Total assembled break-reads into contig/s\">
    ##INFO=<ID=KDIV,Number=1,Type=Float,Description=\"Kmer diversity of assembled break-reads\">
    ##INFO=<ID=PE,Number=1,Type=Integer,Description=\"Total discordant paired-end reads\">
    ##INFO=<ID=PPE,Number=1,Type=Integer,Description=\"Phred-scale discordant p-value\">
    ##INFO=<ID=RRD,Number=1,Type=Float,Description=\"Read-depth ratio at flanking coordinates\">
    ##INFO=<ID=MADRD,Number=1,Type=Float,Description=\"Median Absoute Deviation of the segmented read-depth\">
    ##INFO=<ID=PRD,Number=1,Type=Integer,Description=\"Phred-scale flanking read-depth difference\">
    ##INFO=<ID=CN,Number=1,Type=Integer,Description=\"Copy number of segment containing breakend\">
    ##INFO=<ID=GC,Number=1,Type=Integer,Description=\"GC content between REF and ALT alleles\">
    ##INFO=<ID=NSNV,Number=1,Type=Float,Description=\"Total SNVs found within the SV interval\">
    ##INFO=<ID=percHomSNV,Number=1,Type=Float,Description=\"Proportion of homozygous SNVs\">
##################################
sub getInfo {
    my $line = shift;
    my @tmp = split (/\t/, $line);

    $tmp[12] = $tmp[12] > 0 ? sprintf "%.2f", $tmp[12] : 0;
    $tmp[13] = $tmp[13] > 0 ? sprintf "%.2f", $tmp[13] : 0;
    $tmp[17] = $tmp[17] > 0 ? sprintf "%.2f", $tmp[17] : 0;
    $tmp[18] = sprintf "%.2f", $tmp[18]/100;
    $tmp[19] = sprintf "%.2f", $tmp[19];
    $tmp[20] = sprintf "%.2f", $tmp[20];

    my %hash = (
        SVTYPE => $tmp[4],
        CNVR   => $tmp[5],
        SVLEN  => $tmp[6],
        MQ     => $tmp[7],
        KDIV   => $tmp[8],
        GC     => $tmp[18],
        BR     => $tmp[9],
        ASBR   => $tmp[10],
        PE     => $tmp[11],
        PPE    => $tmp[14],
        RRD    => $tmp[12],
        MADRD  => $tmp[13],
        PRD    => $tmp[15],
        NSNV   => $tmp[16],
        HSNVR => $tmp[17],
        "5pRRD"  => $tmp[19],
        "3pRRD"  => $tmp[20],
        EV => $tmp[21]
       # inSegDup  => $tmp[21],
    );
    return %hash;
}

return 1;