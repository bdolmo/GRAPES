#!/usr/bin/env perl
package VAF;

use strict;
use Getopt::Long;
use File::Basename; 
use List::Util qw(min max);
use Sort::Key::Natural qw(natsort);

# This module is intended for calculating the B-allele frequencies from each sample (if VCF file is available):
# "The B-Allele Frequency is a normalized measure of the allelic intensity ratio of two alleles (A and B), 
# such that a BAF of 1 or 0 indicates the complete absence of one of the two alleles (e.g. AA or BB), 
# and a BAF of 0.5 indicates the equal presence of both alleles (e.g. AB)."
# Detection of allelic imbalances such as those caused by duplications (e.g. AAB/BBA) or mosaic deletions in the sample. 
# Such imbalances can be identified on a BAF plot by the presence of SNPs at frequencies between 0.5 and 0 or 1.
# For example, the theoretical BAF values of triploid regions (AAA, AAB, ABB or BBB) are 0, 0.33, 0.66 and 1 respectively."

my $minFreebQual= 200;
my $minVarDepth = 30; 
my $varType     = "snp";

sub selectVariants {

    foreach my $sample ( natsort keys %::sampleHash ) {

        if (-e "$::sampleHash{$sample}{SNPVCF}" ) {

            my $filteredVCF = "$::outDir/BAF_DATA/" . basename($::sampleHash{$sample}{SNPVCF});
            $filteredVCF =~s/.vcf/.filtered.vcf/;

            open (VCF, "<", "$::sampleHash{$sample}{SNPVCF}") || die " ERROR: Unable to open $::sampleHash{$sample}{SNPVCF}\n";

            open (FVCF, ">", $filteredVCF) || die " ERROR: Unable to open $filteredVCF\n";

            my $BAF = "$::outDir/BAF_DATA/" . $sample . ".BAF.bed";
            open (BAF, ">", $BAF) || die " ERROR: Unable to open $BAF\n";

            while (my $line=<VCF>) {

                chomp $line;

                # Print header
                if ($line =~/^#/) {
                    print FVCF "$line\n";
                    next;
                }

                my @tmp = split (/\t/, $line);
                my $Qual= $tmp[5];
                ###INFO=<ID=DP4,Number=4,Type=Integer,Description="Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases">
                #chrM	410	.	A	T	10.1993	.	DP=1;SGB=-0.379885;MQ0F=0;ICB=1;HOB=0.5;AC=1;AN=2;DP4=0,0,1,0;MQ=60	GT:PL	0/1:38,3,0
                #chrM	150	.	TCT	CCC	395.542	.	AB=0;ABP=0;AC=2;AF=1;AN=2;AO=13;CIGAR=1X1M1 GT:DP:AD:RO:QR:AO:QA:GL	1/1:12:0,12:0:0:12:444:-28.245,-3.61236,0
                my @info = split (/;/, $tmp[7]);
                my ($DP) = grep ($_=~/^DP=/, @info);
                $DP =~s/DP=//;
                my $Var = '.';
                if (length($tmp[3]) > 1 || length ($tmp[4]) > 1) {
                    $Var = "indel";
                    }
                else {
                    $Var = "snp";
                }

                my $Freq;
                if ($line=~/DP4=/) {
                    my ($DP4) = grep ($_=~/DP4=/, @info);
                    $DP4=~s/DP4=//;
                    #DP4=0,0,1,0
                    #ref-fwd, ref-rev, alt-fwd, alt-rev
                    my @tmpDP = split (/,/, $DP4);
                    my $ref = $tmpDP[0]+$tmpDP[1];
                    my $alt = $tmpDP[2]+$tmpDP[3];
                    $Freq = sprintf "%.3f",$alt/$DP; 
                }
                else {
                    my @details = split (/:/, $tmp[9]);
                    my $AD = $details[2];
                    my @tmpAD = split (/,/, $AD);
                    $Freq = sprintf "%.3f",$tmpAD[1]/$DP;
                }
                # Only selecting variants that pass the filters specified               
                if ($Qual < $minFreebQual) {
                    next;
                }
                if ($DP < $minVarDepth) {
                    next;
                }
                if ($Var ne $varType) {
                    next;
                }
                print FVCF "$line\n";
                my $End = $tmp[1]+1;

                print BAF "$tmp[0]\t$tmp[1]\t$End\t$Freq\n";
            }
            close VCF;

            my $segOnOffData= -e "$::outDir/SEGMENT_DATA/toplot.segmented.$sample.bed" 
            ? "$::outDir/SEGMENT_DATA/toplot.segmented.$sample.bed" : undef;
            my $segOnData   = -e "$::outDir/ON_TARGET/SEGMENT_DATA/toplot.segmented.$sample.bed" 
            ? "$::outDir/ON_TARGET/SEGMENT_DATA/toplot.segmented.$sample.bed" : undef;
            my $segOffData  = -e "$::outDir/OFF_TARGET/SEGMENT_DATA/toplot.segmented.$sample.bed" 
            ? "$::outDir/OFF_TARGET/SEGMENT_DATA/toplot.segmented.$sample.bed" : undef;

            if (defined $segOnOffData) {
                annotateBAF($segOnOffData, $BAF, $sample, "MAIN");
            }
            if (defined $segOnData) {
                annotateBAF($segOnData, $BAF, $sample, "ON_TARGET");
            }
            if (defined $segOffData) {
                annotateBAF($segOffData, $BAF, $sample, "OFF_TARGET");
            }
        }
    }
}

#####################
sub varCallFreebayes {
    my $bamFile = shift;
    my $sample  = shift;

    # Output is uncompressed VCF
    print " INFO: Calling SNVs on sample $sample (freebayes)\n";
    my $cmd = "$::freebayes --min-alternate-count 10 -f $::genome $bamFile > $::input/$sample.vcf";
	print "$cmd\n" if $::verbose;
	system $cmd;
}

#####################
sub varCallSamtools {
    my $bamFile = shift;
    my $sample  = shift;

    # Output is uncompressed VCF
    print " INFO: Calling SNVs on sample $sample (samtools)\n";
    my $cmd = "$::samtools mpileup -uf $::genome $bamFile $::devNullStderr | $::bcftools call -mv -Ov  > $::input/$sample.vcf ";
	print "$cmd\n" if $::verbose;
    system $cmd;
}

#####################
sub annotateBAF {
    
    my $segFile   = shift; # segment file
    my $snpFile   = shift;
    my $sample    = shift;
    my $dirType   = shift;
    my $bafSegmentFile;

    if ($dirType eq 'MAIN') {
        $bafSegmentFile = "$::outDir/BAF_DATA/$sample.baf_segment.bed";
        $::sampleHash{$sample}{BAF} = $bafSegmentFile;
    }
    if ($dirType eq 'ON_TARGET') {
        $bafSegmentFile = "$::outDir/BAF_DATA/$sample.baf_segment.bed";
        $::sampleHash{$sample}{BAF_ONTARGET} = $bafSegmentFile;
    }
    if ($dirType eq 'OFF_TARGET') {
        $bafSegmentFile = "$::outDir/BAF_DATA/$sample.baf_segment.bed";
        $::sampleHash{$sample}{BAF_OFFTARGET} = $bafSegmentFile;
    }

    if (-e $segFile && -e $snpFile) {
        my $cmd = "$::bedtools intersect -a $snpFile -b $segFile -wa -wb > $bafSegmentFile";
        system $cmd;
    }
}



return 1;