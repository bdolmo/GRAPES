#!/usr/bin/env perl

package Annotate;

use strict;
use Getopt::Long;
use File::Basename; 
use List::MoreUtils qw(uniq);

# Purpose of this module is to annotate VCF files with gnomAD database

#Hardcoded default fields 
my @popDefaults = (
    "EVIDENCE",
    "AF",
    "AN",
    "AC",
    "AFR_AC",
    "AFR_AN",
    "AFR_AF",  
    "AMR_AC",
    "AMR_AN",
    "AMR_AF",
    "EAS_AC",
    "EAS_AN",
    "EAS_AF",
    "EUR_AC",
    "EUR_AN",
    "EUR_AF",
    "OTH_AC",
    "OTH_AN",
    "OTH_AF"
);

################
sub doAnnotation {
    
    my $inputVCF = shift;

    # Get a hash with vcf fields from gnomad
    my %headerVCF = returnVcfFields();

    # Annotate with gnomAD
    annotateGnomad($inputVCF, \%headerVCF);
}

################
sub getGenomicCode {
    my $chr    = shift;
    my $start  = shift;
    my $end    = shift;
    my $svtype = shift;
}

################
sub annotateGnomad {

    my $inputVCF = shift;
    my $hashRef  = shift;

    my %headerVCF = %$hashRef;

    ( my $annotVCF = $inputVCF ) =~ s/\.vcf/\.tmp.vcf/;
    open (VCF, ">", $annotVCF) || die " ERROR: Unable to open $annotVCF\n";

    if ( isGzipped($inputVCF) ) {
        open (IN, "$::zcat $inputVCF |") 
            || die "ERROR: Unable to open $inputVCF\n";
    }
    else {
        open (IN, "<", $inputVCF) 
            || die "ERROR: Unable to open $inputVCF\n";
    }

    my %seenVar = ();

    my $flag = 0;
    while (my $line=<IN>) {
        chomp $line;
        if ($line=~/^#/) {
            if ($line=~/##ALT/ && !$flag) {
                $flag = 1;
                foreach my $field (@popDefaults) {
                    
                    my $outField = $headerVCF{$field}{DESCRIPTION};
                    $outField =~s/EVIDENCE/gnomAD_SV_method/;

                    print VCF "$outField\n";
                }
                print VCF "$line\n";
            }
            else {
                print VCF "$line\n";
            }
            next;
        }

        $seenVar{$line}++;
        next if $seenVar{$line} > 1;

        my @tmp = split (/\t/, $line);
        my @info= split (/;/, $tmp[7]);

        my $chr   = $tmp[0];
        my $start = $tmp[1];

        my ($end) = grep ($_=~/^END=/, @info);        
        $end=~s/END=//;

        my ($svtype) = grep ($_=~/^SVTYPE=/, @info);
        ($svtype) =~s/SVTYPE=//;

        my $precision = $info[0];

        my @gnomad 
            = tabixQuery($chr, $start, $end, $svtype, $precision);

        if (@gnomad) {
            foreach my $entry (@gnomad) {
                my @info = split(";", $entry);
                my @gathered = ();

                foreach my $desired (@popDefaults) {
                    my ($x) = grep ($_=~/^$desired=/, @info);
                    push @gathered, $x;
                }
                
                my $annot = join(";", @gathered);
                $annot =~s/EVIDENCE/gnomAD_SV_method/;

                my $extendedVar = join("\t", @tmp[0..7]) . ";" . $annot . join ("\t", @tmp[8..@tmp-1]);

                $seenVar{$extendedVar}++;
                next if $seenVar{$extendedVar} > 1;

                print VCF join("\t", @tmp[0..7]) . ";" . $annot . "\t". $tmp[8] . "\t" . $tmp[9] . "\n";
            }
        }
        else {
            my @tmp = split (/\t/, $line);
            my %popHash = map { $popDefaults[$_] => '.' } 0..$#popDefaults;
            my @popData = map { $_ . '=' . $popHash{$_} } keys %popHash;
            print VCF join("\t", @tmp[0..7]) . ";" . join(";", @popData) ."\t". $tmp[8] . "\t" . $tmp[9] . "\n";
        }
    }
    close IN;
    close VCF;

    # Now erase unnanotated VCF and rename tmp annotated
    unlink($inputVCF);
    rename $annotVCF, $inputVCF;

}

################
sub tabixQuery {

    my $chrQuery   = shift;
    my $startQuery = shift;
    my $endQuery   = shift;
    my $typeQuery  = shift;
    my $precision  = shift;

    # Remove chr if present
    if ($chrQuery =~/^chr/) {
        $chrQuery =~s/chr//;
    }

    my $str = `$::tabix $::gnomADvcf $chrQuery:$startQuery-$endQuery`;   
    chomp $str;
    print "$::tabix $::gnomADvcf $chrQuery:$startQuery-$endQuery\n" 
        if $::verbose;

    my @tmpStr = split(/\n/, $str);

    my @records = ();

    foreach my $entry (@tmpStr) {

        my @tmp = split (/\t/, $entry);
        my @info= split (/;/, $tmp[7]);

        my $chrDB   = $tmp[0];
        my $startDB = $tmp[1];

        my ($endDB) = grep ($_=~/^END=/, @info);        
        $endDB=~s/END=//;

        my ($typeDB) = grep ($_=~/^SVTYPE=/, @info);
        ($typeDB)=~s/SVTYPE=//;

        # Calculate overlap between query and database
        my ($overlapQuery, $overlapDB)= calculateOverlap($chrQuery, 
        $startQuery, $endQuery, $chrDB, $startDB, $endDB);

        my $joinedInfo = join(";", @info);

        # Variant class must be the same
        if ($typeDB eq $typeQuery) {

            # Do not apply reciprocal overlap for IMPRECISE calls
            if ($precision eq 'IMPRECISE') {

                if ($overlapQuery == 1) {
                    push @records, $joinedInfo;
                }
            }
            # But do it for PRECISE calls
            else {
                if ($overlapQuery >= 0.7 && $overlapDB >= 0.7 ) {
                    push @records, $joinedInfo;               
                }
            }
        }
    }
    return @records;
}

################
sub isGzipped {
    my $inputFile = shift;
    if ($inputFile =~/.gz/) {
        return 1;
    }
    else {
        return 0;
    }
}

################
sub returnVcfFields {
    
    my %headerVCF = ();

    my @infoArray = ();
    open (IN, "$::zcat $::gnomADvcf | ") || die " Unable to open $::gnomADvcf\n";
    while (my $line=<IN>) {
        chomp $line;
        if ($line =~/^##INFO/) {
            push @infoArray, $line;
        }
        last if $line !~/^##/;
    }
    close IN;

    my $i = 0;
    foreach my $field (@infoArray) {
        my @info = split(/,/, $field);
        # ##INFO=<ID=AC,Number=A,Type=Integer,Description="Number of non-reference alleles observed (for biallelic sites) or individuals at each copy state (for multiallelic sites).">
        my $tag = $info[0];
        $tag =~s/INFO=<ID=//;
        $tag =~s/#//g;
        $headerVCF{$tag}{DESCRIPTION} = $field;
        $headerVCF{$tag}{IDX} = $i;
        $i++;
    }
    return %headerVCF;
}


######################
sub calculateOverlap {

    my $chr   = shift;
    my $start = shift;
    my $end   = shift;

    my $chrom = shift;
    my $st    = shift;
    my $ed    = shift;

    my $length_A = $end-$start;
    my $length_B = $ed-$st;
    my $olap_A;
    my $olap_B;

    # Minimum reciprocal overlap of 0.9 
    # Case1: coord fits within intron
    #A          ////////////   cnv (start, end)
    #B        ^^^^^^^^^^^^^^^^ intron (st, ed)
    if ($start >= $st && $end <= $ed) {
        $olap_A = 1.00;
        $olap_B = 1.00 - ( ( ($start-$st) + ($ed-$end) ) / $length_B);
    } 
    # Case2: coord left-overlap
    #A      ////////////////   cnv (start, end)
    #B        ^^^^^^^^^^^^^^^^ intron (st, ed)
    if ($start < $st && $end < $ed) {
        $olap_A = 1.00 - ( ( $st-$start ) / $length_A);
        $olap_B = 1.00 - ( ( ($ed-$end) ) / $length_B);
    }
    # Case3: coord right-overlap
    #A           /////////////// cnv (start, end)
    #B        ^^^^^^^^^^^^^^^^   intron (st, ed)
    if ($start > $st && $end > $ed) {
        $olap_A = 1.00 - ( ( $end-$ed ) / $length_A);
        $olap_B = 1.00 - ( ( ($start-$st) ) / $length_B);
    }
    # Case4: coord fits entirely
    #A      ///////////////////// cnv (start, end)
    #B        ^^^^^^^^^^^^^^^^    intron (st, ed)
    if ($start <= $st && $end >= $ed) {
        $olap_A = 1.00 - ( ( ($st-$start) + ($end-$ed) ) / $length_A);
        $olap_B = 1.00;
    }
    return ($olap_A, $olap_B);    
}

return 1;