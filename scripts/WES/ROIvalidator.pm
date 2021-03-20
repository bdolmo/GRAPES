#!/usr/bin/env perl

package ROIvalidator;

use strict;
use warnings;
use File::Basename;

# Here we will validate ROI coordinates as follows:
# ROI file is sorted alphanumerically, if not, dump a warning and sort it automatically
# ROI file must contain 4 fields: chr[tab]start[tab]end[tab]roi
# ROI file and BAM files must match in chromosome convention usage (e.g chr1 or 1)
# Assembled chromosomes, alt contigs, etc must match between ROI and BAM files

sub validate {

    if (!$::bed) {
        print " ERROR: missing BED regions file\n";
        exit; 
    } 
    if (!-e $::bed) {
        print " ERROR: $::bed BED regions file does not exist!\n";
        exit; 
    }

    # Check if ROI bed is sorted
    my $checkSorted = `$::sort -c -V $::bed`;
    chomp $checkSorted;
    if ($checkSorted) {

        print " WARNING: ROI file is not sorted alpha-numerically\n";
 
        #sort it automatically
        print " INFO: Sorting $::bed. Sorted bed will be located at $::outDir\n";
        my $bedName = basename ($::bed);
        my $cmd = "$::sort -V $::bed > $::outDir/$bedName";
        system $cmd;
        $::bed = "$::outDir/$bedName";
    }
    my %ChromLength = ();
    foreach my $bam ( @::bams) {
        # Checking that bam file contains the same chromosome nomenclature that the input ROI file
        my $header = `$::samtools view -H $bam | $::grep '^\@SQ'`;
        chomp $header;
        my @tmpHeader = split(/\n/, $header);
        foreach my $line (@tmpHeader) {
            chomp $line;
            my @tmp = split (/\t/, $line);
            my $chr = $tmp[1];
            my $length = $tmp[2];
            $chr =~s/SN://;
            $length =~s/LN://;
            $ChromLength{$chr} = $length;
        }
    }
    if ( checkROIconsistency(\%ChromLength) ) {
        print " INFO: $::bed ROI file is well formatted!\n";
    }
}

###########################
sub checkROIconsistency {

    my $hash = shift;
    my %ChromLength = %$hash;
    my %ROIchr = ();

    open (BED, "<", $::bed) || die " ERROR: Unable to open $::bed\n";
    my $entry = 1;
    my $isOk  = 1;
    while (my $line =<BED>) {
        chomp $line;
        my @tmp = split (/\t/, $line);
        if (@tmp < 4) {
            print " ERROR: ROI file (entry $entry). $::bed requires 4 fields per record (" . scalar(@tmp) . ") detected\n"; 
            exit;
        }
        if ($tmp[2] < $tmp[1]) {
            print " ERROR: ROI file (entry $entry). End coordinate $tmp[2] must be greater than Start coordinate $tmp[1]\n";
            exit;
        }
        if (!exists $ChromLength{$tmp[0]}) {
            print " ERROR: ROI file (entry $entry). Chromosome $tmp[0] is not found on BAM header. ROI bed and BAM require the same genome assembly\n";
            exit;
        }
        if ($tmp[2] > $ChromLength{$tmp[0]} ) {
            print " ERROR: ROI file (entry $entry). End coordinate $tmp[2] is greater than expected $tmp[0] length (%$ChromLength{$tmp[2]})\n";
            exit;
        }

        $ROIchr{$tmp[0]}++;
        $entry++;
    }
    return $isOk;
}

return 1;