#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use Cwd;
our $dirname = dirname (__FILE__);

my $countFile     = $ARGV[0];
my $outDir        = $ARGV[1];
my $outName       = $ARGV[2];
my $chromo        = $ARGV[3];

if (@ARGV < 3) {
    print " Usage: $0 <cov_file> <out_dir> <out_name>\n";
    exit;  
}
my $cmd;

$chromo = "" if !$chromo;

# Compressing count info to ascii-based coverage level
$cmd = "perl $dirname/convertCount2Coverage.pl $countFile $outDir/$outName $chromo";
system $cmd;

$cmd = "perl $dirname/totalCoverageByGC.pl $outDir/$outName $dirname/../../db/RCP/GCbuckets.hg19.gz $outDir $outName $chromo";
system $cmd;

$cmd = "perl $dirname/normalizeCoverage.pl $outDir/$outName $dirname/../../db/RCP/Illumina.rcp.gz $outDir $outName $chromo";
system $cmd;

$cmd = "perl $dirname/segmentCoverage.pl $outDir/$outName";
system $cmd;