#!/usr/bin/env perl

package CNV;

use strict;
use warnings;
use File::Basename;
use Cwd;
use List::Util qw(sum);
use Sort::Key::Natural qw(natsort);

sub call {

    my $countFile  = shift;
    my $outDir     = shift;
    my $outName    = shift;
    my $chromo     = shift;

    $chromo = "" if !$chromo;
    my $cmd = "perl $::executeRCP $countFile $outDir $outName $chromo";
    #system $cmd;

    open (IN, "$::zcat $outDir/$outName.seg.gz |") || die " ERROR: Unable to open $outDir/$outName.seg.gz\n";
    open (OUT, ">", "$outDir/CNV.tmp.txt") || die " ERROR: Unable to open $outDir/CNV.tmp.txt\n";
    while (my $line=<IN>) {
        chomp $line;
        next if $line =~/^#/;
        my ($chr, $start, $end, $copyNumber, $kbSize) = split (/\t/, $line);
        my $svtype;
        if ($copyNumber < 2) {
            $svtype = "DEL";
        }
        if ($copyNumber > 2) {
            $svtype = "DUP";
        }
        my $svlen = $end-$start;
        if ($svtype) {
            print OUT "$chr\t$start\t$end\tIMPRECISE;SVTYPE=$svtype;CN=$copyNumber\n";
        }
    }
    close IN;
    close OUT;

    # Extract mappability
    $cmd = "";
    my $outWithMapp = "$outDir/$outName.MAP.CNV.bed";
    $cmd .= "$::bedtools intersect -a $outDir/$outName.counts_window.bed.gz -b $outDir/CNV.tmp.txt -wo";
    $cmd .= "| $::awk '{ size=\$9-\$8; marginal=100*(\$NF*\$6)/size; print \$7\"\t\"\$8\"\t\"\$9\"\t\"\$10\"\t\"\$6\"\t\"\$11\"\t\"marginal}'";
    #$cmd .= "| $::awk '{ size=\$9-\$8; marginal=100*(\$NF*\$6)/size; print \$7\"\t\"\$8\"\t\"\$9\"\t\"\$10\"\t\"\$6\"\t\"\$11\"\t\"marginal}' > $outWithMapp";

    $cmd .= "| perl $::extractMapp | $::sort -V  > $outWithMapp";

    system $cmd if !-e $outWithMapp;

    my $outFiltered = filterCNV($outWithMapp);

    calculateMAD($outFiltered, "$outDir/$outName.norm.gz");

    #unlink("$outDir/CNV.tmp.txt", $outWithMapp, $outFiltered);

}
############################
sub filterCNV {
    my $infile  = shift;
    my $outfile = $infile;
    $outfile =~s/.MAP/.TMP/;

    my $cmd = "$::bedtools intersect -a $infile -b $::centromeres -v > $outfile";
    system $cmd;
    return $outfile;
}

############################
sub calculateMAD {
    my $cnvFile  = shift;
    my $normFile = shift;

    my $outfile = $cnvFile;
    $outfile =~s/.TMP//;

     # Hash of Cnv
    my %HoCNV = ();
    open (IN, "$::bedtools intersect -a $normFile -b $cnvFile -wao | $::grep -vP \'\t-1\t\' |");
    while (my $line=<IN>) {
        chomp $line;
        my @tmp = split (/\t/, $line);
        my $cnv = join ("\t", @tmp[4..7]);
        my $ratio;
        if ($tmp[3] > 0) {
            $ratio = $tmp[3]/100;
        }
        else {
            $ratio = 0;
        }
        push @{$HoCNV{$cnv}},$ratio;
    }
    close IN;

    open (OUT, ">", $outfile) || die " ERROR: Unable to open $outfile\n";
    foreach my $cnv (natsort keys %HoCNV) {
        my $MAD = MAD( @{$HoCNV{$cnv}} );
        my $RDratio = mean( @{$HoCNV{$cnv}} );
        print OUT "$cnv;RDmad=$MAD;RDratio=$RDratio\n";
    }
    close OUT;
}

############################
sub mean {
    return 0 if @_ == 0;

    return sprintf "%.2f", sum(@_)/@_;
}
############################
sub MAD {
	# Calculate Median absolute Deviation

	my @array = @_;

	my $median = Median(@array);
	my @substractMedian = ();

	foreach my $value (@array) {
		my $x = abs($value-$median);
		push @substractMedian, $x;
	}

	my $MAD = abs( Median (@substractMedian) );
	return sprintf "%.2f",$MAD;
}
############################
sub Median
{
    my @vals = sort {$a <=> $b} @_;
    my $len = @vals;
    if($len%2) #odd?
    {
        return $vals[int($len/2)];
    }
    else #even
    {
        return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
    }
}


return 1;