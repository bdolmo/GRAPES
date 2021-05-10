#!/usr/bin/env perl

package Breakpoint;

use strict;
use Getopt::Long;
use File::Basename;
use File::Copy;

sub AnalyzeWES {

  my $bam     = shift;
  my $genome  = shift;
  my $outdir  = shift;
  my $outname = shift;
  my $nDiscordants = shift;
  my $nSds    = shift;
  my $nBreaks = shift;

  my $rawSvCalls      = "$outdir/ON_TARGET/$outname/$outname.tmp.rawcalls.bed";
  my $ontargetSvCalls = "$outdir/ON_TARGET/$outname/$outname.breakpoints.bed";

	# Calling SV through grapes_sv binary
  my $cmd = "$::grapesBreak -b $bam --find-small on -j $::maxSizeSV -g $genome -o $outdir -n $outname -t $::threads --wes on -c $::minDiscordants -s $::minDiscordantsSD -r $::minBreakReads -e $::centromeres $::devNull";
	print"$cmd\n" if $::verbose;
  system $cmd if !-e $ontargetSvCalls;

	move "$outdir/$outname.tmp.rawcalls.bed", $rawSvCalls if -e "$outdir/$outname.tmp.rawcalls.bed";

	# Merging SV breakpoints that point to the same variant
	my $mergedSVs = mergeSV::Merge( $rawSvCalls );

	# Selecting SV calls that overlap with ROI file
	my $str =
	 `$::sort -V $mergedSVs | $::uniq | $::bedtools intersect -a stdin -b $::bed -wa -wb`;
	chomp $str;

	my @tmpBreaks = split (/\n/, $str);
	my %seen = ();

	open (OUT, ">", $ontargetSvCalls) || die " ERROR: Unable to open $ontargetSvCalls\n";
	foreach my $line (@tmpBreaks) {

		# Skipping imprecise calls
		if ($::filterDiscordantOnly) {
			next if $line =~ 'IMPRECISE';
		}

		my @tmp    = split (/\t/, $line);
		my $segment = join ('\t', @tmp[0..3]);

		$seen{$segment}++;
		next if $seen{$segment} > 1;

		my $length = $tmp[2]-$tmp[1];

		# Skipping SVs smaller than defined min and max sizes
		next if $length < $::minSizeSV;
		next if $length > $::maxSizeSV;

		# getting all ROI targets (that come from the 4th column on the BED file) inside the segment
		my @targets   = grep ($_ =~/$segment/, @tmpBreaks);

		# Total numbe rof rois affected by the SV
		my $nROIs = scalar @targets;

		# Merging affected ROIs if there is a multi-exonic CNV
		my $affectedROIs = getAffectedROIs( $segment, \@targets );

		my $precision   = $tmp[3];
		my $svtype      = $tmp[4];
		my $mapQ        = $tmp[7];
		my $kdiv        = sprintf "%.2f", $tmp[8];
		my $AF          = $tmp[9];
		my $breakReads  = $tmp[10];
		my $assembled   = $tmp[11];
		my $PE          = $tmp[12];
		my $rdRatio     = $tmp[13];
		my $rdMad       = $tmp[14];
		my $isLOh       = $tmp[15];
		my $cumulative  = $tmp[16];
		my $nins        = $tmp[17];

		my $cipos="-10,10";
		my $ciend="-10,10";

		print OUT "$tmp[0]\t$tmp[1]\t$tmp[2]\t$precision;CIPOS=$cipos;CIEND=$ciend;SVTYPE=$svtype;SVLEN=$length;EV=SR;REGIONS=$nROIs;GENE=$affectedROIs;MAPQ=$mapQ;KDIV=$kdiv;AF=$AF;BREAKREADS=$breakReads;ASSEMBLED=$assembled;PE=$PE;RRD=$rdRatio;MADRD=$rdMad;LOH=.;CSDISC=$cumulative;NINS=$nins\n";
	}
	close OUT;

	# Removing temporal files
  removeSvTmpFiles($outdir, $outname);

 }
####################################################################
 sub getAffectedROIs {

	my $segment = shift; # coordinate from the called cnv segment
	my $arrRef  = shift; # Array reference of each entry from the ROI bed

	my @arr = @$arrRef;

	# getting all ROI targets (that come from the 4th column on the BED file) inside the segment
	my @targets   = grep ($_ =~/$segment/, @arr);

	#Defining left-most and right-most exons from the segment if available
	my ($firstEx, $firstGene, $firstExon);
	my ($lastEx, $lastGene, $lastExon);

	# 0-index will be the first exon
	$firstEx     = $targets[0];
	my @tmpFirst = split (/\t/, $firstEx);

	#-1 index is the last exon
	$lastEx  = $targets[-1];
	my @tmpLast = split (/\t/, $lastEx);

	my @arrROIs;
	my $affectedROIs;

	# if ROIs contains the expected format (e.g NM_707070_3_2;DUMMYGENE)
	if ( $tmpFirst[-1] =~/.{1,}_.{1,}_.{1,}_.{1,};.{1,}$/
		&& $tmpLast[-1] =~/.{1,}_.{1,}_.{1,}_.{1,};.{1,}$/ )  {
		$firstGene = ( split /[_;]/, $tmpFirst[-1] )[4];
		$firstExon = ( split /[_;]/, $tmpFirst[-1] )[3];
		$lastGene  = ( split /[_;]/, $tmpLast[-1] )[4];
		$lastExon  = ( split /[_;]/, $tmpLast[-1] )[3];

		# For multi-exonic CNVs displaying first and last exon
		if (@targets > 3) {
			$affectedROIs = "$firstGene\_$firstExon" . "," . "$lastGene\_$lastExon";
		}
		else {
			foreach my $roi (@targets) {
				my @tmp = split (/\t/, $roi);
				push @arrROIs, $tmp[-1];
			}
			$affectedROIs = join (",", @arrROIs );
		}
	}
	# However, if ROI name contains other format we just concatenate its content
	else  {
		# For multi-exonic CNVs displaying first and last exon
		if (@targets > 3) {
			$affectedROIs = $tmpFirst[-1] . "," . $tmpLast[-1];
		}
		else {
			foreach my $roi (@targets) {
				my @tmp = split (/\t/, $roi);
				push @arrROIs, $tmp[-1];
			}
			$affectedROIs = join (",", @arrROIs );
		}
	}
	$affectedROIs =~s/;/_/g;
	return $affectedROIs;
 }

#############################
 sub removeSvTmpFiles {

   my $outdir  = shift;
   my $outname = shift;

   unlink ("$outdir/$outname.FR.bam");
   unlink ("$outdir/$outname.FR.bam.bai");

   unlink ("$outdir/$outname.RF.bam");
   unlink ("$outdir/$outname.RF.bam.bai");

   unlink ("$outdir/$outname.FF.bam");
   unlink ("$outdir/$outname.FF.bam.bai");

   unlink ("$outdir/$outname.RR.bam");
   unlink ("$outdir/$outname.RR.bam.bai");

   unlink ("$outdir/$outname.SR.bam");
   unlink ("$outdir/$outname.SR.bam.bai");

   unlink ("$outdir/$outname.discordantInfo.txt");
   unlink ("$outdir/$outname.fastq");
   unlink ("$outdir/clusters_SR.bed");
   unlink ("$outdir/$outname.tmp.rawcalls.bed");
   #unlink ("$outdir/ON_TARGET/$outname/$outname.tmp.rawcalls.bed");

   unlink( glob "$outdir/$outname.*clusters.bed");
 }

1;
