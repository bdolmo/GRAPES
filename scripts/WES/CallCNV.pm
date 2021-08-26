#!/usr/bin/env perl

package CallCNV;

use strict;
use Getopt::Long;
use File::Basename;
use List::MoreUtils;
use Statistics::Descriptive;
use Parallel::ForkManager;
use List::MoreUtils qw(uniq);
use List::MoreUtils qw(first_index);
use Sort::Key::Natural qw(natsort);
use File::Copy;

 ####################################
 #	        Calling CNVs            #
 ####################################

 our %sampleIndex = ();

 sub call {

  my $inputDir = shift; # on-target, off-target or main outDir
  my $analysis = shift; # Whole-exome or gene-panel
  my $type     = shift; # on-target or off-target or mixed
  my $bedfile  = shift; # BED ROI file for on-target or off-target or mixed
  my $ratioFile= shift; # ratio file for on-target, off-target or mixed

  #Ensuring we have GC and mappability values for every target region being analyzed
  if ( $type eq 'off-target') {
  	Utils::populateMapGcHashFromNormCountsOfftarget();
  }
  elsif ($type eq 'on-target'){
  	Utils::populateMapGcHashFromNormCountsOntarget( $::HoF{MERGED_NORM_COV} );
  }

  # Analyzing gene panel
  if ($analysis eq 'gene-panel') {

    callMultipleExonCNV($inputDir, $analysis, $type, $bedfile, $ratioFile);

	   if ($type eq 'on-target') {
		     callSingleExonCNV($inputDir, $analysis);
  	 }
  }
  # Analyzing exome
  elsif ($analysis eq 'exome') {
  	callMultipleExonCNV($inputDir, $analysis, $type, $bedfile, $ratioFile);
  }

  # Annotate BAF values
#   if ($::doVAF) {
#   	appendVAF($inputDir);
#   }

  dumpFinalCalls($inputDir);

  print " INFO: Adding confidence scores\n";
  scoreCNV($inputDir);

  # Deleting temporary files
  unlink (glob ("$inputDir/*ratios.bed"));
  unlink (glob ("$inputDir/*single.exon.cnv.bed"));
  unlink (glob ("$inputDir/Breaks_and_CNV*"));
  unlink (glob ("$inputDir/cnv_calls.*"));
  unlink (glob ("$inputDir/intersect.*"));
  unlink (glob ("$inputDir/*.exon.cnv.bed"));

 }
####################################################################
sub scoreCNV {

	my $outputDir = shift;

	foreach my $sample (sort keys %::sampleHash) {

		my $finalCalls  = "$outputDir/$sample.CNV.bed";
		my $scoredCalls = "$outputDir/$sample.CNV.score.bed";

		if (-e $finalCalls) {

			my $corr     = $::sampleHash{$sample}{CORRELATION};
			my $nCalls   = `$::cat $finalCalls | $::wc -l`;
			chomp $nCalls;

			my $pctCalls =  ($nCalls/scalar(@::ROIarray))*100;
			my $ratioStd = $::sampleHash{$sample}{ONTARGET_SD_RATIO};

			open (IN, "<", $finalCalls) || die " ERROR: Unable to open $finalCalls\n";
			open (OUT, ">", $scoredCalls) || die " ERROR: Unable to open $scoredCalls\n";
			while (my $line=<IN>){
				chomp $line;
				my @tmp  = split("\t", $line);
				my @info = split(";", $tmp[3]);

				my $precision = $info[0];

				my ($svtype) = grep ($_=~/SVTYPE=/,@info );
				$svtype =~s/SVTYPE=// if $svtype;

				my ($regions) = grep ($_=~/REGIONS=/,@info );
				$regions =~s/REGIONS=// if $regions;

				my ($ratio) = grep ($_=~/RRD=/,@info );
				$ratio =~s/RRD=// if $ratio;

				my ($s2n) = grep ($_=~/SNR=/,@info );
				$s2n =~s/SNR=// if $s2n;

				my ($s2nc) = grep ($_=~/SNRC=/,@info );
				$s2nc =~s/SNRC=// if $s2n;

				my ($cn) = grep ($_=~/CN=/,@info );
				$cn =~s/CN=// if $cn;

				my ($score) = "";
				if ($precision eq 'PRECISE'){
					$score = '.';
				}
				else {
					# Single exon CNV
					if ($regions == 1){
						$score = scoreSingleExon($svtype, $nCalls, $ratioStd, $pctCalls, $corr, $ratio, $s2n, $s2nc, $cn);
					}
					# Multiple exon CNV
					else {
						$score = scoreMultipleExon($svtype, $nCalls, $ratioStd, $pctCalls, $corr, $ratio, $s2n, $s2nc, $cn);
					}
				}
				print OUT "$line;CONF_SCORE=$score\n";
			}
			close IN;
			close OUT;
			unlink $finalCalls;
			rename $scoredCalls, $finalCalls;
		}
	}
}
####################################################################
sub scoreSingleExon {

	my $svtype   = shift;
	my $nRois    = shift;
	my $sampleStd= shift;
	my $pctCalls = shift;
	my $corr     = shift;
	my $ratio    = shift;
	my $s2n      = shift;
	my $s2nc     = shift;
	my $cn       = shift;

	# Sample Level Metrics
	my $stdMetric      = $sampleStd < 0.2 ? 1 : 0;
	my $pctCallsMetric = $pctCalls < 1.5 ? 1 : 0;
	my $corrMetric     = $corr > 0.97 ? 1 : 0;
	my $SLM = 0.45*$stdMetric + 0.45*$pctCallsMetric + 0.1*$corrMetric;

	# Roi LeveL Metrics
	my $expectedRatio = 0.5*$cn;
	my $absDiffRatio  = 0;
	if ($ratio > $expectedRatio) {
		$absDiffRatio = 2.5*(abs($ratio-$expectedRatio));
	}
	else {
		$absDiffRatio =2.5*(abs($expectedRatio-$ratio));
	}
	my $absDiffMetric = 1-$absDiffRatio > 0 ? 1-$absDiffRatio : 0;
	my $signal2noiseMetric = $s2n > 10 ? 1 : 0;
	my $signal2noiseControlsMetric = $s2nc > 10 ? 1 : 0;
	#print"$stdMetric\t$pctCallsMetric\t$corrMetric\t$absDiffMetric\t$signal2noiseMetric\t$signal2noiseControlsMetric\n";

	my $RLM = 0.50*$absDiffMetric + 0.25*$signal2noiseMetric+0.25*$signal2noiseControlsMetric;

	# Score = SLM + RLM
	my $score = 0.4*$SLM + 0.6*$RLM;

}
####################################################################
sub scoreMultipleExon {

	my $svtype   = shift;
	my $nRois    = shift;
	my $sampleStd= shift;
	my $pctCalls = shift;
	my $corr     = shift;
	my $ratio    = shift;
	my $s2n      = shift;
	my $s2nc     = shift;
	my $cn       = shift;

	# Sample Level Metrics
	my $stdMetric      = $sampleStd < 0.2 ? 1 : 0;
	my $pctCallsMetric = $pctCalls < 1.5 ? 1 : 0;
	my $corrMetric     = $corr > 0.97 ? 1 : 0;
	my $SLM = 0.45*$stdMetric + 0.45*$pctCallsMetric + 0.1*$corrMetric;

	# Roi LeveL Metrics
	my $expectedRatio = 0.5*$cn;
	my $absDiffRatio  = 0;
	if ($ratio > $expectedRatio) {
		$absDiffRatio = 2.5*(abs($ratio-$expectedRatio));
	}
	else {
		$absDiffRatio =2.5*(abs($expectedRatio-$ratio));
	}
	my $nRoisMetric   = $nRois > 1 ? 1 : 0;
	my $absDiffMetric = 1-$absDiffRatio > 0 ? 1-$absDiffRatio : 0;
	my $signal2noiseMetric = $s2n > 10 ? 1 : 0;
	my $signal2noiseControlsMetric;

	if ($s2nc eq '.'){
		$signal2noiseControlsMetric = 0;
	}
	else{
		$signal2noiseControlsMetric = $s2nc > 10 ? 1 : 0;
	}
	my $nRoiWeight   = 0.1+($nRois*(0.05));
	my $signalWeight = 1-$nRoiWeight;

	my $RLM = ($nRoiWeight*$nRoisMetric) + ($signalWeight*0.70*$absDiffMetric) + ($signalWeight*0.15*$signal2noiseMetric) + ($signalWeight*0.15*$signal2noiseControlsMetric);
	$RLM = 1 if $RLM > 1;

	# Score = SLM + RLM
	my $score = 0.4*$SLM + 0.6*$RLM;
}
####################################################################
sub dumpFinalCalls {

	my $outputDir = shift;

	foreach my $sample (sort keys %::sampleHash) {

		my $mergedMultipleCalls = "$outputDir/$sample.merged.breaks.multiple.exon.cnv.bed";
		my $mergedAllCalls = "$outputDir/$sample.merged.breaks.single.multiple.exon.cnv.bed";
		my $finalCalls = "$outputDir/$sample.CNV.bed";

		my @tmpFiles;
		if (-e $mergedMultipleCalls) {
			push @tmpFiles, $mergedMultipleCalls;
		}
		if (-e $mergedAllCalls) {
			push @tmpFiles, $mergedAllCalls;
		}
		my $cmd = "$::cat @tmpFiles | $::sort -V | $::uniq > $finalCalls";
		system $cmd if @tmpFiles;
	}
}
####################################################################
sub getCiposCiend {

	my $chr = shift;
	my $pos = shift;
	my $end = shift;
	my $refBed = shift;

	my @ROI = @$refBed;

	# CIPOS tag will be the length difference with the previous ROI
	my ($startIdx) = grep { $ROI[$_] =~ /$pos/ } (0 ..  @ROI-1);
	$startIdx = defined $startIdx ? $startIdx : -1;
	my $previousROI = $startIdx > 0 ? $ROI[$startIdx-1] : $ROI[$startIdx];
	my @tmpPrevious = split(/\t/, $previousROI);
	my $cipos = '.';
	if ($chr eq $tmpPrevious[0]) {
		$cipos = $pos-$tmpPrevious[2];
	}

	# CIEND tag will be the length difference with the next ROI
	my ($endIdx) = grep { $ROI[$_] =~ /$end/ } (0 ..  @ROI-1);
	$endIdx = defined $endIdx ? $endIdx : -1;
	my $nextROI = $endIdx < @ROI-1 ? $ROI[$endIdx+1] : $ROI[$endIdx];
	my @tmpNext = split(/\t/, $nextROI);
	my $ciend = ".";
	if ($chr eq $tmpNext[0]) {
		$ciend = $tmpNext[1]-$end;
	}
	return $cipos, $ciend;
}

####################################################################
sub getAffectedROIs {

	# As input: a CNV segment coordinate and an array reference with all overlapping ROI names
	# Ouput: collapsed ROI info to fill the GENE VCF field

	my $segment = shift; # coordinate from the called cnv segment
	my $arrRef  = shift; # Array reference of each entry from the ROI bed

	my @arr = @$arrRef;

	# getting all ROI targets (that come from the 4th column on the BED file) inside the segment
	my @targets   = grep ($_ =~/$segment/, @arr);

	# Defining left-most and right-most exons from the segment if available
	my ($firstEx, $firstGene, $firstExon);
	my ($lastEx, $lastGene, $lastExon);

	# 0-index will be the first exon
	my $firstEx   = $targets[0];
	my @tmpFirst  = split (/\t/, $firstEx);

	#-1 index is the last exon
	my $lastEx  = $targets[-1];
	my @tmpLast = split (/\t/, $lastEx);

	my @arrROIs;
	my $affectedROIs;

	# if ROIs contains this format (NM_707070_3_2;PKP2)
	if ( $tmpFirst[3] =~/.{1,}_.{1,}_.{1,}_.{1,};.{1,}$/
		&& $tmpLast[3] =~/.{1,}_.{1,}_.{1,}_.{1,};.{1,}$/ )  {
		$firstGene = ( split /[_;]/, $tmpFirst[3] )[4];
		$firstExon = ( split /[_;]/, $tmpFirst[3] )[3];
		$lastGene  = ( split /[_;]/, $tmpLast[3] )[4];
		$lastExon  = ( split /[_;]/, $tmpLast[3] )[3];

		# For multi-exonic CNVs displaying first and last exon
		if (@targets > 3) {
			$affectedROIs = "$firstGene\_$firstExon" . "," . "$lastGene\_$lastExon";
		}
		else {
			foreach my $roi (@targets) {
				my @tmp = split (/\t/, $roi);
				push @arrROIs, $tmp[3];
			}
			$affectedROIs = join (",", @arrROIs );
		}
	}
	# However, if ROI name contains other format we just concatenate its content
	else  {
		# For multi-exonic CNVs displaying first and last exon
		if (@targets > 3) {
			$affectedROIs = $tmpFirst[3] . "," . $tmpLast[3];
		}
		else {
			foreach my $roi (@targets) {
				my @tmp = split (/\t/, $roi);
				push @arrROIs, $tmp[3];
			}
			$affectedROIs = join (",", @arrROIs );
		}
	}
	$affectedROIs =~s/;/_/g;
	return $affectedROIs;
 }

 ###################################
 #	Map sample to an index column  #
 ###################################

 sub getSampleIdxFromRatioHeader {

	my $ratioFile = shift;
	my $strHeaderRatios = `$::zcat $ratioFile | $::head -1`;
	chomp $strHeaderRatios;
	my @tmpHeaderRatios = split (/\t/, $strHeaderRatios );

	my %HashPosAll = ();
	for (my $i = 6; $i < @tmpHeaderRatios-1; $i+=2) {

		$HashPosAll{$i} = $tmpHeaderRatios[$i];
	}
	return %HashPosAll;
 }

 ####################################
 #	Map sample to an index column  #
 ####################################

 sub getSampleIdxFromNormCoverage {

	my $normCoverage = shift;
	my $strHeaderCov = `$::zcat $normCoverage | $::head -1`;
	chomp $strHeaderCov;
	my @tmpHeaderCov = split (/\t/, $strHeaderCov );

	my %sampleIndex = ();
	for (my $i = 5; $i < @tmpHeaderCov; $i++) {
		$sampleIndex{$tmpHeaderCov[$i]}=$i;
		$sampleIndex{$i}=$tmpHeaderCov[$i];
 	}
	return %sampleIndex;
 }

 ###############################################
 #	         Analyze segmented CNVs            #
 ###############################################
 sub callMultipleExonCNV {

  my $inputDir = shift;
  my $analysis = shift; # gene-panel or exome
  my $type     = shift; # on-target or off-target or mixed
  my $bedfile  = shift;
	my $ratioFile= shift;

  # Here we map every sample to its corresponding array index from ratio file
	my %HashPosAll = ();
	if ($type eq 'on-target') {
		%HashPosAll  = getSampleIdxFromRatioHeader($ratioFile);
	}

	if ($analysis eq 'gene-panel' && $type eq 'on-target') {

		if (-e "$inputDir/$::outName.NormalizedCoverage.bed") {
			Utils::compressFile("$inputDir/$::outName.NormalizedCoverage.bed");
		}
		%sampleIndex = getSampleIdxFromNormCoverage("$inputDir/$::outName.NormalizedCoverage.bed.gz");
	}

	# Here we will work always with uncompressed ratios
	if ($ratioFile =~/.gz/) {
		Utils::decompressFile($ratioFile) if -s $ratioFile;
		$ratioFile =~s/.gz//;
	}

  # Obtaining segmented files
  my @segmentedFiles = glob ("$inputDir/SEGMENT_DATA/segmented*");
  @segmentedFiles    = grep (!-z $_, @segmentedFiles);

	if (!@segmentedFiles) {
		print " WARNING: Skipping $type CNV calling. No segmented files were detected\n";
		return 1;
	}

    print " INFO: Calling multi-exon CNVs\n";

	# Get segmented CNVs (multi-exon)
  foreach my $file (@segmentedFiles) {

		my $sampName = basename($file);
		$sampName =~s/segmented.//;
		$sampName =~s/.bed//;
		$sampName=~s/.on_off//;

		my $outBed    = "$inputDir/$sampName.CNV.bed";
		my $breakFile = "$inputDir/$sampName/$sampName.breakpoints.bed";
		my $cmd;

		if (!$::sampleHash{$sampName}{ONTARGET_SD_RATIO} ||
			 $::sampleHash{$sampName}{ONTARGET_SD_RATIO} > 0.2 ) {

			$::sampleHash{$sampName}{ONTARGET_QC_PASS} = "NO";
      print " INFO: Skpping CNV calling from on-target data on sample $sampName\n";

			$cmd = "$::mergeCnvBreaks $breakFile $outBed | $::sort -V | $::uniq > $outBed";
			print "$cmd\n" if $::verbose;
			system $cmd;
			next;
		}
		else {
			$::sampleHash{$sampName}{ONTARGET_QC_PASS} = "YES";
		}
		# Avoid calling CNVs when off-target data is too noisy
		if ($type eq 'off-target') {

			if ($::sampleHash{$sampName}{PERFORM_OFFTARGET} eq 'NO') {

				$::sampleHash{$sampName}{OFFTARGET_QC_PASS} = "NO";
        print " INFO: Skpping CNV calling from off-target data on sample $sampName\n";

				$cmd = "$::mergeCnvBreaks $breakFile $outBed | $::sort -V | $::uniq > $outBed";
				print "$cmd\n" if $::verbose;
				system $cmd;
				next;
			}
			else {
				$::sampleHash{$sampName}{OFFTARGET_QC_PASS} = "YES";
			}
		}
		if ($type eq 'mixed'){
			if ( $::sampleHash{$sampName}{ONOFF_SD_RATIO}> 0.2) {
				$cmd = "$::mergeCnvBreaks $breakFile $outBed | $::sort -V | $::uniq > $outBed";
				print "$cmd\n" if $::verbose;
				system $cmd;
				next;
				$::sampleHash{$sampName}{OFFTARGET_QC_PASS} = "NO";
			}
			else {
				$::sampleHash{$sampName}{OFFTARGET_QC_PASS} = "YES";
			}
		}

		if ($type eq 'off-target'|| $type eq 'mixed') {
			$ratioFile = "$inputDir/$sampName.ratios.txt.gz";
		}

		# Temporal file for intersected segments within the targeted regions
		my $intersect = basename($file);
		$intersect =~s/segmented./intersect.segmentation./;
		$intersect = "$inputDir/$intersect";

		if ($type eq 'on-target' && $analysis ne 'exome') {
			# Removing header for intersecting purposes
			$cmd = "$::zcat $inputDir/$::outName.norm_bylib_coverage.bed.gz";
			$cmd.= "| $::tail -n +2 > $inputDir/$::outName.norm_bylib_coverage_noheader.bed";
			system $cmd;
		}

		# Selecting candidate segments that do not overlap centromeres and that pass the minimunm calling thresholds
		$cmd = "$::bedtools intersect -a $file -b $::centromeres -v | ";
		$cmd.= "$::awk '{ if ((\$5 < $::upperDelCutoff && \$5 > $::lowerDelCutoff) || (\$5 > $::lowerDupCutoff) ";
		$cmd.= "|| (\$5 > 0 && \$5 < 0.12)) { print \$0 }}' | $::sort -V > $intersect";
		system $cmd if !-e $intersect;

		# File for raw multiple-exon cnvs
		my $rawMultipleExonCnv = "$inputDir/$sampName.raw.multiple.exon.cnv.bed";

		# File for merged multiple-exon cnvs. Ad-hoc solution because sometimes DNAcopy does not merge adjacent ROIs
		my $mergedMultipleExonCnv = "$inputDir/$sampName.merged.multiple.exon.cnv.bed";

		# File for merged multiple-exon cnvs + breakpoint calls.
		my $multipleCnvAndBreaks = "$inputDir/$sampName.merged.breaks.multiple.exon.cnv.bed";

		# File for single-exon cnvs
		my $singleExonCnv = "$inputDir/$sampName.raw.single.exon.cnv.bed";

		my $ratios        = $inputDir . "/". $sampName . ".ratios.txt.gz";
		my $singleTmp     = $inputDir . "/". $sampName . ".tmp.single.exon.cnv.bed";

		open (OUTFILE, ">", $rawMultipleExonCnv) || die " ERROR: Unable to open $rawMultipleExonCnv\n";

		# if no segmented CNVs found
		if (-z $intersect && -z $ratios ) {

			# Select single ROIs that pass the minimum calling thresholds
			my $cmd = "$::zcat $ratios | $::awk '{ if ((\$5 < $::upperDelCutoff && \$5 > $::lowerDelCutoff) || (\$5 > $::lowerDupCutoff)";
			$cmd .= "|| (\$5 > 0 && \$5 < 0.12)) { print \$0 }}' | $::sort -V >> $singleExonCnv";
			print "$cmd\n" if $::verbose;
			system($cmd);
		}
		# but if segmented multi-exon CNVs
		else {
			# We might have also single ROIs that pass the minimum calling thresholds
			$cmd = "$::zcat $ratios | $::awk '{ if ((\$5 < $::upperDelCutoff && \$5 > $::lowerDelCutoff)";
			$cmd.= " || (\$5 > $::lowerDupCutoff) || (\$5 > 0 && \$5 < 0.12)) { print \$0 }}' |";
			$cmd.= "$::bedtools intersect -a stdin -b $intersect -v |" if -s $intersect;
			$cmd.= "$::sort -V >> $singleTmp";
			print "$cmd\n" if $::verbose;
			system($cmd);

			# Get ROI level information for the segmented CNV
			my $str = `$::bedtools intersect -a $ratios -b $intersect -wo`;
			my @tmp = split (/\n/, $str);
			my %seen = ();
			foreach my $line (@tmp) {

				my @tmpLine = split (/\t/, $line);

				# These are the coordinates of the CNV
				my $segment = join ("\t", @tmpLine[6..8]);

				# Skipping already visited coordinates
				$seen{$segment}++;
				next if $seen{$segment} > 1;

				# Coordinates.
				# badEnd: in some cases DNAcopy outputs weird end coordinates. We will fix them afeterwards.
				my ($chr, $start, $badEnd) = split (/\t/, $segment);

				# GEt all ROIs that belong to the CNV
				my @targets   = grep ($_ =~/$segment/, @tmp);

				# Get a summarized tag (e.g first ROI and last ROI)
				# instead of outputting all ROI names
				my $affectedROIs = getAffectedROIs( $segment, \@targets );

				# Getting corrected End position
				my $lastEx  = $targets[-1];
				my @tmpLast = split (/\t/, $lastEx);
				my $End     = $tmpLast[2];

				# Number of overlapping ROIs
				my $nexons  = scalar @targets;

				# Define CNV type by hard thresholds
				my $cnvType  = $tmpLine[10] > 0.71 ? 'DUP' : 'DEL';

				# Assigning integer copy number
				my $copyNumber = int ($tmpLine[10]*2+0.5);

				my %HoZ = ();
				my @ArrGC = ();
				my @ArrMAP= ();
				my @signal2noise = ();
				my @ArrMeanZscore =  ();
				my @ArrRatiosSample = ();
				my @ArrRatiosControl = ();
				my @ArrRatiosOfftarget = ();

				open (TMP, ">", "$inputDir/PLOT_DATA/$sampName.$chr.$start.tmp.txt")
				|| die " ERROR: Unable to open $inputDir/PLOT_DATA/$sampName.$chr.$start.tmp.txt\n";

				# Then, for each ROI from the segmented CNV
				foreach my $target ( @targets ){

					my @tmp = split (/\t/, $target);
					my $coordinate  = join ("\t", @tmp[0..2]);
					my $coordPlusRoi= join ("\t", @tmp[0..3]);
					my $ratioSample = $tmp[4];

					my $str;
					if ($type eq 'on-target') {
						$str = `$::grep -P \'^$coordinate\' $ratioFile`;
					}
					else {
						$str = `$::zgrep -e \'$coordinate\' $ratioFile`;
					}
					print "$::grep -P ^$coordinate $ratioFile\n" if $::verbose;

					my @tmpStr = split (/\t/, $str);

					# Get mappability
					my $map = int($::ExonFeatures{$coordinate}{MAP});

					# Get GC content
					my $gc  = int($::ExonFeatures{$coordinate}{GC});

					# And store them
					push @ArrMAP, $map;
					push @ArrGC, $gc;

					my $meanControlRatio;
					my $sdControls;
					my $sdSample;
					my $zscoreSample;

					# Get copy ratios, std.dev and z-scores
					if ($type eq 'on-target') {
						for (my $i = 6; $i < @tmpStr; $i+=2) {

							if ($HashPosAll{$i} ne $sampName) {
								push @ArrRatiosControl, $tmpStr[$i];
							}
							else {
								push @ArrRatiosSample, $tmpStr[$i];
							}
						}
						$meanControlRatio =	$::sampleHash{$sampName}{ONTARGET_MEAN_RATIO};
						$sdSample         = $::sampleHash{$sampName}{ONTARGET_SD_RATIO};
						$sdSample = 0.01 if $sdSample == 0;
						$zscoreSample  = sprintf "%.3f", (($ratioSample - $meanControlRatio)/$sdSample);

						push @ArrMeanZscore, $zscoreSample;
					}
					elsif ($type eq 'off-target') {
						$meanControlRatio =	$::sampleHash{$sampName}{OFFTARGET_MEAN_RATIO};
						$sdSample         = $::sampleHash{$sampName}{OFFTARGET_SD_RATIO};
						$sdSample = 0.01 if $sdSample == 0;
						$zscoreSample     = sprintf "%.3f", (($ratioSample - $meanControlRatio)/$sdSample);
						push @ArrMeanZscore, $zscoreSample;
						push @ArrRatiosSample, $ratioSample;

						print TMP "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t$sampName\t$ratioSample\t$zscoreSample\t$sampName\n";
					}
					elsif ($type eq 'mixed') {
						$meanControlRatio =	$::sampleHash{$sampName}{ONOFF_MEAN_RATIO};
						$sdSample         = $::sampleHash{$sampName}{ONOFF_SD_RATIO};
						$sdSample = 0.01 if $sdSample == 0;
						$zscoreSample     = sprintf "%.3f", (($ratioSample - $meanControlRatio)/$sdSample);
						push @ArrMeanZscore, $zscoreSample;
						push @ArrRatiosSample, $ratioSample;

						print TMP "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t$sampName\t$ratioSample\t$zscoreSample\t$sampName\n";
					}
					# Calculating Z-score for controls
					for (my $i = 6; $i < @tmpStr; $i+=2) {
						if ($HashPosAll{$i} eq $sampName) {
							# COORDINATE	# SAMPLE	# ZSCORE
							$HoZ{$coordPlusRoi}{$HashPosAll{$i}}{CLASS}= $sampName;
							$HoZ{$coordPlusRoi}{$HashPosAll{$i}}{RATIO}= $tmpStr[$i];
							$HoZ{$coordPlusRoi}{$HashPosAll{$i}}{ZSCORE}= $zscoreSample;
						}
						else {
							my $zscoreControl  = sprintf "%.3f", (($tmpStr[$i] - $meanControlRatio)/$sdSample);
							$HoZ{$coordPlusRoi}{$HashPosAll{$i}}{CLASS}= "control";
							$HoZ{$coordPlusRoi}{$HashPosAll{$i}}{RATIO}= $meanControlRatio;
							$HoZ{$coordPlusRoi}{$HashPosAll{$i}}{ZSCORE}= $zscoreControl;
						}
					}
					push @ArrRatiosOfftarget, $tmp[4] if $type eq 'mixed';
					push @ArrRatiosOfftarget, $tmp[4] if $type eq 'off-target';
					push @signal2noise, $tmp[5] if $type eq 'on-target';
				}

				if ($type eq "on-target") {
					foreach my $coord ( natsort keys %HoZ) {
						foreach my $sample ( natsort keys  %{$HoZ{$coord}}) {
							print TMP "$coord\t$sample\t$HoZ{$coord}{$sample}{RATIO}\t$HoZ{$coord}{$sample}{ZSCORE}\t$HoZ{$coord}{$sample}{CLASS}\n";
						}
					}
				}
				%HoZ = ();
				close TMP;

				# Calulation of quality metrics and other annotations
				# Median absolute deviation of the segmented ratios
				my $MAD  = Utils::MAD(@ArrRatiosSample);

				# Mean mappability of the segment
				my $meanMap = Utils::meanArray(@ArrMAP);

				# Mean GC of the segment
				my $meanGC = Utils::meanArray(@ArrGC);

				my $meanSignal2Noise     = "";
				my $signal2noiseControls = "";
				my $meanRatioOfftarget   = "";
				if ($type eq 'on-target') {
					$meanSignal2Noise = Utils::meanArray(@signal2noise);
					$signal2noiseControls = Utils::signal2noise(@ArrRatiosControl);
				}
				else {
					($meanRatioOfftarget, $meanSignal2Noise) = Utils::signal2noise(@ArrRatiosOfftarget);
					$MAD = Utils::MAD(@ArrRatiosOfftarget);
					$signal2noiseControls = ".";
				}
				# Mean z-score of the segmented ratios
				my $meanZscore =
					sprintf "%.3f", (($tmpLine[10]-$::sampleHash{$sampName}{ONTARGET_MEAN_RATIO})/$::sampleHash{$sampName}{ONTARGET_SD_RATIO});

				my ($cipos, $ciend) = getCiposCiend($chr, $start, $End, \@::ROIarray);
				my $size = $End-$start;

				next if abs($meanZscore) < $::minZscore;
				next if $meanSignal2Noise < 5;
				print OUTFILE "$chr\t$start\t$End\tIMPRECISE;CIPOS=$cipos;CIEND=$ciend;SVTYPE=$cnvType;SVLEN=$size;GC=$meanGC;MAP=$meanMap;";
				print OUTFILE "GENE=$affectedROIs;REGIONS=$nexons;RRD=$tmpLine[10];MADRD=$MAD;CN=$copyNumber;SNR=$meanSignal2Noise;";
				print OUTFILE "SNRC=$signal2noiseControls;ZSCORE=$meanZscore\n";
			}

			$sampName = basename($sampName);
			$sampName =~s/.on_off//;
			$sampName =~s/;//;

			if ( $type eq 'mixed') {
				$breakFile = "$inputDir/ON_TARGET/$sampName/$sampName.breakpoints.bed";

				$cmd = "$::cat $::outDir/OFF_TARGET/$sampName.offtarget.bed $::bed | $::sort -V > $::outDir/$sampName.mixed.bed";
				print "$cmd\n" if $::verbose;
				system $cmd;

				$bedfile = "$::outDir/$sampName.mixed.bed";
			}
			if ( $type eq 'off-target') {
				$breakFile = "$inputDir/../ON_TARGET/$sampName/$sampName.breakpoints.bed";
				$bedfile   = "$::outDir/OFF_TARGET/$sampName.offtarget.bed";
			}

			my $cmd = "$::mergeSegments $rawMultipleExonCnv $bedfile > $mergedMultipleExonCnv";
			print "$cmd\n" if $::verbose;
			system $cmd;

			$cmd = "$::mergeCnvBreaks $breakFile $mergedMultipleExonCnv | $::sort -V | $::uniq > $multipleCnvAndBreaks";
			print "$cmd\n" if $::verbose;
			system $cmd;

			# Plotting stuff to be deleted
			open (MERGED, "<", $mergedMultipleExonCnv ) || die "ERROR: Unable to open $mergedMultipleExonCnv\n";
			while (my $line=<MERGED>) {

				chomp $line;

				my @tmp   = split (/\t/, $line);
				my $chr   = $tmp[0];
				my $start = $tmp[1];
				my $end   = $tmp[2];
				my $coordinates = "$chr:$start-$end";
				my $size  = $end-$start;
				my $roi_info= $tmp[3];
				my @info  = split (/;/, $roi_info);

				my ( $regions ) = grep ($_=~/REGIONS=/, @info);
				$regions =~s/REGIONS=//;

				my ( $cnvtype ) = grep ($_=~/SVTYPE=/, @info);
				$cnvtype =~s/SVTYPE=//;

				my $candidateTxt        = "$inputDir/PLOT_DATA/$sampName.$chr.$start.$end.txt";
				my $candidateBed        = "$inputDir/PLOT_DATA/$sampName.$chr.$start.$end.bed";
				my $candidateRatios     = "$inputDir/PLOT_DATA/$sampName.$chr.$start.$end.ratios.bed";
				my $candidateRatiosFinal= "$inputDir/PLOT_DATA/$sampName.$chr.$start.$end.final_ratios.bed";

				# Bed temporal file used for intersection
				open (TXT, ">", $candidateTxt);
				print TXT "$chr\t$start\t$end\n";
				close TXT;

				if ($::plotLargeCNV) {

					# Plotting large segmented CNVs for gene-panel and all exome CNVs
					if ( $regions > 15 || $analysis eq 'exome' || $type eq 'off-target' || $type eq 'mixed') {

						if (!-e "$inputDir/$sampName.$chr.$start.$end.png") {

							print " INFO: Plotting $regions total regions from sample $sampName => $chr:$start-$end\n";

							$cmd = "$::bedtools intersect -a $inputDir/PLOT_DATA/$sampName.$chr.$start.tmp.txt -b $candidateTxt -wa > $candidateRatios";
							print "$cmd\n" if $::verbose;
							system($cmd);

							#Plot::plotLongSegment("$inputDir/PLOT_DATA/$sampName.$chr.$start.$end.tmp.txt",
							#	$inputDir, $sampName, $bedfile, $chr, $start, $end);

							Utils::loadLongPlot2Json($candidateRatios, $coordinates, $sampName);

							# Move ratio file to sample directory
							move $candidateRatios, "$inputDir/PLOT_DATA/";
							Utils::compressFile($candidateRatios);

							unlink ("$inputDir/PLOT_DATA/$sampName.$chr.$start.tmp.txt") if -e "$inputDir/PLOT_DATA/$sampName.$chr.$start.tmp.txt";
							unlink ("$inputDir/$sampName.$chr.$start.$end.tmp.txt") if -e "$inputDir/$sampName.$chr.$start.$end.tmp.txt";
						}
					}
					# Plotting small segmented CNVs for gene-panel analysis
					if ( ( $regions > 1 && $regions < 15) && ($analysis eq 'gene-panel') && ($type eq 'on-target')) {
						print " INFO: Plotting $regions total regions from sample $sampName => $chr:$start-$end\n";

						$cmd = "$::bedtools intersect -a $inputDir/$::outName.norm_bylib_coverage_noheader.bed -b $candidateTxt > $candidateBed";
						print "$cmd\n" if $::verbose;
						system $cmd;

						open (C, "<", $candidateBed) || die " ERROR: Unable to open $candidateBed\n";
						open (OUT, ">", $candidateRatios) || die " ERROR: Unable to open $candidateRatios\n";
						while (my $line=<C>) {
							chomp $line;

							my @tmp = split (/\t/, $line);
							my $index = $sampleIndex{$sampName};
							chomp $index;

							my $sample_value = $tmp[$index];
							my @control_values;

							for (my $i=5;$i<@tmp;$i++){
								next if $i == $index;
								push @control_values, $tmp[$i];
							}
							my $control_mean = Utils::meanArray(@control_values);
							if ($cnvtype eq 'DUP') {
								print OUT "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t$sample_value\t1_$sampName\n";
								print OUT "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t$control_mean\t2_Controls\n";
							}
							else {
								print OUT "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t$sample_value\t2_$sampName\n";
								print OUT "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t$control_mean\t1_Controls\n";
							}
						}
						close OUT;
						close C;

						Utils::loadShortPLot2Json($candidateRatios, $coordinates, $sampName);

						fillGaps($candidateRatios, $inputDir, $sampName, $cnvtype, $chr, $start, $end,);
						unlink ($candidateRatios);

						rename $candidateRatiosFinal, $candidateRatios;
						#Plot::plotShortSegment($candidateRatios, $inputDir, $sampName, $bedfile, $chr, $start, $end, $cnvtype);

						Utils::compressFile($candidateRatios);
						move "$candidateRatios.gz", "$inputDir/PLOT_DATA/";
					}
			    }
				#unlink ($candidateRatiosFinal);
				#unlink ($candidateBed);
				#unlink ($candidateRatios);
				#unlink ($candidateTxt);
		    }
			close MERGED;

		if (!-e $singleExonCnv) {
			rename $singleTmp, $singleExonCnv;
		}
		else {
			$cmd = "$::cat $singleTmp $singleExonCnv | $::bedtools intersect -a stdin -b $rawMultipleExonCnv -v | $::sort -V | $::uniq > $::outDir/merged.tmp.single.txt";
			system $cmd;

			print "$cmd\n" if $::verbose;
			unlink ($singleExonCnv);

			$cmd = "mv $::outDir/merged.tmp.single.txt $singleExonCnv";
			system $cmd;
			print "$cmd\n" if $::verbose;

			unlink("$::outDir/merged.tmp.single.txt");
		}

		unlink (glob ("$inputDir/intersect*"));
		unlink (glob ("$inputDir/*tmp.txt"));
		unlink ($singleTmp);

	}
	if ($type eq 'off-target'|| $type eq 'mixed') {
		Utils::compressFile($ratioFile) if -s $ratioFile;
	}
  }

  Utils::compressFile($ratioFile) if -s $ratioFile;
}

 ################################################
 #	        Analysing single-exon CNVs          #
 ################################################

 sub callSingleExonCNV {

  my $outputDir = shift;
  my $analysis  = shift;

  my @singleCallFiles = glob ("$outputDir/*raw.single.exon.cnv.bed");
  my $ratiosFile      = $::HoF{NORM_COVERAGE_ON} . ".gz";
  my $ungzRatioFile   = $::HoF{NORM_COVERAGE_ON};

  if (!-e $ratiosFile && !-e $ungzRatioFile) {
	  print " WARNING: non-existent $ratiosFile. Skipping single-exon CNV analysis\n";
	  return 0;
  }

  Utils::decompressFile($ratiosFile) if -e $ratiosFile;

  # Selecting sample names from header
  my $str = `$::cat $ungzRatioFile | $::head -1`;
  chomp $str;
  my @tmpStr  = split (/\t/, $str);
  my @samples = @tmpStr[5..@tmpStr-1];

  # Emptying reference id file
  foreach my $sample ( sort keys %::sampleHash ) {
	@{$::sampleHash{$sample}{REF_ID}} = ();
  }
  # Select columns that correspond to the sample baseline;
  foreach my $sample( @samples ) {
	if ($::doCaseControl) {
		next if $::sampleHash{$sample}{CONTROL};
	}
	my @idx = ();
    my $n = 0;

	foreach my $element ( @{$::sampleHash{$sample}{REFERENCE}} ) {

		next if $element eq $sample;
		if ($element eq 'none') {
 			@{$::sampleHash{$sample}{REF_ID}} = ();
		}
		else {
			$n++;
			if ($n > $::maxSampleSizeCluster) {
				$n = 0;
				last;
			}
			for (my $index = 0; $index <@samples; $index++) {
				if ($element eq $samples[$index]) {
					push @{$::sampleHash{$sample}{REF_ID}}, $index;
				}
			}
		}
	}
  }

  open (MASTER_FILE, ">", "$outputDir/master_single_calls.bed");
  foreach my $file (@singleCallFiles) {
	my $sampleName = basename($file);
	$sampleName =~s/.raw.single.exon.cnv.bed//;
	open (IN,"<", $file) || die " ERROR: Unable to open $file\n";
	while (my $line=<IN>){
		chomp $line;
		my @tmp = split("\t", $line);
		print MASTER_FILE "$line\t$sampleName\n";
	}
	close IN;
	unlink($file);
  }
  close MASTER_FILE;

  print " INFO: Calling single-exon CNVs\n";

  my $cmd = "$::Rscript $::analyzeSingleExonR --raw_calls $outputDir/master_single_calls.bed --coverage_data $ungzRatioFile ";
  $cmd.= " --references $::HoF{REFERENCES} --min_zscore $::minZscore --min_s2n 5 --lower_del_cutoff $::lowerDelCutoff ";
  $cmd.= " --upper_del_cutoff $::upperDelCutoff --lower_dup_cutoff $::lowerDupCutoff --output_dir $outputDir";
  system($cmd);

  my @tmpRdata = glob("$outputDir/*single.exon.cnv.bed");
  if (!@tmpRdata) {
	  print " ERROR: Missing single exon tmp data\n";
  }

  foreach my $sample( sort keys %::sampleHash ) {

	my $tmpRcalls = "$outputDir/$sample.single.exon.cnv.bed";
	if (!-e $tmpRcalls) {
		next;
	}

	my $toMergeSegs          = "$outputDir/$sample.tomerge.single.exon.cnv.bed";
	my $mergedAdjacentRois   = "$outputDir/$sample.merged.multiple.exon.cnv.bed";
	my $mergedSingleSegments = "$outputDir/$sample.merged.single.multiple.exon.cnv.bed";
	my $mergedCnvBreakpoints = "$outputDir/$sample.merged.breaks.single.multiple.exon.cnv.bed";
	my $breakpointCalls      = "$outputDir/$sample/$sample.breakpoints.bed";

	if (-e $mergedAdjacentRois ) {
		$cmd = "$::cat $mergedAdjacentRois > $toMergeSegs";
		system $cmd;
	}

  open (IN, "<", $tmpRcalls) || die " ERROR: Unable to open $tmpRcalls\n";
	open (OUT, ">>", $toMergeSegs) || die " ERROR: Unable to open $toMergeSegs\n";
  while(my $line=<IN>) {
    chomp $line;
		my @tmp = split("\t", $line);

		my $chr   = $tmp[0];
		my $start = $tmp[1];
		my $end   = $tmp[2];

		# Dumping results to BED
		my $genename = $tmp[4];
		$genename =~s/;/,/;

		# Get GC and Mappability
		my $GC  = int ($::ExonFeatures{"$chr\t$start\t$end"}{GC});
		my $MAP = int ($::ExonFeatures{"$chr\t$start\t$end"}{MAP});

		# Calculate confidence interval for imprecise events
		my ($cipos, $ciend) = getCiposCiend($chr, $start, $end, \@::ROIarray);
		my $cnv_type        = $tmp[3];
		my $size            = $end-$start;
		my $testSampleRatio = $tmp[5];
		my $s2nTest         = $tmp[6];
		my $s2nControls     = $tmp[7];
		my $zscore          = $tmp[8];

		#chrX	31140035	31140047	DEL	NM_004006_78_79;DMD	0.653179190751445	20.3063211047918	13.3728451863809	-4.613666036998
		print OUT "$chr\t$start\t$end\tIMPRECISE;CIPOS=$cipos;CIEND=$ciend;SVTYPE=$cnv_type;SVLEN=$size;GENE=$genename;GC=$GC;";
		print OUT "MAP=$MAP;REGIONS=1;RRD=$testSampleRatio;MADRD=.;SNR=$s2nTest;SNRC=$s2nControls;ZSCORE=$zscore\n";
	}
	close IN;
	close OUT;
	unlink($tmpRcalls);

	if (-e $toMergeSegs) {

		# Re-segment single exon CNVs that may belong to longer CNVs
		$cmd = "perl $::mergeSegments $toMergeSegs $::bed > $mergedSingleSegments";
		system $cmd;
		print "$cmd\n" if $::verbose;

		$mergedAdjacentRois = $mergedSingleSegments;
	}

	# Merge CNV and Breakpoint predictions if they point to the same variant
	$cmd = "perl $::mergeCnvBreaks $breakpointCalls $mergedAdjacentRois ";
	$cmd.= "| $::sort -V | $::uniq > $mergedCnvBreakpoints";
	system $cmd;
	print "$cmd\n" if $::verbose;

	unlink($toMergeSegs);
  }
  Utils::compressFile($ungzRatioFile);
  return;

  foreach my $file (@singleCallFiles) {

    my $test_sample = basename($file);
	$test_sample =~s/.single.exon.cnv.bed//;

	my $toMergeSegs          = "$outputDir/TMP.tomerge.segments.$test_sample.bed";
	my $mergedAdjacentRois   = "$outputDir/cnv_calls.merged.$test_sample.bed";
	my $mergedSingleSegments = "$outputDir/merged.single_and_segments.$test_sample.bed";
	my $mergedCnvBreakpoints = "$outputDir/Breaks_and_CNV.$test_sample.bed";

	my %seen = ();
	open (IN, "<", $file) || die " ERROR: Cannot open $file\n";
	while (my $line=<IN>) {
		chomp $line;

		my @tmp     = split (/\t/, $line);
		my $chr     = $tmp[0];
		my $start   = $tmp[1];
		my $end     = $tmp[2];
		my $size    = $tmp[2]-$tmp[1];
		my $pattern = $tmp[3];

		my $coordinates = "$chr:$start-$end";

		$seen{"$chr\t$start\t$end"}++;
		next if $seen{"$chr\t$start\t$end"} > 1;

		my $str = `LC_ALL=C $::grep '$pattern' $ungzRatioFile`;
		chomp $str;

		print "LC_ALL=C $::grep '$pattern' $ungzRatioFile\n" if $::verbose;

		# Arrays for ratios and signal-2-noise
		my @arrOfTest = ();
		my @arrOfCont = ();
		my @arrOfRatios = ();

		# Now reading base by base
		my @tmpStr = split (/\n/, $str);
		my @matrix = ();

		foreach my $ntd ( @tmpStr ) {

			my @tmp = split (/\t/, $ntd);
			my $coord = join ("\t", @tmp[0..3]);

			my $ntdStart = $tmp[1];
			my $ntdEnd   = $tmp[2];
			if ($ntdStart < $start || $ntdEnd > $end) {
				next;
			}

			my $j = 0;
			for (my $i = 5; $i < @tmp; $i++) {
				my @tmpRefs;
				my $sample = @samples[$j];
				$j++;
				if (  $sample ne $test_sample && !grep( /$sample/, @{$::sampleHash{$test_sample}{REFERENCE}} ) ) {
					next;
				}
				#exit;
				foreach my $idx ( @{$::sampleHash{$sample}{REF_ID}} ) {
					push @tmpRefs, $tmp[5+$idx];
				}

				if (@tmpRefs) {
					my ($medianBaseline, $s2n) = Utils::signal2noise(@tmpRefs);
					my $ratio =  $medianBaseline > 0 ? $tmp[$i]/$medianBaseline : "0.00";
					my $group;
					if ($sample eq $test_sample) {
						push @arrOfTest, $ratio;
						$group = $test_sample;
					}
					else {
						push @arrOfCont, $ratio;
						$group = $sample;
					}
					push @arrOfRatios, $ratio;
					push @matrix, "$coord\t$sample\t$ratio\t$group\n";
				}
			}
		}

		# Check mean signal-2-noise of control samples
		my ($testSampleRatio, $s2nTest) = Utils::signal2noise(@arrOfTest);
		my $sd   = Utils::std(@arrOfCont);
		my $MAD  = Utils::MAD(@arrOfRatios);

   	    $sd = 0.01 if $sd == 0;
		my ($mean, $s2nControls)  = Utils::signal2noise(@arrOfCont);

		if (scalar @arrOfCont <= 1) {
			$s2nControls = 9;
		}

		# Calculating z-score
		#z = (x – μ) / σ
		my $zscore = sprintf "%.3f", (($testSampleRatio - $mean)/$sd);

		# Check mean ratio of test sample
		my $tmpZscore = $zscore < 0 ? $zscore*-1 : $zscore;

		if ($::verbose) {
			print "$pattern=>s2nTest($s2nTest)\ts2nControls($s2nControls)\tZ-score($tmpZscore)\tSampleRatio=>$testSampleRatio\n";
		}

		if (-e $mergedAdjacentRois ) {
			$cmd = "$::cat $mergedAdjacentRois > $toMergeSegs";
			system $cmd;
		}

		# cnv_calls.merged contains segmented CNVs and single-exon CNVs
		open (TMP, ">>", $toMergeSegs) || die " ERROR: Cannot open $toMergeSegs\n";
		open (BED, ">>", $mergedAdjacentRois) || die " ERROR: Cannot open $mergedAdjacentRois\n";

		# Check if single CNV has enough quality to be reported as a candidate
		my $do_output = 0;
		my $cnv_type = '.';

		if ($s2nControls > 5 && $s2nTest > 5 && abs($tmpZscore) >= $::minZscore &&
		   ($testSampleRatio < $::sampleHash{$test_sample}{LOWER_LIMIT} ||
		   $testSampleRatio > $::sampleHash{$test_sample}{UPPER_LIMIT})){

			if ($testSampleRatio >= $::lowerDelCutoff && $testSampleRatio < $::upperDelCutoff) {
				$cnv_type = "DEL";
				$do_output = 1;
			}
			elsif ($testSampleRatio >= 0.00 && $testSampleRatio < 0.12) {
				$cnv_type = "DEL";
				$do_output = 1;
			}
			elsif ($testSampleRatio > $::lowerDupCutoff) {
				$cnv_type = "DUP";
				$do_output = 1;
			}
			else {
				next;
			}
		}
		if ($do_output) {

			print " INFO: Plotting $pattern on sample $test_sample\n";
			my $baseRatios =
				"$outputDir/$test_sample.$tmp[0].$tmp[1].$tmp[2].ratios.bed";

			# Writing ratios to an output file
			open (BASERATIOS, ">", $baseRatios) || die " ERROR: Cannot open $baseRatios\n";
			foreach my $line (@matrix) {
				chomp $line;
				print BASERATIOS "$line\n";
			}
			close BASERATIOS;

			# Load ratios to Json
			Utils::loadSingleExon2Json($baseRatios, $coordinates, $test_sample);
			if ($::plotSingleExon) {

				# Plotting single exon CNV
				Plot::plotSingleExon2($outputDir, $baseRatios, $pattern, $test_sample );

				#Plot::plotSingleExon($outputDir, $baseRatios, $pattern, $test_sample );
				move $baseRatios, "$outputDir/PLOT_DATA/";
				Utils::compressFile($baseRatios);
			}

			# Dumping results to BED
			my $genename = $tmp[3];
			$genename =~s/;/,/;

			# Get GC and Mappability
			my $GC  = int ($::ExonFeatures{"$tmp[0]\t$tmp[1]\t$tmp[2]"}{GC});
			my $MAP = int ($::ExonFeatures{"$tmp[0]\t$tmp[1]\t$tmp[2]"}{MAP});

			# Calculate confidence interval for imprecise events
			my ($cipos, $ciend) = getCiposCiend($tmp[0], $tmp[1], $tmp[2], \@::ROIarray);

			print TMP "$tmp[0]\t$tmp[1]\t$tmp[2]\tIMPRECISE;CIPOS=$cipos;CIEND=$ciend;SVTYPE=$cnv_type;SVLEN=$size;GENE=$genename;GC=$GC;MAP=$MAP;REGIONS=1;RRD=$testSampleRatio;MADRD=$MAD;SNR=$s2nTest;SNRC=$s2nControls;ZSCORE=$zscore\n";
			print BED "$tmp[0]\t$tmp[1]\t$tmp[2]\tIMPRECISE;CIPOS=$cipos;CIEND=$ciend;SVTYPE=$cnv_type;SVLEN=$size;GENE=$genename;GC=$GC;MAP=$MAP;REGIONS=1;RRD=$testSampleRatio;MADRD=$MAD;SNR=$s2nTest;SNRC=$s2nControls;ZSCORE=$zscore\n";
		}

		close BED;
		close TMP;
		@matrix = ();
	}
	close IN;

	if (-e $toMergeSegs) {

		$cmd = "perl $::mergeSegments $toMergeSegs $::bed > $mergedSingleSegments";
		system $cmd;

		# Catching additional segments that were not properly segmented before
		writePlotDataFromSegments($mergedSingleSegments, $test_sample, $outputDir);

		$mergedAdjacentRois = $mergedSingleSegments;
		#unlink ($toMergeSegs);
	}

	# Merge CNV and Breakpoint predictions if they point to the same variant
	$cmd = "perl $::mergeCnvBreaks $outputDir/$test_sample/$test_sample.breakpoints.bed $mergedAdjacentRois ";
	$cmd.= "| $::sort -V | $::uniq > $mergedCnvBreakpoints";
	system $cmd;
	print "$cmd\n" if $::verbose;

	#unlink ($mergedSingleSegments);
 }
 #unlink("$outputDir/$::outName.norm_bylib_coverage_noheader.bed");
 Utils::compressFile($ungzRatioFile);
}

#####################################
sub writePlotDataFromSegments {

	my $infile   = shift;
	my $sample   = shift;
	my $inputDir = shift;

	my $sampName = basename($sample);
	# Removing header for intersection purposes

	open (IN, "<", $infile) || die " ERROR: Unable to open $infile\n";
	while (my $line=<IN>) {

		chomp $line;

		my @tmp = split (/\t/, $line);
		my @info = split (/;/, $tmp[3]);

		my ($regions) = grep ($_=~/REGIONS=/, @info) ? grep ($_=~/REGIONS=/, @info) : '.';
		$regions =~s/REGIONS=//;

		my ($svtype) = grep ($_=~/SVTYPE=/, @info) ? grep ($_=~/SVTYPE=/, @info) : '.';
		$svtype =~s/SVTYPE=//;

		if ($regions > 1  && $regions < 15 && !-e "$inputDir/PLOT_DATA/$sampName.$tmp[0].$tmp[1].$tmp[2].regions.bed.gz") {
			my $bed = "$inputDir/PLOT_DATA/$tmp[0].$tmp[1].$tmp[2].bed";

			open (BED, ">", $bed);
			print BED "$tmp[0]\t$tmp[1]\t$tmp[2]\n";
			close BED;

			my $str = `$::bedtools intersect -a $inputDir/$::outName.norm_bylib_coverage_noheader.bed -b $bed`;
			print "$::bedtools intersect -a $inputDir/$::outName.norm_bylib_coverage_noheader.bed -b $bed" if $::verbose;
			chomp $str;
			unlink($bed);

			my $candidateRatios = "$inputDir/PLOT_DATA/$sampName.$tmp[0].$tmp[1].$tmp[2].ratios.bed";

			open (OUT, ">", $candidateRatios) || die " ERROR: Unable to open $candidateRatios\n";
			my @tmpStr = split (/\n/, $str);
			foreach my $line (@tmpStr) {
				chomp $line;

				my @tmp = split (/\t/, $line);
				my $index = $sampleIndex{$sampName};
				chomp $index;

				my $sample_value = $tmp[$index];
				my @control_values;

				for (my $i=5;$i<@tmp;$i++){
					next if $i == $index;
					#print "$line\n";
					push @control_values, $tmp[$i];
				}
				my $control_mean = Utils::meanArray(@control_values);
				if ($svtype eq 'DUP') {
					print OUT "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t$sample_value\t1_$sampName\n";
					print OUT "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t$control_mean\t2_Controls\n";
				}
				else {
					print OUT "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t$sample_value\t2_$sampName\n";
					print OUT "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t$control_mean\t1_Controls\n";
				}
			}
			close OUT;

			my $cmd = "$::gzip -f $candidateRatios" if -e $candidateRatios;
			system $cmd;
		}

	}
	close IN;
}
#####################################
# Add single CNV when segmented CNV has a single region
# This is for plotting purposes
sub addSingleCNV {

	my $singleCnvFile = shift;
	my $chr   = shift;
	my $start = shift;
	my $end   = shift;
	my $exon  = shift;
	my $ratio = shift;

	open (TMP, ">","$::outDir/$chr.$start.$end.add.txt");
	print TMP "$chr\t$start\t$end\t$exon\t$ratio\t$ratio\n";
	close TMP;

	my $cmd = "$::cat $singleCnvFile $::outDir/$chr.$start.$end.add.txt | $::sort -V | $::uniq > $::outDir/$chr.$start.$end.added.txt";
	system $cmd;

	unlink($singleCnvFile);
	$cmd = "mv $::outDir/$chr.$start.$end.added.txt $singleCnvFile";
	system $cmd;

	unlink("$::outDir/$chr.$start.$end.add.txt", "$::outDir/$chr.$start.$end.added.txt");

}

#####################################
sub fillGaps {

    my $input      = shift;
    my $inputDir   = shift;
    my $sampleName = shift;
	my $cnvtype    = shift;
	my $chr        = shift;
	my $start      = shift;
	my $end        = shift;

    my $outfile = $input;
    $outfile =~s/.ratios.bed/.final_ratios.bed/;
    my $noFirstLine = $input;
    $noFirstLine =~s/.bed/.nofirstline.bed/;

    my $ini_coord = `$::head -2 $input | $::tail -1`;
    chomp $ini_coord;

    my @tmpCoord = split (/\t/, $ini_coord);
    my $chr   = $tmpCoord[0];
    my $start = $tmpCoord[1];
    my $numberSamples = scalar @tmpCoord -4;

    my $end_coord = `$::tail -1 $input`;
    chomp $end_coord;

    my @tmpEnd = split (/\t/, $end_coord);
    my $end = $tmpEnd[2];

    my $tmpStart;
    my $tmpFinal;

    if ($start > $end) {
        $tmpStart = $end;
        $tmpFinal = $start;
        $start = $tmpStart;
        $end = $tmpFinal;
    }
    my $entire_bed = "$inputDir/$sampleName.$chr.$start.$end.entire_bed.bed";
    open (BED, ">", $entire_bed) || die " ERROR: Unable to open $entire_bed\n";
    for (my $i = $start; $i <= $end; $i++) {
        my $sum = $i+1;
        print BED "$chr\t$i\t$sum\n";
    }
    close BED;

    `$::bedtools intersect -a $entire_bed -b $input -v > $inputDir/$sampleName.$chr.$start.$end.no_coverage.bed`;
    `$::cat $inputDir/$sampleName.$chr.$start.$end.no_coverage.bed $input | $::uniq | $::sort -V > $inputDir/$sampleName.$chr.$start.$end.merged.tmp.bed`;

    open ( IN, "<", "$inputDir/$sampleName.$chr.$start.$end.merged.tmp.bed");
    open (OUT, ">", $outfile);

    while (my $line=<IN>) {
        chomp $line;
        my @tmp = split (/\t/, $line);
        if (@tmp > 3) {
            print OUT "$line\n";
        }
        else {
			if ($cnvtype eq 'DEL') {
				print OUT "$tmp[0]\t$tmp[1]\t$tmp[2]\tND\t0\t1_Controls\n";
				print OUT "$tmp[0]\t$tmp[1]\t$tmp[2]\tND\t0\t2_$sampleName\n";
			}
			else {
				print OUT "$tmp[0]\t$tmp[1]\t$tmp[2]\tND\t0\t2_Controls\n";
				print OUT "$tmp[0]\t$tmp[1]\t$tmp[2]\tND\t0\t1_$sampleName\n";
			}
        }
    }
    close OUT;
    unlink ("$inputDir/$sampleName.$chr.$start.$end.no_coverage.bed");
    unlink ("$inputDir/$sampleName.$chr.$start.$end.merged.tmp.bed");
    unlink $entire_bed;
}
###########################
sub getBafSegmentFile {

	my $dirClass = shift;
	my $sample   = shift;

	my $bafSegmentFile;
	if ($dirClass =~/ON_TARGET/) {
		$bafSegmentFile = $::sampleHash{$sample}{BAF_ONTARGET};
	}
	elsif ($dirClass =~/OFF_TARGET/) {
		$bafSegmentFile = $::sampleHash{$sample}{BAF_OFFTARGET};
	}
	else {
		$bafSegmentFile = $::sampleHash{$sample}{BAF};
	}
	return $bafSegmentFile;
}
###########################
# Input: single BED file with sample CNV calls. Note that input BED is the final file before VCF creation
# Output: single BED file with BAF-annotated CNV calls

sub appendVAF {

	my $outputDir = shift;

	#print "Appending BAF $outputDir\n";
	foreach my $sample ( natsort keys %::sampleHash ) {

		my $inputBed  = "$outputDir/$sample.CNV.bed";
		next if !-e $inputBed;

		my $outputBed = $inputBed;
		$outputBed =~s/.bed/.BAF.bed/;
		my $bafSegmentFile = getBafSegmentFile($outputDir, $sample);

		my %regionDict = ();
		if (-s $inputBed && -s $bafSegmentFile ) {

			my $str = `$::bedtools intersect -a $inputBed -b $bafSegmentFile -wa -wb`;
			chomp $str;
			my @tmpStr = split (/\n/, $str);

			foreach my $line (@tmpStr) {
				my @tmp = split (/\t/, $line);
				my $region = join ("\t", @tmp[0..3]);
				push( @{ $regionDict { $region } }, $tmp[7]);
			}
		}

		# Calculating median BAF value per CNV
		open (IN, "<", $inputBed)   || die " ERROR: Unable to open $inputBed\n";
		open (OUT, ">", $outputBed) || die " ERROR: Unable to open $outputBed\n";
		while (my $variant=<IN>) {
			chomp $variant;
			my $medianBAF = '.';
			if ($regionDict{$variant}) {
				$medianBAF = Utils::medianRefVar($regionDict{$variant});
			}
			# replace BAF tag and append new value
			$variant =~s/;BAF=\.//;
			$variant.= ";BAF=$medianBAF";
			print OUT "$variant\n";
		}
		close OUT;
		close IN;

		# Removing input file, and renaming output with the same name
		unlink $inputBed;
		rename $outputBed, $inputBed;
	}
}

 #################   DEPRECATED  #####################
 #	  Filtering CNV predictions when:		         #
 #	- Intron is small (less than 200 bp)		     #
 #	- Intron has enough coverage ( >30X)		     #
 #	- Intron does not have any breakpoint support	 #
 #####################################################

 sub filterCNVsByBreakpoint {

   #print " INFO: Filtering CNVs\n";

   foreach my $sample (	sort keys %::sampleHash ) {

	my $cnv_calls = "$::outDir/ON_TARGET/Breaks_and_CNV.$sample.bed";
	my $breakpoint_calls = "$::outDir/ON_TARGET/$sample/$sample.breakpoints.bed";

	next if !-e $breakpoint_calls;

	my $filtered_calls = "$::outDir/ON_TARGET/FILTERED_CNV_CALLS.txt";
	open (FILTERED_CNV, ">>", $filtered_calls) || die " ERROR: Unable to open $filtered_calls\n";

	if ( !-e "$::outDir/ON_TARGET/Breaks_and_CNV.$sample.bed") {
		open (CLEAN_VCF, ">", "$::outDir/ON_TARGET/$sample.CNV.bed") || die " ERROR: Unable to open $::outDir/ON_TARGET/$sample.CNV.bed\n";
		next;
	}

	open (CLEAN_VCF, ">", "$::outDir/ON_TARGET/$sample.CNV.bed") || die " ERROR: Unable to open $::outDir/ON_TARGET/$sample.CNV.bed\n";
	open (CNV, "$::cat $cnv_calls | $::sort -V |") || die " ERROR: Unable to open $cnv_calls\n";

	while (my $line=<CNV>) {
		chomp $line;

   		my @tmp  = split (/\t/, $line);
		my @info = split (/;/, $tmp[3]);
		my ($regions) = grep ($_=~/REGIONS=/, @info);

		if ($regions) {
			$regions =~s/REGIONS=//;

			# Dumping results ifs segmented call
			if ($regions > 1) {
				print CLEAN_VCF "$line\n";
			}
			# If single exon call we still check the presence of short intronic regions with high mappability and coverage
			else {
				my @tmp = split (/\t/, $line);
				my $coordinate = `$::grep '$tmp[0]\t$tmp[1]\t$tmp[2]' $::bed`;
				chomp $coordinate;

				if ($::ExonFeatures{$coordinate}{PREVIOUS_INTRON}) {

					my $intron = $::ExonFeatures{$coordinate}{PREVIOUS_INTRON};
					my $previous_map = $::IntronFeatures{$intron}{MAP};
					my $previous_gc  = $::IntronFeatures{$intron}{GC};
					my $previous_cov = $::IntronFeatures{$intron}{$sample};

					if ($previous_map >= 90 && $previous_cov >= 20) {
						# Checking if breakpoint is present
						open (TMP, ">", "$::outDir/ON_TARGET/call.tmp.txt") || die " ERROR: Unable to open $::outDir/ON_TARGET/call.tmp.txt\n";
						my $modCall = "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]";
						print TMP "$modCall\n";

						if (-z "$::outDir/ON_TARGET/$sample/$sample.breakpoints.bed") {
							print FILTERED_CNV "$line\tABSENCE_OF_BREAKPOINT\t$intron;MAP=$previous_map;GC=$previous_gc;MEAN_COV=$previous_cov\n";
							next;
						}

						my $str = `$::bedtools intersect -a $::outDir/ON_TARGET/call.tmp.txt -b $::outDir/ON_TARGET/$sample/$sample.breakpoints.bed -f 0.5 -r -wo`;
						chomp $str;

						if ($str) {
							print CLEAN_VCF "$str\n";
						}
						else{
							print FILTERED_CNV "$line\tABSENCE_OF_BREAKPOINT\t$intron;MAP=$previous_map;GC=$previous_gc;MEAN_COV=$previous_cov\n";
						}
						unlink ("$::outDir/ON_TARGET/call.tmp.txt");
						next;
					}
				}
				elsif ($::ExonFeatures{$coordinate}{NEXT_INTRON}) {

					my $intron   = $::ExonFeatures{$coordinate}{NEXT_INTRON};
					my $next_map = $::IntronFeatures{$intron}{MAP};
					my $next_gc  = $::IntronFeatures{$intron}{GC};
					my $next_cov = $::IntronFeatures{$intron}{$sample};

					if ($next_map >= 90 && $next_cov >= 30) {
						# Checking if breakpoint is present
						open (TMP, ">", "$::outDir/ON_TARGET/call.tmp.txt") || die " ERROR: Unable to open $::outDir/ON_TARGET/call.tmp.txt\n";
						my $modCall = "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]";
						print TMP "$modCall\n";

						if (-z "$::outDir/ON_TARGET/$sample/$sample.breakpoints.bed") {
							print FILTERED_CNV "$line\tABSENCE_OF_BREAKPOINT\t$intron;MAP=$next_map;GC=$next_gc;MEAN_COV=$next_cov\n";
							next;
						}

						my $str = `$::bedtools intersect -a $::outDir/ON_TARGET/call.tmp.txt -b $::outDir/ON_TARGET/$sample/$sample.breakpoints.bed -f 0.5 -r -wo`;
						chomp $str;
						if ($str) {
							print CLEAN_VCF "$str\n";
						}
						else{
							print FILTERED_CNV "$line\tABSENCE_OF_BREAKPOINT\t$intron;MAP=$next_map;GC=$next_gc;MEAN_COV=$next_cov\n";
						}
						unlink ("$::outDir/ON_TARGET/call.tmp.txt");
						next;
					}
				}
				else {
					print CLEAN_VCF "$line\n";
				}
			}
		}
	}
	}
 }


1;
