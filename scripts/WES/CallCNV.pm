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

  my $inputDir = shift; # on-target, off-target or outDir
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

    next if !-e $finalCalls;
		# if (-e $finalCalls) {

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
						$score = scoreSingleExon($svtype, $nCalls, $ratioStd, $pctCalls,
              $corr, $ratio, $s2n, $s2nc, $cn, $regions);
					}
					# Multiple exon CNV
					else {
						$score = scoreMultipleExon($svtype, $nCalls, $ratioStd, $pctCalls,
              $corr, $ratio, $s2n, $s2nc, $cn, $regions);
					}
				}
				print OUT "$line;CONF_SCORE=$score\n";
			}
			close IN;
			close OUT;
			unlink $finalCalls;
			rename $scoredCalls, $finalCalls;
		# }
	}
}

####################################################################
sub scoreSingleExon {

	my $svtype   = shift;
	my $nCalls   = shift;
	my $sampleStd= shift;
	my $pctCalls = shift;
	my $corr     = shift;
	my $ratio    = shift;
	my $s2n      = shift;
	my $s2nc     = shift;
	my $cn       = shift;
  my $nRois    = shift;

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
  return $score;
}
####################################################################
sub scoreMultipleExon {

	my $svtype   = shift;
	my $nCalls   = shift;
	my $sampleStd= shift;
	my $pctCalls = shift;
	my $corr     = shift;
	my $ratio    = shift;
	my $s2n      = shift;
	my $s2nc     = shift;
	my $cn       = shift;
  my $nRois   = shift;

	# Sample Level Metrics
	my $stdMetric      = $sampleStd < 0.2 ? 1 : 0;
	my $pctCallsMetric = $pctCalls < 1.5 ? 1 : 0;
	my $corrMetric     = $corr > 0.97 ? 1 : 0;
	my $SLM = 0.45*$stdMetric + 0.45*$pctCallsMetric + 0.1*$corrMetric;

	# Roi LeveL Metrics
  # Difference between observed ratio and expected ratio
	my $expectedRatio = 0.5*$cn;
	my $absDiffRatio  = 0;
	if ($ratio > $expectedRatio) {
		$absDiffRatio = 2.5*(abs($ratio-$expectedRatio));
	}
	else {
		$absDiffRatio =2.5*(abs($expectedRatio-$ratio));
	}
	my $absDiffMetric = 1-$absDiffRatio > 0 ? 1-$absDiffRatio : 0;

  # Signal to noise
	my $signal2noiseMetric = $s2n > 10 ? 1 : 0; # For the analyzed sample
	my $signal2noiseControlsMetric; # For the controls
	if ($s2nc eq '.'){
		$signal2noiseControlsMetric = 0;
	}
	else{
		$signal2noiseControlsMetric = $s2nc > 10 ? 1 : 0;
	}

  # Nbmber of targets
  my $nRoisMetric   = $nRois > 1 ? 1 : 0;
	my $nRoiWeight   = 0.2+($nRois*(0.1));
	my $signalWeight = 1-$nRoiWeight;

	my $RLM = ($nRoiWeight*$nRoisMetric) + ($signalWeight*0.70*$absDiffMetric)
    + ($signalWeight*0.15*$signal2noiseMetric) + ($signalWeight*0.15*$signal2noiseControlsMetric);

	$RLM = 1 if $RLM > 1;

	# Score = SLM + RLM
	my $score = 0.4*$SLM + 0.6*$RLM;

  return $score;

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
  #chr6	118880084	118880243	NM_002667_1_2;PLN	39.375000	100	1.314	8.78	chr6	118880084	118880243	2	1.314	159

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
			if ( $::sampleHash{$sampName}{ONOFF_SD_RATIO} > 0.2) {
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
			my $cmd = "$::zcat $ratios | $::awk '{ if ((\$5 < $::upperDelCutoff && \$5 > $::lowerDelCutoff) " .
      "|| (\$5 > $::lowerDupCutoff) || (\$5 > 0 && \$5 < 0.12)) { print \$0 }}' " .
      "| $::sort -V >> $singleExonCnv";
			print "$cmd\n" if $::verbose;
			system($cmd);
		}
		# but if segmented multi-exon CNVs
		else {
			# We might have also single ROIs that pass the minimum calling thresholds
			$cmd = "$::zcat $ratios | $::awk '{ if ((\$7 < $::upperDelCutoff && \$7 > $::lowerDelCutoff)";
			$cmd.= " || (\$7 > $::lowerDupCutoff) || (\$7 > 0 && \$7 < 0.12)) { print \$0 }}' |";
			$cmd.= "$::bedtools intersect -a stdin -b $intersect -v |" if -s $intersect;
			$cmd.= "$::sort -V >> $singleTmp";
			print "$cmd\n" if $::verbose;
			system($cmd);

			# Get ROI level information for the segmented CNV
			my $str = `$::bedtools intersect -a $ratios -b $intersect -wo`;
			my @tmp = split (/\n/, $str);
			my %seen = ();

			foreach my $line (@tmp) {

        #chr6	118880084	118880243	NM_002667_1_2;PLN	39.375000	100	1.314	8.78	chr6	118880084	118880243	2	1.314	159
				my @tmpLine = split (/\t/, $line);

				# These are the coordinates of the CNV
				my $segment = join ("\t", @tmpLine[8..10]);

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
				my $cnvType  = $tmpLine[12] > $::upperDelCutoff ? 'DUP' : 'DEL';

				# Assigning integer copy number
				my $copyNumber = int (($tmpLine[12]*2)+0.5);

				my @ArrGC = ();
				my @ArrMAP= ();
				my @signal2noise = ();
				my @ArrMeanZscore =  ();
				my @ArrRatiosSample = ();
				my @ArrRatiosControl = ();
				my @ArrRatiosOfftarget = ();

				# Then, for each ROI from the segmented CNV
				foreach my $target ( @targets ){

          #chr6	118880084	118880243	NM_002667_1_2;PLN	39.375000	100	1.314	8.78	chr6	118880084	118880243	2	1.314	159
					my @tmp = split (/\t/, $target);
					my $coordinate  = join ("\t", @tmp[0..2]);
					my $coordPlusRoi= join ("\t", @tmp[0..3]);
					my $ratioSample = $tmp[6];

					my $str = $type eq 'on-target' ? `$::grep -P \'^$coordinate\' $ratioFile` :
            `$::zgrep -e \'$coordinate\' $ratioFile`;
          chomp $str;
					my @tmpStr = split (/\t/, $str);

					# Get mappability
          my $gc  = $tmp[4];
          push @ArrGC, $gc;

          # Get GC content
          my $map = $tmp[5];
					push @ArrMAP, $map;

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
          }
          elsif ($type eq 'mixed') {
            $meanControlRatio =	$::sampleHash{$sampName}{ONOFF_MEAN_RATIO};
            $sdSample         = $::sampleHash{$sampName}{ONOFF_SD_RATIO};
          }
          $sdSample = 0.01 if $sdSample == 0;
          $zscoreSample     = sprintf "%.3f", (($ratioSample - $meanControlRatio)/$sdSample);
          push @ArrMeanZscore, $zscoreSample;
          push @ArrRatiosSample, $ratioSample;

          if ($type eq 'off-target' or $type eq 'mixed') {
            push @ArrRatiosOfftarget, $tmp[6]
          }
          else {
            push @signal2noise, $tmp[7];
          }
				}

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
					sprintf "%.3f", (($tmpLine[12]-$::sampleHash{$sampName}{ONTARGET_MEAN_RATIO})/$::sampleHash{$sampName}{ONTARGET_SD_RATIO});

				my ($cipos, $ciend) = getCiposCiend($chr, $start, $End, \@::ROIarray);
				my $size = $End-$start;

				next if abs($meanZscore) < $::minZscore;
				next if $meanSignal2Noise < 5;

        my %InfoFields = (
          CIPOS => $cipos,
          CIEND => $ciend,
          SVTYPE=> $cnvType,
          SVLEN => $size,
          GC    =>$meanGC,
          MAP   =>$meanMap,
          GENE  =>$affectedROIs,
          REGIONS => $nexons,
          RRD => $tmpLine[12],
          MADRD => $MAD,
          CN => $copyNumber,
          SNR => $meanSignal2Noise,
          SNRC => $signal2noiseControls,
          ZSCORE => $meanZscore
        );

        my @vcfData = map { $_ . '=' . $InfoFields{$_} } keys %InfoFields;
        print OUTFILE "$chr\t$start\t$End\tIMPRECISE;" . join(";", @vcfData) . "\n";
			}

			$sampName = basename($sampName);
			$sampName =~s/.on_off//;
			$sampName =~s/;//;

			if ( $type eq 'mixed') {
				$breakFile = "$inputDir/ON_TARGET/$sampName/$sampName.breakpoints.bed";

				$cmd = "$::cat $::HoF{GLOBAL_OFFTARGET_BED} $::bed | $::sort -V > $::outDir/$sampName.mixed.bed";
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
  #exit;
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
	  print " ERROR: Missing $ratiosFile. Skipping single-exon CNV analysis\n";
	  return 0;
  }

  Utils::decompressFile($ratiosFile) if -e $ratiosFile;

  my $masterSingleCnv = $outputDir . "/" . "master_single_calls.bed";

  open (MASTER_FILE, ">", $masterSingleCnv) || die " ERROR: Unable to open $masterSingleCnv\n";
  foreach my $file (@singleCallFiles) {
  	my $sampleName = basename($file);
  	$sampleName =~s/.raw.single.exon.cnv.bed//;
  	open (IN,"<", $file) || die " ERROR: Unable to open $file\n";
  	while (my $line=<IN>){
  		chomp $line;
  		print MASTER_FILE "$line\t$sampleName\n";
  	}
  	close IN;
  	unlink($file);
  }
  close MASTER_FILE;

  print " INFO: Calling single-exon CNVs\n";

  my $cmd = "$::Rscript $::analyzeSingleExonR " .
  " --raw_calls $masterSingleCnv " .
  " --coverage_data $ungzRatioFile " .
  " --references $::HoF{REFERENCES}" .
  " --min_zscore $::minZscore" .
  " --min_s2n 5".
  " --lower_del_cutoff $::lowerDelCutoff " .
  " --upper_del_cutoff $::upperDelCutoff " .
  " --lower_dup_cutoff $::lowerDupCutoff " .
  " --output_dir $outputDir ";
  print " INFO: $cmd\n";
  system($cmd);

  my @tmpRdata = glob("$outputDir/*single.exon.cnv.bed");
  if (!@tmpRdata) {
	  print " ERROR: Missing single exon tmp data\n";
  }

  foreach my $sample( sort keys %::sampleHash ) {

  	my $tmpSinglecalls = $outputDir . "/" . "$sample.single.exon.cnv.bed";
    next if !-e $tmpSinglecalls;

  	my $toMergeSegs          = "$outputDir/$sample.tomerge.single.exon.cnv.bed";
  	my $mergedAdjacentRois   = "$outputDir/$sample.merged.multiple.exon.cnv.bed";
  	my $mergedSingleSegments = "$outputDir/$sample.merged.single.multiple.exon.cnv.bed";
  	my $mergedCnvBreakpoints = "$outputDir/$sample.merged.breaks.single.multiple.exon.cnv.bed";
  	my $breakpointCalls      = "$outputDir/$sample/$sample.breakpoints.bed";

  	if (-e $mergedAdjacentRois ) {
  		$cmd = "$::cat $mergedAdjacentRois > $toMergeSegs";
  		system $cmd;
  	}

    open (IN, "<", $tmpSinglecalls) || die " ERROR: Unable to open $tmpSinglecalls\n";
  	open (OUT, ">>", $toMergeSegs)  || die " ERROR: Unable to open $toMergeSegs\n";
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
      my $copyNumber      = int (($testSampleRatio)*2+0.5);

      my %InfoFields = (
        CIPOS  => $cipos,
        CIEND  => $ciend,
        SVTYPE => $cnv_type,
        SVLEN  => $size,
        GC     => $GC,
        MAP    => $MAP,
        GENE   => $genename,
        REGIONS=> 1,
        RRD    => $testSampleRatio,
        MADRD  => ".",
        CN     => $copyNumber,
        SNR    => $s2nTest,
        SNRC   => $s2nControls,
        ZSCORE => $zscore
      );

      my @vcfData = map { $_ . '=' . $InfoFields{$_} } keys %InfoFields;
      print OUT "$chr\t$start\t$end\tIMPRECISE;" . join(";", @vcfData) . "\n";
  	}
  	close IN;
  	close OUT;
  	unlink($tmpSinglecalls);

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

1;
