#!/usr/bin/env perl

package Offtarget;
use strict;
use Getopt::Long;
use File::Basename;
use List::Util qw(min max);
use Parallel::ForkManager;
use Sort::Key::Natural qw(natsort);

##########################
sub macs2PeakCalling {
	my $inputBam   = shift;
	my $sampleName = shift;
	my $outputDir  = shift;

	my $orfanPeaks = $outputDir . "/" . "$sampleName.orfanPeaks";

  my $cmd = "$::macs2 callpeak -t $inputBam -n $orfanPeaks $::devNull";
	system $cmd;
}

##########################
sub appendOfftargetMetrics2Summary {

	my $offtargetTmp = shift;

	open (IN, "<", $offtargetTmp) || die " ERROR: Unable to open $offtargetTmp\n";

	while (my $line=<IN>) {
		chomp $line;
		my ($sample, $counts, $window) = split (/\t/, $line);
		$::sampleHash{$sample}{READSOFFTARGET} = $counts;
		$::sampleHash{$sample}{OFFTARGETBIN} = $window;
		if ($counts > $::minOfftargetReads) {
			$::sampleHash{$sample}{PERFORM_OFFTARGET} = 'yes';
		}
		else {
			$::sampleHash{$sample}{PERFORM_OFFTARGET} = 'no';
		}
	}
	close IN;

	my $head = `$::head -1 $::outDir/summary_metrics.log`; chomp $head;
	if ($head =~/OFFTARGET_READS/) {
		return 0;
	}

	if (-e "$::outDir/summary_metrics.log") {

    ( my $tmpSummary = "$::outDir/summary_metrics.log" ) =~ s/\.log/\.tmp.log/;

		open ( IN, "<", "$::outDir/summary_metrics.log") || die " ERROR: Unable to open $::outDir/summary_metrics.log\n";
		open ( TMP, ">", $tmpSummary ) || die " ERROR: Unable to open $tmpSummary\n";
		my $nLine = 0;
		while (my $line =<IN> ){

			chomp $line;
			my @tmp = split (/\t/, $line);
			my $sample = $tmp[0];
			if ($nLine == 0) {
				print TMP "$line\tOFFTARGET_READS\tOFFTARGET_BIN_SIZE\n";
			}
			else {
				my $counts = $::sampleHash{$sample}{READSOFFTARGET};
				print TMP "$line\t$::sampleHash{$sample}{READSOFFTARGET}\t$::sampleHash{$sample}{OFFTARGETBIN}\n";
			}
			$nLine++;
		}
		close IN;
		close TMP;

		unlink ("$::outDir/summary_metrics.log");
		rename $tmpSummary, "$::outDir/summary_metrics.log";
	}
}

##########################
sub preProcess {

 my $bed         = shift;
 my $outDir      = shift;
 my $offtargetDir= shift;

 # Creating file with regions to be excluded (padding coding exons to get rid of residual coverage from these regions & centromere concatenation)
 my $cmd = "$::awk '{if (\$2<\$3) {print \$1\"\t\"\$2-1000\"\t\"\$3+1000\"\t\"\$4}".
 " else if (\$2>\$3) {print \$1\"\t\"\$2+1000\"\t\"\$3-1000\"\t\"\$4}}' $bed > $outDir/OFF_TARGET/ontarget_paded_400bp.bed";
 system($cmd);

 $cmd = "$::cat $outDir/OFF_TARGET/ontarget_paded_400bp.bed $::centromeres $::genomePatches ";
 if (!$::hasChr) {
	 $cmd.= " | perl -e -p \"s/chr//\" ";
 }
 $cmd.= " | $::awk '\$1 !~ /_/ {print \$0}' | $::sort -V -u > $outDir/OFF_TARGET/ontarget_paded_400bp_centromers_patches.bed";
 system($cmd);

 # Complementary file generation for the peak detection step
 if (!$::hasChr) {
	$cmd = "cat $::chrFile | perl -e -p \"s/chr//\" | $::bedtools complement -i $outDir/OFF_TARGET/ontarget_paded_400bp_centromers_patches.bed -g stdin ";
 	$cmd.= " | $::awk '{ print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$3-\$2}' "
	. "| $::awk '{if(\$4>=1) print \$1\"\t\"\$2\"\t\"\$3\"\tComplement\"}' > $outDir/OFF_TARGET/offtarget_unprocessed.bed";
 }
 else {
 	$cmd = "$::bedtools complement -i $outDir/OFF_TARGET/ontarget_paded_400bp_centromers_patches.bed -g $::chrFile ";
	$cmd.= " | $::awk '{ print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$3-\$2}' | $::awk '{if(\$4>=1) print \$1\"\t\"\$2\"\t\"\$3\"\tComplement\"}' > $outDir/OFF_TARGET/offtarget_unprocessed.bed";
 }
 system($cmd);

}

##########################
sub removePeaks {

 my $bed         = shift;
 my $outDir      = shift;
 my $offtargetDir= shift;

 if (-e "$offtargetDir/$::outName.MergedNarrowPeaks.bed" ) {
	 return 1;
 }

 foreach my $bam ( @::bams ) {

	my $pid = $::pm -> start() and next;
	my $sample         = basename($bam);
	my $offTargetbam   = $sample;
	my $offTargetPeaks = $sample;

	$sample =~s/.bam//;
	$offTargetbam=~s/.bam//;
	$offTargetPeaks=~s/.bam//;

	if (!-e "$offtargetDir/$sample.narrowPeak.format.bed") {

		# Get an sliced BAM from off-target regions
		print " INFO: Removing off-target peaks on $sample\n";
		my $offtargetBam = $offtargetDir . "/" . "$sample.offTarget.bam";
	 	my $cmd ="$::samtools view -b -L $offtargetDir/offtarget_unprocessed.bed $bam > $offtargetBam";
	 	print "$cmd\n" if $::verbose;
		system($cmd) if !-e $offtargetBam;

		# Index sliced BAM file
		Utils::indexBam( $offtargetBam );

		# Identify off-target peaks using MACS2
		macs2PeakCalling($offtargetBam, $sample, $offtargetDir);

		my $narrowPeaks    = $offtargetDir . "/" . "$sample.orfanPeaks_peaks.narrowPeak";
		my $narrowPeaksFmt = $offtargetDir . "/" . "$sample.narrowPeak.format.bed";

		# Formatting narrowPeak files
		if (-e $narrowPeaks) {
			$cmd = "$::cat $narrowPeaks | $::awk '{print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$4\";\"\$5}' > $narrowPeaksFmt";
			system($cmd);
		}
		else {
			print " WARNING: $sample does not contain off-target peaks\n";
			$::sampleHash{$sample}{READSOFFTARGET} = 0;
		}
	}
	$::pm->finish;
 }
 $::pm->wait_all_children;
 unlink (glob ("$offtargetDir/*summits.bed"));
 unlink (glob ("$offtargetDir/*peaks.xls"));
 unlink (glob ("$offtargetDir/*Peaks_model.r"));
}

##########################
sub createBins {

 my $bed         = shift;
 my $outDir      = shift;
 my $offtargetDir= shift;

 # Creating merged & padded offtarget peak file for all the samples
 my @narrowPeaksFiles = glob ("$offtargetDir/*narrowPeak.format.bed");

 if (!@narrowPeaksFiles) {
	 print "mierda\n";
	foreach my $sample ( sort keys %::sampleHash) {
		$::sampleHash{$sample}{READSOFFTARGET} = 0;
		$::sampleHash{$sample}{PERFORM_OFFTARGET} = 'no';
	}
	return 1;
 }
 my $mergedNarrowPeaksTmp = $offtargetDir . "/" . "$::outName.MergedNarrowPeaks.tmp.bed";

 # Adding 400bp padding coordinates at each peak identified
 my $cmd = "$::cat @narrowPeaksFiles | $::sort -V -u " .
 " | $::awk '{if (\$2<\$3) {print \$1\"\t\"\$2-400\"\t\"\$3+400\"\t\"\"mergedPaddedNarrowPeaks\"}" .
 " else if (\$2>\$3) {print \$1\"\t\"\$2+400\"\t\"\$3-400\"\t\"\$4\"\t\"\"mergedPaddedNarrowPeaks\"}}' - ";
 $cmd.= " | $::grep -v '\_'> $mergedNarrowPeaksTmp";
 system($cmd);

 my $mergedNarrowPeak = $offtargetDir . "/" . "$::outName.MergedNarrowPeaks.bed";
 # Deleting inconsistent coordinate entries when poiting out of the chromosome length
 open (IN, "<", $mergedNarrowPeaksTmp) || die " ERROR: Unable to open $mergedNarrowPeaksTmp\n";
 open (OUT, ">", $mergedNarrowPeak) || die " ERROR: Unable to open $mergedNarrowPeak\n";
 while (my $line=<IN>) {
	chomp $line;
	my @tmp = split (/\t/, $line);
	my $chr   = $tmp[0];
	my $start = $tmp[1];
	my $end   = $tmp[2];
	if ($end > $::ChromosomeLengths{$chr}) {
		$end = $::ChromosomeLengths{$chr};
	}
	print OUT "$chr\t$start\t$end\t$tmp[3]\n";
 }
 close IN;
 close OUT;
 unlink $mergedNarrowPeaksTmp;

 # Merging patched regions and exclusion centromere regions into a single file
 $cmd  = "$::cat $mergedNarrowPeak $::genomePatches $::centromeres";
 if (!$::hasChr)  {
	 $cmd.= " | perl -e -p \"s/chr//\" ";
 }
 $cmd .= " | $::sort -V -u | $::awk '\$1 !~ /_/ {print \$0}' - > $offtargetDir/offtarget_centromeres_patches_peaks.bed";
 print "$cmd\n" if $::verbose;
 system($cmd);

 # Obtaining the complement of the previous exclusion file
 if (!$::hasChr)  {
	 $cmd = "$::cat $::chrFile | perl -e -p \"s/chr//\" ";
	 $cmd.= "$::bedtools complement -i $offtargetDir/offtarget_centromeres_patches_peaks.bed -g stdin";
 }
 else {
	 $cmd = "$::bedtools complement -i $offtargetDir/offtarget_centromeres_patches_peaks.bed -g $::chrFile";
 }
 $cmd.= " | $::awk '{ print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$3-\$2}'";
 $cmd.= " | $::awk '{if(\$4>=1) print \$1\"\t\"\$2\"\t\"\$3}' - > $offtargetDir/offtarget_centromeres_patches_peaks_negative.bed";
 print "$cmd\n" if $::verbose;
 system($cmd);

 my $offtargetInfoTmp = $::outDir . "/" . "offtarget_info.tmp.txt";
 open (TMP, ">", $offtargetInfoTmp) || die " ERROR: Unable to open $offtargetInfoTmp\n";
 foreach my $sample ( sort keys %::sampleHash) {

    my $bam = $offtargetDir . "/" . "$sample.offTarget.bam";
		if (-e $bam ) {
			my $counts = Utils::getTotalCountsBam($bam,
				"$offtargetDir/offtarget_centromeres_patches_peaks_negative.bed" );

			# Defnining offtarget window size based on this empirically derived rule.
			$::sampleHash{$sample}{OFFTARGETBIN}
				= $counts > 0 ? int ( (4000000 *150000)/$counts ) : 10e6;

			if ($::sampleHash{$sample}{OFFTARGETBIN} < 50000) {
				$::sampleHash{$sample}{OFFTARGETBIN} = 50000;
			}
			$::sampleHash{$sample}{READSOFFTARGET} = $counts;
			print TMP "$sample\t$::sampleHash{$sample}{READSOFFTARGET}\t$::sampleHash{$sample}{OFFTARGETBIN}\n";

			my $countsHuman = Utils::number2human($counts);
			my $binSizeHuman= Utils::number2human($::sampleHash{$sample}{OFFTARGETBIN});
			print " INFO: $sample off-target reads: $countsHuman\toff-target bin size: $binSizeHuman\n";

			# my $sample_offtarget_bed
			# 	= createPseudowindows($offtargetDir, $sample, $::sampleHash{$sample}{OFFTARGETBIN});
		}
		else {
			print " INFO: Skipping off-target analysis for sample $sample\n";
			$::sampleHash{$sample}{READSOFFTARGET}    = 0;
			$::sampleHash{$sample}{PERFORM_OFFTARGET} = 'no';
		}
 }

 # Calculate median bin size from all samples
 my @binArray = ();
 for my $sample (sort keys %::sampleHash) {
	 if ($::sampleHash{$sample}{OFFTARGETBIN}) {
		 push @binArray, $::sampleHash{$sample}{OFFTARGETBIN};
	 }
 }
 my $medianBinSize = Utils::medianArray(@binArray);
 if ($medianBinSize) {
	 print " INFO: Median off-target bin size: " . Utils::number2human($medianBinSize) . "\n";
 }

 # Create a global offtarget bed file
 if ($medianBinSize) {
	 	$::HoF{GLOBAL_OFFTARGET_BED} = createPseudowindows($offtargetDir, "global", $medianBinSize)
 }

 #$::pm->wait_all_children;
 close TMP;

 # Append sample offtarget info to main summary file
 appendOfftargetMetrics2Summary("$::outDir/offtarget_info.tmp.txt");

 unlink($offtargetInfoTmp);
 unlink("$offtargetDir/ontarget_paded_400bp.bed");
 unlink("$offtargetDir/offtarget_unprocessed.bed");
 #unlink("$outDir/offTarget.pseudowindowed.bed");
 #unlink("$outDir/offtarget_temp_with_centromeres.bed");
 #unlink("$outDir/offtarget_centromeres_patches_peaks.bed");
 unlink("$offtargetDir/ontarget_paded_400bp_centromers_patches.bed");
 #unlink("$outDir/offtarget_centromeres_patches_peaks_negative.bed");
}

#########################
sub estimateWindowSize {

	my $outDir  = shift;
	my $sample  = shift;
	my $binSize = shift;
	my $bam     = shift;

	my $cmd = "$::bedtools bamtobed -i $bam | $:cut -f 1,2,3 | $::sort -V > $outDir/$sample.offtarget.coordinates.bed";
	system $cmd if !-e "$outDir/$sample.offtarget.coordinates.bed";

	my $off_bed = createPseudowindows($outDir, $sample, $binSize);
	my $ini_std = caculateRatioStd($off_bed, "$outDir/$sample.offtarget.coordinates.bed");

	if ($ini_std > 0.2) {
		my $maxBinSize = $binSize + (50000*5);
		my $count = 0;
		for(my $bin=$binSize; $maxBinSize; $bin+=50000) {
			my $off_bed = createPseudowindows($outDir, $sample, $bin);
			my $bin_std = caculateRatioStd($off_bed, "$outDir/$sample.offtarget.coordinates.bed");
			$count++;
			#print "$sample\t$bin\t$bin_std\n";
			last if $count == 10;

		}
	}

	# Now we should get counts and produce bin-level ratios
}

#########################
sub caculateRatioStd {

	my $off_bed   = shift;
	my $coord_bed = shift;

	my %joinedWindows = ();
	my %seen = ();
	my $min = 10e20;
	my $max = 0;
	open (IN, "$::bedtools intersect -a $off_bed -b $coord_bed -c -sorted |");
	while(my $line=<IN>) {
		chomp $line;
		my @tmp = split("\t", $line);
		my @info = split (";", $tmp[3]);
		my $window_name = $info[0];

		if (!exists $seen{$window_name} ) {
			$min = 10e20;
			$max = 0;
		}
		$seen{$window_name}++;
		$joinedWindows{$window_name}{CHR} = $tmp[0];
		if ($tmp[1] < $min) {
			$joinedWindows{$window_name}{START} = $tmp[1];
			$min = $tmp[1];
		}
		if ($tmp[2] > $max) {
			$joinedWindows{$window_name}{END} = $tmp[2];
		}
		my $length = $tmp[2]-$tmp[1];
		$joinedWindows{$window_name}{COUNTS}+=$tmp[4];
	}
	close IN;

	#my $offCounts = "$outDir/$sample.offtarget.counts.bed";
	my @arrCounts = ();
	#open (OUT, ">", $offCounts) || die " ERROR: unable to open $offCounts\n";
	foreach my $region (natsort keys %joinedWindows) {
		#print OUT "$joinedWindows{$region}{CHR}\t$joinedWindows{$region}{START}\t$joinedWindows{$region}{END}\t$region\t$joinedWindows{$region}{COUNTS}\n";
		push @arrCounts,$joinedWindows{$region}{COUNTS};
	}
	my $medianCounts = Utils::medianArray(@arrCounts);
	#close OUT;
	my @arrRatios = ();
	foreach my $region (natsort keys %joinedWindows) {
		my $ratio = $joinedWindows{$region}{COUNTS}/$medianCounts;
		push @arrRatios,$ratio;
	}
	my $sdRatios = Utils::MAD(@arrRatios);
	return $sdRatios;
}

##########################
sub createPseudowindows {

 #Goal here is to create a sample-specific off-target BED regions file
 #Author: Jesús Matés Ramírez, PhD.
 my $outDir  = shift;
 my $sample  = shift;
 my $binsize = shift;

 my $off_bed = "$outDir/$sample.offtarget.bed";
 if (-e $off_bed) {
 	#return 1;
 }

 my $inputpw  = "$outDir/offtarget_centromeres_patches_peaks_negative.bed";
 my $outputpw = "$outDir/$sample.offTarget.pseudowindowed.bed";

 open (INPW, "<", $inputpw) or die "cannot open $inputpw";
 open (OUTPW, ">", $outputpw) or die "cannot open $outputpw";

 my $sum     =  0;
 my $presum  =  0;
 my $lc      =  0;
 my $winc    =  1;
 my $prechr  = '';
 my $cntpiece=  1;
 my $newend  = '';
 my $newst   = '';
 my $prend   = '';
 my $prest   = '';

 while(my $line=<INPW>) {
	chomp $line;
	$lc++;
	my @tmp = split /\t/,$line;
	my $chr = $tmp[0];
	my $st  = $tmp[1];
	my $end = $tmp[2];

	next if(($st eq '0' && $end eq '0'));

	if($lc==1){ $prechr=$chr; }

	my $len=$end-$st;
	$sum+=$len;

	if($chr ne $prechr){
		$prechr=$chr;
		$winc++;
		$cntpiece=1;
		$sum=$len;
		$presum=0;
	}

	if(($presum==$sum)&&($sum==$binsize)){
		$prend=$st;
		$presum=0;
		next;
	}

	if($sum>$binsize){
		my $diff=$binsize-$presum; #if($st==62158568){ print "$sum>$binsize / $diff=$binsize-$presum \n"; }
		my $iters=int(($sum/$binsize) + 0.5);
		for(my $i=0; $i<$iters; $i++){
			if($i==0){
				$newst=$st;
				$newend=$newst+$diff;
				$sum=$binsize;

				next if $newst >= $::ChromosomeLengths{$chr};

				if ($newend >= $::ChromosomeLengths{$chr}) {
					$newend = $::ChromosomeLengths{$chr};
				}

				print OUTPW "$chr\t$newst\t$newend\tpdwindow_$winc;piece_$cntpiece\n";

				$cntpiece=1;
				$winc++;
				$prend=$newend;
				$prest=$newst;
				$presum=$sum;

				my $nlen=$end-$prend;
				if($nlen<=$binsize){
					$sum=$nlen;
					$newst=$prend+1;
					$newend=$end;

					if ($newend > $::ChromosomeLengths{$chr}) {
						$newend = $::ChromosomeLengths{$chr};
					}

					print OUTPW "$chr\t$newst\t$newend\tpdwindow_$winc;piece_$cntpiece\n";

					$cntpiece++;
					$presum=$sum;
					last;
				}
				else{
					$sum=$binsize;
					$newst=$prend+1;
					$newend=$newst+$binsize;

					if ($newend > $::ChromosomeLengths{$chr}) {
						$newend = $::ChromosomeLengths{$chr};
					}

					print OUTPW "$chr\t$newst\t$newend\tpdwindow_$winc;piece_$cntpiece\n";

					$cntpiece=1;
					$winc++;
					$prend=$newend;
					$prest=$newst;
					next;
				}
			}
			else{
				my $nlen=$end-$prend;
				if($nlen<=$binsize){
					$sum=$nlen;
					$newst=$prend+1;
					$newend=$end;
					next if $newst >= $::ChromosomeLengths{$chr};
					if ($newend > $::ChromosomeLengths{$chr}) {
						$newend = $::ChromosomeLengths{$chr};
					}

					print OUTPW "$chr\t$newst\t$newend\tpdwindow_$winc;piece_$cntpiece\n";

					$cntpiece++;
					$prend=$newend;
					$prest=$newst;
					$presum=$sum;
					last;
				}
				else{
					$newst=$prend+1;
					$newend=$newst+$binsize;
					$sum=$binsize;
					next if $newst >= $::ChromosomeLengths{$chr};
					if ($newend > $::ChromosomeLengths{$chr}) {
						$newend = $::ChromosomeLengths{$chr};
					}

					print OUTPW "$chr\t$newst\t$newend\tpdwindow_$winc;piece_$cntpiece\n";

					$cntpiece=1;
					$winc++;
					$prend=$newend;
					$prest=$newst;
					next;
				}
			}
		}
		next;
	}
	else{
		if ($end > $::ChromosomeLengths{$chr}) {
			$end = $::ChromosomeLengths{$chr};
		}
		print OUTPW "$chr\t$st\t$end\tpdwindow_$winc;piece_$cntpiece\n";

		if($sum==$binsize){
			$cntpiece=1;
			$winc++;
			$prend=$end;
			$prest=$st;
			$presum=$sum;
			next;
		}
		else{
			$cntpiece++;
			$prend=$end;
			$prest=$st;
			$presum=$sum;
			next;
		}
	}
	$presum=$sum;
	if($presum==$binsize){
		$presum=0;
	}
 }
 close INPW;
 close OUTPW;

 my $cmd = "$::awk '{if( \$3>\$2) print \$0}' $outDir/$sample.offTarget.pseudowindowed.bed | ";
 $cmd .= "$::awk '{if (\$1 !~/_/) { print \$0 }}' | $::awk '{if(\$3!=\$2) print \$0}' | ";
 $cmd .= "$::bedtools intersect -a stdin -b $::centromeres $::genomePatches -v | ";
 $cmd .= "$::sort -V >  $outDir/$sample.offtarget_temp_with_centromeres.bed";
 system($cmd);

 $cmd = "$::awk '{if( \$4 ~ /centromere/) { print \$1\"\t\"\$2-1000000\"\t\"\$3+1000000\"\t\"\$4 }}' $::centromeres |";
 $cmd.="$::bedtools intersect -a $outDir/$sample.offtarget_temp_with_centromeres.bed -b stdin -v | $::sort -V > $off_bed";
 system $cmd;

 unlink("$outDir/$sample.offTarget.pseudowindowed.bed");
 unlink("$outDir/$sample.offtarget_temp_with_centromeres.bed");
 return $off_bed;
}

1;
