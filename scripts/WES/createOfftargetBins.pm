#!/usr/bin/env perl

package createOfftargetBins;
use strict;
use Getopt::Long;
use File::Basename; 
use List::Util qw(min max);
use Parallel::ForkManager;

##########################

sub appendOfftargetMetrics2Summary {

	my $offtargetTmp = shift;
	open (IN, "<", $offtargetTmp) || die " ERROR: Unable to open $offtargetTmp\n";
	while (my $line=<IN>) {
		chomp $line;
		my ($sample, $counts, $window) = split (/\t/, $line);
		$::sampleHash{$sample}{READSOFFTARGET} = $counts;
		$::sampleHash{$sample}{OFFTARGETBIN} = $window;
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
 my $cmd = "$::awk '{if (\$2<\$3) {print \$1\"\t\"\$2-1000\"\t\"\$3+1000\"\t\"\$4} else if (\$2>\$3) {print \$1\"\t\"\$2+1000\"\t\"\$3-1000\"\t\"\$4}}' $bed > $outDir/OFF_TARGET/ontarget_paded_400bp.bed";
 system($cmd);

 $cmd = "$::cat $outDir/OFF_TARGET/ontarget_paded_400bp.bed $::centromeres $::genomePatches | $::awk '\$1 !~ /_/ {print \$0}' | $::sort -V -u > $outDir/OFF_TARGET/ontarget_paded_400bp_centromers_patches.bed";
 system($cmd);

 # Complementary file generation for the peak detection step
 $cmd = "$::bedtools complement -i $outDir/OFF_TARGET/ontarget_paded_400bp_centromers_patches.bed -g $::chrFile | $::awk '{ print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$3-\$2}' | $::awk '{if(\$4>=1) print \$1\"\t\"\$2\"\t\"\$3\"\tComplement\"}' > $outDir/OFF_TARGET/offtarget_unprocessed.bed";
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

 # Offtarget peak detection
 my $count = 0; 

 foreach my $bam ( @::bams ) {

	my $pid = $::pm -> start() and next;
	$count++;

	my $sample         = basename($bam);
	my $offTargetbam   = $sample;
	my $offTargetPeaks = $sample;

	$sample =~s/.bam//;
	$offTargetbam=~s/.bam//;
	$offTargetPeaks=~s/.bam//;

	if (!-e "$offtargetDir/$sample.narrowPeak.format.bed") {

		print " INFO: Processing $sample bam file\n";
	 	my $cmd="$::samtools view -b -L $offtargetDir/offtarget_unprocessed.bed $bam > $offtargetDir/$sample.offTarget.bam";
		system($cmd);

	 	$cmd="$::samtools index $offtargetDir/$sample.offTarget.bam";
		system($cmd);

		print " INFO: Offtarget peak calling for $sample.offTarget.bam file\n";
	 	$cmd = "$::macs2 callpeak -t $offtargetDir/$sample.offTarget.bam -n $offtargetDir/$sample.orfanPeaks > /dev/null 2>&1";
		system($cmd);

		if (-e "$offtargetDir/$sample.orfanPeaks_peaks.narrowPeak") {
			$cmd = "$::cat $offtargetDir/$sample.orfanPeaks_peaks.narrowPeak | $::awk '{print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$4\";\"\$5}' > $offtargetDir/$sample.narrowPeak.format.bed";
			system($cmd);
		}
		else {
			print " WARNING: $sample does not contain off-target peaks\n";
		}
	}
	$::pm->finish;
 }
 $::pm->wait_all_children;
 unlink (glob ("$offtargetDir/*summits.bed"));
 unlink (glob ("$offtargetDir/*peaks.xls"));
 #unlink (glob ("$offtargetDir/*narrowPeak.format.bed"));
 unlink (glob ("$offtargetDir/*Peaks_model.r"));

}

##########################
sub createBins {

 my $bed         = shift; 
 my $outDir      = shift;
 my $offtargetDir= shift;

 # Creating merged & padded offtarget peak file for all the samples
 my @narrowPeaksFiles = glob ("$offtargetDir/*narrowPeak.format.bed");

 # Returning if no peaks have been detected
 if (!@narrowPeaksFiles) {
	 return 0;
 }

 # Adding some padding coordinates at each peak identified
 my $cmd = "$::cat $offtargetDir/*narrowPeak.format.bed | $::sort -V -u | $::awk '{if (\$2<\$3) {print \$1\"\t\"\$2-400\"\t\"\$3+400\"\t\"\"mergedPaddedNarrowPeaks\"} else if (\$2>\$3) {print \$1\"\t\"\$2+400\"\t\"\$3-400\"\t\"\$4\"\t\"\"mergedPaddedNarrowPeaks\"}}' - | $::grep -v '\_'> $offtargetDir/$::outName.MergedNarrowPeaks.tmp.bed";  
 system($cmd);

 # Deleting inconsistent coordinate entries when poiting out of the chromosome length
 open (IN, "<", "$offtargetDir/$::outName.MergedNarrowPeaks.tmp.bed") || die " ERROR: Unable to open $offtargetDir/$::outName.MergedNarrowPeaks.tmp.bed\n";
 open (OUT, ">", "$offtargetDir/$::outName.MergedNarrowPeaks.bed") || die " ERROR: Unable to open $offtargetDir/$::outName.MergedNarrowPeaks.bed\n";
 while (my $line=<IN>) {
	chomp $line;
	my @tmp = split (/\t/, $line);
	if ($tmp[2] > $::ChromosomeLengths{$tmp[0]}) {
		$tmp[2] = $::ChromosomeLengths{$tmp[0]};
	}
	print OUT "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\n";
 }
 close IN;
 close OUT;
 unlink ("$offtargetDir/$::outName.MergedNarrowPeaks.tmp.bed");

 # Merging patched regions and exclusion centromere regions into a single file
 $cmd = "$::cat $offtargetDir/$::outName.MergedNarrowPeaks.bed $::genomePatches $::centromeres | $::sort -V -u | $::awk '\$1 !~ /_/ {print \$0}' - > $offtargetDir/offtarget_centromeres_patches_peaks.bed";
 system($cmd);

 # Obtaining the complement of the previous exclusion file
 $cmd = "$::bedtools complement -i $offtargetDir/offtarget_centromeres_patches_peaks.bed -g $::chrFile | $::awk '{ print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$3-\$2}' | $::awk '{if(\$4>=1) print \$1\"\t\"\$2\"\t\"\$3}' - > $offtargetDir/offtarget_centromeres_patches_peaks_negative.bed"; 
 system($cmd);

 open (TMP, ">", "$::outDir/offtarget_info.tmp.txt") || die " ERROR: Unable to open $::outDir/offtarget_info.tmp.txt\n";
 my @offbams = glob("$offtargetDir/*bam");
 foreach my $bam ( @offbams ) {
	my $pid = $::pm -> start() and next;

	my $counts = `$::samtools view -c $bam -L $offtargetDir/offtarget_centromeres_patches_peaks_negative.bed`;
	chomp $counts;

	my $sample = basename($bam);
	$sample =~s/.offTarget.bam//;
	
	# Defnining offtarget window size based on this empirically derived rule. 
	# May be worth to add a more advanced method: https://academic.oup.com/bioinformatics/article/30/13/1823/2422194
	$::sampleHash{$sample}{OFFTARGETBIN}   = $counts > 0 ? int ( (4000000 *150000)/$counts ) : 10e6; 
	$::sampleHash{$sample}{READSOFFTARGET} = $counts;

	print TMP "$sample\t$::sampleHash{$sample}{OFFTARGETBIN}\t$::sampleHash{$sample}{READSOFFTARGET}\n";

	my $countsHuman = Utils::number2human($counts);
	my $binSizeHuman= Utils::number2human($::sampleHash{$sample}{OFFTARGETBIN});

	print " INFO: $sample off-target reads: $countsHuman\toff-target bin size: $binSizeHuman\n";

	my $sample_offtarget_bed = "$offtargetDir/$sample.offtarget.bed";
	
	createPseudowindows($outDir, $sample, $::sampleHash{$sample}{OFFTARGETBIN}, $sample_offtarget_bed);

	$::pm->finish;
 }
 $::pm->wait_all_children;
 close TMP;

 # Append sample offtarget info to main summary file
 appendOfftargetMetrics2Summary("$::outDir/offtarget_info.tmp.txt");

 unlink("$::outDir/offtarget_info.tmp.txt");
 unlink("$offtargetDir/ontarget_paded_400bp.bed");
 unlink("$offtargetDir/offtarget_unprocessed.bed");
 #unlink("$outDir/offTarget.pseudowindowed.bed");
 #unlink("$outDir/offtarget_temp_with_centromeres.bed");
 #unlink("$outDir/offtarget_centromeres_patches_peaks.bed");
 unlink("$offtargetDir/ontarget_paded_400bp_centromers_patches.bed");
 #unlink("$outDir/offtarget_centromeres_patches_peaks_negative.bed");
}

##########################
sub createPseudowindows {

 #Goal here is to create a sample-specific off-target BED regions file

 my $outDir  = shift;
 my $sample  = shift; 
 my $binsize = shift;
 my $off_bed = shift;

 if (-e $off_bed) {
 	return 1;
 } 

 my $inputpw  = "$outDir/OFF_TARGET/offtarget_centromeres_patches_peaks_negative.bed";
 my $outputpw = "$outDir/OFF_TARGET/$sample.offTarget.pseudowindowed.bed";

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

 my $cmd = "$::awk '{if( \$3>\$2) print \$0}' $outDir/OFF_TARGET/$sample.offTarget.pseudowindowed.bed | $::awk '{if (\$1 !~/_/) { print \$0 }}' | $::awk '{if(\$3!=\$2) print \$0}' | $::bedtools intersect -a stdin -b $::centromeres $::genomePatches -v | $::sort -V >  $outDir/OFF_TARGET/$sample.offtarget_temp_with_centromeres.bed";
 system($cmd); 
 
 $cmd = "$::awk '{if( \$4 ~ /centromere/) { print \$1\"\t\"\$2-1000000\"\t\"\$3+1000000\"\t\"\$4 }}' $::centromeres | $::bedtools intersect -a $outDir/OFF_TARGET/$sample.offtarget_temp_with_centromeres.bed -b stdin -v | $::sort -V > $off_bed";
 system $cmd;

 unlink("$outDir/OFF_TARGET/$sample.offTarget.pseudowindowed.bed");
 unlink("$outDir/OFF_TARGET/$sample.offtarget_temp_with_centromeres.bed");

}

1;
