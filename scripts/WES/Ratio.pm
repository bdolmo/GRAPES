#!/usr/bin/env perl


package Ratio;

use strict;
use Getopt::Long;
use File::Basename; 
use List::MoreUtils;
use Statistics::Descriptive;
use Parallel::ForkManager;
use List::MoreUtils qw(uniq);
 
 
 #####################################################################################
 #	  Calculate Copy Number Ratios for the merged file for ON and OFF target         #
 #####################################################################################
 sub calculateRatioOfftarget {

	my $offtargetDir = shift;

	foreach my $sample (sort keys %::sampleHash ) {

		# Skipping ratio calculation when not applicable
		next if exists $::referenceHash{$sample};
		next if $::sampleHash{$sample}{PERFORM_OFFTARGET} eq "no";
		next if !-e "$offtargetDir/$sample.normalized.bed.gz";

		open RATIOS, ">", "$offtargetDir/$sample.ratios.txt";
		open NORMALIZED, "$::zcat $offtargetDir/$sample.normalized.bed.gz |";

		my @ratios = ();
	 	while (my $line=<NORMALIZED>) {			
			 chomp $line;
			#chr1	578745	771659	pdwindow_4;piece_1	0.417414	26.1383443952743	275
			my @tmp = split (/\t/, $line);
			my $chr   = $tmp[0];
			my $start = $tmp[1];
			my $end   = $tmp[2];
			my $info  = $tmp[3];
			my $counts= $tmp[7];
			my $mappability = $tmp[5];
			my $ratio;
			my $coord = "$chr\t$start\t$end";

			next if $::filteredRois{$coord};

			$::ExonFeatures{$coord}{MAP} = $mappability;
			$::ExonFeatures{$coord}{GC} = int ($tmp[4]*100);

			# Skipping low mappability windows.
			next if $mappability < 90;

			# Skipping low mappability windows.
			if ( isChrX($chr) ) {
				$ratio = sprintf "%.3f", $counts/$::sampleHash{$sample}{OFFTARGET_NORMALIZED_X};
			}
			else {
				$ratio = sprintf "%.3f", $counts/$::sampleHash{$sample}{OFFTARGET_NORMALIZED};
			}

  			if ( $chr !~ "Y" ) {
				push @ratios, $ratio;
			}
			print RATIOS "$chr\t$start\t$end\t$info\t$ratio\t5\n";
		}
		close RATIOS;
		$::sampleHash{$sample}{OFFTARGET_MEAN_RATIO}= sprintf "%.3f", Utils::meanArray(@ratios);
		$::sampleHash{$sample}{OFFTARGET_SD_RATIO}  = sprintf "%.3f", Utils::std(@ratios);

		if ($::sampleHash{$sample}{OFFTARGET_SD_RATIO} < $::minOfftargetSD) {
			$::sampleHash{$sample}{PERFORM_OFFTARGET} = 'yes'
		}
		else {
			$::sampleHash{$sample}{PERFORM_OFFTARGET} = 'no';
		}

		# compress file
		Utils::compressFile("$offtargetDir/$sample.ratios.txt");

		print " INFO: $sample mean ratio: $::sampleHash{$sample}{OFFTARGET_MEAN_RATIO} sd: $::sampleHash{$sample}{OFFTARGET_SD_RATIO}\n";
	}
 }

 ##############################
 #	  Merging ON/OFF ratios   #
 ##############################
 sub mergeOnOffRatios {

	# Mergin On/Off target data for combined analysis 
	my $ontargetDir  = shift;
	my $offtargetDir = shift;

	foreach my $sample (sort keys %::sampleHash ) {

		my $ontargetRatios   = "$ontargetDir/$sample.ratios.txt.gz";
		my $offtargetRatios  = "$offtargetDir/$sample.ratios.txt.gz";
		my $OnOfftargetRatios= "$::outDir/$sample.ratios.txt.gz";

		# Concatenating on/off ratios into a unified file, but just adding offtarget data that shows less than 0.2 standard deviation
		if ( -e $ontargetRatios && -e $offtargetRatios ) {
			my $cmd = "$::zcat $ontargetRatios $offtargetRatios | $::sort -V | $::gzip > $OnOfftargetRatios\n";
			system $cmd;
		}
		elsif ( !-e $ontargetRatios && -e $offtargetRatios ){
			my $cmd = "$::zcat $offtargetRatios | $::sort -V | $::gzip > $OnOfftargetRatios\n";
			system $cmd;			
		}
		elsif ( -e $ontargetRatios && !-e $offtargetRatios ) {
			my $cmd = "$::zcat $ontargetRatios | $::sort -V | $::gzip > $OnOfftargetRatios\n";
			system $cmd;				
		}

		if (-s $OnOfftargetRatios) {
			my @ratios = ();
			open (IN, "$::zcat $OnOfftargetRatios |") || die " ERROR: Unable to open $OnOfftargetRatios\n";
			while (my $line =<IN>) {
				chomp $line;
				my @tmp = split (/\t/, $line);
				my $coord = "$tmp[0]\t$tmp[1]\t$tmp[2]";
				next if $::filteredRois{$coord};
				next if scalar @tmp < 5;
				my $ratio = $tmp[4];
				push @ratios, $ratio;
			}
			close IN; 
			$::sampleHash{$sample}{ONOFF_MEAN_RATIO}= Utils::meanArray(@ratios);
			$::sampleHash{$sample}{ONOFF_SD_RATIO}  = Utils::std(@ratios);
			@ratios = ();
		}
	}
 }

 ###################################################
 #	  Calculating On-target Copy Number Ratios     #
 ###################################################

 sub calculateOnTargetRatios {

   my $inputDir= shift;
   my $type    = shift;

   # Ratios will be calculated from MergedNormalizedCounts.bed file
   # that contains data from the pooled samples and the database refs.

   if (!-e $::HoF{MERGED_NORM_COV}) {
	   print " WARNING: skipping copy ratio calculation\n";
	   return 0;
   }

   # Select samples from header
   my $str = $::HoF{MERGED_NORM_COV} =~/.gz/ ? `$::zcat $::HoF{MERGED_NORM_COV} | $::head -1`
   : `$::cat $::HoF{MERGED_NORM_COV} | $::head -1`;
   chomp $str;

   my @tmpStr  = split (/\t/, $str);
   my @samples = @tmpStr[6..@tmpStr-1];

   # Select columns that correspond to the sample baseline
   foreach my $sample ( @samples ) {

	   	# Skipping non-analyzable samples
		next if exists $::referenceHash{$sample};
        if ($::doCaseControl) {
            next if $::sampleHash{$sample}{CONTROL};
        }
	    
		my @idx = ();
		foreach my $element ( @{$::sampleHash{$sample}{REFERENCE}} ) {

			next if $element eq $sample;
			if ($element eq 'none') {
				@{$::sampleHash{$sample}{REF_ID}} = ();
			}
			else {
				my $index = List::MoreUtils::first_index {$_ eq $element} @samples;
				push @{$::sampleHash{$sample}{REF_ID}}, $index;
			}
		}
   }

   open (TOPLOTRATIOS, ">", "$inputDir/toplotratios.bed") || die " ERROR: Cannot open $inputDir/toplotratios.bed\n";
   open (RATIOS, ">", $::HoF{RATIOS_ON}) || die " ERROR: Cannot open $::HoF{RATIOS_ON}\n";
   my @headerArr = ();

   open (IN, "<", $::HoF{MERGED_NORM_COV}) || die " ERROR: Cannot open $::HoF{MERGED_NORM_COV}\n";
   while (my $line=<IN>) {
	chomp $line;
	my @tmp = split (/\t/, $line);

	# Printing header section
	if ($line=~/^chr\tstart/) {
		my @tmpHeader = @tmp[6..@tmp-1];
		my @newArr;
		@headerArr = @tmpHeader;

		foreach my $sample (@tmpHeader) {
			next if !exists $::sampleHash{$sample};
			if ($::doCaseControl) {
				next if $::sampleHash{$sample}{CONTROL};
			}
			if (exists $::sampleHash{$sample}{REF_ID}) {
				next if @{$::sampleHash{$sample}{REF_ID}} < 1;
			}
			else {
				next;
			}
			push @newArr, $sample;
			push @newArr, "$sample\_s2n";
		}
		my @sampleArray = ();
		for my $s (@newArr) {
			if ($s !~/s2n/) {
				if ($::sampleHash{$s}{CONTROL}) {
					next;
				}
			}
			push @sampleArray, $s;
		}
		my $sampnames = join("\t", @sampleArray);
		print RATIOS "chr\tstart\tend\tinfo\tGC\tmap\t$sampnames\n";
		next;
	}
	
	# Body section
	my $coord = "$tmp[0]\t$tmp[1]\t$tmp[2]";

	# Excluding filtered rois if present
	next if $::filteredRois{$coord};

	my $coordinate = join "\t", @tmp[0..3];
	print RATIOS join "\t", @tmp[0..5];

	# Now calculating ratios for every sample. 
	my $j = 0;
	for (my $i = 6; $i < @tmp; $i++) {
		my $sample = @samples[$j];
        if ($::doCaseControl) {
            next if $::sampleHash{$sample}{CONTROL};
        }

		if ( !exists $::sampleHash{$sample} ) {
			$j++;
			next;
		}

		if (exists $::sampleHash{$sample}{REF_ID}) {
			if ( @{$::sampleHash{$sample}{REF_ID}} < 1 ) {
				$j++;
				next;
			}
		}
		else {
			next;
		}
		$j++;

		my @tmpRefs;
		# Calculating a baseline from the median corrected counts from each reference	
		foreach my $idx ( @{$::sampleHash{$sample}{REF_ID}} ) {
			push @tmpRefs, $tmp[6+$idx];
		}
		if (@tmpRefs) {
			# Median and signal to noise calculation
			my $meanBaseline = Utils::medianArray(@tmpRefs);
			my ($medianBaseline, $s2n) = Utils::signal2noise(@tmpRefs);
			$s2n = 5 if @tmpRefs == 1;
			my $ratio =  $meanBaseline > 0 ? sprintf "%.3f", $tmp[$i]/$meanBaseline : "0.00";

			print TOPLOTRATIOS join("\t", @tmp[0..5]) . "\t$sample\t$ratio\n";
			print RATIOS "\t$ratio\t$s2n";

			# Load JSON
			push @{$::analysisJson{$sample}{On_target}{Rois}}, $tmp[3];
			push @{$::analysisJson{$sample}{On_target}{Ratios}}, $ratio;
		}
	}
	print RATIOS "\n";
   }

   my $j = 7;
   my $k = 8;

   # Splitting Ratios file for segmentation and deleting outlier ratios
   foreach my $sample (	sort keys %::sampleHash ) {

	   	# Skipping database references 
		if (exists $::referenceHash{$sample}) {
			next;
		}

	   	# Skipping samples with no available reference
		if ( exists $::sampleHash{$sample}{REF_ID} && @{$::sampleHash{$sample}{REF_ID}} >= 1 ) {

			# Here, skipping header and selecting sample columns
			my $cmd = "$::cat $inputDir/$::outName.Ratios.bed | $::tail -n +2 | $::cut -f 1,2,3,4,$j,$k |"; 
			$cmd .=  "$::awk '{if (\$4 !~/info/) { print \$0 } }' > $inputDir/$sample.ratios.txt";
			system($cmd);

			#Utils::compressFile("$inputDir/$sample.ratios.txt") if !-s "$inputDir/$sample.ratios.txt.gz";
			Utils::compressFile("$inputDir/$sample.ratios.txt");

			# Now calculating mean and std deviation statistics
			my $statistics = `$::zcat $inputDir/$sample.ratios.txt.gz | $::tail -n +2 | $::cut -f 5 | $::awk '{if (\$1<0){\$1=-\$1}else{\$1=\$1} sum+=\$1; sumsq+=\$1*\$1} END {print sum/NR, sqrt(sumsq/NR - (sum/NR)^2)}'`;
			my ($mean, $stdev) = split(/\s/,$statistics);

			$::sampleHash{$sample}{ONTARGET_MEAN_RATIO}  = $mean;
			$::sampleHash{$sample}{ONTARGET_SD_RATIO}    = $stdev;

			my $sdevs = $::sampleHash{$sample}{ONTARGET_SD_RATIO}*$::minZscore;
			my $upper_limit = $::sampleHash{$sample}{ONTARGET_MEAN_RATIO} + $sdevs;
			my $lower_limit = $::sampleHash{$sample}{ONTARGET_MEAN_RATIO} - $sdevs;

			$::sampleHash{$sample}{UPPER_LIMIT} = $upper_limit;
			$::sampleHash{$sample}{LOWER_LIMIT} = $lower_limit;
			print " INFO: $sample\tratio:$mean\tStd.dev:$stdev\tUp_limit:$upper_limit\tlow_limit:$lower_limit\n";

			$j=$j+2;
			$k=$k+2;
		}
   }
      Utils::compressFile($::HoF{RATIOS_ON});

 }

####################
 sub isChrX {

  my $chr = shift;

  if ($chr =~/(chr)?(X|23)/) {
	return 1;
  }
  else {
	return 0;
  }

 }

 1;
