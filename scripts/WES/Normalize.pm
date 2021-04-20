#!/usr/bin/env perl

package Normalize;

use strict;
use Getopt::Long;
use File::Basename; 
use Data::Dumper;
use Parallel::ForkManager;
use Sort::Key::Natural qw(natsort);

##########################
 sub normalizeOntarget {

	my $outputDir    = shift;
	my $analysisMode = shift; # exome or gene-panel
	my $type 	     = shift; # ontarget
    my $hashOfRegion = shift; # Reference of a hash containing roi info

   	# Corrected read counts
	my $outCorrectedCounts = $::HoF{NORM_COUNTS_ON};

	if ($analysisMode eq 'gene-panel') {
		if ( (-s $::HoF{NORM_COUNTS_ON} ||  -s "$::HoF{NORM_COUNTS_ON}.gz") && (-s $::HoF{NORM_COVERAGE_ON} || -s "$::HoF{NORM_COVERAGE_ON}.gz") ) {
			return 1;
		}
	}
	# Read Raw Depth
	readRawDepth($outputDir, $hashOfRegion, $type);

	if ($type eq "offtarget" && !-e"$outputDir/$::outName.ReadCounts.bed.gz") {
		print " WARNING: skipping off-target normalization\n";
		return 0;
	}

	my @tmpOfftarget    = ();
	my @tmpOfftargetChr = ();

   	# Normalization
	foreach my $sample ( natsort keys %::sampleHash ) {

		$::sample{$sample}{MEDIANCOUNTS} = Utils::medianRefVar($::sample{$sample}{ARR_COUNTS});

		if (!$::sample{$sample}{MEDIANCOUNTS}) {
			print " WARNING: Non-existent $sample median counts\n";
		}

		# Calculate binned median counts per gc interval
		foreach my $gc_interval (sort keys %{$::GC{$sample}}) {
			if ($::GC{$sample}{$gc_interval}{ARR_COUNTS}) {
				$::GC{$sample}{$gc_interval}{MEDIAN} = Utils::medianRefVar($::GC{$sample}{$gc_interval}{ARR_COUNTS});
				$::GC{$sample}{$gc_interval}{ARR_COUNTS} = ();
			}
			if ($::GC{$sample}{$gc_interval}{ARR_COUNTS_X}) {
				$::GC{$sample}{$gc_interval}{MEDIAN_X} = Utils::medianRefVar($::GC{$sample}{$gc_interval}{ARR_COUNTS_X});
				$::GC{$sample}{$gc_interval}{ARR_COUNTS_X} = ();
			}
		}
		if ($type eq 'ontarget') {
			# Calculate binned median counts per length interval
			foreach my $length_interval ( sort keys %::ExonLength ) {
				if ($::ExonLength{$length_interval}{ARR_COUNTS}) {
					$::ExonLength{$length_interval}{MEDIAN} = Utils::medianRefVar($::ExonLength{$length_interval}{ARR_COUNTS});
					$::ExonLength{$length_interval}{ARR_COUNTS} = ();
				}
				if ($::ExonLength{$length_interval}{ARR_COUNTS_X}) {
					$::ExonLength{$length_interval}{MEDIAN_X} = Utils::medianRefVar($::ExonLength{$length_interval}{ARR_COUNTS_X});
					$::ExonLength{$length_interval}{ARR_COUNTS_X} = ();
				}
			}
		}

		# Performing GC normalization using rolling median procedure
		foreach my $coordinate( natsort keys %$hashOfRegion ) {

			my $chr = ( split /\t/, $coordinate)[0];

			my $gc_interval = $::ExonFeatures{$coordinate}{GC_INTERVAL};

			# 2. Normalizing scaled-Counts by GC
			if ( ( isChrX($chr) ) && (exists $::GC{$sample}{$gc_interval}{MEDIAN_X} ) ) {
				$hashOfRegion->{$coordinate}->{ $sample }->{GC_NORM} = sprintf "%.3f", ( $hashOfRegion->{$coordinate}->{ $sample }{LIBNORM}*($::sample{$sample}{MEDIANCOUNTS}/$::GC{$sample}{$gc_interval}{MEDIAN_X}));		
			}
			elsif  ( (!isChrX($chr)) && (exists $::GC{$sample}{$gc_interval}{MEDIAN} ) ) {
				$hashOfRegion->{$coordinate}->{ $sample }->{GC_NORM} = sprintf "%.3f", ( $hashOfRegion->{$coordinate}->{ $sample }{LIBNORM}*($::sample{$sample}{MEDIANCOUNTS}/$::GC{$sample}{$gc_interval}{MEDIAN}));
			}
			else {
				$hashOfRegion->{$coordinate}->{ $sample }->{GC_NORM} = $hashOfRegion->{$coordinate}->{ $sample }{LIBNORM};
			}

			if ( isChrX($chr) ) {
				push @tmpOfftargetChr, $hashOfRegion->{$coordinate}->{ $sample }->{GC_NORM} if $type eq 'offtarget';
			}
			else {
				push @tmpOfftarget, $hashOfRegion->{$coordinate}->{ $sample }->{GC_NORM} if $type eq 'offtarget';
			}

			my $length_interval = $::ExonFeatures{$coordinate}{LENGTH_INTERVAL};

			# 3. Normalizing GC-corrected counts by region length
			if ($::ExonFeatures{$coordinate}{LENGTH_INTERVAL} && $type eq 'ontarget') {
				$hashOfRegion->{$coordinate}->{ $sample }->{LENGTH_NORM}
				 = sprintf "%.3f", ( $hashOfRegion->{$coordinate}->{ $sample }{LIBNORM}*($::sample{$sample}{MEDIANCOUNTS}/$::ExonLength{$length_interval}{MEDIAN}));
			}
		}
		if ($type eq 'offtarget') {

			if (@tmpOfftarget) {
				$::sampleHash{$sample}{OFFTARGETMEDIAN} = Utils::medianArray(@tmpOfftarget);
				@tmpOfftarget = ();
			}
			if (@tmpOfftargetChr) {
				$::sampleHash{$sample}{OFFTARGETMEDIAN_X} = Utils::medianArray(@tmpOfftargetChr);
				@tmpOfftargetChr = ();
			}
		}
	}

	# Temporal file for plotting purposes
	my $biasInfo = "$outputDir/bias_info.txt";
	
	open (TMP, ">", $biasInfo) || die " ERROR: Cannot open $biasInfo\n";
	open (CORRECTED, ">", $::HoF{NORM_COUNTS_ON}) || die " ERROR: Cannot open $::HoF{NORM_COUNTS_ON}\n";

	# Grouping header info
	my @tmpSamp;
	my $str;
	foreach my $sample ( natsort keys %::sampleHash ) {

		my $header_sample = $sample;
		$header_sample =~s/.bam//;

		push @tmpSamp, "$header_sample\_counts";
		push @tmpSamp, "$header_sample\_gc_corrected";
		push @tmpSamp, "$header_sample\_length_corrected";

		$str.= "\t$header_sample";
	}
	print TMP "chr\tstart\tend\tinfo\tgc\tmap\t" . join ("\t", @tmpSamp) . "\n";
	print CORRECTED "chr\tstart\tend\tinfo\tgc\tmap$str\n";

	# Writing normalized values
	foreach my $coordinate( natsort keys %$hashOfRegion ) {

		my @tmp = split (/\t/, $coordinate);
		my $chr = $tmp[0];
		my $coord = "$tmp[0]\t$tmp[1]\t$tmp[2]";

		print TMP "$coordinate\t$::ExonFeatures{$coordinate}{GC}\t$::ExonFeatures{$coord}{MAP}";
		print CORRECTED "$coordinate\t$::ExonFeatures{$coordinate}{GC}\t$::ExonFeatures{$coord}{MAP}";

		# Accessing secondary reference
		my $sample_ref =  $hashOfRegion->{$coordinate};

		foreach my $sample ( natsort keys %$sample_ref) {
			next if $sample eq 'GC';

			$sample_ref->{ $sample }->{LENGTH_NORM} = $sample_ref->{ $sample }->{GC_NORM};

			print CORRECTED "\t$sample_ref->{ $sample }->{GC_NORM}";
			print TMP "\t$sample_ref->{$sample}->{LIBNORM}\t$sample_ref->{$sample}->{GC_NORM}\t$sample_ref->{ $sample }->{LENGTH_NORM}";
		}
		print TMP "\n";
		print CORRECTED "\n";
	}

	Utils::compressFile($outCorrectedCounts);

	%::GC = ();
	%::ExonLength = ();

	if ($analysisMode eq 'gene-panel') {
		normalizePerBase($outputDir);
	}
 }

 ##########################
 # Re-normalizing for single-exon calling purpose
 sub normalizePerBase {
  
   my $outputDir = shift;

   print " INFO: Normalizing per base coverage\n";

   # Now for Coverage	
   my $coverageFile        = "$outputDir/$::outName.Coverages.bed.gz";
   my $outNormLibCoverages = "$outputDir/$::outName.norm_bylib_coverage.bed";
   my $normalizedCoverage =  "$outputDir/$::outName.NormalizedCoverage.bed";

   if (!-e $normalizedCoverage . ".gz") {

		open (COV, "$::zcat $coverageFile |") || die " ERROR: Cannot open $coverageFile\n";
		open (NORMLIBCOVS, ">", $outNormLibCoverages) || die " ERROR: Cannot open $outNormLibCoverages\n";

		my $count = 0;
		my @header;
		my $coordinate;
		my @norms;
		while (my $line=<COV>) {
			chomp $line;
			my @tmp = split (/\t/, $line);
			if ($count == 0) {
				@header = @tmp;
				$count++;
				print NORMLIBCOVS "$line\n";
				next;
			}
			else {
				# Filling hash of Ontarget (exons)
				for (my $i = 5; $i <= @tmp-1; $i++) {

					$tmp[$i] = 1 if $tmp[$i] == 0;
					$coordinate = join ("\t", @tmp[0..3]);
			
					# Normalized Coverage
					my $normalized;
					if ( isChrX($tmp[0]) ) {
						$normalized = sprintf "%.3f", ( $tmp[$i]/$::sampleHash{$header[$i]}{MEANCOVX_ONTARGET} );
						push @norms, $normalized;
					}
					else {
						$normalized = sprintf "%.3f", ( $tmp[$i]/$::sampleHash{$header[$i]}{MEANCOV_ONTARGET} );
						push @norms, $normalized;
					}
				}
			}
			print NORMLIBCOVS "$coordinate\t$tmp[4]\t" . join ("\t" , @norms) . "\n";
			@norms = ();
			$count++;
		}
		close COV;
		close NORMLIBCOVS;

   my $cmd = "$::sort -k5,5 $outNormLibCoverages > $outputDir/Sorted_GC_coverages.bed";
   system($cmd) if -e "$outputDir/Sorted_GC_coverages.bed";

   my $count = 0;
   my $n = 7;

   if ( -z "$outputDir/$::outName.NormalizedCoverage.bed.gz" ) {

		my %GC_intervals = ();

		for (my $i = 35;$i<=80; $i++) {
			$GC_intervals{$i}{INI} = undef;
			$GC_intervals{$i}{FI} = undef;
		}

		my $j = 35;
		my $count = 1;
		my @tmparray = ();
		open (IN, "<", "$outputDir/Sorted_GC_coverages.bed");
		while (my $line =<IN>){
			chomp $line;
			if ($line =~/start/) {
				$count++; next;
			}
			my @tmp = split (/\t/, $line);
			my $gc = $tmp[4];

			if ($gc < 35 || $gc > 80) {
				$count++;
				next;
			}
			else {
				if ($gc >= $j && $gc <= $j+1) {
					push @tmparray, $count;
				}
				elsif ($gc > $j+1) {
					my $ini = $tmparray[0];
					my $end = $tmparray[-1];

					if ($ini && $end) {
						my $tmpStr = `$::head -n $end $outputDir/Sorted_GC_coverages.bed | $::tail -n \$(($end-$ini+1))`;
						chomp $tmpStr;

						my @tmpLines = split (/\n/, $tmpStr);
						my $i = 5;
						foreach my $sample (natsort keys %::sampleHash) {
							my @array = ();
							foreach my $mline (@tmpLines) {
								#chr1	237821243	237821244	NM_001035_53_54;RYR2	35.000000	0.478	0.450	0.417	0.557	0.430
								my @tmp = split (/\t/, $mline);
								push @array, $tmp[$i];
							}
							my $sum = 0;
							foreach my $element (@array) {
								$sum+=$element;
							}
							my $mean = $sum/scalar(@array);
							$::GC{$sample}{$i}{MEDIAN} = $mean;
							$i++;
						}
						$j++;
						@tmparray = ();
					}
				}
			}
		}
		$count++;
	}
   
	#unlink( "$outputDir/Sorted_GC_coverages.bed" );

	# Applying normalization
	open (IN, "<", $outNormLibCoverages) || die " ERROR: Cannot open $outNormLibCoverages\n";
	open (OUT, ">", $normalizedCoverage) || die " ERROR: Cannot open $normalizedCoverage\n";
	while (my $line =<IN>) {
		chomp $line;
		if ($line =~/^chr\tstart/) {
			print OUT "$line\n";
			next;
		}
		my @tmp = split (/\t/, $line);
		my $n = 5;
		my $gc_interval;
		for (my $i = 10;$i<=90; $i++) {
			if ($tmp[4] > $i-1 && $tmp[4] <=$i) {
				$gc_interval = $i;
				last;
			}
		}		
		#print "$line\t$gc_interval\n";

		my $str;

		foreach my $sample (natsort keys %::sampleHash) {
			my $gc_norm;
			if (defined $gc_interval) {
				if (exists $::GC{$sample}{$gc_interval}{MEDIAN}) {
					$gc_norm = $tmp[$n];
				}
				else {
					$gc_norm = $tmp[$n];
				}
			}
			else {
				$gc_norm = $tmp[$n];
			}
			$str.= "\t$gc_norm";
			$n++;
		}
		print OUT  join("\t", @tmp[0..4]) . $str . "\n";	
	   }
	   Utils::compressFile($normalizedCoverage);
	   Utils::compressFile($outNormLibCoverages);
   }
 }

 #####################
 sub readMappability {

	my $mappabilityFile = shift;

	 open (MAP, "<", $mappabilityFile) || die " ERROR: Cannot open $mappabilityFile\n";
	 while (my $line=<MAP>){
		chomp $line;
		my @tmp = split (/\t/, $line);
	    my $coordinate = join "\t", @tmp[0..2];
		my $map = $tmp[-1] > 100 ? 100 : $tmp[-1];
		$::ExonFeatures{$coordinate}{MAP} = $map if !exists $::ExonFeatures{$coordinate}{MAP};
	 }
	 close MAP;
 }

 ############################################################################################################
 # Here we read the file with all counts, and we normalize them by library size (ALL_norm_bylib_counts.bed) #
 # Additionally, we save every counts per sample per GC interval					  					    #
 ############################################################################################################

 sub readRawDepth {
	
	my $outputDir = shift;
	my $hashOfRegion = shift;
	my $type = shift;

	# Input Raw read counts and coverages
 	my $countFile;	
	my $mappabilityFile;

	if ($type eq 'ontarget') {					
	 	$countFile       = "$outputDir/$::outName.ReadCounts.bed.gz";	
		$mappabilityFile = "$outputDir/mappability_ontarget.bed";	
	}								
	if ($type eq 'offtarget') {	
	 	$countFile       = "$outputDir/$::outName.ReadCounts.joint.bed.gz";
		$mappabilityFile = "$outputDir/mappability_offtarget.bed";
	}			

	# add mappability to %::ExonFeatures 
	readMappability($mappabilityFile);
		
	# Output-> normalized counts/covs by library size
	#my $outNormLibCounts    = "$outputDir/$::outName.norm_bylib_counts.bed";

	# library normalized counts file
	#open (NORMLIBCOUNTS, ">", $outNormLibCounts) || die " ERROR: Cannot open $outNormLibCounts\n";

	# Regions not meeting a minimum creteria of GC content, mappability and length are filtered and annotated in this file
	open (FILTERED, ">", "$outputDir/filtered_regions.txt") || die " ERROR: Cannot open $outputDir/filtered_regions.txt\n"; 
	print FILTERED "chr\tstart\tend\tinfo\tavg_coverage\t\%gc\t\%mappability\n";

	if (!-e $countFile) {
		print " ERROR: $countFile has not been created\n"; exit;
	}
	if (-z $countFile) {
		print " ERROR: $countFile is empty\n"; exit;
	}

	# opening raw counts input file (ALL_ReadCounts.[joint]bed)
	open (COUNT, "$::zcat $countFile |") || die " ERROR: Cannot open $countFile\n";

	my $flag = 0;
	my $count = 0;
	my $coordinate;
	my @header;
	my @norms = ();

	while (my $line=<COUNT>) {

		chomp $line;
		my @tmp = split (/\t/, $line);

		# Printing header
		if ($count == 0) {
			@header = @tmp;
			$count++;
			next;
		}
		else {

			my $chr    = $tmp[0];
			my $gc_int = int ($tmp[4]) +1;
			my $length = $tmp[2] - $tmp[1] +10;

			# Filling hash of sample's normalized coverage
			for (my $i = 5; $i <= @tmp-1; $i++) {

				# Ad-hoc solution for setting 1 read when there's no coverage
				$tmp[$i] = 1 if $tmp[$i] == 0;

				$coordinate = join ("\t", @tmp[0..3]);

				my $length     = $tmp[2] - $tmp[1];
				my $sample     = $header[$i];
				my $raw_counts = $tmp[$i];
				my $gc         = $tmp[4];
				my $map        = $tmp[5];

				my @tmp2  = split (/\t/, $coordinate);
				my $coord = "$tmp2[0]\t$tmp2[1]\t$tmp2[2]";

				# Filtering
				my $avg_coverage = $raw_counts * (100/$length);

				$::ExonFeatures{$coordinate}{GC} = $gc;
				$::ExonFeatures{$coord}{GC} = $gc;
				$::ExonFeatures{$coord}{MAP} = $map if $type eq 'offtarget';

				# Region GC
				$hashOfRegion->{$coordinate}->{ $sample }->{GC} = $gc;

				# Filtering criteria
				if ($type eq 'ontarget') {
					#if ($::ExonFeatures{$coord}{MAP} < 40 || $length < 10 ) {
					if ($length < 10) {
						$::filteredRois{$coord}++;
						print FILTERED "$coordinate\t$avg_coverage\t$gc\t$::ExonFeatures{$coord}{MAP}\tSMALL_SIZE\t$length\n";
						$flag = 1;
						delete ($hashOfRegion->{$coordinate});
						last;
					}
					if ($::ExonFeatures{$coord}{MAP} < 50 ) {
						$::filteredRois{$coord}++;
						print FILTERED "$coordinate\t$avg_coverage\t$gc\t$::ExonFeatures{$coord}{MAP}\tLOW_MAPPABILITY\t$length\n";
						$flag = 1;
						delete ($hashOfRegion->{$coordinate});
						last;
					}
				}
				if ($type eq 'offtarget') {												
					if ($::ExonFeatures{$coord}{MAP} < 50 ) {
						$::filteredRois{$coord}++;			
						print FILTERED "$coordinate\t$avg_coverage\t$gc\t$::ExonFeatures{$coord}{MAP}\tLOW_MAPPABILITY\t$length\n";			
						$flag = 1;												
						delete ($hashOfRegion->{$coordinate});									
						last;													
					}														
				}															
					
				# Raw Counts will be stored inside this global hash
				$hashOfRegion->{$coordinate}->{ $sample }{RAW} = $raw_counts;

				# Region Length
				$::ExonFeatures{$coordinate}{LENGTH} =  $length;
				$::ExonFeatures{$coord}{LENGTH} =  $length;

				#  Normalize Counts by Region length
				my $perBaseCount = $type eq 'ontarget' ? sprintf "%.3f",  10e6*($raw_counts/$length) : $raw_counts;

				#  Normalize Counts by library size
				if ($type eq 'ontarget') {
					if ( isChrX($chr) ) {
						$hashOfRegion->{$coordinate}->{ $sample }->{LIBNORM} =  sprintf "%.3f", 10e2*( $perBaseCount/$::sampleHash{$sample}{TOTALREADS_ONTARGET_X});
					}
					else {
						$hashOfRegion->{$coordinate}->{ $sample }->{LIBNORM} =  sprintf "%.3f", 10e2*( $perBaseCount/$::sampleHash{$sample}{READSONTARGET});
					}
				}
				if ($type eq 'offtarget') {
					if ( isChrX($chr) ) {
						$hashOfRegion->{$coordinate}->{ $sample }->{LIBNORM} =  sprintf "%.6f", 10e6*( $perBaseCount/$::sampleHash{$sample}{TOTALREADS_OFFTARGET_X});
					}
					else {
						$hashOfRegion->{$coordinate}->{ $sample }->{LIBNORM} =  sprintf "%.6f", 10e6*( $perBaseCount/$::sampleHash{$sample}{READSOFFTARGET});
					}
				}

				my $libnormalized = $hashOfRegion->{$coordinate}->{ $sample }->{LIBNORM};

				if ($chr !~ 'Y' || $chr ne 'chr23')  {
					push @norms, $libnormalized;
					push @{$::sample{$sample}{ARR_COUNTS}}, $libnormalized;
				}

				if ( $::GC{$sample}{$gc_int} ) {
					if ( isChrX($chr) ) {
						push @{$::GC{$sample}{$gc_int}{ARR_COUNTS_X}}, $libnormalized;
						$::ExonFeatures{$coordinate}{GC_INTERVAL} = $gc_int;
					}
					else {
						push @{$::GC{$sample}{$gc_int}{ARR_COUNTS}}, $libnormalized;
						$::ExonFeatures{$coordinate}{GC_INTERVAL} = $gc_int;
					}
				}
				if ( $::ExonLength{$length} ) {
					if ( isChrX($chr) ) {
						push @{$::GC{$sample}{$gc_int}{ARR_COUNTS_X}}, $libnormalized;
						$::ExonFeatures{$coordinate}{GC_INTERVAL} = $gc_int;
					}	
					else {
						push @{$::ExonLength{$length}{ARR_COUNTS}}, $libnormalized;
						$::ExonFeatures{$coordinate}{LENGTH_INTERVAL} = $length;
					}
				}
			}
		}
		if ($flag == 0) {
			#print NORMLIBCOUNTS "$coordinate\t$tmp[4]\t" . join ("\t" , @norms) . "\n";
		}
		$flag = 0;
		@norms = ();
		$count++;
	}
	close COUNT;
	#close NORMLIBCOUNTS;
 }

 sub readOfftargetFile {

	 my $countFile = shift;
	 my $sample    = shift;
	 my $GC_offtarget = ();

	 my @tmpFile    = ();
	 my @addCounts  = ();
	 my @addCountsX = ();

	 open (FILE, "$::zcat $countFile |") || die " ERROR: Unable to open $countFile\n";
	 while (my $line=<FILE>) {

	 	chomp $line;
		push @tmpFile, $line;

		my @tmp = split (/\t/, $line);
		my $chr = $tmp[0];
		my $start  = $tmp[1];
		my $end    = $tmp[2];
		my $info   = $tmp[3];
		my $counts = $tmp[6];
		my $gc_int = int ($tmp[4]*100) +1;
		my $coordinate = "$chr\t$start\t$end\t$info";

		if ( isChrX($chr) ) {
			push @addCountsX, $counts;
		}
		else {
			push @addCounts, $counts;
		}

		if ( $::GC{$sample}{$gc_int} ) {
			if ( isChrX($chr) ) {
				push @{$::GC{$sample}{$gc_int}{ARR_COUNTS_X}}, $counts;
			} 
			else {
				push @{$::GC{$sample}{$gc_int}{ARR_COUNTS}}, $counts;
			}
		}
	}
	close FILE;

	$::sampleHash{$sample}{OFFTARGETMEDIAN} = Utils::medianArray(@addCounts);
	$::sampleHash{$sample}{OFFTARGETMEDIANX} = Utils::medianArray(@addCountsX);

	return @tmpFile;
 }

##########################
 sub joinBins {
	
	my $inputFile = shift;
	my $outputDir = shift;

	my %HoP = ();
	my %Regions = ();

	my $name = basename($inputFile);
	$name =~s/.normalized.bed/.normalized.joined.bed/;
	my $output = "$outputDir/$name";

	open (IN, "<", $inputFile) || die " ERROR: unable to open $inputFile\n";

	while (my $line=<IN>) {
		chomp $line;
		my @tmp = split (/\t/, $line);
		my @info = split (";", $tmp[3]);

		$Regions{"$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t$tmp[4]\t$tmp[5]\t$tmp[6]"} = $info[0];
		$HoP{$info[0]}{COUNTS}+=$tmp[7];
		$HoP{$info[0]}{MAP}+=$tmp[6];
	}
	close IN;

	open (OUT, ">", $output) || die " ERROR: unable to open $output\n";
	foreach my $region (natsort keys %Regions) {
		print OUT "$region\t$HoP{$Regions{$region}}{COUNTS}\n";
	}
	close OUT;
	unlink ($inputFile);
 }

##########################
 sub normalizeOfftarget {

	my $offtargetDir = shift;

	foreach my $sample (natsort keys %::sampleHash ) {

		my @normCounts;
		my @normCountsX;

		# Re-initializing GC hash
		for (my $i = 10;$i<=90; $i++) {
			$::GC{$sample}{$i}{ARR_COUNTS}     = undef;
			$::GC{$sample}{$i}{MEDIAN_COUNTS}  = undef;
			$::GC{$sample}{$i}{ARR_COUNTS_X}   = undef;
			$::GC{$sample}{$i}{MEDIAN_COUNTS_X}= undef;
		}	

		my $countFile = "$offtargetDir/$sample.offtarget_joined_counts.bed.gz";

		if (!-e $countFile || -z $countFile) {
			next;
		}

		my @tmpFile   = readOfftargetFile($countFile, $sample);

		foreach my $gc_interval (natsort keys %{$::GC{$sample}}) {
			if ($::GC{$sample}{$gc_interval}{ARR_COUNTS} && $::GC{$sample}{$gc_interval}{ARR_COUNTS} > 50 ) {
				$::GC{$sample}{$gc_interval}{MEDIAN} = Utils::medianRefVar($::GC{$sample}{$gc_interval}{ARR_COUNTS});
				$::GC{$sample}{$gc_interval}{ARR_COUNTS} = ();
			}
			if ($::GC{$sample}{$gc_interval}{ARR_COUNTS_X} && $::GC{$sample}{$gc_interval}{ARR_COUNTS_X} > 50) {
				$::GC{$sample}{$gc_interval}{MEDIAN_X} = Utils::medianRefVar($::GC{$sample}{$gc_interval}{ARR_COUNTS_X});
				$::GC{$sample}{$gc_interval}{ARR_COUNTS_X} = ();
			}
		}

		open OUT, ">", "$offtargetDir/$sample.normalized.bed";
		foreach my $line (@tmpFile) {
			chomp $line;
			#chr1	578745	771659	pdwindow_4;piece_1	0.417414	26.1383443952743	275
			my @tmp = split (/\t/, $line);
			my $chr    = $tmp[0];
			my $start  = $tmp[1];
			my $end    = $tmp[2];
			my $info   = $tmp[3];
			my $counts = $tmp[6];
			my $gc_int = int ($tmp[4]*100) +1;
			my $normGC;
			if ( $::GC{$sample}{$gc_int} ) {
				if ( isChrX($chr) ) {
					if ($::GC{$sample}{$gc_int}{MEDIAN_X}) {
						$normGC = sprintf "%.3F", ( $counts*($::sampleHash{$sample}{OFFTARGETMEDIANX}/$::GC{$sample}{$gc_int}{MEDIAN_X}));
						push @normCountsX, $normGC;
					}
					else {
						next;
					}
				}
				else {
					if ($::GC{$sample}{$gc_int}{MEDIAN}) {
						$normGC = sprintf "%.3F", ( $counts*($::sampleHash{$sample}{OFFTARGETMEDIAN}/$::GC{$sample}{$gc_int}{MEDIAN}));
					}
					else {
						$normGC = $counts;
					}
					push @normCounts, $normGC;
				}
				print OUT "$chr\t$start\t$end\t$info\t$tmp[4]\t$tmp[5]\t$counts\t$normGC\n";
			}
		}
		close OUT;
		$::sampleHash{$sample}{OFFTARGET_NORMALIZED}   = Utils::medianArray(@normCounts)  if @normCountsX > 0;
		$::sampleHash{$sample}{OFFTARGET_NORMALIZED_X} = Utils::medianArray(@normCountsX) if @normCountsX > 0;

		@normCounts = ();
		@normCountsX = ();
          
        joinBins("$offtargetDir/$sample.normalized.bed", $offtargetDir);

		unlink "$offtargetDir/$sample.normalized.bed.gz";

		rename "$offtargetDir/$sample.normalized.joined.bed", "$offtargetDir/$sample.normalized.bed";
	    Utils::compressFile("$offtargetDir/$sample.normalized.bed");
	}
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
