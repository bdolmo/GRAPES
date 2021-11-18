#!/usr/bin/env perl
package Extract;
use strict;
use Getopt::Long;
use File::Basename;
use List::Util qw(min max);
use Sort::Key::Natural qw(natsort);

###############################
# sub convertMosdepth2Grapes {
# 	my $sample = shift;
# 	# Output files fro mosdepth:
# 	#*.mosdepth.global.dist.txt
# 	#*.mosdepth.region.dist.txt
# 	#*.mosdepth.summary.txt
# 	#*.regions.bed.gz
# 	#*.regions.bed.gz.csi

# 	my $cmd = "$::gunzip -c $::offtargetDir/$sample.offtarget.regions.bed.gz > $::offtargetDir/$sample.offtarget_counts.tmp.bed";
# 	print "$cmd\n" if $::verbose;
# 	system $cmd;

# 	unlink("$::offtargetDir/$sample.offtarget.mosdepth.global.dist.txt");
# 	unlink("$::offtargetDir/$sample.offtarget.mosdepth.region.dist.txt");
# 	unlink("$::offtargetDir/$sample.offtarget.mosdepth.summary.txt");
# }

###############################
sub annotateGC {
  my $inputBed  = shift;
  my $outputBed = shift;

  open (OUT, ">", $outputBed) || die " ERROR: Unable to open $outputBed\n";
  open IN, "$::bedtools nuc -fi $::genome -bed $inputBed |";
  while (my $line=<IN>)  {
    chomp $line;
    next if $line=~/^#/;
    my @tmp = split("\t", $line);
    my $gc  = sprintf "%.3f",(100*$tmp[5]);
    print OUT join("\t", @tmp[0..3]) . "\t" . $gc . "\n";
  }
  close IN;
  close OUT;

  return $outputBed;
}

###############################
sub addOfftargetRegionName {

   my $offtargetBed = shift;
   my $megadepthBed = shift;

   # Read offtarget bed and save coordinate names in a hash
   my %HoRegion = ();
   open (BED, "<", $offtargetBed) || die " ERROR: Unable to open offtarget bed $offtargetBed\n";
   while (my $line=<BED>) {
     chomp $line;
     my ($chr, $start, $end, $info) = (split "\t", $line);
     my $coordinate = $chr . "\t" . $start . "\t" . $end;
     $HoRegion{$coordinate} = $info;
   }
   close BED;

   # Now recover region name and rewrite megadepth output
   my $tmpBed = $megadepthBed;
   $tmpBed=~s/.bed/.tmp.bed/;
   open (TMP, ">", $tmpBed) || die " ERROR: Unable to open $tmpBed\n";
   open (IN, "<", $megadepthBed) || die " ERROR: Unable to open $megadepthBed\n";

   while (my $line=<IN>) {
     chomp $line;
     my @tmp = split("\t", $line);
     my $coordinate = "$tmp[0]\t$tmp[1]\t$tmp[2]";
     my $regionName = $HoRegion{$coordinate};
     print TMP "$tmp[0]\t$tmp[1]\t$tmp[2]\t$regionName\t$tmp[3]\n";
   }
   close IN;
   close TMP;

   unlink $megadepthBed;
   rename $tmpBed, $megadepthBed;

 }

###############################
sub extractOfftargetCounts {

	my $bamDir = shift;
	my @offtargetBams= glob("$bamDir/*bam");

  my $offtargetGc = $::HoF{GLOBAL_OFFTARGET_BED};
  $offtargetGc =~s/.bed/.gc.bed/;

  my $offtargetGcMap  = $offtargetGc;
  $offtargetGcMap =~s/.bed/.map.bed/;

  # Annotate GC content
  if (!-e $offtargetGc)  {
    $offtargetGc = annotateGC( $::HoF{GLOBAL_OFFTARGET_BED}, $offtargetGc );
  }
  # Annotate mappability
  if (!-e $offtargetGcMap)  {
    $offtargetGcMap = annotateMappabilityOfftarget( $offtargetGc,
      $offtargetGcMap );
  }

	foreach my $bam (@offtargetBams) {

		my $sample = basename($bam);
		$sample =~s/.offTarget.bam//;

    # offtarget regions bed
    #my $offtargetBed          = "$bamDir/$sample.offtarget.bed";
    my $offtargetBed          = $::HoF{GLOBAL_OFFTARGET_BED};
    # offtarget regions bed with counts at 4th column
    my $offtargetTmpCounts    = "$bamDir/$sample.offtarget_counts.tmp.bed";

    # joined smaller regions that belong to the parent window
    my $offtargetJoinedCountsGz = "$bamDir/$sample.offtarget_joined_counts.bed.gz";
    my $offtargetJoinedCounts   = "$bamDir/$sample.offtarget_joined_counts.bed";

		# Skipping off-target analysis on samples with less than specified off-target reads
		if ($::sampleHash{$sample}{READSOFFTARGET} < $::minOfftargetReads) {
			print " INFO: Skipping off-target analysis on $sample (low num: $::sampleHash{$sample}{READSOFFTARGET})\n";
			#$::pm->finish;
			#next;
		}
		else {
			print " INFO: Extracting offtarget counts of sample $sample\n";
		}

		if (!-s $offtargetJoinedCountsGz) {

			if ( -s $offtargetBed ) {

        if (!-e $offtargetTmpCounts) {
          my $cmd = "$::offtargetExtractor $bam $::genome $offtargetGcMap > $offtargetTmpCounts";
          print "$cmd\n";
          system $cmd;
        }

				# Joining single pieces that come from the same window.
				# This is because we split each defined window based on surrounding peak/ROIs
        if (!-e $offtargetJoinedCounts) {
				  joinBins( $offtargetTmpCounts, $offtargetJoinedCounts );
        }

        # Deleting temoporary files
				unlink("$bamDir/$sample.offtarget.regions.bed.gz");
				unlink("$bamDir/$sample.offtarget.regions.bed.gz.csi");
				# unlink($offtargetCountsGc);
				#unlink($offtargetTmpCounts);
			}
		}
 	}
  #exit;
}

###############################
sub joinBins {

	my $inputBed  = shift;
	my $outputBed = shift;

  my %Regions = ();
	my $min = 10e20;
	my $max = 0;
	my %seen = ();

	open (IN, "<", $inputBed) || die " ERROR: unable to open $inputBed\n";
	while (my $line=<IN>) {
		chomp $line;
    #chr1	843529	1084543	pdwindow_1;piece_16	61.829	96.083	58
		my @tmp = split (/\t/, $line);
		my @info = split (";", $tmp[3]);
		my $windowName = $info[0];

		if (!exists $seen{$windowName} ) {
			$min = 10e20;
			$max = 0;
		}

		$seen{$windowName}++;
		$Regions{$windowName}{CHR} = $tmp[0];

		my $length = $tmp[2]-$tmp[1];
		$Regions{$windowName}{COUNTS}+=$tmp[-1];
		$Regions{$windowName}{MAP}+=($tmp[5]*$length);
		$Regions{$windowName}{GC}+=($tmp[4]*$length);
		$Regions{$windowName}{EFFECTIVE_LENGTH}+=$length;

		if ($tmp[1] < $min) {
			$Regions{$windowName}{START} = $tmp[1];
			$min = $tmp[1];
		}
		if ($tmp[2] > $max) {
			$Regions{$windowName}{END} = $tmp[2];
		}
	}
	close IN;
	open (OUT, ">", $outputBed) || die " ERROR: unable to open $outputBed\n";
	foreach my $region (natsort keys %Regions) {

		my $length = $Regions{$region}{END}-$Regions{$region}{START};
		my $map = sprintf "%.3f", $Regions{$region}{MAP}/$Regions{$region}{EFFECTIVE_LENGTH};
		my $gc  = sprintf "%.3f", $Regions{$region}{GC}/$Regions{$region}{EFFECTIVE_LENGTH};

		print OUT "$Regions{$region}{CHR}\t$Regions{$region}{START}\t$Regions{$region}{END}\t"
		. "$region\t$gc\t$map\t$Regions{$region}{COUNTS}\n";
	}
	close OUT;

	#unlink ($countFile);
	Utils::compressFile($outputBed);
 }

###############################
sub annotateMappabilityOfftarget {

  my $inputBed  = shift;
  my $outputBed = shift;

  my $outputDir = dirname($outputBed);
  my $sample    = basename($inputBed);

  my $mappability = $outputDir . "/" . "mappability.offtarget.bed";

  if (!-e $inputBed) {
	  print " WARNING: non existent $inputBed. Skipping mappability extraction\n";
	  return 0;
  }
  if (-z $inputBed ) {
	  print " WARNING: empty $inputBed. Skipping mappability extraction\n";
	  return 0;
  }
  my $cmd;
  if (!-s $outputBed) {
  	# Intersecting mappability track with the targeted regions
    if ($::hasChr) {
      $cmd = "$::bedtools intersect -a $::mappTrack -b $inputBed -wo | ";
    }
    else {
      $cmd = "perl -e -p \"s/chr//\" $inputBed | $::bedtools intersect -a $::mappTrack -b stdin -wo -sorted | ";
    }

    #chr1	8500026	8500030	0.333333	chr1	8483908	8521275	1497	0.442824	4
    $cmd .= "$::awk '{ size=\$7-\$6; marginal=100*(\$NF*\$4)/size; print \$5\"\t\"\$6\"\t\"\$7\"\t\"\$8\"\t\"\$4\"\t\"\$9\"\t\"marginal}'";
  	$cmd .= "| perl $::getMapOfftarget | $::sort -V  > $mappability";
  	print "$cmd\n" if $::verbose;
  	system($cmd);

  	# File editing
  	$cmd = "$::cut -f 5 $mappability > $outputDir/only_mapp.txt";
    system $cmd;

  	$cmd = "$::paste $inputBed $outputDir/only_mapp.txt > $outputBed";
    system $cmd;

  	open (MAP, "<", $outputBed) || die " ERROR: Cannot open $outputBed\n";
  	while (my $line=<MAP>){
  		chomp $line;
  		my @tmp = split (/\t/, $line);
  		my $coordinate = join "\t", @tmp[0..2];
  		my $map = $tmp[-1] > 100 ? 100 : $tmp[-1];
  		$::ExonFeatures{$coordinate}{MAP} = $map;
  	}
  	close MAP;
  	unlink("$outputDir/only_mapp.txt");
  }
  else {
	   print " INFO: skipping Mappability extraction. File was already created\n";
  }
  return $outputBed;
}

###############################
sub getMappability {

  my $bedfile      = shift;
  my $type         = shift;
  my $outputDir    = shift;
  my $hashOfRegion = shift;

  my $outfile = "$outputDir/mappability.bed";

  if ($type eq 'ontarget') {
	   $outfile = "$outputDir/mappability_ontarget.bed";
  }
  if ($type eq'offtarget') {
	   $outfile = "$outputDir/mappability_offtarget.bed";
  }

  if (!-e $bedfile) {
	  print " WARNING: non existent $bedfile. Skipping mappability extraction\n";
	  return 0;
  }
  if (-z $bedfile ) {
	  print " WARNING: empty $bedfile. Skipping mappability extraction\n";
	  return 0;
  }

  if (!-e $outfile || -z $outfile) {

  	# Intersecting mappability track with the targeted regions
    my $cmd;
    if ($::hasChr) {
      $cmd .= "$::bedtools intersect -a $::mappTrack -b $bedfile -wo | ";
    }
    else {
      $cmd .= "perl -e -p \"s/chr//\" $bedfile | $::bedtools intersect -a $::mappTrack -b stdin -wo | ";
    }
    $cmd .= "$::sort -V | ";
    $cmd .= "$::awk '{ size=\$7-\$6; marginal=100*(\$NF*\$4)/size; print \$5\"\t\"\$6\"\t\"\$7\"\t\"\$8\"\t\"\$4\"\t\"\$9\"\t\"marginal}'";

  	if ($type eq 'offtarget') {
        $cmd .= "| perl $::getMapOfftarget | $::sort -V  > $outfile";
  	}
  	if ($type eq 'ontarget') {
        $cmd .= " | perl $::getMapOntarget | $::sort -V  > $outfile";
  	}
  	if ($::verbose) {
  		print "$cmd\n";
  	}
  	system($cmd);
 }
 else {
	print " INFO: skipping Mappability extraction. File was already created\n";
 }

 if ($type eq 'ontarget' || $type eq 'offtarget') {
	 my @exon_len;
	 open (MAP, "<", $outfile) || die " ERROR: Cannot open $outfile\n";
	 while (my $line=<MAP>){
		chomp $line;
		my @tmp = split (/\t/, $line);
		my $coordinate = join "\t", @tmp[0..2];
		#my $map = $tmp[-1] > 100 ? 100 : $tmp[-1];
    my $map = $tmp[-1];

		$::ExonFeatures{$coordinate}{MAP} = $map;

		my $length = $tmp[2]-$tmp[1];
		push @exon_len, $length;
	 }
	 close MAP;

	 # Extract exon length
	 $::median_exon_length = Utils::medianArray(@exon_len);
	 my $smallest_exon = min (@exon_len);
	 my $biggest_exon  = max (@exon_len);

	 for (my $i = $smallest_exon; $i < $biggest_exon; $i=$i+10) {
		$::ExonLength{$i}{ARR_COUNTS} = undef;
	 }
 }

}

1;
