#!/usr/bin/env perl
package Extract;
use strict;
use Getopt::Long;
use File::Basename; 
use List::Util qw(min max);
use Sort::Key::Natural qw(natsort);

###############################
sub convertMosdepth2Grapes {
	my $sample = shift;
	# Output files fro mosdepth:
	#*.mosdepth.global.dist.txt
	#*.mosdepth.region.dist.txt
	#*.mosdepth.summary.txt
	#*.regions.bed.gz
	#*.regions.bed.gz.csi

	my $cmd = "$::gunzip -c $::offtargetDir/$sample.offtarget.regions.bed.gz > $::offtargetDir/$sample.offtarget_counts.tmp.bed";
	print "$cmd\n" if $::verbose;
	system $cmd;

	unlink("$::offtargetDir/$sample.offtarget.mosdepth.global.dist.txt");
	unlink("$::offtargetDir/$sample.offtarget.mosdepth.region.dist.txt");
	unlink("$::offtargetDir/$sample.offtarget.mosdepth.summary.txt");
}
###############################
 sub extractOfftargetCounts {

	my $bamDir = shift;
	
	my @offtargetBams= glob("$bamDir/*bam");

	foreach my $bam (@offtargetBams) {
		
		my $pid = $::pm -> start() and next;

		my $sample = basename($bam);
		$sample =~s/.offTarget.bam//;

		# Skipping off-target analysis on samples with less than specified off-target reads
		if ($::sampleHash{$sample}{READSOFFTARGET} < $::minOfftargetReads) {
			print " WARNING: Skipping off-target analysis on $sample (low num: $::sampleHash{$sample}{READSOFFTARGET})\n";
			$::pm->finish;
			#next;
		}
		else {
			print " INFO: Extracting offtarget counts of sample $sample\n";
		}

		if (!-s "$bamDir/$sample.offtarget_joined_counts.bed.gz") {

			my $off_bed = "$::outDir/OFF_TARGET/$sample.offtarget.bed";

			if ( -s $off_bed ) {

				# extract mean coverage per region using Mosdepth
				#my $cmd = "$::mosdepth -t $::threads -n -x --by $off_bed $bamDir/$sample.offtarget $bam";
				#print "$cmd\n" if $::verbose;
				#system $cmd;

				# Make mosdepth output files compatible with downstream processes
				#convertMosdepth2Grapes($sample);

				# Extract read counts across all off-target regions
				my $cmd = "$::offtargetExtractor $bam $::genome $off_bed > $bamDir/$sample.offtarget_counts.tmp.bed";
     			system $cmd;

				# Get GC content, this is slow for very large windows. Use precomputed file? Tabix?
				$cmd = "$::bedtools nuc -fi $::genome -bed $bamDir/$sample.offtarget_counts.tmp.bed |";
				$cmd.= " $::cut -f 1,2,3,4,5,7 |";
				$cmd.= " $::awk '{ print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$4\"\t\"\$6\"\t\"\$5 }' > $bamDir/$sample.offtarget_counts_nomap.bed";
				print "$cmd\n" if $::verbose;
				system $cmd if !-s "$bamDir/$sample.offtarget_counts_nomap.bed";

				# Remove header
				$cmd = "$::sed -i '1d' $bamDir/$sample.offtarget_counts_nomap.bed";
				system $cmd;

				# Annotate mappability of every window. Again, thi slows down the process. Use mean MAPQ value per region instead?
				getMappabilityOfftarget("$bamDir/$sample.offtarget_counts_nomap.bed", $bamDir);

				# Joining single pieces that come from the same window. 
				# This is because we split each defined window based on surrounding peak/ROIs
				joinBins("$bamDir/$sample.offtarget_counts.bed", $bamDir);

				unlink("$bamDir/$sample.offtarget.regions.bed.gz");
				unlink("$bamDir/$sample.offtarget.regions.bed.gz.csi");
				unlink("$bamDir/$sample.offtarget_counts_nomap.bed");
				unlink("$bamDir/$sample.offtarget_counts.tmp.bed");

			}
		}
		#else {
	#		print " INFO: Skipping Offtarget count extraction of sample $sample\n";
		#}
		$::pm->finish;
 	}
	$::pm->wait_all_children;	
 }

###############################
 sub joinBins {
	
	my $countFile = shift;
	my $outputDir = shift;

	my %Regions = ();

	my $name = basename($countFile);
	$name =~s/.offtarget_counts.bed/.offtarget_joined_counts.bed/;
	my $output = "$outputDir/$name";

	my $min = 10e20;
	my $max = 0;
	my %seen = ();

	open (IN, "<", $countFile) || die " ERROR: unable to open $countFile\n";
	while (my $line=<IN>) {
		chomp $line;
		my @tmp = split (/\t/, $line);
		
		my @info = split (";", $tmp[3]);

		my $window_name = $info[0];

		if (!exists $seen{$window_name} ) {
			$min = 10e20;
			$max = 0;
		}

		$seen{$window_name}++;
		$Regions{$window_name}{CHR} = $tmp[0];

		my $length = $tmp[2]-$tmp[1];
		$Regions{$window_name}{COUNTS}+=$tmp[5];
		$Regions{$window_name}{MAP}+=($tmp[6]*$length);
		$Regions{$window_name}{GC}+=($tmp[4]*$length);
		$Regions{$window_name}{EFFECTIVE_LENGTH}+=$length;

		if ($tmp[1] < $min) {
			$Regions{$window_name}{START} = $tmp[1];
			$min = $tmp[1];
		}
		if ($tmp[2] > $max) {
			$Regions{$window_name}{END} = $tmp[2];
		}
	}
	close IN;

	open (OUT, ">", $output) || die " ERROR: unable to open $output\n";
	foreach my $region (natsort keys %Regions) {

		my $length = $Regions{$region}{END}-$Regions{$region}{START};
		my $map = sprintf "%.3f", $Regions{$region}{MAP}/$Regions{$region}{EFFECTIVE_LENGTH};
		my $gc  = sprintf "%.3f", $Regions{$region}{GC}/$Regions{$region}{EFFECTIVE_LENGTH};
		print OUT "$Regions{$region}{CHR}\t$Regions{$region}{START}\t$Regions{$region}{END}\t$region\t$gc\t$map\t$Regions{$region}{COUNTS}\n";
	}
	close OUT;
	
	#unlink ($countFile);
	Utils::compressFile($output);

 }

###############################
sub getMappabilityOfftarget {

  my $countsFile= shift;
  my $outputDir = shift;

  my $sample = basename($countsFile);
  $sample =~s/.offtarget_counts_nomap.bed//;

  my $mappability = "$outputDir/$sample.offtarget.mappability.bed";
  my $outfile     = "$outputDir/$sample.offtarget_counts.bed";
    
  if (!-e $countsFile) {
	  print " WARNING: non existent $countsFile. Skipping mappability extraction\n";
	  return 0;
  }
  if (-z $countsFile ) {
	  print " WARNING: empty $countsFile. Skipping mappability extraction\n";
	  return 0;
  }

  if (!-s $outfile) {

	# Intersecting mappability track with the targeted regions
	my $cmd .= "$::bedtools intersect -a $::mappTrack -b $countsFile -wo | ";
    $cmd .= "$::awk '{ size=\$7-\$6; marginal=100*(\$NF*\$4)/size;  print \$5\"\t\"\$6\"\t\"\$7\"\t\"\$8\"\t\"\$4\"\t\"\$9\"\t\"marginal}'";
	$cmd .= "| perl $::getMapOfftarget | $::sort -V  > $mappability";
	print "$cmd\n" if $::verbose;
	system($cmd);

	# File editing
	`$::cut -f 5 $mappability > $outputDir/only_mapp.$sample.txt`;
	`$::paste $countsFile $outputDir/only_mapp.$sample.txt > $outfile`;

	open (MAP, "<", $outfile) || die " ERROR: Cannot open $outfile\n";
	while (my $line=<MAP>){
		chomp $line;
		my @tmp = split (/\t/, $line);
		my $coordinate = join "\t", @tmp[0..2];
		my $map = $tmp[-1] > 100 ? 100 : $tmp[-1];
		$::ExonFeatures{$coordinate}{MAP} = $map;
	}
	close MAP;
	unlink("$outputDir/only_mapp.$sample.txt");
  }
  else {
	print " INFO: skipping Mappability extraction. File was already created\n";
  }
}

###############################
 sub getMappability {

  my $bedfile   = shift;
  my $type      = shift;
  my $outputDir = shift;
  my $hashOfRegion = shift;

  my $outfile;
  if ($type eq 'intronic') {
	$outfile = "$outputDir/mappability_intronic.bed";
  }
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
	my $cmd;
	# Intersecting mappability track with the targeted regions
	$cmd .= "$::bedtools intersect -a $::mappTrack -b $bedfile -wo | ";
	if ($type eq 'offtarget') {
	    $cmd .= "$::sort -V | ";
	    $cmd .= "$::awk '{ size=\$7-\$6; marginal=100*(\$NF*\$4)/size;  print \$5\"\t\"\$6\"\t\"\$7\"\t\"\$8\"\t\"\$4\"\t\"\$9\"\t\"marginal}' | perl $::getMapOfftarget | sort -V  > $outfile";
	}
	if ($type eq 'ontarget' || $type eq 'intronic') {
	    $cmd .= "$::sort -V | ";
	    $cmd .= "$::awk '{ size=\$7-\$6; marginal=100*(\$NF*\$4)/size;  print \$5\"\t\"\$6\"\t\"\$7\"\t\"\$8\"\t\"\$4\"\t\"\$9\"\t\"marginal}' | perl $::getMapOntarget | sort -V  > $outfile";
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
		my $map = $tmp[-1] > 100 ? 100 : $tmp[-1];

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
