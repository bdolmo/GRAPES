#!/usr/bin/env perl

package Intron;
use strict;
use Getopt::Long;
use File::Basename; 
use List::Util qw(min max);

 sub Analyze {

	my $str = `$::cat $::bed | $::sort -V`;

	if (!$str) {
		print " ERROR: $::bed file could not be sorted properly\n"; exit;
	}

	my @ExonArray = split (/\n/, $str);
	my $index = 0;
	foreach my $exon (@ExonArray) {
		#print "$exon\n";
		$::ExonFeatures{$exon}{MAP} = undef;
		my @current = split (/\t/, $exon);
		if ($index == 0) {
			my $previousExon = "none";
			my $nextExon     = $ExonArray[$index+1];
			my @tmpNext      = split (/\t/, $nextExon);
			my $nextLength   = $tmpNext[1] - $current[2];
			
			if ($tmpNext[0] eq $current[0] && $nextLength <= 200 && $nextLength > 0) {

				$::ExonFeatures{$exon}{NEXT_INTRON} = "$tmpNext[0]\t$current[2]\t$tmpNext[1]\tINTRON_$current[2];$tmpNext[1]";
				$::IntronFeatures{$::ExonFeatures{$exon}{NEXT_INTRON}}{MAP} = undef;
				$::IntronFeatures{$::ExonFeatures{$exon}{NEXT_INTRON}}{GC} = undef;
				$::IntronFeatures{$::ExonFeatures{$exon}{NEXT_INTRON}}{MEAN_COV} = undef;
			}
		}
		elsif ($index > 0 && $index < @ExonArray-1) {

			my $previousExon = $ExonArray[$index-1];
			my $nextExon     = $ExonArray[$index+1];

			# info from the previous exon
			my @tmpPrevious = split (/\t/, $previousExon);
			my $prevLength  = $current[1] - $tmpPrevious[2];
			
			# info from the next exon
			my @tmpNext     = split (/\t/, $nextExon);
			my $nextLength  = $tmpNext[1] - $current[2];

			if ($tmpPrevious[0] eq $current[0] && $prevLength <= 200  && $prevLength > 0) {

				$::ExonFeatures{$exon}{PREVIOUS_INTRON} = "$tmpPrevious[0]\t$tmpPrevious[2]\t$current[1]\tINTRON_$tmpPrevious[2];$current[1]";
				$::IntronFeatures{$::ExonFeatures{$exon}{PREVIOUS_INTRON}}{MAP} = undef;
				$::IntronFeatures{$::ExonFeatures{$exon}{PREVIOUS_INTRON}}{GC} = undef;
				$::IntronFeatures{$::ExonFeatures{$exon}{PREVIOUS_INTRON}}{MEAN_COV} = undef;
			}

			if ($tmpNext[0] eq $current[0] && $nextLength <= 200 && $nextLength > 0) {
				$::ExonFeatures{$exon}{NEXT_INTRON} = "$tmpNext[0]\t$current[2]\t$tmpNext[1]\tINTRON_$current[2];$tmpNext[1]";
				$::IntronFeatures{$::ExonFeatures{$exon}{NEXT_INTRON}}{MAP} = undef;
				$::IntronFeatures{$::ExonFeatures{$exon}{NEXT_INTRON}}{GC} = undef;
				$::IntronFeatures{$::ExonFeatures{$exon}{NEXT_INTRON}}{MEAN_COV} = undef;
			}
		}
		elsif ( $index == @ExonArray-1 ) {
			my $previousExon = $ExonArray[$index-1];
			my $nextExon     = "none";
			my @tmpPrevious = split (/\t/, $previousExon);
			my $prevLength  = $current[1] - $tmpPrevious[2];

			if ($tmpPrevious[0] eq $current[0] && $prevLength <= 200 && $prevLength > 0) {
				$::ExonFeatures{$exon}{PREVIOUS_INTRON} = "$tmpPrevious[0]\t$tmpPrevious[2]\t$current[1]\tINTRON_$tmpPrevious[2];$current[1]";
				$::IntronFeatures{$::ExonFeatures{$exon}{PREVIOUS_INTRON}}{MAP} = undef;
				$::IntronFeatures{$::ExonFeatures{$exon}{PREVIOUS_INTRON}}{GC} = undef;
				$::IntronFeatures{$::ExonFeatures{$exon}{PREVIOUS_INTRON}}{MEAN_COV} = undef;
			}
		}
		$index++;
	}

	open (INTRON, ">", "$::outDir/short_introns.tmp.bed") || die " ERROR: unable to open $::outDir/short_introns.tmp.bed\n";
	foreach my $intron ( sort { $a cmp $b } keys %::IntronFeatures) {
		print INTRON "$intron\n";
	}
	close INTRON;

	my $cmd = "$::sort -V $::outDir/short_introns.tmp.bed > $::outDir/SHORT_INTRONS/short_introns.bed";
	system $cmd;
	
	unlink ("$::outDir/short_introns.tmp.bed");

	my $intronDir = "$::outDir/SHORT_INTRONS";

	if (!-e "$intronDir/SHORT_INTRONS.Coverages.bed") {

		if ($::input) {
			$cmd = "$::targetDepth -i $::input -o $intronDir -n SHORT_INTRONS -g $::genome -b $::outDir/SHORT_INTRONS/short_introns.bed -c -d -t 4";
			system ($cmd);
		}
		if ($::control && $::test) {
			$cmd = "$::targetDepth -i $::control,$::test -n SHORT_INTRONS -o $intronDir -g $::genome -b $::outDir/SHORT_INTRONS/short_introns.bed -c -d -t 4";
			system ($cmd);
		}
	}

	# Extracting mappability
 	print " INFO: Extracting Intron Mappability\n";
 	Extract::getMappability("$::outDir/SHORT_INTRONS/short_introns.bed", "intronic", $intronDir, \%::IntronFeatures);
	getValidIntrons();

 }

 sub getValidIntrons {

  my $str = `cat $::outDir/SHORT_INTRONS/mappability_intronic.bed`;
  chomp $str;
 
  my %tmpMap = ();
  my @tmpStr = split (/\n/, $str);
  foreach my $c (@tmpStr) {
	my @tmp = split (/\t/, $c);
	my $coordinate = join ("\t", @tmp[0..3]);
	my $map = $tmp[4];
	$::tmpMap{$tmp[3]}{MAP} = $map;
  }
  
  # Selecting introns that accomplish:
  # Minimum mean coverage of 30
  # Minimum mappability of 90%

  # Reading coverage file
  my @header;
  my $count = 0;
  my $flag = 0;
  my $exonName_A;
  my $chr_A;
  my $exonName_B;
  my $chr_B;
  my %coverageBySample = ();
  my $meanCov;

  my $countLines = `$::cat $::outDir/SHORT_INTRONS/SHORT_INTRONS.Coverages.bed | $::wc -l | $::cut -f1`; 
  chomp $countLines;


  open (OUT_TMP, ">", "$::outDir/SHORT_INTRONS/intron_summary.txt");	
  print OUT_TMP "chr\tstart\tend\tID\t";
  open (COV, "<", "$::outDir/SHORT_INTRONS/SHORT_INTRONS.Coverages.bed") || die " ERROR: Unable to open $::outDir/SHORT_INTRONS/SHORT_INTRONS.Coverages.bed\n";
  while (my $line=<COV>) {
	chomp $line;

	my @tmp = split (/\t/, $line);
	if ($count == 0) {
		print OUT_TMP "$line\n";
		@header = @tmp;
		$count++;
		next;
	}
	else {
		if ($flag == 0) {
			$chr_A = $tmp[0];
			$exonName_A = $tmp[3];
			$flag = 1;	
		}
		if ($flag == 1) {
			$exonName_B = $tmp[3];
			$chr_B = $tmp[0];
			my $mappability = $::tmpMap{$exonName_B}{MAP};

			if ($exonName_A eq $exonName_B) {
				for (my $i = 5; $i <= @tmp-1; $i++) {
					push( @{ $coverageBySample { $header[$i] } }, $tmp[$i]); 
				}	
			}
			if ($exonName_A ne $exonName_B || $countLines-1) {
				my @tmporal = split (/[_;]/, $exonName_A);
                my $tmpcoord = "$chr_A\t$tmporal[1]\t$tmporal[2]\t$exonName_A";
                print OUT_TMP "$tmpcoord";

				foreach my $sample ( sort keys %coverageBySample ) {
					$meanCov = Utils::meanHash($coverageBySample{$sample});
					my @tmp2 = split (/[_;]/, $exonName_A);
					my $coordinate = "$chr_A\t$tmp2[1]\t$tmp2[2]\t$exonName_A";
						
					$::IntronFeatures{$coordinate}{MAP} = $mappability;
					$::IntronFeatures{$coordinate}{GC}  = $tmp[3];
					$::IntronFeatures{$coordinate}{$sample} = $meanCov;

					print OUT_TMP "\t$meanCov";
				}
                print OUT_TMP "\n";
				$exonName_A = $tmp[3];
				$chr_A = $tmp[0];
				%coverageBySample = ();
			}
		 }
	}
	$count++;
  }
 }
 

return 1;
