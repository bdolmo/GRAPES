#!/usr/bin/perl
package LOH;
use strict;
use Getopt::Long;
use File::Basename; 
use List::Util qw(min max);

# Rules for accepting variants:

# 1. Avoid indels 
# 2. Minimum coverage accepted of 30X
# 3. Minimum MAPQ of 50
# 4. Minimum mean base quality of 25

# Methods for identifying LOH events

sub varCall {

	my $bamFile = shift;
	my $genome  = shift;
	my $coordinate = shift;
	my $meanCoverage = shift;
	my $SVTYPE  = shift;

	# This array will be populated with found SNV along with their genotypes
	my @ArrSNV = ();

	# Initiallizing hash of frequencies
	my %Freqs = (
		"A" => 0,
		"C" => 0,
		"T" => 0,
		"G" => 0
	);

	open (MPILEUP , "$::samtools mpileup -r $coordinate $bamFile -q 30 -C 50 -d 180 -f $genome $devNullStderr |");
	while (my $line=<MPILEUP>) {
		chomp $line;
		#chr2	155829801	g	32	.$..,,..,..,.,..,,,...,,...,,..,,	CDDJJFFIJJJJJJJIHHGGJIJGHIHHDJEE

        if ($line =~/\-|\+/) {
         # print "$line\n"; exit;
        }
		my @tmp = split (/\t/, $line);
		my $depth = $tmp[3];

		# Skipping positions with less than 15 reads
		next if $depth < 10;
        $tmp[4] = uc ($tmp[4]);

		my @ntds = split (//, $tmp[4]);

		my $i = 0;
		my @idxs = ();
        my $flag_indel = 0;
		foreach my $nuc (@ntds) {

            if ($nuc =~/\-|\+/) {
                $flag_indel = 1;
            }
            if ($flag_indel == 1) {
                if ($nuc !~/[0-9]+/ && $nuc !~/[ACTG]/ && $nuc !~/\-|\+/) {
                    $flag_indel = 0;
                }
            }
            if ($flag_indel == 0) {
                if ( exists $Freqs{$nuc}) {
                    push @idxs, $i;
                    $Freqs{$nuc}++;
                }
          		$i++;
            }
		}
		my $max = 0;
		my $alternative;
		foreach my $base (sort keys %Freqs) {
			if ($Freqs{$base} > $max) {
				$max = $Freqs{$base};
				$alternative = $base;
			}
		}
        #print "$max\n";
        %Freqs = (
		    "A" => 0,
		    "C" => 0,
		    "T" => 0,
		    "G" => 0
	    );
		my $freq;
		if ($depth > 0) {
			$freq = sprintf "%.2f", $max/$i;
		}
		# Simple threshold genotype
		my $genotype;
		if ($freq >= 0.12 && $freq < 0.8) {
			$genotype = "HET";
                   # print "$tmp[0]\t$tmp[1]\t$tmp[2]\t$alternative\t$genotype\t$freq\n";
		}
		elsif ($freq > 0.8) {
			$genotype = "HOM";
                  #  print "$tmp[0]\t$tmp[1]\t$tmp[2]\t$alternative\t$genotype\t$freq\n";
		}

    	# Using hard-threshold similarly than in VarDict
		if ($freq >= 0.12) {
			push @ArrSNV, $genotype;
		}
	}
	close MPILEUP;

	my %Concordance = (
		"HET" => 0,
		"HOM" => 0
	);
	foreach my $genotype (@ArrSNV) {
		$Concordance{$genotype}++;
	}

	my $HMZ_RATE = 0;
	my $nSNV = scalar @ArrSNV;

	if ($nSNV > 0) {
		$HMZ_RATE = sprintf "%.2f", $Concordance{"HOM"}/scalar(@ArrSNV);
	}
	if ($HMZ_RATE >= 0.8) {
		print " INFO: DEL\t$coordinate\tSNV=$nSNV\tHOM_rate=$HMZ_RATE\n";
		
	}
	return $nSNV, $HMZ_RATE;

}

return 1;