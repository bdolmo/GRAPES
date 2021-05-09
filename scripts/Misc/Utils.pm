#!/usr/bin/perl

package Utils;
use strict;
use File::Basename;
use Statistics::Descriptive;
use Sort::Key::Natural qw(natsort);
use Scalar::Util qw(looks_like_number);


##############################################
sub loadSingleExon2Json {

	my $dataFile    = shift;
	my $coordinates = shift;
	my $sample      = shift;

	$::analysisJson{$sample}{Calls}{$coordinates}{Plot_type} = "Single_exon";
	open (IN, "<", $dataFile) || die " ERROR: Unable to open $dataFile\n";
	while (my $line=<IN>) {
		chomp $line;
		my @tmp = split (/\t/, $line);

		my $sampleName = $tmp[4];
		push @{$::analysisJson{$sample}{Calls}{$coordinates}{$sampleName}{Ratios}}, $tmp[5];
	}
	close IN;
}


##############################################
sub loadShortPLot2Json {

	my $dataFile    = shift;
	my $coordinates = shift;
	my $sample      = shift;
	$::analysisJson{$sample}{Calls}{$coordinates}{Plot_type} = "Short";
	open (IN, "<", $dataFile) || die " ERROR: Unable to open $dataFile\n";
	while (my $line=<IN>) {
		chomp $line;
		my @tmp = split (/\t/, $line);
		#print OUT "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t$sample_value\t2_$sampName\n";
		my @sampleTmp = split("_", $tmp[5]);
		my $sampleName = $sampleTmp[1];
		push @{$::analysisJson{$sample}{Calls}{$coordinates}{$sampleName}{Ratios}}, $tmp[4];
	}
	close IN;
}

##############################################
sub loadLongPlot2Json {

	my $dataFile    = shift;
	my $coordinates = shift;
	my $sample      = shift;
	$::analysisJson{$sample}{Calls}{$coordinates}{Plot_type} = "Long";
	open (IN, "<", $dataFile) || die " ERROR: Unable to open $dataFile\n";
	while (my $line=<IN>) {
		chomp $line;
		my @tmp = split (/\t/, $line);
		my $sampleName = $tmp[4];
		#chr2	30142807	30143575	chr2:30142808-30143575	20B3334.rmdup	1.78528	-0.271	control
		#chr9	136291319	136291320	NM_139025.4_5_6;ADAMTS13	RB18972_9999999.rmdup	0.903041825095057	RB18972_9999999.rmdup
		push @{$::analysisJson{$sample}{Calls}{$coordinates}{$sampleName}{Ratios}}, $tmp[5];
	}
	close IN;
}

##############################################
sub readCountsFromBamIdx {
	my $sample = shift;

	my $counts = `$::samtools idxstats $sample | $::awk -F \'\t\' \'{s+=\$3+\$4}END{print s}'`;
	chomp $counts;
	return $counts;
}

##############################################
sub populateMapGcHashFromNormCountsOfftarget {

	my @offtargetNormalized = glob ("$::offtargetDir/*.normalized.bed.gz");
	foreach my $file (@offtargetNormalized) {
		open (IN, "$::zcat $file |") || die " ERROR: Unable to open $file\n";
		while (my $line=<IN>) {
			chomp $line;
			my @tmp = split (/\t/, $line);
			#chr1	701420	1260606	pdwindow_2	0.577	88.577	1.54	0.254			my @tmp = split (/\t/, $line);
			my $coordinate = join "\t", @tmp[0..2];
			my $gc  = $tmp[4];
			my $map = $tmp[5];
			$::ExonFeatures{$coordinate}{GC}  = $gc;
			$::ExonFeatures{$coordinate}{MAP} = $map;
		}
		close IN;
	}
}

#############################################
sub populateMapGcHashFromNormCountsOntarget {

	my $infile = shift;
	if ($infile =~/.gz/) {
		open (IN, "$::zcat $infile |") || die " ERROR: Unable to open $infile\n";
	}
	else {
		open (IN, "<", $infile ) || die " ERROR: Unable to open $infile\n";
	}
	while (my $line=<IN>) {
		chomp $line;
		#chr1	11906065	11906071	NM_006172_2_3;NPPA	42	100	53	49	59	63	63	54
		my @tmp = split (/\t/, $line);
		my $coordinate = join "\t", @tmp[0..2];
		my $gc  = $tmp[4];
		my $map = $tmp[5];
		$::ExonFeatures{$coordinate}{GC}  = $gc;
		$::ExonFeatures{$coordinate}{MAP} = $map;
	}
	close IN;

}
###################
sub compressFileBgzip {

	my $infile = shift;
	return if $infile =~/.gz/;

	my $cmd = "$::bgzip -c $infile > $infile.gz";
	system $cmd if -s $infile;
}
###################
sub tabixIndex {
	my $infile   = shift;
	my $filetype = shift;

	if ($infile !~/.gz$/) {
		print " ERROR: Could not index $infile. $infile is not gzipped\n";
	}

	my $cmd;
	if ($filetype eq 'vcf') {
		$cmd = "$::tabix -f -p vcf $infile";
	}
	elsif ($filetype eq 'bed') {
		$cmd = "$::tabix -f -p bed $infile";
	}
	else  {
		print " ERROR: Could not index $infile. filetype $filetype is not supported. Supported formats are: bed, vcf\n";
	}
	system $cmd if -s $infile;
}

###################
# forcing -f compression/decompression
sub decompressFile {

	my $infile = shift;
	return if $infile !~/.gz/;

	my $cmd = "$::gunzip -f $infile";
	system $cmd if -s $infile;
}
###################
sub compressFile {

	my $infile = shift;
	return if $infile =~/.gz/;

	my $cmd = "$::gzip -f $infile";
	system $cmd if -s $infile;
}
###################
# This function converts a genomic coordinate (e.g chr1:1000-1100) to an absolute position
# Taken from gendicall. Original author is Manuel Rueda, PhD
sub number2human {
    my $number = shift;
       if(looks_like_number($number)) {
        $number >= 1000000 ? sprintf( "%0.2fM", $number / 1000000 )
      : $number >= 1000    ? sprintf( "%0.0fK", $number / 1000 )
      :                      $number;
       }
       else
       {
            $number;
       }
}

###################
# This function converts a genomic coordinate (e.g chr1:1000-1100) to an absolute position
sub chrCoordinate2Absolute {
	my $chr = shift;
	my $genomicCoordinate = shift;

	my $chrLen = $::ChromosomeLengths{$chr};
	my $absolutePosition = 0;
	foreach my $chromosome ( natsort keys %::ChromosomeLengths) {
		if ($chr eq $chromosome) {
			$absolutePosition+=$genomicCoordinate;
			last;
		}
		$absolutePosition+=$::ChromosomeLengths{$chromosome};
	}
	return $absolutePosition;
}
###################
# Calculate Median absolute Deviation. This is more robust to outliers than standard deviation
sub MAD {

	my @array = @_;

	my $median = medianArray(@array);
	my @substractMedian = ();

	foreach my $value (@array) {
		my $x = abs($value-$median);
		push @substractMedian, $x;
	}

	my $MAD = abs( medianArray (@substractMedian) );
	return sprintf "%.2f",$MAD;
}

###################
# Median calculation using a reference as input param
sub medianRefVar {

 # Using reference
 my $data = shift;

 # Derreferencing
 my @data = @$data;
 my $stat = Statistics::Descriptive::Full->new();
 $stat -> add_data(@data);
 my $median = $stat->median();
 return $median;
}

###################
# Median calculation using a copied array values as input param
sub medianArray {
 my @data = @_;
 my $stat = Statistics::Descriptive::Full->new();
 $stat -> add_data(@data);
 my $mean = $stat->median();
 $mean = 0 if !$mean;
 return sprintf "%.3f",$mean;
}

###################
# Mean calculation using a copied array values as input param
sub meanArray {
 my @data = @_;
 my $stat = Statistics::Descriptive::Full->new();
 $stat -> add_data(@data);
 my $mean = $stat->mean();
 return sprintf "%.3f",$mean;
}

###################
sub meanHash {
 # Using reference
 my $data = shift;
 # Derreferencing
 my @data = @$data;
 my $stat = Statistics::Descriptive::Full->new();
 $stat -> add_data(@data);
 my $mean = $stat->mean();
 return $mean;
}

###################
sub TrimmedMean {

 my @array = @_;
 my $stat = Statistics::Descriptive::Full->new();
 $stat -> add_data(@array);
 my $tm  =  $stat->trimmed_mean(.10);
 return $tm;
}

###################
sub std {

 my @array = @_;
 my $stat = Statistics::Descriptive::Full->new();
 $stat -> add_data(@array);
 my $sd = $stat->standard_deviation();
 $sd = 0 if !$sd;
 return $sd;
}

###################
sub signal2noise {

 my @array =@_;
 my $median = Utils::medianArray(@array);
 my $sd =  Utils::std(@array);

 my $signal2noise = 0;
 if ($sd > 0 ) {
 	$signal2noise = sprintf "%.2f", $median/$sd;
 }
 return ($median, $signal2noise);
}

1;
