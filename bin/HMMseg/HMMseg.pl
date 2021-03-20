#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use Statistics::Descriptive;
use Config;
use Sort::Key::Natural qw(natsort);

our $dirname = dirname (__FILE__);

my $java = `which java`;
chomp $java;

my $input       = $ARGV[0];
my $binSize     = $ARGV[1];
my $outfile     = $ARGV[2];
my $sampleName  = $ARGV[3];
my $outputDir   = $ARGV[4];
my $centromeres = $ARGV[5];
my $bamFile     = $ARGV[6];
my $genome      = $ARGV[7];
my $meanCoverage= $ARGV[8];

if (@ARGV < 9 ) {
	print "\nUsage: perl $0 <copy_ratios.bed> <bin_size> <output_file> <sample_name> <output_dir> <centromers> <blacklist.bed> <bam_file> <genome_reference> <mean_cov>\n\n";exit;
}

my $hmmseg = "$dirname/HMMSeg.jar";
my $extractMapp = "$dirname/offtargetMap.pl";

our $grep;
our $sort;
if ($Config{osname} =~/darwin/) {
	$grep  = `which egrep`; chomp $grep;
	$sort  = `which gsort`; chomp $sort;
}
else {
	$grep  = `which grep`; chomp $grep;
	$sort  = `which sort`; chomp $sort;
}

our $cat  = `which cat`;
chomp $cat;
if (!$cat) {
	print "ERROR: cat command not found\n";
		exit;
}

our $wc  = `which wc`;
chomp $wc;
if (!$wc) {
	print "ERROR: wc command not found\n";
	exit;
}

our $awk  = `which awk`;
chomp $awk;
if (!$awk) {
	print "ERROR: awk command not found\n";
	exit;	
}

our $cut  = `which cut`;
chomp $cut;
if (!$cut) {
	print "ERROR: cut command not found\n";
	exit;
}

my $head  = `which head`;
chomp $head;
if (!$head) {
	print "ERROR: head command not found\n";
	exit;
}

my $bedtools = `which bedtools`;
chomp $bedtools;

if (!$bedtools) {
    print " ERROR: BEDtools was not found on PATH\n";
	exit;
}
else {
    my $bedtools_version = `$bedtools --version | $cut -d \" \" -f 2 | $awk '{ split (\$0, a, \".\");  print a[2] }'`;
    chomp $bedtools_version;
	if ($bedtools_version < 19) { 
    	print " ERROR: BEDtools version $bedtools_version is too old. Please, consider installing the newest version\n"; exit;
    }
}

if (!$input || !-e $input) {
	print " ERROR: please introduce an input file\n";
}

my $checkChr = `$head -1 $input`;
chomp $checkChr;

my @chrs = ();
my $mappTrack = "$dirname/../../mappability/wgEncodeCrgMapabilityAlign100mer.chr.bedgraph.gz";
my $blacklist   = "$dirname/../../db/QDNAseq_blacklisted_hg19.chr.bed";

if ($checkChr =~/^chr/) {
	@chrs = qw( chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY);
}
else {
	@chrs = qw( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y);
	$mappTrack = "$dirname/../../mappability/wgEncodeCrgMapabilityAlign100mer.nochr.bedgraph.gz";
	$blacklist   = "$dirname/../../db/QDNAseq_blacklisted_hg19.nochr.bed";
}

if (!-e "$outputDir/tosegment.$sampleName" ) {

	# Formatting input file
	open (IN, "<", $input) || die " Unable to open $input\n";
	open (OUT, ">", "$outputDir/tosegment.$sampleName") || die " Unable to open $outputDir/tosegment.$sampleName\n";

	while (my $line=<IN>){
		chomp $line;
		my @tmp = split (/\t/, $line);

		next if $line =~/\t0\t0\t0/;

		# Skipping low mappability bins
		if ($tmp[3] < 10) {
			next;
		}
		print OUT "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t$tmp[6]\n";
	}
	close IN;
	close OUT;
}

createModel();
segmentRatios("$outputDir/tosegment.$sampleName", $sampleName);

unlink("$outputDir/tosegment.$sampleName");
unlink ("$outputDir/regions.list");
unlink ("$outputDir/model.hmm");

my $outWithMapp = $outfile;
$outWithMapp=~s/.CNV.bed/.MAP.CNV.bed/;

print " INFO: Annotating mappability\n";

# Getting mappability
my $cmd .= "$bedtools intersect -a $mappTrack -b $outfile -wo";
$cmd .= "| awk '{ size=\$7-\$6; marginal=100*(\$NF*\$4)/size;  print \$5\"\t\"\$6\"\t\"\$7\"\t\"\$8\"\t\"\$4\"\t\"\$9\"\t\"marginal}'";
$cmd .= "| perl $extractMapp | $sort -V  > $outWithMapp";
system $cmd;
 
unlink $outfile;

filterLowQualCNVs($outWithMapp);

##########################
sub filterLowQualCNVs {

	my $infile = shift;

	my $cmd = "$bedtools intersect -a $infile -b $centromeres -v | $bedtools intersect -a $input -b stdin -wa -wb > $outputDir/before_filtering.tmp.txt";
	system $cmd;
	
	my %Calls = ();

	open (IN, "<", "$outputDir/before_filtering.tmp.txt") || die " ERROR: Unable to open $outputDir/before_filtering.tmp.txt";
	while (my $line=<IN>) {
		chomp $line;

		my @tmp = split (/\t/, $line);
		my @info = split (/;/, $tmp[10]);

		my ($mappability) = grep ($_=~/MAP=/, @info);
		$mappability =~s/MAP=//;

		my $ratio = $tmp[6];
	    my $mapq = $tmp[4];

		my $coordinate = "$tmp[7]\t$tmp[8]\t$tmp[9]\t$tmp[10]";

		push @ { $Calls{$coordinate}{RATIOS} }, $ratio;
		push @ { $Calls{$coordinate}{MAP} }, $mapq;
	}
	close IN;

	my $outfile = $outWithMapp;
	$outfile=~s/.MAP//;

	open (OUT, ">", $outfile) || die " ERROR: Unable to open $outfile\n";
	foreach my $call (natsort keys %Calls) {

		my $mean_ratio = getMean( @{$Calls{$call}{RATIOS} } );
		my $mapq = getMean( @{$Calls{$call}{MAP} } );

		my $MAD = MAD( @{$Calls{$call}{RATIOS} } );

		next if $mapq < 20;

		my @tmp = split ("\t", $call);
		my @info= split (";", $tmp[3]);
		my ($svtype) = grep ($_=~/SVTYPE=/, @info);
		$svtype =~s/SVTYPE=//;

		my $length = $tmp[2]-$tmp[1];
		next if $length < 1000;
		print OUT "$tmp[0]\t$tmp[1]\t$tmp[2]\tIMPRECISE\t$svtype\t$mapq\t$mean_ratio\t$MAD\n"; 
	}
	close OUT;

	unlink ("$outputDir/before_filtering.tmp.txt");
	unlink $infile;
}

# Adding dummy coordinates at genomic gaps due to incorrect execution of Hmmseg
###############################################################################
sub fakeSpaces {

	my $infile = shift;
	my $outfile = $infile;

	$outfile =~s/.tmp//;

	my $flag = 0;
	my @tmp_A;	
	open (IN, "<", $infile) || die " ERROR: Unable to open $infile\n";
	open (OUT, ">", $outfile) || die " ERROR: Unable to open $outfile\n";
	while (my $line=<IN>) {
		chomp $line;

		if ($flag == 0) {
			@tmp_A = split (/\t/, $line);
			print OUT "$line\n";
			$flag = 1;
			next;
		}
		elsif ($flag == 1) {
			my @tmp_B = split (/\t/, $line);
			my $difference = $tmp_B[1]-$tmp_A[2];
			if ($tmp_A[0] ne $tmp_B[0]) {
				print OUT "$line\n";
			}
			if ($tmp_A[0] eq $tmp_B[0] && $difference <= $binSize ) {
				print OUT "$line\n";
			}
			if ($tmp_A[0] eq $tmp_B[0] && $difference > $binSize ) {
				for (my $i = $tmp_A[2]+1; $i<=$tmp_B[2]; $i=$i+$binSize+1) {
					my $start = $i;
					my $end = $i +$binSize;
					print OUT "$tmp_A[0]\t$start\t$end\t50\t1\n";
				}
			}
			@tmp_A = split (/\t/, $line);
		}
	}
	close IN;
	close OUT;
}


####################################
sub segmentRatios {

	my $input = shift;
	my $sampleName = shift;
	foreach my $chromosome (@chrs) {
		my $cmd = "$grep -P '^$chromosome\t' $input | $bedtools intersect -a stdin -b $centromeres -v | $bedtools intersect -a stdin -b $blacklist -v > $outputDir/$chromosome.tosegment.tmp.txt";
		system $cmd;

		fakeSpaces("$outputDir/$chromosome.tosegment.tmp.txt");
		unlink ("$outputDir/$chromosome.tosegment.tmp.txt");

		if (!-z "$outputDir/$chromosome.tosegment.txt") {
			
		    print " INFO: segmenting $chromosome\n";			
			my $dir_input = dirname ("$outputDir/$chromosome.tosegment.txt");

			`echo $chromosome.tosegment.txt > $outputDir/regions.list`;

			#$cmd = "$java -jar $hmmseg --model $outputDir/model.hmm --input-bed $outputDir/regions.list";

			$cmd = "$java -jar $hmmseg --model $outputDir/model.hmm --input-bed $outputDir/regions.list >/dev/null 2>&1";
			system $cmd;

			if (!-e "$outputDir/$chromosome.tosegment.txt.wig") {
				unlink ("$outputDir/$chromosome.tosegment.txt");
				print " INFO: skipping $chromosome\n";
				next;
			}

			# Creating human-readable output
			open (IN, "<", "$outputDir/$chromosome.tosegment.txt.wig") || die "Unable to open $outputDir/$chromosome.tosegment.txt.wig\n";
			open (TMP, ">", "$outputDir/tmp.CNV.bed");
			open (CHR, ">", "$outputDir/$chromosome.plot.txt") || die "Unable to open $outputDir/$chromosome.plot.txt\n";

			my $flag = 0;
			while (my $line=<IN>) {
				chomp $line;
				if ($line =~/viterbi_segmentation/) {
					$flag = 1;
					next;
				}
				if ($flag == 1) {
					my ($chr, $start, $end, $state) = split (/\t/, $line);
					my $type;
					my $cn;
					if ($state == 0 ) {
						$cn = 0;
						$type = "DEL";
					}
					elsif ($state == 1) {
						$cn = 1;
						$type = "DEL";
					}	
					elsif ($state == 2) {
						$type = "normal";
						$cn = 2;
					}
					elsif ($state == 3) {
						$cn = 3;
						$type = "DUP";
					}
					elsif ($state == 4) {
						$cn = 4;
						$type = "DUP";
					}
					my $size = $end-$start;
					next if $size == $binSize;

					if ($type eq 'DEL' || $type eq 'DUP') {
						print TMP "$chr\t$start\t$end\tIMPRECISE;SVTYPE=$type;SVLEN=$size;CN=$cn\n";
					}
					my $copy_ratio = $cn/2;
					print CHR "$chr\t$start\t$end\t$copy_ratio\t$cn\n";
				}
			}
			`$bedtools intersect -a $outputDir/tmp.CNV.bed -b $centromeres -v >> $outfile`;

			unlink ("$outputDir/$chromosome.tosegment.tmp.txt");
			unlink ("$outputDir/$chromosome.tosegment.txt.wig");
			unlink ("$outputDir/tmp.CNV.bed");

			# Plotting whole-chromosome ratios and segmentation
			#plotChromosomeCNVs( "$outputDir/$chromosome.tosegment.txt", "$outputDir/$chromosome.plot.txt", $chromosome, $sampleName );
			unlink ("$outputDir/$chromosome.tosegment.txt");
			unlink ("$outputDir/$chromosome.plot.txt");
		}
		else {
			print " INFO: skipping $chromosome\n";
			unlink ("$outputDir/$chromosome.tosegment.txt");
			unlink ("$outputDir/$chromosome.tosegment.tmp.txt");
		}
		unlink ("$outputDir/$chromosome.tosegment.txt");
		unlink ("$outputDir/$chromosome.tosegment.tmp.txt");
	}
}

###################
sub Median {
 my @data = @_;
 my $stat = Statistics::Descriptive::Full->new();
 $stat -> add_data(@data);
 my $mean = $stat->median();
 return sprintf "%.3f",$mean;
}

##########################
sub MAD {
	# Calculate Median absolute Deviation

	my @array = @_;

	my $median = Median(@array);
	my @substractMedian = ();

	foreach my $value (@array) {
		my $x = abs($value-$median);
		push @substractMedian, $x;
	}

	my $MAD = abs( Median (@substractMedian) );
	return sprintf "%.2f",$MAD;
}
##########################
sub plotChromosomeCNVs {

	my $a   = shift;
	my $b   = shift;
	my $chr = shift;
	my $name= shift;

	`$bedtools intersect -a $a -b $b -wa -wb | awk '{print \$1"\t"\$2"\t"\$3"\t"\$4"\t"\$5"\t"\$9"\t"\$10}' > $outputDir/toplot.txt`;

	open (PLOT, ">", "$outputDir/toplot.R") || die "Unable to open $outputDir/toplot.txt\n";
 	print PLOT "library(ggplot2)\n";
	print PLOT "format_by_2 <- function(x){ x/2 }\n";
	print PLOT "png(\"$outputDir/$name\_$chr.png\", height = 1000, width =3000, res=350)\n";
	print PLOT "mydata<-read.table(file=\"$outputDir/toplot.txt\", sep=\"\t\", header=FALSE) \n";
	print PLOT "attach(mydata)\n";
	print PLOT "megabase<-V2/1000000\n";
	print PLOT "logSeg <- V7\n";
	print PLOT "myplot <- ggplot(mydata, aes(megabase, V5, colour=factor(V7))) + geom_point(aes(fill=factor(V7), shape = '.', alpha = 0.5 ), size = 0.1, alpha = 0.2)+ scale_colour_manual(values=c(\"#FF4C4C\", \"#FF4C4C\", \"#7E7E7E\", \"#007F00\")) + geom_point(aes(y=V6),shape =\".\", size=0.5, colour=\"red\") + theme_bw() + theme(panel.border = element_rect(colour = \"black\", size = 1, fill=NA)) + xlab(\"Genomic Position (mb)\") + theme(axis.text.x = element_text(size = 12, hjust = 1)) + theme(axis.text.y = element_text(size = 12, hjust = 1)) + ylab(\"Copy ratio\") + theme(axis.title = element_text(family = \"helvetica\", color=\"black\", face=\"bold\", size=14)) + ggtitle(\" $name - $chr copy number profile\") + theme(plot.title = element_text(\"helvetica\", face=\"bold\", \"black\", size=15)) + theme(legend.title = element_text(colour=\"black\", size=14, face=\"bold\")) + theme(legend.position=\"none\")\n";
	print PLOT "myplot + ylim(0,4)\n";
	print PLOT "dev.off()\n";
	`Rscript $outputDir/toplot.R >/dev/null 2>&1`;

	unlink ("$outputDir/toplot.txt", "$outputDir/toplot.R");
}

##########################
sub createModel {

	open (MODEL, ">", "$outputDir/model.hmm") || die "Unable to open $outputDir/model.hmm\n";
	print MODEL "HMMSEG 3
State Start
Emissions none
Transitions 0.0025 0.0025 0.99 0.0025 0.0025
State 0
Emissions  0.1 0.0001
Transitions 0.99 0.0025 0.0025 0.0025 0.0025
State 1
Emissions  0.5 0.02
Transitions 0.0025 0.99 0.0025 0.0025 0.0025
State 2
Emissions  1.0 0.03
Transitions 0.0025 0.0025 0.99 0.0025 0.0025
State 3
Emissions  1.5 0.04
Transitions 0.0025 0.0025 0.0025 0.99 0.0025
State 4
Emissions  2 0.06
Transitions 0.0025 0.0025 0.0025 0.0025 0.99\n";
close MODEL;

}

###################
sub getMean {
 my @data = @_;
 my $stat = Statistics::Descriptive::Full->new();
 $stat -> add_data(@data);
 my $mean = $stat->mean();
 return sprintf "%.3f",$mean;
}

###################
sub getSNR {

 my @array =@_;
 my $mean = getMean(@array);
 my $sd =  getStd(@array);

 my $signal2noise = 0;
 if ($sd > 0 ) {
 	$signal2noise = sprintf "%.3f", $mean/$sd;
 }
 return ($signal2noise);
}

###################
sub getStd {

 my @array = @_;
 my $stat = Statistics::Descriptive::Full->new();
 $stat -> add_data(@array);
 my $sd = $stat->standard_deviation();
 return $sd;
}