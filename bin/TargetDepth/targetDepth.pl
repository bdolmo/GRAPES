#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use Parallel::ForkManager;
use Cwd qw(cwd abs_path);
use Getopt::Long;
our $dirname = dirname (__FILE__);

# Check for TargetDepth binary execution
# Do not move this script from the actual directory
my $targetDepthExe = ( -x "$dirname/TargetDepth" )
? "$dirname/TargetDepth"
: die " ERROR: Cannot execute TargetDepth\n";

# Here we will check if samtools is set on Path to get total read numbers
my $samtools = `which samtools`;
chomp $samtools;

my $cat = `which cat`;
chomp $cat;

my $cut = `which cut`;
chomp $cut;

my $paste = `which paste`;
chomp $paste;

my $sed = `which sed`;
chomp $sed;

my $rm = `which rm`;
chomp $rm;

if (!$samtools) {
	print " WARNING: SAMtools was not set on PATH, which is not mandatory, but required to get total number of reads if needed\n";
}

my $input;
my $outdir;
my $outName;
my $report_counts;
my $report_coverage;
my $bed;
my $version;
my $genome;
my $verbose;
# Default cpu's is 1. Using more than 4 threads does not improve performance due to IO Disk constraints.
my $threads = 1;

# Displaying usage panel
Help () if (@ARGV<1 or !GetOptions(
	'i=s' =>\$input,
	'o=s' =>\$outdir,
	'n=s' =>\$outName,
	'c' =>\$report_counts,
	'd' =>\$report_coverage,
	'b=s' =>\$bed,
	'g=s' =>\$genome,
	't=i' =>\$threads,
	'v' =>\$version,
	'verbose' =>\$verbose,
	)
);

 if ($version) {
	# Following https://en.wikipedia.org/wiki/Software_versioning
	print " targetDepth.pl v1.0\n";
	exit;
 }

 # Now checking input bams
 if (!$input) {
	print " ERROR: No input file was provided. You can do it like this:
			option1: introduce input directory. All bam files will be analyzed
			option2: Introduce comma separated bams(BAM1,BAM2..) 
			option3: Specify bam paths in a text file (one path per line)\n";
	exit;
 }
 if (!$genome) {
	print " ERROR: No genome in FASTA format was introduced\n";
	exit;
 }
 if (!-e $genome) {
	print " ERROR: $genome genome does not exist\n";
	exit;
 }
 if (!$bed) {
	print " ERROR: No BED file introduced\n";
	exit;
 }
 if (!-e $bed) {
	print " ERROR: $bed genome does not exist\n";
	exit;
 }
 
 #Checking BED consistency
 open (BED, "<$bed") || die " ERROR: Unable to open $bed\n";
 while (my $line=<BED>) {
	my @tmp = split (/\t/, $line);
	if (@tmp < 4) {
		print " ERROR: bed file requires a fourth field containing region information (e.g gene name)\n";
		exit;
	}
 }
 close BED;
 
 my $depthOption;
 if ($report_counts && $report_coverage) {
	$depthOption = "-c -d";
 }
 elsif (!$report_counts && $report_coverage) {
	$depthOption = "-d";
 }
 elsif ($report_counts && !$report_coverage) {
	$depthOption = "-c";
 }

 my $date = localtime();

 if ($verbose) {
 	print " CMD: $0 -i $input -g $genome -b $bed -o $outdir $depthOption \n\n";
 }

 # Adding multi BAM input option separated by comma
 my @bams;
 my @tmp = split ("," , $input);
 if (@tmp > 1 ) {
	print " INFO: multiple input samples recognized\n";
	foreach my $sample (@tmp) {
		print " INFO: sample $sample\n";
		push @bams, $sample;
	}
 }
 if ( -d $input ) {
	print " INFO: Analyzing all bam files from $input directory\n";
	@bams = glob "$input/*bam";
 }
 else {
	if (-B $input) {
		print " INFO: Binary Alignment Map $input\n";
	}
	elsif (-f $input) {
		print " INFO: Checking list of BAM files\n";
		open (IN, "<", $input) || die " ERROR: Cannot open $input\n";

		while (my $sample=<IN>) {
			chomp $sample;
			if (!-e $sample || !-B $sample) {
				print " ERROR: $sample is not BAM file\n";
				exit;
			}
			push @bams, $sample;
		}
		close IN;
	}
 }

 # Now do the job
 if (!-e $outdir) {
	mkdir $outdir;
 }
 my $pm = Parallel::ForkManager->new($threads);

 my $index = 0;
 my $count = 0;
 my @batch;
 foreach my $bam (@bams) {

	$count++;
	$index++;
	push @batch, $bam;
	if ($count == $threads || $index == scalar @bams) {
		 foreach my $bam (@batch) {
			print " INFO: Extracting $bam read depth\n";
			my $pid = $pm -> start() and next; 
			if (!-e "$bam.bai" && $samtools) {
				print " INFO: Indexing $bam\n";
				my $cmd = "$samtools index $bam\n";
				system($cmd);
			}
			my $cmd = "$targetDepthExe -i $bam $depthOption -g $genome -b $bed -o $outdir >/dev/null 2>&1";

			#my $cmd = "$targetDepthExe -i $bam $depthOption -g $genome -b $bed -o $outdir";
			if ($verbose) {
				print " INFO: CMD: $cmd\n";
			}
			system($cmd);
			print " INFO: done with $bam\n";
		   	$pm->finish;
		 }
		 $pm->wait_all_children;
		 @batch = ();
		 $count = 0;
	}
 } 
 # Editing summary metrics file
 my $summary    = "$outdir/summary_metrics.log";
 my $newSummary = "$outdir/summary_metrics.tmp.log";

 if (!-e $summary || -z $summary ) {
	 print " ERROR: missing $summary file\n";

 }

 if (-e $summary ) {
	my $capture = `$cat $summary`;
	chomp $capture;
	open (IN, ">", $newSummary) || die " ERROR: Cannot open $newSummary\n";	
	if ($samtools) {
		print IN "SAMPLE\tTOTAL_READS\tREADS_ON_TARGET\tREADS_CHRX\t\%ROI\tMEAN_COVERAGE\tMEAN_COUNTS\tMEAN_ISIZE\tSD_ISIZE\tMEAN_COVERAGE_X\tMEAN_COUNTS_X\n";
	}
	else {
		print IN "SAMPLE\tREADS_ON_TARGET\tREADS_CHRX\tMEAN_COVERAGE\tMEAN_COUNTS\tMEAN_ISIZE\tSD_ISIZE\tMEAN_COVERAGE_X\tMEAN_COUNTS_X\n";
	}
	my @tmpStr = split (/\n/, $capture);
	foreach my $line (@tmpStr) {

		my ($sample, $reads_on_target, $reads_X, $mean_coverage, $mean_counts, $mean_isize, $sd_isize, $mean_coverage_X, $mean_counts_X) = split (/\t/, $line);
		my $outsample = $sample;
		$outsample =~s/.bam//;
		if ($samtools) {
			my @samplePaths = grep ($_=~/$sample/, @bams);

       		my $total_reads = `$samtools idxstats $samplePaths[0] | awk '{i+=\$3} END {print i}'`;
			chomp $total_reads;

			my $roi = sprintf "%.2f", 100* ($reads_on_target/$total_reads);
			print IN "$outsample\t$total_reads\t$reads_on_target\t$reads_X\t$roi\t$mean_coverage\t$mean_counts\t$mean_isize\t$sd_isize\t$mean_coverage_X\t$mean_counts_X\n";
		}
		else {
			print IN "$line\n";
		}
	}
	close IN;
	#`$rm $summary`;
	rename $newSummary, $summary ;
 }
 

 # Merging temporal files 
 my @countFiles    = glob ("$outdir/*_counts.bed");
 my @coverageFiles = glob ("$outdir/*_coverage.bed");
 my @isizeFiles = (); # Now this is deprecated
 #my @isizeFiles    = glob ("$outdir/*_isizes.bed");

 my $masterCounts = "$outdir/$outName.ReadCounts.bed";
 my $masterCovs   = "$outdir/$outName.Coverages.bed";
 #my $masterIsizes = "$outdir/$outName.InsertSizes.bed";

 # Getting coordinates
 
 if($report_counts) {
 	`$cat $countFiles[0] | $cut -f 1,2,3,4,5 > $outdir/coords1.txt`;
 	#`$cat $isizeFiles[0] | $cut -f 1,2,3,4,5 > $outdir/coords2.txt`;
 }
 if($report_coverage) {
 	`$cat $coverageFiles[0] | $cut -f 1,2,3,4,5 > $outdir/coords3.txt`;
 }

 # Inserting fields
 my @tmps;
 my @samp;
 my $str;
 my $header;
 if ($report_counts) {
	 foreach my $ct (@countFiles) {
		my $sample = basename($ct);
		$sample=~s/_counts.bed//;
		$sample=~s/.bam//;
		push @samp, $sample;	
		`$cat $ct | $cut -f 6 > $ct.tmp`;
		push @tmps, "$ct.tmp";
	 }

	 $str = join ("\t", @samp);
	 $header = "chr\tstart\tend\texon\t\%GC\t$str";
	 `$paste $outdir/coords1.txt @tmps > $masterCounts `;
	 `$sed -i '1i\'"$header" $masterCounts`; 
 
	 @tmps = ();
	 @samp = ();
	 $str = "";
	# foreach my $is (@isizeFiles) {
	#	my $sample = basename($is);
	#	$sample=~s/_isizes.bed//;
	#	$sample=~s/.bam//;
	#	push @samp, "$sample\_Mean_Isize";
	#	push @samp, "$sample\_SD_Isize";
	#	`$cat $is | $cut -f 6,7 > $is.tmp`;
	#	push @tmps, "$is.tmp";
	# }
	# $str = join ("\t", @samp);
	# $header = "chr\tstart\tend\texon\t\%GC\t$str";
	# `$paste $outdir/coords2.txt @tmps > $masterIsizes `;
	# `$sed -i '1i\'"$header" $masterIsizes`; 
 }
 if ($report_coverage) {
	 @tmps = ();
	 @samp = ();
	 $str = "";
	 foreach my $cv (@coverageFiles) {
		my $sample = basename($cv);
		$sample=~s/_coverage.bed//;
		$sample=~s/.bam//;
		push @samp, $sample;
		`$cat $cv | $cut -f 6 > $cv.tmp`;
		push @tmps, "$cv.tmp";
	 }
	 $str = join ("\t", @samp);
	 $header = "chr\tstart\tend\texon\t\%GC\t$str";
	 `$paste $outdir/coords3.txt @tmps > $masterCovs `;
	 `$sed -i '1i\'"$header" $masterCovs`;
	 `$rm $outdir/coords3.txt`; 
 	 `$rm $outdir/*_coverage.bed`;
 }

  # Deleting temporal files
 `$rm $outdir/*.tmp`;
 `$rm $outdir/coords1.txt`;
 `$rm $outdir/*_counts.bed`;
 `$rm $outdir/*_isizes.bed`;


sub Help {

 	print("
Description: targetDepth.pl
Purpose: Extract read depth metrics from a list of files
Version 0.1
Usage:  perl targetDepth.pl <OPTIONS>
options:
	-i	STRING	input BAM/s. Choose between:
			option1: introduce input directory. All bam files will be analyzed
			option2: Introduce comma-separated bams (BAM1,BAM2..) 
			option3: Specify bam paths in a text file (one path per line)
	-o	STRING	Output directory. Default is '.'
	-c	BOOL	Report exon counts
	-d	BOOL	Report per base coverage
	-b	STRING	BED regions file
	-g	STRING	Genome file in FASTA format
	-t	INT	Number of threads
	-v	BOOL	Display version
\n");
exit;

}


