#!/usr/bin/env perl
package Segment;
use strict;
use File::Basename; 
use Statistics::Descriptive;
use Parallel::ForkManager;


 sub CBS {

   my $inputDir = shift;
   my $type     = shift;
   my $analysis = shift;
 
   print " INFO: Segmenting ratios (CBS)\n";

   my $j = 0;
   my $i = 6;

   my @ratios;
   foreach my $sample ( sort keys %::sampleHash) {

	   next if exists $::referenceHash{$sample};

	   if ($type eq 'offtarget') {
	   	next if $::sampleHash{$sample}{OFFTARGET_SD_RATIO} > 0.2;
	   }

	   push @ratios, "$inputDir/$sample.ratios.txt.gz" if -s "$inputDir/$sample.ratios.txt.gz";
   }

   foreach my $file (@ratios) {

	 my $pid = $::pm -> start() and next;

	 my $inputTmp = $file;
	 $inputTmp=~s/.ratios.txt.gz/.ratios.tmp.txt/;

	 my $cmd = "$::zcat $file | $::awk \'{ gsub(\"chr\", \"\", \$0); print \$0 }\' > $inputTmp";
	 system $cmd;

	 my $checkNumColumns = `$::head -1 $inputTmp | perl -ne \'\@tmp=split(\"\t\", \$_); print scalar \@tmp;\'`;
	 chomp $checkNumColumns;

	 if ($checkNumColumns < 5) {
		$::pm->finish;
		next;
	 }
	
	 my $sample = basename($file);
	 $sample    =~s/.ratios.txt.gz//;
     my $segmentRscript = "$inputDir/$sample.CBS.R";
	 my $tmpSegment     = "$inputDir/$sample.tosegment_ratios.txt";

	 open (OUT, ">", $tmpSegment);
	 open (IN , "<", $inputTmp ) || die " Unable to open $inputTmp\n";
	 while (my $line =<IN>) {
		chomp $line;
		my @tmp = split (/\t/, $line);

		next if scalar @tmp < 5;

		$tmp[0] = renameChrForCBS($tmp[0]);

		print OUT "$tmp[0]\t$tmp[1]\t$tmp[1]\t$tmp[3]\t$tmp[4]\n";
		print OUT "$tmp[0]\t$tmp[2]\t$tmp[2]\t$tmp[3]\t$tmp[4]\n";
	 }
	 close IN;
	 close OUT;

	 if (!-e "$inputDir/SEGMENT_DATA/segmented.$sample.bed") {

		 open (R, ">", $segmentRscript) || die "Unable to open $segmentRscript\n";	
		 print R "options(scipen = 999)\n";

		 if ($type eq 'ontarget' && $analysis eq 'gene-panel') {
			
			 print R "library(PSCBS)\n";
			 print R "mydata<-read.table(\"$tmpSegment\", sep =\"\t\",header=FALSE)\n";
			 print R "attach(mydata)\n";
			 print R "mydata<-mydata[, c(\"V1\", \"V2\",\"V5\")]\n";
			 print R "names(mydata)<-c(\"chromosome\", \"x\", \"y\")\n";
			 print R "gaps <- findLargeGaps(mydata, minLength = 1e+06)\n";
			 print R "knownSegments <- gapsToSegments(gaps)\n";
			 print R "fit <- segmentByCBS(mydata, avg= \"median\", knownSegments = knownSegments, min.width=2, alpha = 0.000001, seed = 48879, verbose = -10, joinSegments=FALSE)\n";
			 print R "pathname <- writeSegments(fit, name=\"$inputDir/SEGMENT_DATA/PSCBS.$sample.bed\", simplify=TRUE)\n";
		 }
		 # For off-target data or exomes
		 else {
			 print R "library(DNAcopy)\n";
			 print R "cn <- read.table(\"$tmpSegment\", header=F)\n"; 
			 print R "CNA.object <-CNA( genomdat = cn[,5], chrom = cn[,1], maploc = cn[,2], data.type = 'logratio')\n";
			 print R "CNA.smoothed <- smooth.CNA(CNA.object)\n";
			 print R "segs <- segment(CNA.object, verbose=0, min.width=2, alpha = 0.0001 )\n";
			 print R "segs2 = segs\$output\n";
			 print R "write.table(segs2[,2:6], file=\"$inputDir/SEGMENT_DATA/segmented.$sample.tmp.bed\", row.names=F, col.names=F, quote=F, sep=\"\t\")\n";
		 }
		 close (R);

		 $cmd = "$::Rscript $segmentRscript $::devNull";
		 system ($cmd);
		 
		 if ($type eq 'ontarget' && $analysis eq 'gene-panel') {
			
			if (-e "$inputDir/SEGMENT_DATA/PSCBS.$sample.bed.tsv") {
				$cmd = "$::cat $inputDir/SEGMENT_DATA/PSCBS.$sample.bed.tsv | $::grep -v '#' | $::grep -v 'sampleName'| $::grep -v '\tNA'| $::awk '{ print \"chr\"\$2\"\t\"\$3\"\t\"\$4\"\t\"\$5\"\t\"\$6  }'| $::awk '{ if (\$1 ~ /23/) { print \"chrX\"\"\t\"\$2\"\t\"\$3\"\t\"\$4\"\t\"\$5} else { print \$0} }' |   $::awk '{ if (\$1 ~ /^24/) { print \"chrY\"\"\t\"\$2\"\t\"\$3\"\t\"\$4\"\t\"\$5} else { print \$0} }'  > $inputDir/SEGMENT_DATA/segmented.$sample.bed";
				system $cmd;

				$cmd = "$::cat $inputDir/SEGMENT_DATA/PSCBS.$sample.bed.tsv | $::grep -v '#' | $::grep -v 'sampleName'| $::grep -v '\tNA' | awk '{ print \"chr\"\$2\"\t\"\$3\"\t\"\$4\"\t\"\$5\"\t\"\$6  }' | $::awk '{ if (\$1 ~ /^23/) { print \"chrX\"\"\t\"\$2\"\t\"\$3\"\t\"\$4\"\t\"\$5} else { print \$0} }' |   $::awk '{ if (\$1 ~ /^24/) { print \"chrY\"\"\t\"\$2\"\t\"\$3\"\t\"\$4\"\t\"\$5} else { print \$0} }' | $::bedtools intersect -a $file -b stdin -wo | awk '{ print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$4\"\t\"\$5\"\t\"\$11}' > $inputDir/SEGMENT_DATA/toplot.segmented.$sample.bed"; #t
				system $cmd;
			}
		 }
		 else { 
			if (-e "$inputDir/SEGMENT_DATA/segmented.$sample.tmp.bed") {

				# Add chr at the beggining of the file
				$cmd = "$::cat $inputDir/SEGMENT_DATA/segmented.$sample.tmp.bed | perl -ne \'if (\$_=~/^23/) { \$_=~s/23/chrX/; print \$_ }elsif (\$_=~/^24/) { \$_=~s/24/chrY/; print \$_ }else { print \"chr\$_\"} \' > $inputDir/SEGMENT_DATA/segmented.$sample.bed";
				system $cmd;

				# Remove temporary file
				unlink("$inputDir/SEGMENT_DATA/segmented.$sample.tmp.bed");

			 	$cmd = "$::cat $inputDir/SEGMENT_DATA/segmented.$sample.bed | $::bedtools intersect -a $file -b stdin -wo | awk '{ print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$4\"\t\"\$5\"\t\"\$11}' > $inputDir/SEGMENT_DATA/toplot.segmented.$sample.bed";
			}
		 }
		 system($cmd);
		 $j++;
		 $i++;
	 }
	 unlink ($tmpSegment);
	 unlink($segmentRscript);
	 unlink($inputTmp);
	 $::pm->finish;
  }
  $::pm->wait_all_children;

 }

###########################
 sub renameChrForCBS {
	# editting Chr format to
	# compatibilize with DNAcopy

	my $old_chr = shift;
	my $new_chr = $old_chr;
	if ($old_chr eq 'chrX') {
		$new_chr = "23";
	}
	if ($old_chr eq 'X') {
		$new_chr = "23";
	}
	if ($old_chr eq 'chrY') {
		$new_chr = "24";
	}
	if ($old_chr eq 'Y') {
		$new_chr = "24";
	}
	return $new_chr;
 }

return 1;
