#!/usr/bin/env perl

package BuildReference;

use strict;
use Getopt::Long;
use File::Basename; 
use List::MoreUtils qw(uniq);
use Sort::Key::Natural qw(natsort);
use Statistics::Descriptive;

########################
 sub clusterBatches {

  my $outputDir = shift;
  my %nodes = ();
  my %seen  = ();
  
  open (IN, "<", $::HoF{CORRELATIONS_ON}) || die " ERROR: Cannot open $::HoF{CORRELATIONS_ON}\n";
  my $count = 0;
  my @samples = ();
  while (my $line =<IN>) {

	chomp $line;
	$line=~s/"//g;
	my @tmp = split (/\t/, $line);
	if ($count == 0) {
		@samples = @tmp;
		$count++;
		next;
	}
	my $j = 0;

	for (my $i = 1; $i < @tmp; $i++) {
		# skipping seen nodes

		# skipping self correlation
		if ($tmp[0] eq $samples[$j]) {
			$j++;
			next;
		}
		# Creating nodes
	 	$nodes{$tmp[0]}{$samples[$j]} = $tmp[$i];
	 	$nodes{$samples[$j]}{$tmp[0]} = $tmp[$i];
		$j++;
	}
	$seen{$tmp[0]}++;
	$count++;
 }
 close IN;

 open (LOG, ">", $::HoF{REFERENCES}) || die " ERROR: Cannot open $::HoF{REFERENCES}\n";
 print LOG "N\tSAMPLE\tREFERENCE_SET\tMEAN_CORRELATION\n";
 my @cluster = ();
 $count = 0;
 my $sum = 0;
 my $mean_corr;
 my $sampleCount = 0;
 foreach my $sample1 (@samples) {
	$count++;
	my %sampleCorr = ();
	foreach my $sample2 (@samples) {
		next if $sample1 eq $sample2;
		if ($nodes{$sample1}{$sample2} >= $::minCorrelation) {
			my $corr = $nodes{$sample1}{$sample2};
			$sampleCorr{$corr} = $sample2;
		}
	}
	# Selecting the 15th most well correlated samples
	foreach my $corr (reverse sort keys %sampleCorr) {

		$sampleCount++;

		if ($::input) {
			last if $sampleCount > $::maxSampleSizeCluster;
		}
		push @{$::sampleHash{$sample1}{REFERENCE}}, $sampleCorr{$corr};
		$sum+=$corr;
		my $corrdec = sprintf "%.3f", ($corr); 
		push @cluster, "$sampleCorr{$corr}($corrdec)";
	}
	$sampleCount = 0;
	my $reference;
	if (@cluster == 0) {
		$mean_corr = 1;	
		$reference = "none";
		push @{$::sampleHash{$sample1}{REFERENCE}}, $reference;
	}
	else {
		my $n = scalar @cluster > 0 ? scalar@cluster : 1;
		$mean_corr = sprintf "%.3f", ($sum/$n);
		$reference = join (",", @cluster);
	}

	$::sampleHash{$sample1}{CORRELATION} = $mean_corr;

	print LOG "$count\t$sample1\t" . join ("\t", @cluster) . "\t" . $mean_corr . "\n";
	@cluster = ();
	$sum = 0;
 }
}

######################## 
#Subroutines for runREF
 sub getSampleIdxFromRatioHeader {

	my $ratioFile = shift;
	my $strHeaderRatios = $ratioFile =~/.gz/ ? `$::zcat $ratioFile | $::head -1`
	: `$::cat $ratioFile | $::head -1`;
	chomp $strHeaderRatios;
	my @tmpHeaderRatios = split (/\t/, $strHeaderRatios );

	my %HashPosAll = ();
	for (my $i = 6; $i < @tmpHeaderRatios; $i++) {
		$HashPosAll{$tmpHeaderRatios[$i]} = $i;
	} 
	return %HashPosAll;
 }

########################
 sub getCorrAndPlot {

	my $outputDir = shift;

	my $libraries = "library(ggplot2)\nlibrary(RColorBrewer)\nlibrary(corrplot)\nlibrary(gplots)\nlibrary(grid)\nlibrary(gridExtra)\nlibrary(reshape2)\n";
	my $myFile = qq("$outputDir/$::outName.merged.bed\");
	open (R, ">", "$outputDir/plotCorrelationMatrix.R") || die " ERROR: Cannot open plotCorrelationMatrix.R\n";
    print R "$libraries\n";
	print R "mydata<-read.table(file=$myFile ,sep ='\t', check.names = FALSE, header=TRUE)\n";
	print R "attach(mydata)\n";
	#print R "names(mydata) <- sname<-gsub(\"_gc_corrected\",\"\",names(mydata))\n";
	print R "corrected_seq <- as.matrix(mydata[seq (7, ncol(mydata))])\n";

	# Calculating Pearson correlation
	print R "cor_corrected<-cor(corrected_seq, method=\"pearson\")\n";
	print R "scale_cor_corrected <- scale(cor_corrected)\n";
	my $greyStr;
	for (my $i = 8; $i <=87; $i++) {
		$greyStr .= qq(\"grey$i\",);
	}
	print R "melt_corrected <- melt(cor_corrected)\n";
	print R "par(mar=c(7,4,4,2)+0.1) \n";
	print R "Height <- ncol(cor_corrected)*34\n";
	print R "png(\"$outputDir/heatmap_corrected_depth.png\", res = 300, height=Height, width=Height)\n";

	# Minimum correlation value will be the lower bound of the bar legend
	print R "minlegend1 <- min(cor_corrected)\n";
	print R "melt_corrected <- melt(cor_corrected)\n";
	print R "minlegend2 <- min(cor_corrected)\n";
	print R "my_palette <- colorRampPalette(c(\"darkblue\", \"#196aff\",\"#6f6fff\", \"#327aff\", \"white\", \"#ffdfc0\", \"#ffccb3\", \"#CC483A\", \"#790000\"))(n = 150)\n";

	if ($::plotClusters) {
		# Creating correlation matrix
		print R " if (min(cor_corrected) > 0.90) { heatmap.2(cor_corrected, col=my_palette, trace=\"none\", sepcolor=\"black\",symm=TRUE, margins =c(9,9), main=\"Corrected read counts\",cexRow=0.001, cexCol=0.001, breaks = seq(0.9, 1, length.out = 151)) } else { heatmap.2(cor_corrected, col=my_palette, sepcolor=\"black\", trace=\"none\",margins =c(12,9),breaks = seq(min(cor_corrected), 1, length.out = 151)) }\n";
	}
	print R "dev.off()\n";
	print R "write.table(cor_corrected, \"$outputDir/$::outName.cor_corrected.txt\", sep=\"\\t\")\n";
	close R;
	`$::Rscript $outputDir/plotCorrelationMatrix.R`;

	#`$::Rscript $outputDir/plotCorrelationMatrix.R $::devNull`;
	#unlink("$outputDir/bias_info.tmp.txt");
 }

########################
 sub joinNormalizedCounts {

	my @normFiles = glob ("$::outDir/ON_TARGET/*NormalizedCounts.bed.gz");
	my $refDB = "$::outDir/ON_TARGET/References.sqlite.txt.gz";

	# creating filtered regions file without header
	my $noheaderROI = "filtered_regions.bed";
	$noheaderROI =~s/.bed/.noheader.bed/;
	my $cmd = "$::tail -n +2 $::outDir/ON_TARGET/filtered_regions.txt > $::outDir/ON_TARGET/$noheaderROI";
	system $cmd;	

    # Remove filtered regions from initial bedfile
	my $filteredROI = basename($::bed);
	$filteredROI =~s/.bed/.filtered.bed/;
	$cmd = "$::bedtools intersect -a $::bed -b $::outDir/ON_TARGET/$noheaderROI -v > $::outDir/ON_TARGET/$filteredROI";
	system $cmd;
	
	# Goal here is to merge all normalized counts files into a single master file
	# 1. Check files that have the same coordinates than the input bed file
	my @filesToMerge = ();

	foreach my $file (@normFiles) {
		push @filesToMerge, $file;
	}
	push @filesToMerge, $refDB if -e $refDB;

	# Here we merge all files
	my $n = 1;
	my $headerFirst;
	my $masterFile = "$::outDir/ON_TARGET/$::outName.merged.tmp.bed";
	open (OUT, ">", $masterFile);

	foreach my $file ( @filesToMerge) {

		my $name = basename($file);

		if ($n == 1) {
			$headerFirst = `$::zcat $file | $::head -1`; 
			chomp $headerFirst;
			print OUT "$headerFirst\n";
			close OUT;

			my $cmd = "$::zcat $file | $::tail -n +2 > $file.noheader.txt";
			system $cmd;	

			#`$::cut -f 1,2,3,4,5,6 --complement $file > $::outDir/$name.tmp.txt`;
			`$::cat $masterFile $file.noheader.txt > $::outDir/ON_TARGET/$name.tmp.master.txt`;
			$masterFile = "$::outDir/ON_TARGET/$name.tmp.master.txt"; 
			unlink ("$file.noheader.txt");
		}
		else {
			`$::zcat $file | $::cut -f 1,2,3,4,5,6 --complement > $::outDir/ON_TARGET/$name.tmp.txt`;
			`$::paste $masterFile $::outDir/ON_TARGET/$name.tmp.txt > $::outDir/ON_TARGET/$name.tmp.master.txt`;
			$masterFile = "$::outDir/ON_TARGET/$name.tmp.master.txt"; 
			unlink ("$file.noheader.txt");
		}
		$n++;
	}
	rename $masterFile, "$::outDir/ON_TARGET/$::outName.merged.bed";
	unlink( glob ("$::outDir/ON_TARGET/*tmp.master.txt"));
	unlink( glob ("$::outDir/ON_TARGET/*tmp.txt"));
	unlink( glob ("$::outDir/ON_TARGET/*tmp.bed"));

	# Perform correlation matrix
	getCorrAndPlot("$::outDir/ON_TARGET");

	# Get all references
	createReferences();
 }

########################
 sub createReferences {

   # Goal: group samples by high pairwise correlation
   # Goal: Dump all possible references up to a mean correlation greater than 0.95 or defined

  my %nodes = ();
  my %seen  = ();  
  my $count = 0;
  my @samples = ();
  open (IN, "<", "$::HoF{CORRELATIONS_ON}") || die " ERROR: Cannot open $::HoF{CORRELATIONS_ON}\n";

  while (my $line =<IN>) {
	chomp $line;

	$line=~s/"//g;
	my @tmp = split (/\t/, $line);

	# Sample header
	if ($count == 0) {
		@samples = @tmp;
		$count++;
		next;
	}
	my $j = 0;
	for (my $i = 1; $i < @tmp; $i++) {

	 	$nodes{$tmp[0]}{$samples[$j]} = $tmp[$i];
	 	$nodes{$samples[$j]}{$tmp[0]} = $tmp[$i];
		$j++;
	}
	$seen{$tmp[0]}++;
	$count++;
 }
 close IN;

 my %visited;
 my $refN = 0;
 my %Refs = ();
 foreach my $sample1 (@samples) {
	#print "$sample1\n";

	if (exists $visited{$sample1}) {
		next if $visited{$sample1} > 0;
	}

    if ($::doCaseControl) {
        next if $::sampleHash{$sample1}{CONTROL};
    }

	# Clustering
	foreach my $sample2 (@samples) {
        if ($::doCaseControl) {
            next if $::sampleHash{$sample2}{CASE};
        }
		if ($nodes{$sample1}{$sample2} >= $::minCorrelation) {
			my $corr = $nodes{$sample1}{$sample2};
			#print "$refN => $sample2\n";
			push @{ $Refs{$refN} }, $sample2;
			$visited{$sample2}++;
		}
	}
	$refN++;
 }


 if ($::doPooled) {

	my %hashOfSample = getSampleIdxFromRatioHeader( $::HoF{MERGED_NORM_COV} );

	foreach my $refN (natsort keys %Refs) {

		# Here we will construct every reference based on the samples with high correlation
		open (REF, ">", "$::outDir/ON_TARGET/$::outName.$refN.ref");
		open (MERGED , "<", $::HoF{MERGED_NORM_COV} ) || die " ERROR: Unable to open $::HoF{MERGED_NORM_COV}\n";
		my @idxs = ();
		foreach my $sample (@{$Refs{$refN} } ) {
			push @idxs, $hashOfSample{$sample};
		}
		my $n = 0;
		while (my $line =<MERGED>) {
			chomp $line;
			my @tmp = split (/\t/, $line);
			print REF join("\t", @tmp[0..5]);
			if ($n == 0) {
				print REF "\tref_$refN" . "\n";
			}
			else {
				my @normArr = ();
				foreach my $idx (@idxs) {
					push @normArr, $tmp[$idx];
				}
				my $medianNormCov = scalar @normArr > 1 ? medianArray(@normArr) : $normArr[0];
				print REF "\t$medianNormCov\n";
			}
			$n++;
		}
		close MERGED;
		close REF;
	}
 }
}

#####################
sub medianArray {
 my @data = @_;
 my $stat = Statistics::Descriptive::Full->new();
 $stat -> add_data(@data);
 my $median = $stat->median();
 return sprintf "%.2f", $median;
}
1;
