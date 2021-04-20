#!/usr/bin/env perl

package Plot;

use strict;
use Getopt::Long;
use File::Basename; 
use WES::BuildReference;
use List::MoreUtils qw(uniq);
use Chart::Gnuplot;
 
########################
sub plotSingleExon2 {

   my $outputDir = shift;
   my $inputData = shift;
   my $exonName  = shift;
   my $sampleName= shift;

   # Data
   my @x;
   my @y;
   my @z;
   my %HoCase = ();
   my %HoControls = ();

   my $pos = 0;
   open (IN, "<", $inputData) || die " ERROR: Unable to open $inputData\n";
   while (my $line=<IN>) {
	    chomp $line;
		my ($chr, $start, $end, $exon, $sample, $ratio, $group) = split (/\t/, $line);
		push @x, $start;
		push @y, $ratio;
		if ($group eq $sampleName) {
			push( @{$HoCase{$sample} }, $ratio);
   		}
		else {
			push( @{$HoControls{$sample} }, $ratio);
   		}
		$pos++;
	}

	# Defining X ticks to be shown
	my $minx = $x[0];
	my $midx = $x[int(scalar@x/2)];
	my $maxx = $x[-1];

	# Create chart object and specify the properties of the chart
	my $chart = Chart::Gnuplot->new(
		terminal => 'pngcairo size 1000, 700 noenhanced font \'Verdana,11\'',
		output => "$outputDir/$exonName.png",
		yrange => [0, 2],
   		grid   => "on",
		title  => {
			text => "$exonName",
			font => "arial, 25",
		},
 		bmargin => 5,
 		tmargin => 5,	
 		lmargin => 5,	
 		rmargin => 5,	
		xlabel  => {
			text => "Coordinates (Mb)",
			font => "arial, 20",
		},
		ylabel  => {
			text => "Ratio",
			font => "arial, 20",
		},
		xtics  => {
			labels => [$minx, $midx, $maxx],
			font  => "arial, 15",
			mirror => 'off',
		},
		ytics  => {
			font  => "arial, 15",
			mirror => 'off',
		},
		border => {
			sides    => "bottom, left",
			width    => 2,
		},
		samples => 10000
    );

	# Add a line
	$chart->line(
		from  => "$minx, $::lowerDupCutoff",
		to    => "$maxx, $::lowerDupCutoff",
		width => 1,
    	color => "blue",
    	linetype => "longdash",
	);
	# Add a line
	$chart->line(
		from  => "$minx, $::upperDelCutoff",
		to    => "$maxx, $::upperDelCutoff",
		width => 1,
    	color => "red",
    	linetype => "longdash",
	);

	@x = uniq(@x);
	my @datasets;
	# Create dataset object and specify the properties of the dataset
	push @datasets, Chart::Gnuplot::DataSet->new(
		xdata => \@x,
		ydata => \@{$HoCase{$sampleName}},
		title => $sampleName,
		style => "lines",
   		width => 5,
		color => "#cc0000"
	);

	my %hash = ();
	my$i = 0;
	foreach my $control ( keys %HoControls) {
		push @datasets, Chart::Gnuplot::DataSet->new(
			xdata => \@x,
			ydata => \@{$HoControls{$control}},
			title => $control,
			style => "lines",
		    width => 5,
			color => "grey"
		);
    	$chart->plot2d(@datasets);
	}

}

########################
 sub plotSingleExon {

   my $outputDir = shift;
   my $inputData = shift;
   my $exonName  = shift;
   my $sampleName= shift;
   #my $libraries = "library(ggplot2)\nlibrary(ggsignif)\nlibrary(grid)\nlibrary(gridExtra)\n";
   my $panel_border = "theme(panel.border = element_rect(colour = \"black\", fill=NA, size=1))";

   my $colourLine = "";
   my @order = ("control", $sampleName);
   @order = sort @order;
   if ($order[0] eq "control") {
	$colourLine = "scale_fill_manual(values=c(\"darkgrey\", \"#ffb2b2\")) +  scale_colour_manual(values=c(\"darkgrey\", \"red\"))";
   }
   else {
	$colourLine = "scale_fill_manual(values=c(\"#ffb2b2\",\"darkgrey\")) + scale_colour_manual(values=c(\"red\", \"darkgrey\"))";
   }
   my $colourgroup =  "group.colors <- c(control = \"darkgrey\", '$sampleName' = \"#ffb2b2\")\n";

   open (R, ">", "$outputDir/plotSingleExon.R") || die " ERROR: Cannot open $outputDir/plotSingleExon.R\n";
   #print R "$libraries\n";
   print R "mydata<-read.table(file = \"$inputData\", sep =\"\t\", check.names = FALSE, header = FALSE)\n";
   print R "attach(mydata)\n";
   print R "group.colors <- c(control = \"darkgrey\", 'sampleName' = \"#ffb2b2\")\n";
   print R "samplename <- paste(\"$sampleName.$exonName\", \".png\", sep=\"\")\n";
   print R "filename <- paste(\"$outputDir/\", samplename, sep=\"\")\n";
   print R "png(filename, res =237, width = 1400, height=1400)\n";
   print R "Class<-V7\n";
   print R "plot(x = V2, y = V6, type=\"l\", bg = colors[ unclass(V5) ], cex = 3, pch=21)\n";

   #print R "myplot<-ggplot(mydata, aes(x=V2, y=V6, group=V5)) + geom_line(aes(colour=Class))";
   #print R "+ scale_colour_manual(values=group.colors) + theme_bw() + ylim(0,2)";
   #print R "+ geom_hline(yintercept = 0.6, colour=\"red\", linetype=\"dotted\", size=0.3)";
   #print R "+ geom_hline(yintercept = 1.4, colour=\"dark green\", linetype=\"dotted\", size=0.3)";
   #print R "+ theme(plot.title = element_text(hjust = 0.5))+ xlab(\"genomic position\") + ylab(\"Ratio\")";
   #print R "+ $panel_border + ggtitle(\"$exonName\") + theme(plot.title = element_text(hjust = 0.5)) + theme(legend.title=element_blank())";
   #print R "\n";
   #print R "density_plot <- ggplot(mydata, aes(x=V7, y=V6, fill=Class)) + geom_boxplot(colour=\"black\")";
   #print R "+ geom_signif(comparisons = list(c(\"control\", \"$sampleName\")), map_signif_level = TRUE)";
   #print R "+ scale_fill_manual(values=group.colors) + xlab(\"Class\") + ylab(\"Ratio\") +  theme_classic()";
   #print R "+ ggtitle(\"Significance test\")";
   #print R "\n";

   #print R "grid.arrange(myplot, density_plot, nrow = 2 )\n";
   print R "dev.off()\n";
   `$::Rscript $outputDir/plotSingleExon.R`;
   unlink ("$outputDir/plotSingleExon.R");
 }

 ########################
 sub plotScatter {

	my $outputDir   = shift;
	my $analysisType= shift;
	my $on_or_off   = shift;

	my @segmented = glob ("$outputDir/SEGMENT_DATA/toplot.segmented.*");
	my $libraries = "library(ggplot2)\nlibrary(plyr)\nlibrary(grid)\nlibrary(gridExtra)\nlibrary(gtools)\n";
	my $shape;
	my $alpha;
	my $size;
	if ($analysisType eq 'exome' && $on_or_off eq 'on-target') {
		$shape = "";
		$size  = "size=0.05,";
		$alpha = "alpha=0.2,";
	}
	elsif ($analysisType eq 'exome' && $on_or_off eq 'off-target') {
		$shape = "";
		$size  = "size=0.3,";
		$alpha = "alpha=0.5,";
	}	
	else {
		$size  = "size=0.4,";
		$shape = "";
		$alpha = "";
	}

	foreach my $file (@segmented) {
		my $samplename = basename($file);
		$samplename =~s/toplot.segmented.//;
		$samplename =~s/.bed//;
		
		if (!-e "$outputDir/$samplename.genomewide.png") {

			print" INFO: Plotting sample $samplename\n";

			open (R, ">", "$outputDir/$samplename.plotScatter.R") || die " ERROR: Cannot open $outputDir/$samplename.plotScatter.R\n";
			print R "$libraries\n";
			print R "mydata<-read.table(file= \"$file\", sep =\"\t\",check.names = FALSE, header=FALSE)\n";
			print R "attach(mydata)\n";
			print R "chrsort <- mixedsort(levels(factor(V1)))\n";
			print R "mydata\$V1<-factor(mydata\$V1, levels=chrsort)\n";
			print R "mydata<-mydata[order(mydata\$V1),]\n";
			print R "xval<-seq(1,length(mydata\$V1))\n";
			print R "mycount<-count(mydata\$V1)\n";
			print R "attach(mycount)\n";
			print R "mysum<-cumsum(freq)\n";
			print R "mysum<-append(mysum,0, after=0)\n";
			print R "mysum<-mysum[-length(mysum)]\n";
			print R "png(\"$outputDir/$samplename.genomewide.png\", res=137, width=1600, height=600)\n";
			print R "myplot<-ggplot(mydata,aes(x=xval, y =V5), group=V1) + geom_point($size $shape $alpha aes(colour=V1))";
			print R "+ xlab(\"Ratio\") + geom_point(aes(y=V6), $shape, size=0.2, colour=\"red\")";
			print R "+ scale_x_continuous(\"chromosome\",breaks = mysum, labels=chrsort) + theme_bw()";
			print R "+ ylim(0,3) + theme(panel.grid.major.x = element_line(colour = \"black\"),panel.grid.minor.y=element_blank(), panel.grid.major.y = element_blank(),  panel.background = element_rect(colour = \"black\", size=1))";
			print R "+ theme(axis.text.x = element_text(size = 12, hjust = 1, angle=45))";
			print R "+ theme(axis.text.y = element_text(size = 12, hjust = 1))";
			print R "+ scale_colour_manual(values=c(\"#4f6b76\", \"#b2bbc0\", \"#4f6b76\", \"#b2bbc0\",\"#4f6b76\", \"#b2bbc0\",\"#4f6b76\", \"#b2bbc0\",\"#4f6b76\", \"#b2bbc0\",\"#4f6b76\", \"#b2bbc0\",\"#4f6b76\", \"#b2bbc0\",\"#4f6b76\", \"#b2bbc0\",\"#4f6b76\", \"#b2bbc0\",\"#4f6b76\", \"#b2bbc0\", \"#4f6b76\", \"#b2bbc0\", \"#4f6b76\", \"#b2bbc0\", \"#4f6b76\", \"#b2bbc0\", \"#4f6b76\", \"#b2bbc0\")) + theme(legend.position=\"none\")";
			print R "+ geom_hline(yintercept = 1.4, colour=\"blue\", linetype=\"dashed\", size=0.3)";
			print R "+ geom_hline(yintercept = 0.6, colour=\"red\", linetype=\"dashed\", size=0.3) +ylab(\"Ratios\")";
			print R "+ ggtitle(\"$samplename\") + theme(plot.title = element_text(hjust = 0.5))\n";
			print R "myplot\n";
			print R "dev.off()\n";

			my $cmd="$::Rscript $outputDir/$samplename.plotScatter.R $::devNull";
			system $cmd;
			#unlink("$outputDir/$samplename.plotScatter.R");
		}
	}
 }

#############################
  sub plotCorrelationMatrix {

	my $outputDir = shift;

	my $countFile = "$outputDir/$::outName.NormalizedCounts.bed.gz";
	Utils::decompressFile($countFile);
	$countFile =~s/.gz//;

	#if (!-e "$outputDir/heatmap_corrected_depth.pdf") {

		my $libraries = "library(ggplot2)\nlibrary(RColorBrewer)\nlibrary(corrplot)\nlibrary(gplots)\nlibrary(grid)\nlibrary(gridExtra)\nlibrary(reshape2)\n";
		#my $myFile = qq("$outputDir/bias_info.txt\");

		open (R, ">", "$outputDir/plotCorrelationMatrix.R") || die " ERROR: Cannot open plotCorrelationMatrix.R\n";
		print R "$libraries\n";
		
		print R "mydata<-read.table(file=\"$countFile\", sep =\'\t\', check.names = FALSE, header=TRUE)\n";
		print R "png(\"$outputDir/heatmap_corrected_depth.png\", res = 300, height=2800, width=2800)\n";
		print R "attach(mydata)\n";
		print R "names(mydata) <- sname<-gsub(\"_gc_corrected\",\"\",names(mydata))\n";

		print R "scaled_seq <- as.matrix(mydata[ seq (7, ncol(mydata))])\n";
		print R "corrected_seq <- as.matrix(mydata[seq (7, ncol(mydata))])\n";

		# Calculating Pearson correlation
		print R "cor_scaled<-cor(scaled_seq, method=\"pearson\")\n";
		print R "cor_corrected<-cor(corrected_seq, method=\"pearson\")\n";

		print R "scale_cor_corrected <- scale(cor_corrected)\n";
		print R "scale_cor_scaled <- scale(cor_scaled)\n";
		my $greyStr;
		for (my $i = 8; $i <=87; $i++) {
			$greyStr .= qq(\"grey$i\",);
		}
		print R "melt_scaled <- melt(cor_scaled)\n";
		print R "par(mar=c(7,4,4,2)+0.1) \n";

		# Minimum correlation value will be the lower bound of the bar legend
		print R "minlegend1 <- min(cor_scaled)\n";
		print R "melt_corrected <- melt(cor_scaled)\n";
		print R "minlegend2 <- min(cor_corrected)\n";
		print R "my_palette <- colorRampPalette(c(\"darkblue\", \"#196aff\",\"#6f6fff\", \"#327aff\", \"white\", \"#ffdfc0\", \"#ffccb3\", \"#CC483A\", \"#790000\"))(n = 150)\n";

		# Creating correlation matrix
		print R " if (min(cor_corrected) > 0.90) { heatmap.2(cor_corrected, col=my_palette, trace=\"none\", symm=TRUE, margins =c(9,9),sepcolor=\"black\", main=\"Corrected read counts\", breaks = seq(0.9, 1, length.out = 151)) } else { heatmap.2(cor_corrected, col=my_palette, sepcolor=\"black\", trace=\"none\",margins =c(12,9),breaks = seq(min(cor_corrected), 1, length.out = 151)) }\n";
	
		print R "dev.off()\n";
		#print R "write.table(cor_scaled, \"$outputDir/cor_scaled.txt\", sep=\"\\t\")\n";
		print R "write.table(cor_corrected, \"$outputDir/$::outName.cor_corrected.txt\", sep=\"\\t\")\n";
		close R;

		`$::Rscript $outputDir/plotCorrelationMatrix.R $::devNull`;
		unlink("$outputDir/bias_info.tmp.txt") if -s "$outputDir/bias_info.tmp.txt";
		unlink("$outputDir/plotCorrelationMatrix.R");
		Utils::compressFile($countFile);

	#}
 }

########################
sub plotShortSegment {

   my $infile   = shift;
   my $inputDir = shift;
   my $sampName = shift;
   my $bed      = shift;
   my $chr      = shift;
   my $start    = shift;
   my $end      = shift;
   my $cnvtype  = shift;

   my $outname = "$sampName.$chr\_$start\_$end.png";
   if (-e $outname )  {
	return 1;
   }
   my $panel_border  = "theme(panel.border = element_rect(colour = \"black\", fill=NA, size=1))";
   my $remove_x_axis = "theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())";
   my $remove_y_axis = "theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())";
   my $remove_legend_title = "theme(legend.title=element_blank())";
   my $tmpExons = `$::bedtools intersect -a $bed -b $infile | $::cut -f 4 | $::uniq | $::grep -v 'info'`;

   # GGbio is creating weird results for big segmented regions, so we decided to use plain ggplot2
   chomp $tmpExons;
   my @uniqExons = split (/\n/, $tmpExons);

   my $append_rect;
   my $append_seg;
   my @append_all = ();
   my @arraySeg   = ();
   my @arrayWrite = ();
   my $counter = 0;
   my $writeExon;
   my $firstAnnotate = "false";

   foreach my $exon (@uniqExons) {
        next if $exon eq 'ND';
	     my $geneA = ( split /;/, $exon)[1];
		my $str = `egrep '$exon' $bed`;
		chomp $str;
		my @tmpStr= split (/\n/, $str);
		my @first = split (/\t/, $tmpStr[0]);
		my @last  = split (/\t/, $tmpStr[0]);
		my ($start, $end);	

		my @tmpGene = split (/[;_]/, $exon);
		my $exonNumber = $tmpGene[3] ? $tmpGene[3] : $counter;

		# Will always be sorted?
		# If we have 5' -> 3'
		if ($first[1] < $last[2]) {	
			$start = $first[1];
			$end   = $last[2];
		}
		else {
			$start = $last[2];	
			$end   = $first[1];
		}
		my $geneB;
		# Appending segment line if previous exon comes from the same gene
	        if ($counter >= 0) {
			my $previous_exon = $uniqExons[$counter-1];
            	$geneB = ( split /;/, $previous_exon)[1];
			if ($geneA eq $geneB) {
				my $str = `egrep '$previous_exon' $bed`;
				chomp $str;
				my @tmpStr = split (/\n/, $str);
				my @first = split (/\t/, $tmpStr[0]);
				my @last  = split (/\t/, $tmpStr[0]);
				my $prev_start = $first[1];
				my $prev_end   = $last[2];	
				my $intron_start = $prev_end+1;
				my $intron_end   = $start-1;
               	$writeExon .= "+ annotate(\"text\",x=($start+$end)/2, size=2, y = 1, label = \"$exonNumber\")";
				push @arrayWrite, "+ annotate(\"text\",x=($start+$end)/2, size=2, y = 0.75, label = \"$exonNumber\")";
				$append_seg .= " + geom_segment(data=mydata, aes(x=$prev_end, y = 1.5, xend = $intron_end , yend = 1.5),fill=\"#00004C\", colour=\"#00004C\", size = 0.1)";
				push @arraySeg,  " + geom_segment(data=mydata, aes(x=$prev_end, y = 1.5, xend = $intron_end , yend = 1.5),fill=\"#00004C\", colour=\"#00004C\", size = 0.1)";
            }
            if ($geneB ne $geneA) {
                $firstAnnotate = "false";
            }	
            if ($exonNumber && $firstAnnotate eq "false") {

				my $tmplabel .= "+ annotate(\"text\",x=($start+$end)/2, size=2, y = 2, label = \"$geneA\")";
				push @arrayWrite, $tmplabel;
                $firstAnnotate = "true";
            }
	    $writeExon = join (" ", uniq (@arrayWrite) );
	    my $unique_seg = join (" ", uniq (@arraySeg) );
	    push @append_all, $append_rect;
	    $append_rect .= "+ geom_rect(data=mydata, aes(xmin = $start, xmax = $end, ymin = 1.3, ymax =1.7), fill=\"#00004C\", colour=\"#00004C\") $writeExon $unique_seg";
        }
        $counter++;
   }

   open (R, ">", "$inputDir/$sampName.$chr.$start.$end.R");
   print R "library(ggplot2)\n";
   print R "library(egg)\n";
   print R "library(grid)\n";
   print R "mydata<-read.table(file=\"$infile\", sep = \"\t\", check.names = FALSE, header =FALSE)\n";
   print R "attach(mydata)\n";
   print R "xmin <-min(V2)\n";
   print R "xmax <- max(V2)\n";

   if ($cnvtype eq 'DUP') {
   	print R "group.colors <- c('2_Controls' = \"darkgrey\", '1_$sampName' = \"#a50f15\")\n";
   }
   else {
   	print R "group.colors <- c('1_Controls' = \"darkgrey\", '2_$sampName' = \"#a50f15\")\n";
   }
   print R "png(\"$inputDir/$outname\", res = 320, width = 2700, height = 1200)\n";
   print R "Ratios<-V5\n";
   print R "myplot<-ggplot(mydata, aes(x=V2, y=V5, fill=V6)) + xlab(\"Coordinates\") + geom_area(position = \"identity\",size =0.15, color= \"black\") + scale_fill_manual(values=group.colors) +  scale_colour_manual(values=group.colors) + theme(panel.border = element_rect(colour = \"black\", fill=NA, size=1))\n";
   print R "trackplot <- ggplot() $append_rect  + ylim(0.5, 2.5)  + ylab(\"Exons\") + ggtitle(\"$chr:$start-$end\") + theme(plot.title = element_text(hjust = 0.5)) + geom_text(size=6, check_overlap =TRUE) +$remove_x_axis +$remove_y_axis + $remove_legend_title\n";


   print R "trackplot <- trackplot + xlim(xmin,xmax) + theme_minimal() + theme(panel.background = element_rect(colour = \"grey\"))+ theme(panel.border = element_rect(colour = \"black\", fill=NA, size=1)) + $remove_x_axis + $remove_y_axis\n";
   print R "myplot<- myplot + theme_bw() + theme(panel.background = element_rect(colour = \"white\"))+ theme(panel.border = element_rect(colour = \"black\", fill=NA, size=1)) + ylab(\"Normalized coverage\") + theme(legend.title=element_blank())\n";
   print R "ggarrange( trackplot, myplot, ncol = 1, heights = c(1, 4))\n";
   print R "dev.off()\n";
   #`$::Rscript $inputDir/$sampName.R`;
   
   `$::Rscript $inputDir/$sampName.$chr.$start.$end.R $::devNull`;

	unlink("$inputDir/$sampName.$chr.$start.$end.R");
}

########################
sub plotLongSegment {

   my $infile   = shift;
   my $inputDir = shift;
   my $sampName = shift;
   my $bed      = shift;
   my $chr      = shift;
   my $start    = shift;
   my $end      = shift;

   my $panel_border  =  "theme(panel.border = element_rect(colour = \"black\", fill=NA, size=1))";
   my $remove_x_axis =  "theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())";
   my $remove_y_axis =  "theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())";
   my $title = "$chr:$start-$end";
   my $remove_legend_title = "theme(legend.title=element_blank())";

   open (R, ">", "$inputDir/$sampName.$chr.$start.$end.R");
   print R "library(ggplot2)\n";
   print R "mydata<-read.table(file=\"$infile\", sep = \"\t\", check.names = FALSE, header =FALSE)\n";
   print R "attach(mydata)\n";
   print R "group.colors <- c(Controls = \"darkgrey\", '$sampName' = \"#ffa500\")\n";
   print R "Ratios<-V6\n";
   print R "max_ratio <-max(V6)\n";
   print R "ZSCORE<-V7\n";
   print R "Coordinates<-V2\n";
   print R "Class<-V8\n";
   print R "myplot<-ggplot(mydata, aes(x=Coordinates, y=Ratios)) + geom_hline(yintercept=0.71, color=\"red\", linetype=7, size = 0.7) + geom_hline(yintercept=1.25, color=\"dark blue\", linetype=7, size = 0.7) +  ggtitle(\"$title\") + geom_line(aes(group=Class), linetype=2) + geom_point(aes(fill=ZSCORE, size=0.01), colour=\"black\", pch=21)  +scale_fill_gradient2()  + ylim(0, max_ratio+0.2) +guides(size = FALSE)\n";
   print R "png(\"$inputDir/$sampName.$chr.$start.$end.png\", res = 250, width = 2000, height = 1000)\n";
   #print R "myplot<-ggplot(mydata, aes(x=Coordinates, y=Ratios)) +ggtitle(\"$title\") + geom_point(aes(shape=Class, fill="red" ,color =\"black\"), size =3)  + scale_colour_gradient2() + ylim(0, max_ratio+0.2)\n";
   print R "myplot + $panel_border\n";
   print R "dev.off()\n";

   `$::Rscript $inputDir/$sampName.$chr.$start.$end.R $::devNull`;
   unlink("$inputDir/$sampName.$chr.$start.$end.R");
}


1;
