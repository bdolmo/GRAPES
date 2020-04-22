#!/usr/bin/env perl


use warnings;
use strict;
use diagnostics;
use threads;
use Data::Dumper;
use PDF::API2;
use PDF::Table;
use File::Basename;

my $dirname = dirname (__FILE__);

################ Process input #1 RUN20170210_147.metrics.log

my $run_name    = $ARGV[0];
my $metricsFile = $ARGV[1];
my $run_dir     = $ARGV[2];
my $pngs        = $ARGV[3];
my $outDir      = $ARGV[4];
my $xipname     = $ARGV[5];

our %samples = ();

if (@ARGV < 4) {
	print "Usage: perl grapes_report <RUN_NAME> <METRICS_FILE.log> <PATHS_VCFs.log> <PATH_PNG> <OUTPUT_DIR> <BED_NAME>\n";
	exit;
}

my $in      = $metricsFile;
my $runname = $in; 

my $outpdf  = "$outDir/$run_name.GRAPES.report.pdf";


################ Process input #2 RUN20170210_147.samplepath.log

my @paths = glob ("$run_dir/FINAL_CALL*");
foreach my $vcf (@paths) { print "$vcf\n"; }

################

my $i=1;
my $numsamples= `wc -l $in`; chomp $numsamples;
$numsamples   = (split /\s+/, $numsamples)[0];

my @nocalls   = ();

my $pdf = PDF::API2->new(-file => "$outpdf"); # Create a blank PDF file
$pdf->mediabox('A4'); # Set the PDF file size
my $font = $pdf->corefont('Helvetica'); # Add a built-in font to the PDF

coverpage($run_name, $numsamples, $xipname, $pdf); # Makes cover page

my $last_y = runpage($in, $run_name, $pdf, $numsamples); # Makes run page

callpages($pdf); # Generates 1 page for each filled report

nocallpages($run_name, $i, $pdf); # Last page summarizes samples with no SV detected

$pdf->save();
$pdf->end;

pagenumbers($pdf,$outpdf); # Enumerate report pages

###############################################################
sub coverpage {
	my ($run_name,$numsamples,$xipname,$pdf) = @_;
	my ($sec,$min,$hour,$day,$month,$year)   = getdate();
	my $coverpage = $pdf->page(); # Add a blank page

my $job = <<"END";
$run_name
END

my $job2 = <<"END";
Structural Variant Report
END

my $inforun = <<"END";
# Total samples: $numsamples
  Gene panel: $xipname
    GRAPES v.0.9.2
     $day/$month/$year
END


my $outNumSamples = <<"END";
Samples analyzed: $numsamples
END

my $genePanel = <<"END";
Gene panel: $xipname
END

my $grapesVersion = <<"END";
GRAPES v0.9.2
END

my $outDate = <<"END";
$day/$month/$year
END


my $advice = <<'END';
Condifentiality Notice
This communication, including attachments, is for the exclusive use of addressees and may contain
confidential, privileged information or communications. If you are not the intended recipient, any
use, copying, disclosure, dissemination or distribution is strictly prohibited. If you are not the
intended recipient, please notify the sender immediately by return e-mail, delete this communication
and destroy all copies.
END

	# Adding RUN name
   	my $content = $coverpage->text();
	$content->translate(300, 720);
	$content->font($pdf->corefont('Helvetica'), 35);
	$content->text_center($job);
	
	# Adding Additional run info
	$content->translate(300, 640);
	$content->font($pdf->corefont('Helvetica'), 25);
	$content->text_center($job2);

	# Adding Additional run info
	$content->translate(300, 580);
	$content->font($pdf->corefont('Helvetica'), 15);
	$content->text_center($outNumSamples);

	# Adding Additional run info
	$content->translate(300, 560);
	$content->font($pdf->corefont('Helvetica'), 15);
	$content->text_center($genePanel);

	# Adding Additional run info
	$content->translate(300, 540);
	$content->font($pdf->corefont('Helvetica'), 15);
	$content->text_center($grapesVersion);

	# Adding Additional run info
	$content->translate(300, 520);
	$content->font($pdf->corefont('Helvetica'), 15);
	$content->text_center($outDate);

	super_textlabel($coverpage, 75,  130,  $pdf->corefont('Helvetica'), 10, $advice, 0,  10);

	my $logogr="$pngs/grapes_new_logo2.png";
	my $image = $pdf->image_png($logogr);
	my $gfx = $coverpage->gfx;
	$gfx->image($image, 75, 156, 0.11);


}

###############################################################
sub svColour {

 my $line = shift;
 my $color;
 # For deletions
 if ($line =~/DEL/) {
	$color = "#FFCCCC";
 }
 if ($line =~/DUP/) {
	$color = "#b2efb2";
 }
 if ($line =~/INV/) {
	$color = "#CCF9FB";
 }
 return $color;
}

###############################################################
sub runpage{

	my $run_file = shift;
	my $runname  = shift;
	my $pdf      = shift;
	my $nSamples = shift;

	my $runpage = $pdf->page(0);


	my $runtext = $runpage->text(); # Add some text to the page
	my $bold = $pdf->corefont('Helvetica-Bold');
	$runtext->font($bold, 14);
	$runtext->translate(300, 780); #(1_540 hor / 1_830 vert)
	$runtext->text_center(" Sequencing summary for $runname");

	#$runtext->font($bold, 12);
	#$runtext->translate(70, 780); #(1_540 hor / 1_830 vert)
	#$runtext->text("Table 1.  Sequencing summary for $runname");
	my $cell_props = [];

	my @metrics = ();
	$metrics[0][0] = "Sample";
	$metrics[0][1] = "Total_Reads";
	$metrics[0][2] = "ROI_Reads";
	$metrics[0][3] = "\%ROI";
	$metrics[0][4] = "Mean Coverage";
	$metrics[0][5] = "Mean Counts";
	$metrics[0][6] = "Mean insert size";
	$metrics[0][7] = "Std insert size";
	$metrics[0][8] = "Correlation";

	my $altitud = 780;

	my $j = 1;
	open (IN, "<", "$in") || die "cannot open file $in";
	while (my $line =<IN>){
		$altitud = $altitud-25;
		if ($j == 1) {
			$j++;
			next;
		}
		chomp $line;
		$line=~s/\t\D+=/\t/g;

		my @tmp = (split /\t+/, $line);

		$samples{$tmp[0]}{SAMPLE} = $tmp[0];

		$samples{$tmp[0]}{CORR}   = `grep -P '^[0-9]+\t$tmp[0]\t' $run_dir/custom_references.info`;

		chomp $samples{$tmp[0]}{CORR};
		$samples{$tmp[0]}{CORR} = ( split /\t/, $samples{$tmp[0]}{CORR})[-1];
		#print "$tmp[0]\t$samples{$tmp[0]}{CORR}\n";
		my @newTmp = ();
		for (my $k = 0; $k < @tmp; $k++) {
			next if ($k == 3  || $k == 9 || $k == 10);
			push @newTmp, $tmp[$k];
		}
		for (my $k = 0; $k < @newTmp+1; $k++) {
			my $background = "white";
			if ($k == 8) {
				$newTmp[8] = $samples{$newTmp[0]}{CORR};
				if ( $tmp[8] < 0.95) {
				       $background = "#FFCCCC";
				}
			}
			my $background_color;
			my $is_even = $j % 2 == 0;
			if ($is_even) {
				$background_color = "white";

			} 
			else {
				$background_color =  "#e4e4e4";

			}
			$cell_props->[$j][$k] = {
				#Row j cell k
				font_color  => 'black',
				justify	=> 'center',
				background_color => $background_color,
			};
		}
		push @metrics, [@newTmp];
		$j++;
	}

	#@metrics = sort @metrics;
	my ($final_page, $number_of_pages, $final_y) = createSummary($runname, $runpage, \@metrics, $cell_props );


	if ($final_y < 470) {

		my $corplot="$run_dir/heatmap_corrected_depth.png";
		my $image = $pdf->image_png($corplot);
		my $gfx = $runpage->gfx;
		$gfx->image($image, 70, $final_y-450, 0.15);
	}
	else {
		my $newpage = $pdf->page(0);
		my $runtext = $newpage->text(); # Add some text to the page
		my $bold = $pdf->corefont('Helvetica-Bold');
		$runtext->font($bold, 14);
		$runtext->translate(300, 780); #(1_540 hor / 1_830 vert)
		$runtext->text_center("Coverage correlation");
		my $corplot="$run_dir/heatmap_corrected_depth.png";
		my $image = $pdf->image_png($corplot);
		my $gfx = $newpage->gfx;
		$gfx->image($image, 70, 325, 0.15);
	}


	($final_page, $number_of_pages, $final_y) = printFilteredROIs($runpage, $final_y);
	return ($final_y);
}





###############################################################
sub printFilteredROIs {

	my $merda = shift;
	my $altitud = shift;

	$altitud = $altitud-25;
	my $cell_props = [];

	my $runpage = $pdf->page(0);

	my $runtext = $runpage->text(); # Add some text to the page
	my $bold = $pdf->corefont('Helvetica-Bold');
	$runtext->font($bold, 12);
	$runtext->translate(300, 780); #(1_540 hor / 1_830 vert)
	$runtext->text_center(" Excluded ROIs due to short size (<10bp) or low Mappability (<50%)");

	my @ROI = ();
	$ROI[0][0] = "CHR";
	$ROI[0][1] = "START";
	$ROI[0][2] = "END";
	$ROI[0][3] = "ROI";
	$ROI[0][4] = "\%GC";
	$ROI[0][5] = "\%MAP";
	$ROI[0][6] = "FILTER";
	$ROI[0][7] = "SIZE";

	my $j = 0;
	open (IN, "<", "$run_dir/filtered_regions.txt") || die "ERROR: Unable to open $run_dir/filtered_regions.txt\n";
	while (my $line=<IN>) {
		chomp $line;
		next if $line =~/^chr\tstart/;
		my @tmp = split (/\t/, $line);
		my @newTmp = ();
		for (my $i = 0; $i < @tmp; $i++) {
			next if $i == 4;
			push @newTmp, $tmp[$i];
		}	
		for (my $k = 0; $k < @newTmp; $k++) {
			my $background_color;
			my $is_even = $j % 2 == 0;
			if ($is_even) {
				$background_color = "white";
			} 
			else {
				$background_color =  "#e4e4e4";
			}
			$cell_props->[$j][$k] = {
				font_color       => 'black',
				justify	         => 'center',
				background_color => $background_color,
			}
		}
		push @ROI, [@newTmp];
		$j++;	
	}

	my $hdr_props =
		{
			font	   => $pdf->corefont("Helvetica-Bold", -encoding => "utf8"),
			font_size  => 8,
			font_color => '#FFFFFF',
			bg_color   => '#686868',
			repeat	   => 1,
			justify    => 'center'
		};

	my $col_props = [
		{
		    min_w => 20, # Minimum column width of 100
		    max_w => 40, # Maximum column width of 150
		    justify => 'center', # Right text alignment
		    font => $pdf->corefont("Helvetica", -encoding => "latin1"),
		    font_size => 8,
		    font_color=> 'black',
		    background_color_odd  => "white",
		    background_color_even => "#e4e4e4",
		},
	    ];
	my $inici = $altitud - 15;
	my $left_edge_of_table = 65;
	my $metricstable = new PDF::Table; # build the table layout


	 my ($final_page, $number_of_pages, $final_y) = $metricstable->table(
		# required params
		$pdf,
		$runpage,
		\@ROI,
		x => $left_edge_of_table,
		w => 460,
		start_y => 765,
		start_h => 350,

		# optional params
		font_size => 8,
		font => $pdf->corefont("Helvetica", -encoding => "latin1"),
		header_props => $hdr_props,
		column_props => $col_props,
		cell_props   => $cell_props,
		next_y  => 750,
		next_h  => 500,
		padding => 5,
		padding_right => 5,
		background_color_odd  => "white",
		background_color_even => "#e4e4e4",
		border_color => 'white',
	);
	return ($final_page, $number_of_pages, $final_y);
}


###############################################################
sub createSummary {

	my $title = shift;
        my $runpage=shift;
	my $array  =shift;
	my $cell_props = shift;
        
	my @metrics = @$array;

	my $hdr_props =  {
		font	   => $pdf->corefont("Helvetica-Bold", -encoding => "utf8"),
		font_size  => 8,
		font_color => '#FFFFFF',
		bg_color   => '#686868',
		repeat	   => 1,
		justify    => 'center'
	};

	my $col_props = [
		{
		    min_w => 40, # Minimum column width of 100
		    max_w => 40, # Maximum column width of 150
		    justify => 'center', # Right text alignment
		    font => $pdf->corefont("Helvetica", -encoding => "latin1"),
		    font_size => 8,
		    font_color=> 'black',
		    #background_color => 'white',
		    background_color_odd  => "white",
		    background_color_even => "#e4e4e4",
		},
	    ];
	my $metricstable = new PDF::Table; # build the table layout
	my $left_edge_of_table = 65;
	my ($final_page, $number_of_pages, $final_y) = $metricstable->table(
		# required params
		$pdf,
		$runpage,
		\@metrics,
		x => $left_edge_of_table,
		w => 460,
		start_y => 765,
		start_h => 350,

		# optional params
		font_size => 8,
		font => $pdf->corefont("Helvetica", -encoding => "latin1"),
		header_props => $hdr_props,
		column_props => $col_props,
		cell_props	=> $cell_props,
		next_y  => 750,
		next_h  => 500,
		padding => 5,
		padding_right => 5,
		background_color_odd  => "white",
		background_color_even => "#e4e4e4",
		border_color => 'white',
		);

	return ($final_page, $number_of_pages, $final_y);
}

###############################################################
sub callpages{
	my ($pdf) = @_;	
        my $img_count = 0;
	foreach my $way (@paths){

		my $nCalls = `wc -l $way`;
		chomp $nCalls;
		$nCalls = (split /\s/, $nCalls)[0];

		my @calls = ();
		$calls[0][0] = "COORDINATES";
		$calls[0][1] = "SIZE";
		$calls[0][2] = "TYPE";
		$calls[0][3] = "GT";
		$calls[0][4] = "GENE";
		$calls[0][5] = "EXON";
		$calls[0][6] = "RATIO";
		$calls[0][7] = "ZSCORE";		
		$calls[0][8] = "\%GC";
		$calls[0][9] = "\%MAP";
		$calls[0][10] = "BREAKS";	
		#$calls[0][8] = "DGV";
		#$calls[0][9] = "ExAC";
		
		#my @calls=(["Coordinate", "Breakpoint", "SVtype / CN", "Genes", "Regions", "grID"]);		
		my $sample= basename($way); 
		$sample=~s/FINAL_CALLS.//;
		$sample=~s/.bed//;

		my $cell_props = [];
		my $j = 1;

		if (!-z $way){ #file is not empty
			my $callpage = $pdf->page(0);
			$i++;
			open (REP, "<", "$way") or die "cannot open file $way";
			while (my $line=<REP>){
				chomp $line;
			
				my @tmp = split /\t/,$line;
				my $loc ="$tmp[0]:$tmp[1]-$tmp[2]";

				# Precise or imprecise
				my $class = ( split /;/ , $tmp[3])[0] eq 'PRECISE' ? 'YES' : 'NO';

				my @info = split (/;/, $tmp[3]);

				my ($size) = grep ($_ =~/SIZE=/, @info);
				$size =~s/SIZE=// if $size;
				if (!$size) {
					$size = $tmp[2]-$tmp[1];
				}
			
				my ($svtype) = grep ($_ =~/SVTYPE=/, @info);
				$svtype =~s/SVTYPE=//;
			
				my ($CN) = grep ($_ =~/^CN=/, @info);
				$CN =~s/CN=// if $CN;
				$CN = '.' if !$CN;

				my ($zscore) = grep ($_ =~/^ZSCORE=/, @info);
				$zscore =~s/ZSCORE=// if $zscore;
				$zscore = '.' if !$zscore;

				my ($breakQual) = grep ($_ =~/^breakQual=/, @info);
				$breakQual =~s/breakQual=// if $breakQual;
				$breakQual = '.' if !$breakQual;

				my $GT = '.';
				if ($CN ne '.') {
					if ($CN == 1 || $CN == 3) {
						$GT = "1/0";
					}
				}

				my ($CNV_RATIO) = grep ($_ =~/RATIO=/, @info);
				$CNV_RATIO =~s/RATIO=// if $CNV_RATIO;

				if ( $CNV_RATIO) {
					$GT = $CNV_RATIO < 0.8 ? '1/0' : '1/1'; 
				}
				my ($RATIO) = grep ($_ =~/^SOFT_RATIO=/, @info);
				$RATIO =~s/SOFT_RATIO=// if $RATIO;
				if ($GT eq '.' && $RATIO) {
					$GT = $RATIO < 0.8 ? '1/0' : '1/1'; 
				}
	
				my ($AB) = grep ($_ =~/^AB=/, @info);
				$AB =~s/AB=// if $AB;
				if ($GT eq '.' && $AB) {
					$GT = $AB < 0.8 ? '1/0' : '1/1';  
				}

				my ($genes)  = grep ($_ =~/GENE=/, @info);	
				$genes =~s/GENE=// if $genes;
				$genes =~s/GENE=// if $genes;
				$genes = '.' if !$genes;
				my $gene_name = ".";
				my $exon      = ".";
				if ($genes) {
					if ($genes =~/5prime/ && $genes =~/3prime/) {
						$genes=~s/3prime//;
						$genes=~s/5prime//;
						my @tmpGenes = split (/[\_,]/, $genes);
						$gene_name = "$tmpGenes[0],$tmpGenes[2]";	
						$exon = "$tmpGenes[1],$tmpGenes[3]";							
					}
					else  {
						my @tmpGene = split (/[\t\_\,]/, $genes);
						$exon = $tmpGene[3];
						$gene_name = $tmpGene[4];
					}
				}		

				my ($GC) = grep ($_ =~/GC=/, @info);
				$GC =~s/GC=// if $GC;
				$GC = '.' if !$GC;

				my ($MAP) = grep ($_ =~/MAP=/, @info);
				$MAP =~s/MAP=// if $MAP;
				$MAP = '.' if !$MAP;

				my ($grID) = grep ($_ =~/grID/, @info);
				$grID =~s/grID=// if $grID;
				$grID = '.' if !$grID;

				my ($DGV) = grep ($_ =~/DGV_FREQ/, @info);
				$DGV =~s/DGV_FREQ=// if $DGV;
				$DGV = '.' if !$DGV;

				my ($BREAKREADS) = grep ($_ =~/BREAKREADS/, @info);
				$BREAKREADS =~s/BREAKREADS=// if $BREAKREADS;
				$BREAKREADS = '.' if !$BREAKREADS;

				my ($ASSEMBLED) = grep ($_ =~/ASSEMBLED/, @info);
				$ASSEMBLED =~s/ASSEMBLED=// if $ASSEMBLED;
				$ASSEMBLED = '.' if !$ASSEMBLED;

				if ($ASSEMBLED && $BREAKREADS && $ASSEMBLED ne '.' && $BREAKREADS ne '.') {
					$ASSEMBLED = "$ASSEMBLED/$BREAKREADS";
				}
		
				my @exacs = grep ($_=~/ExAC/, @info);

				my $ExAC = ".";
				push @calls, [$loc, $size, $svtype, $GT, $gene_name, $exon, $CNV_RATIO, $zscore, $GC, $MAP, $ASSEMBLED];

				# Filling cell specific colours
				for (my $k = 0; $k <12; $k++) {
					my $background = "white";
					if ($GC ne '.') {
						$GC =~s/\%//;
						if ( $GC < 45 || $GC > 60) {
							$background = "#FFFF99";
						}
					}
					if ($MAP ne '.') {
						$MAP =~s/\%//;
						if ( $MAP < 100) {
							$background = "#FFFF99";
						}
					}
					my $colorFont = "black";
					if ($k == 2) {
						if ($svtype eq 'DEL') {
						    #$background = "#FFCCCC";
						    $colorFont = "#e50000";
						}
						if ($svtype eq 'DUP') {
						    #$background = "#b2efb2";
						    $colorFont = "#007f00";
						}
						if ($svtype eq 'INV') {
						    #$background = "#CCF9FB";
						    $colorFont = "#00898e";
						}
						if ($svtype eq 'INV') {
						    #$background = "#CCF9FB";
						    $colorFont = "#00898e";
						}
					}

					$cell_props->[$j][$k] = {
						#Row j cell k
					        background_color_odd  => "white",
					        background_color_even => "#e4e4e4",

						#background_color => $background,
						font_size => "7",
						font_color       => $colorFont,
						justify	=> 'center',
					};
				}

				#my $normimage = $pdf->image_png("$pngs/results.png");
				#my $gfx3 = $callpage->gfx;
				#$gfx3->image($normimage, 0, 750, 0.1);

				#$calltext->font($font, 14);
				#$calltext->translate(30, 820); #(1_540 hor / 1_830 vert)
				#$calltext->text("Results");

				my $calltext = $callpage->text(); # Add some text to the page
				$calltext->font($font, 24);
				$calltext->translate(300, 760); #(1_540 hor / 1_830 vert)
				$calltext->text_center("$sample");

				$calltext->font($font, 14);
				$calltext->translate(300, 720); #(1_540 hor / 1_830 vert)
				$calltext->text_center("Bias normalization");

				$calltext->font($font, 14);
				$calltext->translate(300, 380); #(1_540 hor / 1_830 vert)
				$calltext->text_center("Genome-wide CNV segmentation");


				my ($biases) = glob ("$run_dir/$sample.rawdepth_biases.png");
				if  ( -e $biases) {

					chomp $biases;
					my $page;

					my $normimage = $pdf->image_png($biases);
					my $gfx3 = $callpage->gfx;
					$gfx3->image($normimage, 94, 430, 0.195);
				}


				# Checking if we have an image for the predicted variant
				my ($genome_wide) = glob ("$run_dir/$sample.genomewide.png");
				if ($genome_wide) {

					chomp $genome_wide;
					my $page;

					my $normimage = $pdf->image_png($genome_wide);
					my $gfx3 = $callpage->gfx;
					$gfx3->image($normimage, 97, 220, 0.245);

				}

			   $j++;
			}
	
			my $newpage = $pdf->page(0);
			my $nexttext = $newpage->text(); # Add some text to the page
			$nexttext->font($font, 12);
			$nexttext->translate(80, 760); #(1_540 hor / 1_830 vert)
			$nexttext->text("$sample - Reported variants");

			my $hdr_props =
				{
				font	   => $pdf->corefont("Helvetica", -encoding => "utf8"),
				font_size  => 7,
				justify	=> 'center',
				font_color => '#FFFFFF',
				bg_color   => '#686868',
				repeat	   => 1,
				justify    => 'center'
				};

			my $metricstable = new PDF::Table; # build the table layout

			$metricstable->table(
				# required params
				$pdf,
				$newpage,
				\@calls,
				x => 80,
				w => 430,
				start_y => 750,
				start_h => 350,

				# optional params
				font_size => 7,
				font => $pdf->corefont("Helvetica", -encoding => "latin1"),
				header_props => $hdr_props,
		next_y  => 735,
		next_h  => 430,
				padding => 4,
				padding_right => 4,
				background_color_odd  => "white",
				background_color_even => "#e4e4e4",
				border_color => 'white',
	       			cell_props => $cell_props,
			);
			
		}
		else{
			push(@nocalls,$sample);
		}
	}
}


###############################################################

sub addSmartImg {

 my $count = shift;
 my $nCalls= shift;

 my ($posX, $posY, $res) = (60, 80, 0.19);

 if ($count > 0 ) {
	$posX = 60;
	$posY = 250;
	$res  = 0.19;
 }
 if ($count == 0 && $nCalls <= 2) {
	$posX = 67; # Horizontal
	$posY = 140; # Vertical
	$res  = 0.19;
 }
 if ($nCalls > 2 && $nCalls < 4) {
	$posX = 67; # Horizontal
	$posY = 130; # Vertical
	$res  = 0.19;
 }

 if ($count == 3 && $nCalls <= 5) {
	$posX = 67; # Horizontal
	$posY = 90; # Vertical
	$res  = 0.19;
 }

 if ($count == 3 && $nCalls <= 6) {
	$posX = 67; # Horizontal
	$posY = 90; # Vertical
	$res  = 0.19;
 }
 if ($count == 0 && $nCalls > 6) {
	$posX = 60;
	$posY = 250;
	$res  = 0.19;
 }
 return ($posX, $posY, $res);
}

###############################################################
sub nocallpages{
	my ($runname,$i,$pdf) = @_;

	my $nocallspage = $pdf->page(0);
	my $loc=775;
	$i++;

	my $nocalltitle = $nocallspage->text(); # Add some text to the page
	$nocalltitle->font($font, 20);
	$nocalltitle->translate(75, $loc); #(1_540 hor / 1_830 vert)
	$nocalltitle->text("Samples with no Structural Variants detected:");
	$loc=$loc-20;

	foreach my $nc (@nocalls){
		$loc=$loc-20;
		my $nocalltext = $nocallspage->text(); # Add some text to the page
		$nocalltext->font($font, 15);
		$nocalltext->translate(75, $loc); #(1_540 hor / 1_830 vert)
		$nocalltext->text("* $nc");
	}
}

###############################################################
sub getdate{
	my ($sec,$min,$hour,$day,$month,$year) = localtime;
	$year += 1900;
	$month += 1;
	return($sec,$min,$hour,$day,$month,$year);
}

#############################################
sub super_textlabel {
    my ($page, $x, $y, $font, $size, $text, $rotate, $lead) = @_;
    my $BIG = 1_000_000;
    $page->gfx()->save();
    my $txt = $page->text();
    $txt->font($font, $size);
    $txt->lead($lead);
    $txt->transform(-translate => [$x, $y], -rotate => $rotate);
    $txt->section($text, $BIG, $BIG);
    $page->gfx()->restore();
}

###############################################################
sub pagenumbers{
    my ($pdf,$outpdf) = @_;
    $pdf = PDF::API2->open("$outpdf");
    my $text_runpage = "Excluded ROIs due to extreme GC or low mappability";
    for my $index (1 .. $pdf->pages) {
	#my $page = $pdf->openpage($index);
	#if ($index == 2) {
		#my $txt  = $page->text;
		#$txt->textlabel(290, 50, $pdf->corefont('Helvetica'), 8, "*Low GC (< 45%), High GC (> 60&) and Low mappability (<100) are shown in yellow.");
	#}
	#else {
		my $page = $pdf->openpage($index);
		my $txt  = $page->text;
		$txt->textlabel(290, 50, $pdf->corefont('Helvetica'), 12, "$index");
	#}
    }
    $pdf->saveas("$outpdf");
    $pdf->end;
}
