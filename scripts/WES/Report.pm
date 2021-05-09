#!/usr/bin/env perl

package Report;

use strict;
use Getopt::Long;
use File::Basename;
use Excel::Writer::XLSX;
use Sort::Key::Natural qw(natsort);

###############################
sub exportCnvXLS  {

  my $workbook  = Excel::Writer::XLSX->new( $::HoF{ALL_CALLS} );
  my $worksheet = $workbook->add_worksheet("Results");
  my $format = $workbook->add_format();
  $format->set_pattern();

  #my $light_green = $workbook->set_custom_color(40, 152, 251, 152);

  my $row = 0;
  $format  = $workbook->add_format();
  $format->set_bold();

  $worksheet->write( $row, 0, 'Chromosome', $format);
  $worksheet->write( $row, 1, 'Start', $format);
  $worksheet->write( $row, 2, 'End', $format);
  $worksheet->write( $row, 3, 'Type', $format);
  $worksheet->write( $row, 4, 'Breakpoint supp.', $format);
  $worksheet->write( $row, 5, 'CN', $format);
  $worksheet->write( $row, 6, 'N Rois', $format);
  $worksheet->write( $row, 7, 'Flanking rois', $format);
  $worksheet->write( $row, 8, 'Genes', $format);
  $worksheet->write( $row, 9, 'Ratio', $format);
  $worksheet->write( $row, 10, 'Score', $format);
  $worksheet->write( $row, 11, 'Sample', $format);

	foreach my $sample ( natsort keys %::sampleHash ) {

    	open(IN, "<", $::sampleHash{$sample}{READY_VCF})
        || die " ERROR: Unable to open $::sampleHash{$sample}{READY_VCF}\n";

      while (my $line=<IN>) {
        chomp $line;
        #Skip VCF header
        next if $line=~/^#/;
        my @tmp = split("\t", $line);
        $row++;

        my $chr = $tmp[0];
        my $pos = $tmp[1];

        my @info = split(";", $tmp[7]);
        my ($end) = grep($_=~/^END=/, @info);
        $end =~s/END=//;
        $end = "." if !$end;

        my ($svtype) = grep($_=~/^SVTYPE=/, @info);
        $svtype =~s/SVTYPE=//;
        $svtype = "." if !$svtype;

        my ($splitReads) = grep($_=~/^BR=/, @info);
        $splitReads =~s/BR=//;
        $splitReads = "." if !$splitReads;

        my ($assembled) = grep($_=~/^ASBR=/, @info);
        $assembled =~s/ASBR=//;
        $assembled = "." if !$assembled;

        my ($pairedEnds) = grep($_=~/^PE=/, @info);
        $pairedEnds =~s/PE=//;
        $pairedEnds = "." if !$pairedEnds;

        my ($CN) = grep($_=~/^CN=/, @info);
        $CN =~s/CN=//;
        $CN = "." if !$CN;

        my ($NRois) = grep($_=~/^REGIONS=/, @info);
        $NRois =~s/REGIONS=//;
        $NRois = "." if !$NRois;

        my ($flankRois) = grep($_=~/^GENE=/, @info);
        $flankRois =~s/GENE=//;
        $flankRois = "." if !$flankRois;

        my ($genes) = grep($_=~/^GENES=/, @info);
        $genes =~s/GENES=//;
        $genes = "." if !$genes;

        my ($ratio) = grep($_=~/^RRD=/, @info);
        $ratio =~s/RRD=//;
        $ratio = "." if !$ratio;

        my ($score) = grep($_=~/^CONF_SCORE=/, @info);
        $score =~s/CONF_SCORE=//;
        $score = "." if !$score;

        $worksheet->write( $row, 0, $chr);
        $worksheet->write( $row, 1, $pos);
        $worksheet->write( $row, 2, $end);
        $worksheet->write( $row, 3, $svtype);

        if ($splitReads ne '.' || $assembled ne '.' || $pairedEnds ne '.') {
          my $results = "Split-reads=$splitReads,Assembled=$assembled,Paired-end=$pairedEnds";
          $worksheet->write( $row, 4, $results);
        }
        else {
          $worksheet->write( $row, 4, ".");
        }

        $worksheet->write( $row, 5, $CN);
        $worksheet->write( $row, 6, $NRois);
        $worksheet->write( $row, 7, $flankRois);
        $worksheet->write( $row, 8, $genes);
        $worksheet->write( $row, 9, $ratio);
        if ($score ne ".") {
          if ($score <= 0.95 && $score > 0.90) {
            my $format1  = $workbook->add_format();
            $format1->set_bg_color("#98FB98");
            $worksheet->write( $row, 10, $score, $format1);
          }
          if ($score <= 0.90 && $score > 0.85) {
            my $format1  = $workbook->add_format();
            $format1->set_bg_color("#ff9999");
            $worksheet->write( $row, 10, $score, $format1);
          }
          if ($score <=  0.85) {
            my $format1  = $workbook->add_format();
            $format1->set_bg_color("#ff3232");
            $worksheet->write( $row, 10, $score, $format1);
          }
          if ($score > 0.95) {
            my $format1  = $workbook->add_format();
            $format1->set_bg_color("#32CD32");
            $worksheet->write( $row, 10, $score, $format1);
          }
        }
        $worksheet->write( $row, 11, $sample);
      }
      close IN;
	}
  $workbook->close();

}

###############################
sub writeRoiQC {

  my $outputDir  = shift;

	my $ratiosFile = $::HoF{RATIOS_ON};
  my %ratioData  = ();

	my $is_gz = 0;
	if (!-e $ratiosFile) {
		$ratiosFile = $ratiosFile . ".gz";
		$is_gz = 1;
	}
	if ($is_gz) {
		open (IN, "$::zcat $ratiosFile |") || die " ERROR: Unable to open $ratiosFile\n";
	}
	else{
		open (IN, "<", $ratiosFile) || die " ERROR: Unable to open $ratiosFile\n";
	}

	my @header = ();
	my $nLine  = 0;
	while (my $line=<IN>) {
		chomp $line;
		my @tmp = split ("\t", $line);
		$nLine++;
		if ($nLine == 1) {
			@header = @tmp;
      next;
		}

    my $chr   = $tmp[0];
    my $start = $tmp[1];
    my $end   = $tmp[2];
    my $exon  = $tmp[3];
    my $roi = join("\t", @tmp[0..3]);
		for (my $i = 6; $i <@tmp; $i+=2){
      my $sampleName = $header[$i];
      $ratioData{$sampleName}{$roi}{RATIO} = $tmp[$i];
      $ratioData{$sampleName}{$roi}{S2N}   = $tmp[$i+1];
      if ($ratioData{$sampleName}{$roi}{S2N} > 5) {
        $::sampleHash{$sampleName}{N_CALL}++;
        $ratioData{$sampleName}{$roi}{EXON_QC} = "PASS";
      }
      else{
        $::sampleHash{$sampleName}{N_NO_CALL}++;
        $ratioData{$sampleName}{$roi}{EXON_QC} = "FAIL";
      }
    }
	}
	close IN;


    my $workbook  = Excel::Writer::XLSX->new( $outputDir . "/" . "roiQC.xlsx" );
    my %HoW = ();
    foreach my $sample (sort keys %::sampleHash) {
        $HoW{$sample} = $workbook->add_worksheet($sample);
        # Write header
        my $row = 0;
        $HoW{$sample}->write( $row, 0, 'Chromosome');
        $HoW{$sample}->write( $row, 1, 'Start');
        $HoW{$sample}->write( $row, 2, 'End');
        $HoW{$sample}->write( $row, 3, 'Exon');
        $HoW{$sample}->write( $row, 4, 'Ratio');
        $HoW{$sample}->write( $row, 5, 'Signal to noise');
        $HoW{$sample}->write( $row, 6, 'QC');
        $row++;
        foreach my $roi (natsort keys %{$ratioData{$sample}}) {
            my ($chr, $start, $end, $exon) = split("\t", $roi);
            $HoW{$sample}->write( $row, 0, $chr);
            $HoW{$sample}->write( $row, 1, $start);
            $HoW{$sample}->write( $row, 2, $end);
            $HoW{$sample}->write( $row, 3, $exon);
            $HoW{$sample}->write( $row, 4, $ratioData{$sample}{$roi}{RATIO});
            $HoW{$sample}->write( $row, 5, $ratioData{$sample}{$roi}{S2N});
            $HoW{$sample}->write( $row, 6, $ratioData{$sample}{$roi}{EXON_QC});
            $row++;
        }

    }
    $workbook->close();
}


###############################
sub writeSampleQC {

    my $outputDir = shift;

    my $sampleQc = $outputDir . "/" . "sampleQC.csv";
    open (SAMPLE_QC, ">", $sampleQc) || die " ERROR: Unable to open $sampleQc\n";
    my @header = (
       "Sample",
       "N Total Reads",
       "N On-target Reads",
       "N Off-target Reads",
       "Mean On-target Ratio",
       "Std On-target Ratio",
       "Mean Off-target Ratio",
       "Std Off-target Ratio",
       "N Calls",
       "N NoCalls",
       "On-target Qual Pass",
       "Off-target Qual Pass");

    print SAMPLE_QC join(",", @header ) . "\n";
    foreach my $sample (sort keys %::sampleHash) {

        $::sampleHash{$sample}{N_NO_CALL} = scalar keys %::filteredRois;
        $::sampleHash{$sample}{N_CALL} = scalar @::ROIarray - $::sampleHash{$sample}{N_NO_CALL};

        my @data = ();
        push @data, $sample;
        push @data, $::sampleHash{$sample}{TOTALREADS} ? $::sampleHash{$sample}{TOTALREADS} : ".";
        push @data, $::sampleHash{$sample}{READSONTARGET} ? $::sampleHash{$sample}{READSONTARGET}: ".";
        push @data, $::sampleHash{$sample}{READSOFFTARGET} ? $::sampleHash{$sample}{READSOFFTARGET} : ".";
        push @data, $::sampleHash{$sample}{ONTARGET_MEAN_RATIO} ? $::sampleHash{$sample}{ONTARGET_MEAN_RATIO} : ".";
        push @data, $::sampleHash{$sample}{ONTARGET_SD_RATIO} ? $::sampleHash{$sample}{ONTARGET_SD_RATIO} : ".";
        push @data, $::sampleHash{$sample}{OFFTARGET_MEAN_RATIO} ? $::sampleHash{$sample}{OFFTARGET_MEAN_RATIO} : ".";
        push @data, $::sampleHash{$sample}{OFFTARGET_SD_RATIO} ? $::sampleHash{$sample}{OFFTARGET_SD_RATIO} : ".";
        push @data, $::sampleHash{$sample}{N_CALL};
        push @data, $::sampleHash{$sample}{N_NO_CALL};
        push @data, $::sampleHash{$sample}{ONTARGET_QC_PASS} ? $::sampleHash{$sample}{ONTARGET_QC_PASS} : ".";
        push @data, $::sampleHash{$sample}{OFFTARGET_QC_PASS} ? $::sampleHash{$sample}{OFFTARGET_QC_PASS} : ".";
        print SAMPLE_QC join(",", @data ) . "\n";
    }
    close SAMPLE_QC;
}

return 1;
