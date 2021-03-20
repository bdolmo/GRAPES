#!/usr/bin/env perl

package QC;

use strict;
use Getopt::Long;
use File::Basename; 
use Excel::Writer::XLSX;
use Sort::Key::Natural qw(natsort);

###############################
sub writeRoiQC {

    my $outputDir  = shift;

	my $ratiosFile = $::HoF{RATIOS_ON};
    my %ratioData = ();

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
	my $nLine = 0;
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