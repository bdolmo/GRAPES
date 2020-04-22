#!/usr/bin/env perl

package mergeSV;

use strict;
use warnings;
use Statistics::Descriptive;

sub Merge {

    #my $bamFile  = shift;
    my $inputBed = shift;

    my $outputBed = $inputBed;
    $outputBed =~s/.tmp.rawcalls.bed/.merged.bed/;

    #my $outputBed= shift;

    my $flag = 0;
    my $count = 0;

    my %HoA = ();
    my %HoB = ();

    my @Starts = ();
    my @Ends   = ();
    my @Mapq   = ();
    my @Kdiv   = ();
    my @Ratios = ();
    my @MAD    = ();
    my @LOH    = ();
    my @Precisions = ();

    my $Breakreads = 0;
    my $Assembled  = 0;
    my $Discordants= 0;
    my $Cumulative = 0;

    my $Chr;
    my $el = 0;
    my %seen = ();

    my $cnvr_A;
    my $cnvr_B;
    my $checkCNVR;

    my $nLines =`$::cat $inputBed | $::wc -l`;
    chomp $nLines;

    open (OUT, ">", $outputBed) || die " ERROR: unable to open $outputBed\n";
    my $num = 0;
    open (IN, "$::sort -V $inputBed | $::uniq |" );
    while (my $line =<IN>) {
        chomp $line;
        $count++;

        my @tmp = split (/\t/, $line);
        my $Size = $tmp[2]-$tmp[1];
        next if $Size < $::minSizeSV;
        next if $Size > $::maxSizeSV;

        if ($flag == 0) {
            %HoA = parseCall ( $line );

            if ($nLines > 1) {
                push @Starts,"$HoA{start}\t$HoA{precision}";
                push @Ends,  "$HoA{end}\t$HoA{precision}";
                push @Mapq, $HoA{mapq};
                push @Kdiv, $HoA{kdiv};
                push @Ratios, $HoA{rdRatio};
                push @MAD, $HoA{rdMad};
                push @LOH, $HoA{lohSup};
                push @Precisions, $HoA{precision};

                $Breakreads+=$HoA{breakreads};
                $Assembled+=$HoA{assembled};
                $Discordants+=$HoA{readPairs};
                if ( $HoA{csdisc} ne '.' && $HoA{csdisc} >0 ) {
                    $Cumulative+=$HoA{csdisc};
                }        
                $flag = 1;
            }
            next if $nLines > 1;
        }
        if ($flag == 1 || $nLines == 1) {
            %HoB = parseCall ( $line );

            my ($is_olap) = reciprocalOverlap( $HoA{chr}, $HoB{chr}, $HoA{start}, $HoA{end}, $HoB{start}, $HoB{end}, $HoA{precision}, $HoB{precision});
            if ($is_olap ) { 
                $checkCNVR = checkCNVR($HoA{svtype}, $HoB{svtype});
            }
            else {
                $checkCNVR = 0;
            }

            if ( ($HoA{chr} eq $HoB{chr} ) && ( $HoA{svtype} eq $HoB{svtype} || $checkCNVR )  &&  ($is_olap)){

                push @Starts,"$HoB{start}\t$HoB{precision}";
                push @Ends,  "$HoB{end}\t$HoB{precision}";
                push @Mapq, $HoB{mapq};
                push @Kdiv, $HoB{kdiv};
                push @Ratios, $HoB{rdRatio};
                push @MAD, $HoB{rdMad};
                push @LOH, $HoB{lohSup};
                push @Precisions, $HoB{precision};

                $Breakreads+=$HoB{breakreads};
                $Assembled+=$HoB{assembled};
                $Discordants+=$HoB{readPairs};
                if ( $HoB{csdisc} ne '.' && $HoB{csdisc} >0 ) {
                    $Cumulative+=$HoB{csdisc};
                }
                $el++;
            }
            if (  ($HoA{chr} ne $HoB{chr}) || ( $HoA{svtype} ne $HoB{svtype} ) || (!$is_olap) || ($nLines == 1) || ($count == $nLines) ) {
                my $Chr      = $HoA{chr};
                my $Start    = updateCoord( @Starts );
                my $End      = updateCoord( @Ends );
                my $Precision= updatePrecision( @Precisions );
                my $SVTYPE   = $HoA{svtype};
                my $LOHsupp  = updateLOH ( @LOH );
                my $ratioRD  = updateRD ( @Ratios );
                my $madRD    = updateRD ( @MAD );
                my $svlen    = $End-$Start;
                my $MAPQ     = int(updateRD( @Mapq));
                my $KDIV     = updateRD( @Kdiv);

                if ( ( $HoA{chr} ne $HoB{chr} )  || ( $HoA{svtype} ne $HoB{svtype} ) || (!$is_olap) ) {
                    $el = 0;
                }
                next if $svlen > $::maxSizeSV;
                next if $svlen < $::minSizeSV;

                my $fallsInCNVR = "no";
                if ($checkCNVR) {
                   $fallsInCNVR = "yes"; 
                }
                if (!$checkCNVR && $cnvr_A) {
                   $fallsInCNVR = "yes"; 
                }
                my $AF = '.';

                my $Size = $End-$Start;
                if ($Size < $::maxSizeSV) {
                    print OUT "$Chr\t$Start\t$End\t$Precision\t$SVTYPE\t$fallsInCNVR\t$svlen\t$MAPQ\t$KDIV\t$AF\t$Breakreads\t$Assembled\t$Discordants\t$ratioRD\t$madRD\t$LOHsupp\t$Cumulative\t$HoA{nins}\n";
                }
                if ( $count == $nLines && !$el) {
                    $svlen = $HoB{end}-$HoB{start};
                    if ($svlen < $::maxSizeSV) {
                        print OUT "$HoB{chr}\t$HoB{start}\t$HoB{end}\t$HoB{precision}\t$HoB{svtype}\t$fallsInCNVR\t$svlen\t$HoB{mapq}\t$HoB{kdiv}\t$AF\t$HoB{breakreads}\t$HoB{assembled}\t$HoB{readPairs}\t$HoB{rdRatio}\t$HoB{rdMad}\t$HoB{lohSup}\t$Cumulative\t$HoB{nins}\n";
                    }
                }
            
                @Starts = (); 
                @Ends = ();
                @Mapq = ();
                @Kdiv = ();
                @Ratios = ();
                @MAD = ();
                @LOH = ();
                @Precisions = ();
                $num = 0;
                push @Starts,"$HoB{start}\t$HoB{precision}";
                push @Ends,  "$HoB{end}\t$HoB{precision}";
                push @Mapq, $HoB{mapq};
                push @Kdiv, $HoB{kdiv};
                push @Ratios, $HoB{rdRatio};
                push @MAD, $HoB{rdMad};
                push @LOH, $HoB{lohSup};
                push @Precisions, $HoB{precision};

                %HoA = %HoB;
                $Breakreads  = $HoB{breakreads};
                $Assembled   = $HoB{assembled};
                $Discordants = $HoB{readPairs};

                if ( $HoB{csdisc} ne '.' ) {
                    $Cumulative = $HoB{csdisc};
                }
                else {
                    $Cumulative = 0;
                }

                $cnvr_A = $checkCNVR;
            }
        }
    }
    close IN;

    return $outputBed;
}

#################################################
sub updateCoord {
    my @array = @_;
    my @I = ();
    my @P = ();
    my $returnPosition;
    foreach my $val (@array) {
        my ($position, $precision) = split (/\t/ , $val);
        push @I, $position if $precision eq 'IMPRECISE';
        push @P, $position if $precision eq 'PRECISE';
    }
    if (@P) {
        $returnPosition = getMostCommon( @P );
    }
    else {
        $returnPosition = getMean( @I );
    }
    return (int $returnPosition);
}
#################################################
sub updateRD {
    my @array = @_;
    my $meanRD;
    my @num = ();
    foreach my $val (@array) {
           push @num, $val if $val ne '.';
    }
    if (@num) {
        $meanRD = getMean (@num);
    }
    else {
        $meanRD = '.';
    }
    return $meanRD;
}
#################################################
sub updateLOH {
    my @array = @_;

    my %hash = ();
    foreach my $val (@array) {
        $hash{$val}++;
    }
    if ( exists $hash{yes} ) {
        return "yes";
    }
    if (exists $hash{no}) {
        return "no";
    }
    if (exists $hash{"."}) {
        return ".";
    }
}

#################################################
sub updatePrecision {
    my @array = @_;

    foreach my $val (@array) {
        if ($val eq 'PRECISE') {
            return "PRECISE";
        }
    }
    return "IMPRECISE";
}

#################################################
sub parseCall {
    
    my $line = shift;

    my @tmp = split (/\t/, $line);
    my @info = split (/;/, $tmp[3]);

    my ($svtype) = grep ($_=~/SVTYPE=/, @info);
	$svtype =~s/SVTYPE=//;	

    my $svlen = $tmp[2]-$tmp[1];

	my ($breakreads) = grep ($_=~/BREAKREADS=/, @info);
	$breakreads =~s/BREAKREADS=//;	

	my ($assembled) = grep ($_=~/ASSEMBLED=/, @info);
	$assembled =~s/ASSEMBLED=//;
	
	my ($kdiv) = grep ($_=~/KDIV=/, @info);
	$kdiv =~s/KDIV=//;

	my ($mapq) = grep ($_=~/MAPQ=/, @info);
	$mapq =~s/MAPQ=//;
    #$mapq = '.' if !$mapq;

	my ($PE) = grep ($_=~/^PE=/, @info);
	$PE =~s/^PE=//;	

	my ($RDratio) = grep ($_=~/RDratio=/, @info);
	$RDratio =~s/RDratio=//;
    $RDratio = '.' if !$RDratio;	

	my ($RDmad) = grep ($_=~/RDmad=/, @info);
	$RDmad =~s/RDmad=//;
    $RDmad = '.' if !$RDmad;

	my ($lohSup) = grep ($_=~/LOHsupp=/, @info);
	$lohSup =~s/LOHsupp=//;
    $lohSup = '.' if !$lohSup;	

    my ($csdisc) = grep ($_=~/CSDISC=/, @info);
	$csdisc =~s/CSDISC=//;
    $csdisc = 0 if !$csdisc;	

    my ($nins) = grep ($_=~/NINS=/, @info);
	$nins =~s/NINS=//;
    $nins = 0 if !$nins;

    my %hash = (
        chr   => $tmp[0],
        start => $tmp[1],
        end   => $tmp[2],
        precision => $info[0],
        svtype => $svtype,
        svlen  => $svlen,
        breakreads => $breakreads,
        assembled => $assembled,
        kdiv => $kdiv,
        mapq => $mapq,
        readPairs => $PE,
        rdRatio => $RDratio,
        rdMad  => $RDmad,
        lohSup => $lohSup,
        csdisc => $csdisc,
        nins => $nins
    );

    return %hash;
}

#################################################

sub returnMinOlap  {

    my $precision_A = shift;
    my $precision_B = shift;

    my $minOlap;
	if ($precision_A eq "IMPRECISE" || $precision_B eq "IMPRECISE") {
		$minOlap = 0.2;
    }
	elsif ($precision_A eq "PRECISE" && $precision_B eq "PRECISE") {
		$minOlap = 0.7;
	}    

    return $minOlap;
}

#################################################

sub getCoverage  {

    my $bamFile = shift;
    my $chr     = shift;
    my $start   = shift;
    my $end     = shift;
    my $precision=shift;

    my $st_A = $start-50;
    my $st_B = $start-49;

    my $ed_A = $start;
    my $ed_B = $start+1;

    if ($precision eq 'IMPRECISE') {
        $ed_A = $start+70;
        $ed_B = $start+71;
    }

    my $coord_upstream   = "$chr:$st_A-$st_B";
    my $coord_downstream = "$chr:$ed_A-$ed_B";

    my $counts_outside = `$::samtools view -c $bamFile $coord_upstream`;
    chomp $counts_outside;

    my $counts_inside = `$::samtools view -c $bamFile $coord_downstream`;
    chomp $counts_inside;

	my $ratio;

	if ($counts_outside == 0) {
		$ratio = 0;
	}
	else {
		$ratio = $counts_inside/$counts_outside;
	}
	return sprintf "%.2f", $ratio;
}

#################################################
sub checkCNVR {
    # Goal here is to intersect DEL and DUP calls with high overlapping rate (>70) that could be Copy-Number Variable Regions (CNVR)
    my $svtype_A = shift;
    my $svtype_B = shift;

    if ($svtype_A eq 'DEL' && $svtype_B eq 'DUP' || $svtype_A eq 'DUP' && $svtype_B eq 'DEL') {
        return 1;
    }
    else {
        return 0;
    }
}

#################################################
sub reciprocalOverlap {

    my $chr_A   = shift;
    my $chr_B   = shift;
    my $start_A = shift;
    my $end_A   = shift;
    my $start_B = shift;
    my $end_B   = shift;
    my $precision_A = shift;
    my $precision_B = shift;    

    my $length_A = $end_A - $start_A;
    my $length_B = $end_B - $start_B;

    my $olap_A;
    my $olap_B;

    if ($chr_A ne $chr_B) {
        return 0;
    }
    if ($start_B > $end_A) {
        $olap_A = 0.00;
        $olap_B = 0.00;
    }
    if ($start_A == $start_B && $end_A == $end_B) {
        $olap_A = 1.00;
        $olap_B = 1.00;
    }
    if ($start_B == $start_A && $end_B > $end_A) {
        #   //////////////
        #   ///////////////
        $olap_A = 1.00; 
		$olap_B = 1.00 - ((($start_B-$start_A) + ($end_B-$end_A) ) / $length_B) ;

    }    
    if ($start_B >= $start_A && $end_B <= $end_A) {
        #   //////////////
        #   //////////
		$olap_A = 1.00 - ((($start_B-$start_A) + ($end_A-$end_B) ) / $length_A) ;
        $olap_B = 1.00;
    }

	#// case3: B is adjacent to A
	if ($start_A < $start_B && $end_B >= $end_A) {

        #///////////
        #    //////////
		$olap_A = 1.00 - ( ($start_B - $start_A )/ $length_A);
		$olap_B = 1.00 - ( ($end_B - $end_A)/ $length_B);
		
		if ($olap_A < 0) { $olap_A = 0; }
		if ($olap_B < 0) { $olap_B = 0; }
	}


    my $minOlap = returnMinOlap($precision_A, $precision_B);

    if ($olap_A >= $minOlap && $olap_B >= $minOlap ) {
        return 1;
    }
    else {
        return 0;
    }
}

###################
sub getMostCommon {
    my @array = @_;

    my %hash = ();
    my $max = 0;
    my $mostCommon;
    foreach my $x (@array) {
        $hash{$x}++;
        if ($max < $hash{$x}) {
            $max = $hash{$x};
            $mostCommon = $x;
        }
    }
    return $mostCommon;
}

###################
sub getMean {
    my @data = @_;
    my $stat = Statistics::Descriptive::Full->new();
    $stat -> add_data(@data);
    my $mean = $stat->mean();
    return sprintf "%.3f",$mean;
}

1;