#!/usr/bin/env perl

package callSV;

use strict;
use warnings;
use Sort::Key::Natural qw(natsort);

my $devNullStderr = '2> /dev/null';

sub call  {

    my $fastq   = shift;
    my $outFile = shift;
    my $minMapQ = shift;
    my $minSize = shift;
    my $maxSize = shift;
    my $numCPU  = shift;
    my $genome  = shift;

    my %readAlign = ();

    my $str = `$::cat $::outDir/$::outName.discordantInfo.txt`;
    chomp $str;

    my @tmpStr = split (/\t/, $str);

    my $delIns = $tmpStr[0];
    my $dupIns = $tmpStr[1];
    my $invIns = $tmpStr[2] + $tmpStr[3];

    my %svIns = (
        DEL => 0,
        DUP => 0,
        INV => 0
    );

    my $cmd = "$::bwa mem -t $numCPU $genome $fastq -B 18 -O 32 -A 2 -w 1000 -M $devNullStderr > $::outDir/$::outName.contigs.sam";
    
    system $cmd if !-e "$::outDir/$::outName.contigs.sam";

    open (OUT, ">", $outFile) || die " ERROR: Unable to open $outFile\n";
    open (IN, "<", "$::outDir/$::outName.contigs.sam" ) || die " ERROR: Unable to open $::outDir/$::outName.contigs.sam\n";
    #open (IN, "<", "$::outDir/test.sam" ) || die " ERROR: Unable to open $::outDir/$::outName.contigs.sam\n";

    while (my $line=<IN>) {
        chomp $line;
        next if $line =~/^@/;
        my @tmp = split (/\t/, $line);
        
        my $readID = $tmp[0];
        my $chr    = $tmp[2];
        my $posStartAlign= $tmp[3];
        my $cigar  = $tmp[5];
        my $mapQual= $tmp[4];

        #my ($SA) = grep ($_=~/SA:Z/, @tmp);
        #$SA =~s/SA:Z://;

        my ($posEndAlign, $softClass, $indelNbases) = cigar2Position($cigar, $posStartAlign);

        my $strand = ($tmp[1] & 0x0010) ? '-' : '+';

        if (!$posEndAlign) {
            next;
          #  print "$chr\t$posStartAlign\t$cigar\n";
        }

        my $suppInfo = "$chr\t$posStartAlign\t$softClass\t$posEndAlign\t$mapQual\t$strand\t$indelNbases";

        push( @{ $readAlign { $readID } }, $suppInfo); 
    }
    close IN;

    foreach my $readID ( natsort keys %readAlign) {
        
        my ($scPos, $scType, $breakReads, $Assembled, $Kdiv) = Qname2Info($readID);

        my @suppArray = natsort (@{$readAlign{$readID}});

        my @mapqArr = ();

        if (@suppArray == 1) { 
            my ($chr, $posStartAlign, $softClass, $posEndAlign, $mapQual, $strand, $indelNbases) = getInfo($suppArray[0]);
            #print "$chr:$posStartAlign\t$posEndAlign\t$softClass\n";

            if (($softClass eq 'DEL' || $softClass eq 'INS') ) {
                my $start = $posEndAlign;
                my $end   = $start+$indelNbases;
                
                print OUT "$chr\t$start\t$end\tPRECISE;SVTYPE=$softClass;BREAKREADS=$breakReads;ASSEMBLED=$Assembled;KDIV=$Kdiv";
                print OUT ";MAPQ=$mapQual;PE=0;CSDISC=0;NINS=0;RDratio=.;RDmad=.;RDsupp=.;LOHsupp=.\n";
                next;
            }
        }
        for (my $i =0; $i < @suppArray-1; $i++) {
            my ($chrA, $posStartAlignA, $softClassA, $posEndAlignA, $mapQualA, $strandA, $indelNbasesA) = getInfo($suppArray[$i]);
            my ($chrB, $posStartAlignB, $softClassB, $posEndAlignB, $mapQualB, $strandB, $indelNbasesB) = getInfo($suppArray[$i+1]);
            my $start;
            my $end;
            my $svtype;
            push @mapqArr, $mapQualA;
            push @mapqArr, $mapQualB;

            if ($chrA ne $chrB) {
                next;
            }
            else {
                if ($strandA eq $strandB) {
                    if ($softClassA eq 'RIGHT') {
                        if ($posStartAlignB > $posEndAlignA) {
                            $start = $posEndAlignA+1;
                            $end   = $posStartAlignB-1;
                            $svtype = "DEL";                    
                        }
                        else {
                            $start = $posStartAlignB;
                            $end   = $posEndAlignA;
                            $svtype = "DUP";                    
                        }
                    }
                    if ($softClassA eq 'LEFT') {
                        if ($posEndAlignB < $posStartAlignA) {
                            $start = $posEndAlignB+1;
                            $end   = $posStartAlignA-1;
                            $svtype = "DEL";
                        }
                        else {
                            $start = $posStartAlignA;
                            $end   = $posEndAlignB;
                            $svtype = "DUP"; 
                        }   
                    }
                }
                else {
                    if ($posEndAlignB > $posEndAlignA) {
                        $start = $posEndAlignA;
                        $end   = $posEndAlignB;
                        $svtype = "INV";
                    }
                    else {
                        $start = $posStartAlignA;
                        $end   = $posStartAlignB;
                        $svtype = "INV";
                    } 
                }
            }
            my $meanMapq = meanArray(@mapqArr);
            if ($svtype) {
                print OUT "$chrA\t$start\t$end\tPRECISE;SVTYPE=$svtype;BREAKREADS=$breakReads;ASSEMBLED=$Assembled;KDIV=$Kdiv";
                print OUT ";MAPQ=$meanMapq;PE=0;CSDISC=0;NINS=$svIns{$svtype};RDratio=.;RDmad=.;RDsupp=.;LOHsupp=.\n";
            }
        }
    }
    close OUT;
}
################

sub meanArray {
    my @array =@_;

    my $sum = 0;
    foreach my $value (@array) {
        $sum+=$value;
    }
    my $mean = $sum/scalar@array;
    return $mean;
}

################
sub Qname2Info {
    my $Qname = shift;

    #r0_chr12_101885983_RIGHT_5_5_0.352518

    my @tmp = split (/_/, $Qname);
    my $scPos      = $tmp[2];
    my $scType     = $tmp[3];
    my $breakReads = $tmp[4];
    my $Assembled  = $tmp[5];
    my $Kdiv       = $tmp[6];

    return ($scPos, $scType, $breakReads, $Assembled, $Kdiv);
}

################
sub getInfo {
    my $info = shift;
    my ($chr, $posStartAlign, $softClass, $posEndAlign, $mapQual, $strand, $indelNbases) = split (/\t/, $info);
    return ($chr, $posStartAlign, $softClass, $posEndAlign, $mapQual, $strand, $indelNbases);
}

################
sub cigar2Position {
    my $cigar = shift;
    my $Pos   = shift;
    #my $SA    = shift;

    #my @tmpSA = split (/[,;]/, $SA);

    my $outPos;
    $cigar =~s/H/S/g;
    my $softClass = "UNK";
    my $indelNbases = 0;

    if ($cigar =~ /^[0-9]+[S]+[0-9]+[M]+[0-9]+S$/) {
        my @tmpCigar = split (/[SM]/, $cigar);
        $outPos = $Pos;

        if ($tmpCigar[2] >= $tmpCigar[0]) {
            $softClass = "RIGHT";
           $outPos = $Pos + $tmpCigar[1];
        }
        else {
            $softClass = "LEFT";
         }
    }
    elsif ($cigar=~ /^[0-9]+[S]+[0-9]+M$/) {
        $softClass = "LEFT";
        my @tmpCigar = split (/[SM]/, $cigar);
        $outPos = $Pos;
    }   
    elsif ($cigar =~/^[0-9]+[M]+[0-9]+S$/) {
        $softClass = "RIGHT";
        my @tmpCigar = split (/[MS]/, $cigar);
        $outPos = $tmpCigar[0] + $Pos-1;
    }
    elsif ($cigar =~/^[0-9]+M[0-9]+D[0-9]+M$/) {
        $softClass = "DEL";
        my @tmpCigar = split (/[MD]/, $cigar);
        $outPos = $tmpCigar[0] + $Pos-1;
        $indelNbases = $tmpCigar[1];
    }
    elsif ($cigar =~/^[0-9]+M[0-9]+I[0-9]+M$/) {
        $softClass = "INS";
        my @tmpCigar = split (/[MI]/, $cigar);
        $outPos = $tmpCigar[0] + $Pos-1;
        $indelNbases = $tmpCigar[1];
    }
    #111S 9M 1D 73M
    elsif ($cigar =~/^[0-9]+[S]+[0-9]+M[0-9]+[DI][0-9]+M$/) {
        $softClass = "LEFT";
        my @tmpCigar = split (/[SM]/, $cigar);
        $outPos = $Pos;
    }
    #M D M S
    elsif ($cigar =~/^[0-9]+[M]+[0-9]+[DI][0-9]+M[0-9]+S$/) {
        $softClass = "RIGHT";
        my @tmpCigar = split (/[MSID]/, $cigar);
        $outPos = $tmpCigar[0] + $tmpCigar[1] + $tmpCigar[2]+ $Pos-1;
    }
    elsif ($cigar =~/^[0-9]+[S]+[0-9]+M[0-9]+[DI][0-9]+M[0-9]+S$/) {
        $softClass = "LEFT";
        my @tmpCigar = split (/[SM]/, $cigar);
        $outPos = $Pos;
    }
    elsif ($cigar =~/^[0-9]+[M]+[0-9]+[DI][0-9]+M[0-9]+S$/) {
        $softClass = "LEFT";
        my @tmpCigar = split (/[SM]/, $cigar);
        $outPos = $Pos;
    }
    #S M I M S
    #chr1	199570742	29M2I6M2I23M109H
    #100H14M2I27M5H


    return $outPos, $softClass, $indelNbases;
}


1;