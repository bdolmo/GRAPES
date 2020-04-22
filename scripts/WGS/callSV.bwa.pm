#!/usr/bin/env perl

package callSV;

use strict;
use warnings;
use Sort::Key::Natural qw(natsort);

my $devNullStderr = '2> /dev/null';


sub call {

    my $fastq   = shift;
    my $outFile = shift;
    my $minMapQ = shift;
    my $minSize = shift;
    my $maxSize = shift;
    my $numCPU  = shift;
    my $genome  = shift;
    
    my $totalFR;
    my $totalRF;
    my $totalFF;
    my $totalRR;

    my %seen = ();

    my $str = `$::cat $::outDir/discordantInfo.txt`;
    chomp $str;

    my @tmpStr = split (/\t/, $str);

    my $delIns = $tmpStr[0];
    my $dupIns = $tmpStr[1];
    my $invIns = $tmpStr[2] + $tmpStr[3];

    my %svIns = (
        DEL => $delIns,
        DUP => $dupIns,
        INV => $invIns
    );

    open (SAM, "$::bwa mem -t $numCPU $genome $fastq -M -U 0 -k 15 -w 500 $devNullStderr |");
    open (OUT, ">", $outFile) || die " ERROR: Unable to open $outFile\n";
    while (my $read =<SAM>) {
        chomp $read;
        next if $read =~/^@/;
        my @tmp = split (/\t/, $read);

        next if $tmp[4] < $minMapQ;

        $seen{$tmp[0]}++;
        if ( $seen{$tmp[0]} > 1) {
            #next;
        }

        my %svCall = bamRecordToSV($read, $minMapQ);

        if (%svCall && $svCall{MapQ} > $minMapQ ) {
            my $svSize = $svCall{End}-$svCall{Start};
            next if $svSize < $minSize;
            next if $svSize > $maxSize;

            print OUT "$svCall{ChrID}\t$svCall{Start}\t$svCall{End}\tPRECISE;SVTYPE=$svCall{Svtype};BREAKREADS=$svCall{BreakReads};ASSEMBLED=$svCall{Assembled};KDIV=$svCall{Kdiv}";
            print OUT ";MAPQ=$svCall{MapQ};PE=0;CSDISC=0;NINS=$svIns{$svCall{Svtype}};RDratio=.;RDmad=.;RDsupp=.;LOHsupp=.\n";
        }
    }
    close SAM;
    close OUT;
    #unlink($fastq);
}

################
sub cigar2Position {
    my $cigar = shift;
    my $Pos   = shift;
    my $outPos;
    $cigar =~s/H/S/;
    if ($cigar=~ /^[0-9]+[S]+[0-9]+[M]*/) {
        my @tmpCigar = split (/[SM]/, $cigar);
        $outPos = $Pos;
    }   
    elsif ($cigar =~/^[0-9]+[M]+[0-9]+[S]*/) {
        my @tmpCigar = split (/[MS]/, $cigar);
        $outPos = $tmpCigar[0] + $Pos-1;
    }

    return $outPos;
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
sub getBreakCoord  {

    my $read = shift;

    my @tmpSplit = split (",", $read);

    my $outChr = $tmpSplit[0];
    my $Pos    = $tmpSplit[1];
    my $strand = $tmpSplit[2];
    my $outStart = cigar2Position($tmpSplit[3], $Pos);
    my $mapq   = $tmpSplit[4];

    return ($outChr, $outStart, $strand, $mapq);
}

################
sub bamRecordToSV {

    my $read = shift;
    my $minMapQ = shift;
    my @tmp  = split (/\t/, $read);

    my @info = split ("_", $tmp[0]);

    my ($SA) = grep ($_ =~/SA:Z/, @tmp);
    $SA=~s/SA:Z:// if $SA;

    my ($scPos, $scType, $breakReads, $Assembled, $Kdiv) = Qname2Info($tmp[0]);

    my $chrP1 = $tmp[2];
    my $posP1 = cigar2Position($tmp[5],$tmp[3]);
  	my $strandP1 = ($tmp[1] & 0x0010) ? '-' : '+';
    my $mapqP1 = $tmp[4];
    $Kdiv = sprintf "%.2f",$Kdiv;

    # Què volem fer?
    # Tracking de tots els possibles reordenaments a través dels split-reads presents a SA:Z
    # Per fer això, hem de tractar cada parell

    my $svtype;
    my ($chrP2, $posP2, $strandP2, $mapqP2);
    my %hash = ();
    if ($SA) {
        my @SplitReads = split (/;/, $SA);
        foreach my $split (@SplitReads) {

            ($chrP2, $posP2, $strandP2, $mapqP2) = getBreakCoord($split);

            #next if $mapqP2 < $minMapQ;

            if ($chrP1 ne $chrP2) {
                #$svtype = "TRA";
            }
            else {
                if ($strandP1 ne $strandP2) {
                 $svtype = "INV";
                }
                else {
                    if ($scType eq 'RIGHT' && $posP2 > $scPos) {
                        $svtype = "DEL";
                    }
                    elsif ($scType eq 'LEFT' && $posP2 < $scPos) {
                        $svtype = "DEL";
                    }
                    elsif ($scType eq 'RIGHT' && $posP2 < $scPos) {
                        $svtype = "DUP";
                    }
                    elsif ($scType eq 'LEFT' && $posP2 > $scPos) {
                        $svtype = "DUP";
                    }                
                }
            }
        }
    }
    else {
        return %hash;
    }

    if (!$svtype) {
        return %hash;
    }

    my $meanMapq = int ($mapqP1+$mapqP2)/2;

    my $tmpStart;
    my $tmpEnd;
    my $Start = $posP1;
    my $End = $posP2;

    my $changeNtd = $svtype ne 'DUP' ? -1 : 0;

    if ($posP1 > $posP2) {
        $tmpStart = $posP2;
        $tmpEnd   = $posP1;
        $Start = $tmpStart;
        $End   = $tmpEnd;
    }

    %hash = (
        Qname => $tmp[0],
        ChrID => $chrP2,
        Start => $Start-$changeNtd,
        End   => $End+$changeNtd,
        Svtype => $svtype,
        BreakReads => $breakReads,
        Assembled => $Assembled,
        Kdiv => $Kdiv,
        MapQ => $meanMapq
    );

    return %hash;
}

1;