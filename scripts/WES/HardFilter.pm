#!/usr/bin/env perl

package HardFilter;

use strict;
use warnings;
use Sort::Key::Natural qw(natsort);

# Annotate on the VCF's FILTER field if the variant passes the quality criteria (PASS) or not (LowQual)
# Todo: if lowqual, annotate why!
my @arrFields = ( "END", "SVTYPE", "CIPOS", "CIEND", "SVLEN", "EV", "REGIONS", "GENE", "MAPQ", "KDIV", "GC", 
"MAP", "BR", "ASBR", "PE", "PPE", "RRD", "MADRD", "CN", "SNR", "ZSCORE", "PRD", "NSNV", "BAF");

#################################
sub filter {

    my $unfilteredVCF = shift; #Input is unfiltered vcf

    my $filteredVCF = $unfilteredVCF;
    $filteredVCF =~s/.vcf/.tmp.vcf/;
    my %seen = ();
    open (IN, "<", $unfilteredVCF) || die " ERROR: Unable to open $unfilteredVCF\n";
    open (OUT, ">", $filteredVCF)  || die " ERROR: Unable to open $filteredVCF\n";
    while (my $line=<IN>) {
        chomp $line;

        if ($line =~/^#/) {
            print OUT "$line\n";
            next;
        }

        my @tmp = split (/\t/, $line);
        my $joined = join ("\t", @tmp[0..4]);  
        $seen{$joined}++;
        
        # Skipping already observed variants
        next if $seen{$joined} > 1;
        my %fields = parse($line);

        my $endPos = $fields{END};
        $endPos =~s/END=//;

        $fields{FILTER} = $tmp[6];
        my $Baf = $fields{BAF};
        $Baf =~s/BAF=//;
        if ($fields{SVTYPE} =~/DEL/) {
            # For deletions we only accept hemizygous variants
            if ($Baf ne "." ) {
                if ( $Baf < 0.85  ) {
                    $fields{FILTER} = "LowQual";
                }
            }
        }
        elsif ($fields{SVTYPE}=~/DUP/) {
            #(AAB, ABB) are 0, 0.33, 0.66 and 1 respectively."
            # For deletions we only accept hemizygous variants
            if ($Baf ne "." ) {
                
                if ( $Baf > 0.4 && $Baf < 0.6  ) {
                    $fields{FILTER} = "LowQual";
                }
            }
        }
        # Check if variant fully overlaps an intron
        my $is_true = annotIntronPseudogene($tmp[0], $tmp[1], $endPos);
        if ($is_true) {
            $fields{FILTER} = "ProcessedPseudogene";
        }

        my @Arr = ();
        foreach my $field (@arrFields) {
            push @Arr, $fields{$field};
        }
       
        $fields{SVTYPE} =~s/SVTYPE=//;
        print OUT "$tmp[0]\t$tmp[1]\t.\tN\t<$fields{SVTYPE}>\t.\t$fields{FILTER}\t$fields{PRECISION};". join(";", @Arr) . "\tGT:CN\t./.:.\n";
        # Todo: filter duplications
    }
    close IN;
    close OUT;

    unlink ($unfilteredVCF);
    rename $filteredVCF, $unfilteredVCF;
}

#################################
sub parse {
    my $line = shift;
    my @t = split (/\t/, $line);

    #chr10	112595625	.	N	<DUP>	.	PASS	IMPRECISE;END=112595736;SVTYPE=DUP;SVLEN=111;EV=RD;REGIONS=1
    my @tmp = split (/;/, $t[7]);
    my %hash = (
        PRECISION => $tmp[0],
        END => "END=" . fetchPattern("END", \@tmp),
        SVTYPE => "SVTYPE=" . fetchPattern("SVTYPE", \@tmp),
        CIPOS  => "CIPOS=" . fetchPattern("CIPOS", \@tmp),
        CIEND  => "CIEND=" . fetchPattern("CIEND", \@tmp),
        SVLEN  => "SVLEN=". fetchPattern("SVLEN", \@tmp),
        EV     => "EV=RD",
        REGIONS=> "REGIONS=". fetchPattern("REGIONS", \@tmp),
        GENE   => "GENE=" . fetchPattern("GENE", \@tmp),
        MAPQ   => "MAPQ=" . fetchPattern("MAP", \@tmp),
        KDIV   => "KDIV=" . fetchPattern("KDIV", \@tmp),
        GC     => "GC=" . fetchPattern("GC", \@tmp),
        MAP    => "MAP=" . fetchPattern("MAP", \@tmp),
        BR     => "BR=" . fetchPattern("BR", \@tmp),
        ASBR   => "ASBR=" . fetchPattern("ASBR", \@tmp),
        PE     => "PE=" . fetchPattern("PE", \@tmp),
        PPE    => "PPE=" . fetchPattern("PPE", \@tmp),
        RRD    => "RRD=" . fetchPattern("RRD", \@tmp),
        CN     => "CN=" . fetchPattern("CN", \@tmp),
        MADRD  => "MADRD=" . fetchPattern("MADRD", \@tmp),
        SNR    => "SNR=" . fetchPattern("SNR", \@tmp),
        ZSCORE => "ZSCORE=" . fetchPattern("ZSCORE", \@tmp),
        PRD    => "PRD=" . fetchPattern("PRD", \@tmp),
        NSNV   => "NSNV=" .fetchPattern("NSNV", \@tmp),
        BAF    => "BAF=" . fetchPattern("BAF", \@tmp)
    );
    return %hash;
}
###########################
sub annotIntronPseudogene {
    # We have noted that most PRECISE SVs called from exome
    # actually come from pseudogenes without intronic sequences
    # We might filter them as LowQual
    # Super fast thanks to tabix index search
    # Retuls TRUE if CNV overlaps an intron

    my $chr  = shift;
    my $start= shift;
    my $end  = shift;

    my $length_A = $end-$start;
    #print "tabix $::intronFile $chr:$start-$end\n";
    #print "$chr\t$start\t$end\n";
    my $found = 0;
    my @query = split ("\n", `$::tabix $::intronFile $chr:$start-$end`);
    foreach my $intron (@query) {
        #print "$intron\n";
        my ($chrom, $st, $ed, $info, $info2, $strand) = split (/\t/, $intron);
        my $length_B = $ed-$st;
        my $olap_A;
        my $olap_B;

        # Minimum reciprocal overlap of 0.9 
        # Case1: coord fits within intron
        #A          ////////////   cnv (start, end)
        #B        ^^^^^^^^^^^^^^^^ intron (st, ed)
        if ($start >= $st && $end <= $ed) {
            $olap_A = 1.00;
            $olap_B = 1.00 - ( ( ($start-$st) + ($ed-$end) ) / $length_B);
        }
        # Case2: coord left-overlap
        #A      ////////////////   cnv (start, end)
        #B        ^^^^^^^^^^^^^^^^ intron (st, ed)
        if ($start < $st && $end < $ed) {
            $olap_A = 1.00 - ( ( $st-$start ) / $length_A);
            $olap_B = 1.00 - ( ( ($ed-$end) ) / $length_B);
            #    print "$start:$end=>$st:$ed\t$olap_A\t$olap_B\n";
        }
        # Case3: coord right-overlap
        #A           /////////////// cnv (start, end)
        #B        ^^^^^^^^^^^^^^^^   intron (st, ed)
        if ($start > $st && $end > $ed) {
            $olap_A = 1.00 - ( ( $end-$ed ) / $length_A);
            $olap_B = 1.00 - ( ( ($start-$st) ) / $length_B);
        }
        # Case4: coord fits entirely
        #A      ///////////////////// cnv (start, end)
        #B        ^^^^^^^^^^^^^^^^    intron (st, ed)
        if ($start <= $st && $end >= $ed) {
            $olap_A = 1.00 - ( ( ($st-$start) + ($end-$ed) ) / $length_A);
            $olap_B = 1.00;
        }
        if ($olap_A >= 0.8 && $olap_B >= 0.8) {
            $found = 1;
            last;
        }
    }  
    return $found;
}

###########################
sub fetchPattern {
    my $pattern = shift;
    my $arr  = shift;
    my @array = @$arr;

    my ($field) = grep ($_ =~/^$pattern=/, @array);
    $field = "$pattern=." if !$field;

    $field =~s/$pattern=//;
    return $field;
}

return 1;