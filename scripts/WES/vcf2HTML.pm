#!/usr/bin/env perl

package vcf2HTML;

use strict;
use warnings;
use diagnostics;
use File::Basename;
use Sort::Key::Natural qw(natsort);
use Scalar::Util qw(looks_like_number);
use List::MoreUtils qw(uniq);

#Goal is to print an static html report

my %varPass   = ();
my %varFilter = ();
my %totalVars = ();
my %annotVars = ();
my %var       = ();
my %vcfInfo   = ();
my $accumulateModal = "";

######################### 
# Bulk function
sub createHTMLreport {

    my $dirName = dirname($::outDir);
    my $runName = basename($::outDir);
    my $HTML    = $dirName . "/$runName/$runName.html";

    # Re-using a global hash that stores sample names, directories, etc.
    foreach my $sample ( sort keys %::sampleHash ) {

        my $vcf = "$::outDir/$sample.CNV.vcf";
        if ( !-e $vcf ) {
            print "ERROR: Missing VCF for sample $sample\n";
            next;
        }
        readVCF($vcf, $sample);
	}

    open (HTML, ">", $HTML) || die " ERROR: Unable to open $HTML\n";
    my $head = printHead("$runName - Structural Variant Results");
    print HTML "$head\n";

    my $navBar = printNavBar();
    print HTML "$navBar\n";

    my $body = printBody();
    print HTML "$body\n";

    my $footer = printFooter();
    print HTML "$footer\n";
    close HTML;
}

##########################
sub getUCSCsnapshot {
    my $chr  = shift;
    my $start= shift;
    my $end  = shift;

    my $url = "<iframe src=\"http://genome.ucsc.edu/cgi-bin/hgRenderTracks?db=hg19&position=$chr%3A$start-$end&hideTracks=1&ncbiRefSeq=dense&guidelines=off&hgt.labelWidth=0\" width=\"1000\" height=\"200\"></iframe>\n";
    return $url;
}
##########################
sub fetchGenes {
    my $chr  = shift;
    my $start= shift;
    my $end  = shift;

    open (BED, ">", "$::outDir/$chr.$start.$end.bed");
    print BED "$chr\t$start\t$end\n";
    close BED;

    my $result = `$::bedtools intersect -a $::geneList -b $::outDir/$chr.$start.$end.bed`;
    chomp $result;

    my %HoG = (); #Hash of Genes
    my @tmpResult = split (/\n/,$result);
    foreach my $line (@tmpResult) {
        my ($chr, $start, $end, $strand, $ensemblg, $name, $geneFeature) = split(/\t/, $line);
        $HoG{"$chr\t$start\t$end"}{NAME}    = $name;
        $HoG{"$chr\t$start\t$end"}{ENSG}    = $ensemblg;
        $HoG{"$chr\t$start\t$end"}{STRAND}  = $ensemblg;
        $HoG{"$chr\t$start\t$end"}{FEATURE} = $ensemblg;
    }

    unlink("$::outDir/$chr.$start.$end.bed");

    return %HoG;
}
########################## Not working
sub returnIdeogram {

  my $chr   = shift;
  my $start = shift;
  my $end   = shift;
  my $id    = shift;
  my $sample= shift;

  my $chromo = $chr;
  $chromo =~s/chr//;

  my $ideogram = "
  <script>
    var config = {
      container: \'#$chr$start$end-ideogram\',
      organism: \'human\',
      chromosome: \'$chromo\',
      orientation: \'horizontal\',
      annotations: [{
        name: \'$chr:$start-$end\',
        chr: \'$chromo\',
        start: $start,
        stop: $end
      }]
    };
    var ideogram = new Ideogram(config);
    </script>
    \n";

  return $ideogram;
}

##########################
sub createMultipleExonPlots2 {

    my $plotFile = shift;
    my $sample   = shift;
    my $id       = shift;
    my $start    = shift;
    my $end      = shift;

    my $multipleExonPlot;
    my @xarr  = ();
    my $count = 0;

    my @sampleArr;
    my @controlArr;
    my %HoROIs      = (); #  Hash of ROIs
    my %seenCoord   = ();
    my @ROIs        = ();
    my @annotations = (); # Saving ROI annotations (4th column from BED) for every point

     if (-e defined $plotFile) {

        open ( P, "$::gunzip -c \"$plotFile\" |") || die "ERROR: cannot open $plotFile\n";
        while (my $line=<P>) {
            chomp $line;
            my @tmp = split (/\t/, $line);

            my $roi = $tmp[3];
            push @ROIs, $roi;

            $seenCoord{"$tmp[0]\t$tmp[1]\t$tmp[2]"}++;

            if (!exists $HoROIs{$roi}) {
                $HoROIs{$roi}{START} = $tmp[1];
            }
            else {
                if ($tmp[2] >= $HoROIs{$roi}{START} ) {
                    $HoROIs{$roi}{END} = $tmp[2];
                }
            }
            if ($line =~/$sample/) {
                push @xarr, "\'$tmp[1]\'";
                push @sampleArr, sprintf "%.3f", $tmp[5];
            }         
            push @annotations, "'$tmp[3]'";
            #chrX	31140035	31140036	NM_004006_78_79;DMD	0.162	2_Controls
            #chrX	31140035    31140036	NM_004006_78_79;DMD	0.266	1_RB23317_9999999.rmdup
            $count++ if $seenCoord{"$tmp[0]\t$tmp[1]\t$tmp[2]"} < 2;
        }
        close P;

        @ROIs = uniq @ROIs;
        #@xarr = uniq @xarr;

        my @arrTickVals = ();
        my @arrTickText = ();

        my $x = join (",", @xarr);
        my $sampleRatios  = join (",", @sampleArr);
        my $sampleAnnotations  = join (",", @annotations);

        my $tickVals = "";
        my $tickText = "";
   
        $multipleExonPlot = "   
        <script>
                var Sample = {
                    name: \'$sample\',
                    x: [$x],
                    y: [$sampleRatios],
                    text: [$sampleAnnotations],
                    marker: {size:12, color: \'#FF4136\', line: { color:\'black\'}},
                    mode: \'lines+markers\',
                    type: \'scatter\',
                    visible:true,
                    responsize: true
                  };

                var data = [Sample];

                var layout$start$end = {
                    margin: {
                        t: 4
                    },
                    shapes: [
                    {
                        type: \'line\',
                        xref: \'paper\',
                        x0: 0,
                        y0: $::upperDelCutoff,
                        x1: 1,
                        y1: $::upperDelCutoff,
                        line:{
                            dash: \'dot\',
                            color: \'red\',
                            width: 1,
                        }
                    },
                    {
                        type: \'line\',
                        xref: \'paper\',
                        x0: 0,
                        y0: $::lowerDupCutoff,
                        x1: 1,
                        y1: $::lowerDupCutoff,
                        line:{
                            dash: \'dot\',
                            color: \'deepskyblue\',
                            width: 1,
                        }
                    }
                  ],
                    legend: {x: 0.4, y: 1.2},
                    xaxis: {
                        //tickvals:[$tickVals],
                        //ticktext:[$tickText],
                        //showgrid: true,
                          zeroline: true,
                          showline: true,
                          mirror: \'ticks\',
                        //gridcolor: \'black\',
                        //gridwidth: 1,
                        //tickangle: 30,
                          linecolor: \'black\',
                          linewidth: 2
                    },
                    yaxis: {
                        showgrid: false,
                        zeroline: true,
                        showline: true,
                        linecolor: \'black\',
                        linewidth: 2,
                        mirror: \'ticks\',

                        //gridcolor: \'black\',
                        //gridwidth: 1,
                        title: {
                            text: \'Copy ratio\',
                            font: {
                                family: \'Helvetica\',
                                size: 16,
                                color: \'black\'
                            }
                        },
                    },
                };
                Plotly.newPlot(\'$id\', data, layout$start$end);
        </script>
            \n";
    }    
    else {
        return 0;
    }
    return $multipleExonPlot;
}   


##########################
sub createMultipleExonPlots1 {

    my $plotFile = shift;
    my $sample   = shift;
    my $id       = shift;
    my $start    = shift;
    my $end      = shift;

    my $multipleExonPlot;
    my @xarr = ();
    my $count  = 0;

    my @sampleArr;
    my @controlArr;
    my %HoROIs = (); #  Hash of ROIs
    my %seenCoord = ();
    my @ROIs = ();

    if (-e $plotFile) {

        open ( P, "$::gunzip -c \"$plotFile\" |") || die "ERROR: cannot open $plotFile\n";
        while (my $line=<P>) {
            chomp $line;
            my @tmp = split (/\t/, $line);
            #print "$tmp[5]\n";
            next if $line =~/ND/;

            my $roi = $tmp[3];
            push @ROIs, $roi;

            $seenCoord{"$tmp[0]\t$tmp[1]\t$tmp[2]"}++;

            if (!exists $HoROIs{$roi}) {
                $HoROIs{$roi}{START} = $count;
            }
            else {
                if ($tmp[2] >= $HoROIs{$roi}{START} ) {
                    $HoROIs{$roi}{END} = $count;
                }
            }
            if ($line =~/$sample/) {
                push @xarr, "\'$tmp[3]\'";
                push @sampleArr, sprintf "%.3f", $tmp[4] 
            }
            
            else {
                push @controlArr, sprintf "%.3f", $tmp[4];
            }
            #chrX	31140035	31140036	NM_004006_78_79;DMD	0.162	2_Controls
            #chrX	31140035    31140036	NM_004006_78_79;DMD	0.266	1_RB23317_9999999.rmdup
            $count++ if $seenCoord{"$tmp[0]\t$tmp[1]\t$tmp[2]"} < 2;
        }
        close P;

        @ROIs = uniq @ROIs;
        #@xarr = uniq @xarr;

        my @arrTickVals = ();
        my @arrTickText = ();

        foreach my $roi (@ROIs) {
            push @arrTickVals, "\'$HoROIs{$roi}{START}\'";
            push @arrTickText, "\'$roi\'";
        }
        # Add end coordinate
        push @arrTickVals,"\'$count\'";
        push @arrTickText,"\'\'";

        my $x = join (",", @xarr);
        my $sampleRatios  = join (",", @sampleArr);
        my $controlRatios = join (",", @controlArr);
        my $tickVals = join (",", @arrTickVals);
        my $tickText = join (",", @arrTickText);

      #  print "$tickVals\n";
      #  print "$tickText\n";
      #  print "$x\n";
      #  print "$sampleRatios\n";
        $multipleExonPlot = "   
        <script>
                var Sample = {
                    name: \'$sample\',
                    marker: {color: \'#FF4136\', line: { color:\'black\'}},
                    x: [$x],
                    y: [$sampleRatios],
                    type: \'box\',
                };

                var Controls = {
                    name:  \'Controls\', 
                    marker: {color: \'grey\', line: { color:\'black\'}},
                    x: [$x],
                    y: [$controlRatios],
                    type: \'box\',
                };

                var data = [Sample, Controls];

                var layout$start$end = {
                    legend: {x: 0.4, y: 1.2},
                    xaxis: {
                        //tickvals:[$tickVals],
                        //ticktext:[$tickText],
                    //    showgrid: true,
                          zeroline: true,
                          showline: true,
                          mirror: \'ticks\',
                    //    gridcolor: \'black\',
                    //    gridwidth: 1,
                    //    tickangle: 30,
                          linecolor: \'black\',
                          linewidth: 2
                    },
                    yaxis: {
                        showgrid: false,
                        zeroline: true,
                        showline: true,
                        linecolor: \'black\',
                        linewidth: 2,
                        mirror: \'ticks\',

                        //gridcolor: \'black\',
                        //gridwidth: 1,
                        title: {
                            text: \'Normalized coverage\',
                            font: {
                                family: \'Helvetica\',
                                size: 16,
                                color: \'black\'
                            }
                        },
                    },
                    boxmode: \'group\',
                };
                Plotly.newPlot(\'$id\', data, layout$start$end);
        </script>
            \n";
    }    
    else {
        return 0;
    }
    return $multipleExonPlot;
}   

##########################
sub createSingleExonPlots {

    my $plotFile = shift;
    my $sample   = shift;
    my $id       = shift;
    #print "$plotFile\n";

    my $singleExonPlot;
    my @xarr = ();
    my $count  = 0;

    my @sampleArr;
    my @controlArr;

    if (defined $plotFile && -e $plotFile) {

        open ( P, "$::gunzip -c \"$plotFile\" |") || die "ERROR: cannot open $plotFile\n";
        while (my $line=<P>) {
            chomp $line;
            my @tmp = split (/\t/, $line);

            #push @xarr, $count;
            if ($line =~/$sample/) {
                push @sampleArr, sprintf "%.3f", $tmp[5];
            }
            else {
                push @controlArr, sprintf "%.3f", $tmp[5];
            }
            #chrX	31366747	31366748	NM_004006_60_61;DMD	RB23419_9999999.rmdup	1.31460674157303	RB23419_9999999.rmdup
            #chrX	31366747	31366748	NM_004006_60_61;DMD	RB23421_9999999.rmdup	0.967177242888403	control
            $count++;
        }
        close P;

        #my $x = join (",", @xarr);
        my $sampleRatios  = join (",", @sampleArr);
        my $controlRatios = join (",", @controlArr);

        $singleExonPlot = "   
        <script>
                var Sample = {
                    name: \'$sample\',
                    marker: {color: \'#FF4136\', line: { color:\'black\'}}, 
                    y: [$sampleRatios],
                    type: \'box\'
                };

                var Controls = {
                    name:  \'Controls\', 
                    marker: {color: \'grey\', line: { color:\'black\'}}, 
                    y: [$controlRatios],
                    type: \'box\'
                };

                var data = [Sample, Controls];

                var layout = {
                    legend: {x: 0.4, y: 1.2},
                    xaxis: {
                        showgrid: true,
                        zeroline: true,
                        showline: true,
                        linecolor: \'black\',
                        linewidth: 2
                    },
                    yaxis: {
                        showgrid: true,
                        zeroline: true,
                        showline: true,
                        linecolor: \'black\',
                        linewidth: 2,
                        title: {
                            text: \'Copy ratio\',
                            font: {
                                family: \'Helvetica\',
                                size: 16,
                                color: \'black\'
                            }
                        },                        
                    }
                };
                Plotly.newPlot(\'$id\', data, layout);
        </script>
            \n";
    }    
    else {
        return 0;
    }
    return $singleExonPlot;
}   

##########################
sub getFilteredRois {

    my $table .= "
        <div id =\"expansio\" style=\"margin:0px;width:100%;\">
        <table id =\"\" class=\"display hover block\" style=\"margin:0px;margin-top:15px;width:100%;padding:0px;\">
            <thead style=\"background-color:white;text-align:left;padding:0px;\">
                <tr style=\"padding:5px;\">
                    <th>Chr</th>
                    <th>Start</th>
                    <th>End</th>
                    <th>Size</th>
                    <th>ROI</th>
                    <th>\%GC</th>
                    <th>\%Mappability</th>
                    <th>Comment</th>
                </tr>
            </thead> 
            <tbody>\n";

    open (FILTERED, "<", $::onTargetFilteredROIs) || die " ERROR: Unable to open $::onTargetFilteredROIs\n";
    while (my $line=<FILTERED>) {
        chomp $line;
        next if $line =~/start/;
        my ($chr, $start, $end, $roi, $meanCov, $gc, $map, $comment) = split (/\t/, $line);
        $meanCov = int($meanCov);
        $gc      = int($gc);
        $map     = int($map);
        my $size = $end-$start;

        $table .= "<tr>";
        $table.="<td style=\"background-color:white;border-bottom: 1px solid black;text-align:left;\">$chr</td>";
        $table.="<td style=\"background-color:white;border-bottom: 1px solid black;text-align:left;\">$start</td>";
        $table.="<td style=\"background-color:white;border-bottom: 1px solid black;text-align:left;\">$end</td>";
        $table.="<td style=\"background-color:white;border-bottom: 1px solid black;text-align:left;\">$size</td>";

        $table.="<td style=\"background-color:white;border-bottom: 1px solid black;text-align:left;\">$roi</td>";
        $table.="<td style=\"background-color:white;border-bottom: 1px solid black;text-align:left;\">$gc</td>";
        $table.="<td style=\"background-color:white;border-bottom: 1px solid black;text-align:left;\">$map</td>";
        $table.="<td style=\"background-color:white;border-bottom: 1px solid black;text-align:left;\">$comment</td>";
        $table .= "</tr>";
    }
    close FILTERED;
    $table .= "</tbody>
    </table>
    </div>\n";

    return $table;
}
##########################
sub getSampleInfo {

    my $numReads =  `$::samtools idxstats $::bam | $::awk -F \'\t\' \'{s+=\$3+\$4}END{print s}\'`; chomp $numReads;
    $numReads = number2human($numReads);

    my $readLength = `$::samtools view $::bam -f 2 -F 512 | $::head -1000 | $::awk '{sum+=length(\$10)} END {print sum/NR}'`;
    chomp $readLength;

    return $numReads, $readLength;
}

##########################
sub number2human {
    my $number = shift;
       if(looks_like_number($number)) { 
        $number >= 1000000 ? sprintf( "%0.2fM", $number / 1000000 )
      : $number >= 1000    ? sprintf( "%0.0fK", $number / 1000 )
      :                      $number;
       }
       else
       {
            $number;
       }
}

##########################
sub getBafData {

    my $bafFile = shift;
    my $dirType = shift;

    my $bafTrace = "";

    if (!-e $bafFile) {
        print "$bafFile\n";
        return $bafTrace
    }

    my @xarr = ();
    my @yarr = ();

    my $xaxisTag;
    my $yaxisTag;
    my $visibility;
    if ($dirType eq "ontarget") {
        $xaxisTag = "x4";
        $yaxisTag = "y4";
        $visibility = "true";
    }
    if ($dirType eq "offtarget") {
        $xaxisTag = "x4";
        $yaxisTag = "y4";
        $visibility = "false";
    }
    if ($dirType eq "main") {
        $xaxisTag = "x4";
        $yaxisTag = "y4";
        $visibility = "false";
    }

    open (IN, "<", $bafFile) || die " ERROR: Unable to open $bafFile\n";
    while (my $line=<IN>) {
        #chr1	21571601	21571602	0.464	chr1	21403001	21636678	pdwindow_89	0.921	1.0128
        chomp $line;
        my @tmp = split (/\t/, $line);
        my $absPosition  = Utils::chrCoordinate2Absolute($tmp[0], $tmp[1]);

        push @xarr, "\'$absPosition\'";
        push @yarr, $tmp[3];
    }
    close IN;

    my $x  = join(',', @xarr);
    my $y  = join(',', @yarr);
 
    #print "$dirType\n";
    #print "$xaxisTag\n";
    #print "$yaxisTag\n";

    $bafTrace = "
    var baftrace$dirType = {
        x: [$x],
        y: [$y],
        mode: \'markers\',
        type: 'scatter',
        visible:$visibility,
        type: \'scattergl\',
        responsize: true,

        xaxis: \'$xaxisTag\',
        yaxis: \'$yaxisTag\',
        
        marker: {
            opacity: 0.5,
            size: 3.5,
            line: {
                color: \'#003333\',
                width: 1
            }
        }, 
    };\n";

    return ($bafTrace);
}

##########################
sub renderGenomeWideScatterPlot {

    my $sample = shift;

    # By default
    my $toPlotData = "$::outDir/SEGMENT_DATA/ON_TARGET/toplot.segmented.$sample.bed";
    if (!-e $toPlotData) {
       $toPlotData = "$::outDir/ON_TARGET/SEGMENT_DATA/toplot.segmented.$sample.bed"; 
    }

    my @xarr = ();
    my @yarr = ();
    my @ontarget = ();
    my @offtarget = ();

    my @segArr = ();

    my @ontargetArr;
    my @offtargetArr;
    my @segArrOntarget  = ();
    my @segArrOfftarget = ();

    my $sampleDivSelector = $sample;
    $sampleDivSelector=~s/\.//g;

    my @allChr;
    my @ontargetChr;
    my @offtargetChr;

    # Declaring arrays for text annotations 
    my @ontargetAnnoArr = ();
    my @offtargetAnnoArr = ();
    my @allAnnoArr = ();

    if (-e $toPlotData) {
        open (IN, "<", $toPlotData);
        my $xval = 0;
      
        while (my $line=<IN>) {
            chomp $line;
            #chr1	26380383	26380455	NM_032588_7_8;TRIM63	1.142	1.048
            my @tmp = split (/\t/, $line);
            my $absPosition  = Utils::chrCoordinate2Absolute($tmp[0], $tmp[1]);

            $xval++;
            #push @xarr, $xval;
            if ($line =~/window/) {
                push @offtarget, $tmp[4];
                #push @offtargetArr,"\'$tmp[3]\'";
                push @offtargetArr,"\'$absPosition\'";
                push @offtargetChr, "\'$tmp[0]\'";
                push @segArrOfftarget, $tmp[5];
                push @offtargetAnnoArr,"\'$tmp[3]\'";
            }
            else {
                push @ontarget, $tmp[4];
                #push @ontargetArr,"\'$tmp[3]\'";
                push @ontargetArr,"\'$absPosition\'";
                push @ontargetChr, "\'$tmp[0]\'";
                push @segArrOntarget, $tmp[5];
                push @ontargetAnnoArr,"\'$tmp[3]\'";

            }
            push @yarr, "\'$tmp[4]\'";
            #push @xarr, "\'$tmp[3]\'";
            push @xarr,"\'$absPosition\'";
            push @allChr,"\'$tmp[0]\'";
            push @segArr, $tmp[5];
            push @allAnnoArr, "\'$tmp[3]\'";

            last if $xval > 5000; # Avoiding huge files
        }
        close IN;
    }
    # Formatting ratio data:
    # For all ratios
    my $x  = join(',', @xarr);
    my $y  = join(',', @yarr);
    my $chrAll = join(',', @allChr);

    # For ontarget ratios
    my $chr_ontarget = join(',', @ontargetChr);
    my $x_ontarget  = join(',',  @ontargetArr);
    my $y_ontarget  = join(',',  @ontarget);

    # For offtarget ratios
    my $chr_offtarget = join(',', @offtargetChr);
    my $x_offtarget   = join(',', @offtargetArr);
    my $y_offtarget   = join(',', @offtarget);

    # Formatting ratio data:
    my $seg           = join(',', @segArr);
    my $segOntarget   = join(',', @segArrOntarget);
    my $segOfftarget  = join(',', @segArrOfftarget);

    # Formatting annotations
    my $ontargetAnnotations = join(',', @ontargetAnnoArr);
    my $offtargetAnnotations = join(',', @offtargetAnnoArr);
    my $allAnnotations = join(',', @allAnnoArr);

    # BAF data
    my $bafTraceAll = "";
    my $bafTraceOntarget = "";
    my $bafTraceOfftarget = "";

    if ($::doVAF) {
        $bafTraceAll      = getBafData("$::outDir/BAF_DATA/$sample.baf_segment.bed", "main");
        $bafTraceOntarget = getBafData("$::outDir/BAF_DATA/$sample.baf_segment.bed", "ontarget");
        $bafTraceOfftarget= getBafData("$::outDir/BAF_DATA/$sample.baf_segment.bed", "offtarget");
    }

    my $yDomain = "0,1";

    my $bafTraces = "";
    if ($bafTraceAll ne "" && $bafTraceOntarget ne "" && $bafTraceOfftarget ne "" ) {
       $bafTraces = ", baftraceontarget, baftraceofftarget, baftracemain";
       $yDomain = "0.5,1";
    }

    my $scatter = "
    <script>
        var updatemenus=[
            {
                buttons: [
                    {
                        args: [{\'visible\': [true, false, true, false, false, false, true, false, false]},\'type\', \'scattergl\'],
                        label: \'on-target\',
                        method: \'update\'
                    },
                    {
                        args: [{\'visible\': [false, true, false, true, false, false, false, true, false]}, \'type\', \'scattergl\'],
                        label:\'off-target\',
                        method:\'update\'
                    },
                    {
                        args: [{\'visible\': [false, false, false, false, true, true, false, false, true]}, \'type\', \'scattergl\'],
                        label:\'All\',
                        method:\'update\'
                    }
                ],
                direction: \'left\',
                pad: {\'r\': 10, \'t\': 10, \'t\': 10},
                showactive: true,
                type: \'buttons\',
                x: 0.1,
                xanchor: \'left\',
                y: 1.2,
                yanchor: \'top\'
            }
        ]

        var onratio = {
            x: [$x_ontarget],
            y: [$y_ontarget],
            mode: \'markers\',
            text: [$ontargetAnnotations],
            marker: {
                opacity: 0.5,
                line: {
                    color: \'##000099\',
                    width: 1
                }
            }, 
            visible:true,
            type: \'scattergl\',
            responsive: true,
        };

        var offratio = {
            x: [$x_offtarget],
            y: [$y_offtarget],
            text: [$offtargetAnnotations],
            mode: \'markers\',
            marker: {
                opacity: 0.5,
                color: \'grey\',
            }, 
            visible:false,
            type: \'scattergl\',
            responsive: true,
        };

        var allratio = {
            x: [$x],
            y: [$y],
            text: [$offtargetAnnotations],
            mode: \'markers\',
            marker: {
                opacity: 0.5,
                color: \'grey\',
            }, 
            visible:false,
            type: \'scattergl\',
            responsive: true,

        };

        var segment = {
            x: [$x],
            y: [$seg],
            mode: \'markers\',
            marker: {
                size: 2.5,
                color: \'red\'
            },
            visible:false,
            type: \'scattergl\',
         };

        var segmentOntarget = {
            x: [$x_ontarget],
            y: [$segOntarget],
            mode: \'markers\',
            marker: {
                size: 2.5,
                color: \'red\'
            }, 
            visible:true,
            type: \'scattergl\',
        };

        var segmentOfftarget = {
            x: [$x_offtarget],
            y: [$segOfftarget],
            mode: \'markers\',
            marker: {
                size: 2.5,
                color: \'red\'
            }, 
            visible:false,
            type: \'scattergl\',
        };
        $bafTraceOntarget
        $bafTraceOfftarget
        $bafTraceAll
        var data = [onratio, offratio, segmentOntarget, segmentOfftarget, allratio, segment $bafTraces];
        var layout$sampleDivSelector = {
            updatemenus: updatemenus,
            showlegend: false,
            hovermode: \'closest\',
            margin: {
                t: 5,
            },
            shapes: [
            {
                type: \'line\',
                xref: \'paper\',
                x0: 0,
                y0: $::upperDelCutoff,
                x1: 1,
                y1: $::upperDelCutoff,
                line:{
                    dash: \'dot\',
                    color: \'red\',
                    width: 1,
                }
            },
            {
                type: \'line\',
                xref: \'paper\',
                x0: 0,
                y0: $::lowerDupCutoff,
                x1: 1,
                y1: $::lowerDupCutoff,
                line:{
                    dash: \'dot\',
                    color: \'deepskyblue\',
                    width: 1,
                }
            }
            ],
            xaxis: {
                showgrid: false,
                zeroline: false,
                showline: true,

                showticklabels: false,
                linewidth: 2,
                domain: [0,1],
            },
            xaxis2: {
                showgrid: false,
                zeroline: false,
                showline: true,

                showticklabels: false,
                linewidth: 2,
                domain: [0,1],
            },
            xaxis3: {
                showgrid: false,
                zeroline: false,
                showline: true,

                showticklabels: false,
                linewidth: 2,
                domain: [0,1],
            },
            xaxis4: {
                showgrid: false,
                zeroline: false,
                showline: true,

                showticklabels: false,
                linewidth: 2,
                domain: [0,1],
                anchor: \'y4\'
            },
            xaxis5: {
                showgrid: false,
                zeroline: false,
                showline: true,
                showticklabels: false,
                linewidth: 2,
                domain: [0,1],
                anchor: \'y5\'
            },
            xaxis6: {
                showgrid: false,
                zeroline: false,
                showline: true,
                showticklabels: false,
                linewidth: 2,
                domain: [0,1],
                anchor: \'y6\'
            },
            yaxis: {
                showgrid: false,
                showline: true,
                automargin: true,
                linewidth: 2,
                mirror: \'ticks\',
                domain: [$yDomain],
                anchor: \'x1\',
                title: {
                    text: \'Copy ratio\',
                    font: {
                        family: \'Helvetica\',
                        size: 16,
                        color: \'black\'
                    }
                },               
               range: [0, 2]
            },
            yaxis2: {
                showgrid: false,
                showline: true,
                automargin: true,
                linewidth: 2,
                mirror: \'ticks\',
                domain: [$yDomain],
                anchor: \'x2\',

                title: {
                    text: \'Copy ratio\',
                    font: {
                        family: \'Helvetica\',
                        size: 16,
                        color: \'black\'
                    }
                },               
               range: [0, 2]
            },
            yaxis3: {
                showgrid: false,
                showline: true,
                automargin: true,
                linewidth: 2,
                mirror: \'ticks\',
                domain: [$yDomain],
                anchor: \'x3\',

                title: {
                    text: \'Copy ratio\',
                    font: {
                        family: \'Helvetica\',
                        size: 16,
                        color: \'black\'
                    }
                },               
               range: [0, 2]
            },
            yaxis4: {
                showgrid: false,
                showline: true,
                automargin: true,
                linewidth: 2,
                mirror: \'ticks\',
                domain: [0,0.4],
                anchor: \"x4\",
                title: {
                    text: \'VAF\',
                    font: {
                        family: \'Helvetica\',
                        size: 16,
                        color: \'black\'
                    }
                },               
               range: [0, 1]
            },
            yaxis5: {
                showgrid: false,
                showline: true,
                automargin: true,
                linewidth: 2,
                mirror: \'ticks\',
                domain: [0,0.4],
                anchor: \"x5\",
                title: {
                    text: \'VAF\',
                    font: {
                        family: \'Helvetica\',
                        size: 16,
                        color: \'black\'
                    }
                },               
               range: [0, 1]
            },
            yaxis6: {
                showgrid: false,
                showline: true,
                automargin: true,
                linewidth: 2,
                mirror: \'ticks\',
                domain: [0,0.4],
                anchor: \"x6\",

                title: {
                    text: \'VAF\',
                    font: {
                        family: \'Helvetica\',
                        size: 16,
                        color: \'black\'
                    }
                },               
               range: [0, 1]
            },
            size: 18,
            color: \'#000000\',
            width: 0.7 * window.innerWidth,
            height: 0.55 * window.innerHeight,
        };
        Plotly.newPlot(\'$sample-scatter\', data, layout$sampleDivSelector, {responsive: true});
    </script>
    \n";

    return $scatter;
}

##########################
sub createSummaryTable {

    my $refs = "$::outDir/ON_TARGET/$::outName.references.txt";
    my %meanCorrs = ();

    open (CORRS, "<", $refs ) || die " ERROR: Unable to open $refs\n";
    while (my $line=<CORRS>) {
        chomp $line;
        my @tmp = split (/\t/, $line);

        next if $line =~/SAMPLE/;
        my $corr = $tmp[-1];
        $meanCorrs{$tmp[1]}{CORR} = $corr;
        if ( scalar @tmp < 5) {
            $meanCorrs{$tmp[1]}{CORR} = "NA";
            $meanCorrs{$tmp[1]}{FLAG} = "<td style=\"background-color:white;border-bottom: 1px solid black;text-align:left;\"><a href=\"#\" title=\"That&apos;s what this widget is\"></a><i class=\"fas fa-exclamation-circle\" style=\"color:#ffae42;\"></i></td>";
            next;
        }
        else {
            if ($corr < $::minCorrelation) {
                $meanCorrs{$tmp[1]}{FLAG} = "<td style=\"background-color:white;border-bottom: 1px solid black;text-align:left;\"><a href=\"#\" title=\"That&apos;s what this widget is\"></a><i class=\"fas fa-exclamation-circle\" style=\"color:#ffae42;\"></i></td>";
            }
            else {
                $meanCorrs{$tmp[1]}{FLAG} = "<td style=\"background-color:white;border-bottom: 1px solid black;text-align:left;\"><i class=\"fas fa-check-circle\" style=\"color:green;\"></i></td>";
            }
        }
    }
    close CORRS;

    my $table .= "
        <table id =\"\" class=\"display hover block\"style=\"margin-top:15px;float:left;width:100%;padding:0px;\">
            <thead style=\"background-color:white;text-align:left;padding:3px;\">
                <tr style=\"padding:5px;\">
                    <th>QC</th>
                    <th>Sample</th>
                    <th>#Reads</th>
                    <th>#Reads on ROI</th>
                    <th>#Reads off-target</th>
                    <th>\%ROI</th>
                    <th>Mean coverage</th>
                    <th>Correlation</th>
                    <th>Mean.isize</th>
                    <th>Std.isize</th>
                </tr>
            </thead> 
            <tbody>\n";

    foreach my $sample (natsort keys %::sampleHash) {
        
        my $numOfftarget = $::sampleHash{$sample}{TOTALREADS_ONTARGET}-$::sampleHash{$sample}{READSONTARGET};
        $numOfftarget = number2human($numOfftarget);
        my $readsOnRoi= number2human($::sampleHash{$sample}{READSONTARGET});
        my $numReads  = number2human ($::sampleHash{$sample}{TOTALREADS_ONTARGET});
        $::sampleHash{$sample}{MEANCOV_ONTARGET} = int($::sampleHash{$sample}{MEANCOV_ONTARGET}) . "X";
        $table .= "
           <tr>
           $meanCorrs{$sample}{FLAG}
            <td style=\"background-color:white;border-bottom: 1px solid black;text-align:left;\">$sample</td>
            <td style=\"background-color:white;border-bottom: 1px solid black;text-align:left;\">$numReads</td>
            <td style=\"background-color:white;border-bottom: 1px solid black;text-align:left;\">$readsOnRoi</td>
            <td style=\"background-color:white;border-bottom: 1px solid black;text-align:left;\">$numOfftarget</td>

            <td style=\"background-color:white;border-bottom: 1px solid black;text-align:left;\">$::sampleHash{$sample}{ROI_ONTARGET}</td>
            <td style=\"background-color:white;border-bottom: 1px solid black;text-align:left;\">$::sampleHash{$sample}{MEANCOV_ONTARGET}</td>
            <td style=\"background-color:white;border-bottom: 1px solid black;text-align:left;\">$meanCorrs{$sample}{CORR}</td>
            <td style=\"background-color:white;border-bottom: 1px solid black;text-align:left;\">$::sampleHash{$sample}{MEANISIZE}</td>
            <td style=\"background-color:white;border-bottom: 1px solid black;text-align:left;\">$::sampleHash{$sample}{SDISIZE}</td>
           </tr>";
    }
    $table .=
    "</tbody>
    </table>\n";

    return $table;
}

#############
sub readVCF {
    my $inputVCF = shift;
    my $sample   = shift;

    open (IN, "<", $inputVCF) || die " ERROR: Unable to open $inputVCF\n";
    while (my $line=<IN>) {
        
        chomp $line;
        next if $line =~/^#/;

        my @tmp = split (/\t/, $line);
        my $sv = $tmp[4];
        $sv =~s/[<|>]//g;

        my @info = split (/;/, $tmp[7]);

        my ($svtype) = grep($_=~/SVTYPE/, @info);
        $svtype=~s/SVTYPE=//;

        my ($End) = grep($_=~/^END/, @info);
        $End=~s/END=//;

        my ($PE) = grep($_=~/^PE/, @info);
        $PE=~s/PE=//;

        my ($BR) = grep($_=~/^BR/, @info);
        $BR=~s/BR=//;

        my ($ASBR) = grep($_=~/^ASBR/, @info);
        $ASBR=~s/ASBR=//;

        #IMPRECISE;END=206944760;SVTYPE=DUP;SVLEN=60;EV=RD;REGIONS=1;GENE=NM_000572.2_1_2,IL10;MAPQ=100;KDIV=.;
        #GC=44;MAP=100;BR=.;ASBR=.;PE=.;PPE=.;RRD=1.359;MADRD=0.06;SNR=35.91;ZSCORE=3.789;PRD=.;NSNV=.;percHomSNV=.	GT:CN	./.:.

        my ($GENE) = grep($_=~/^GENE/, @info);
        $GENE=~s/GENE=//;

        my ($ROIs) = grep($_=~/^REGIONS/, @info);
        $ROIs=~s/REGIONS=//;

        my ($MAPQ) = grep($_=~/^MAPQ/, @info);
        $MAPQ=~s/MAPQ=//;
        $MAPQ = int ($MAPQ) if $MAPQ ne '.';

        my ($GC) = grep($_=~/^GC/, @info);
        $GC=~s/GC=//;
        $GC = int ($GC) if $GC ne '.';

        my ($RDRATIO) = grep($_=~/^RRD/, @info);
        $RDRATIO=~s/RRD=//;

        my ($ZSCORE) = grep($_=~/^ZSCORE/, @info);
        $ZSCORE=~s/ZSCORE=//;

        my ($SNR) = grep($_=~/^SNR/, @info);
        $SNR=~s/SNR=//;

        my $size = number2human($End-$tmp[1]);

        my $filterTag = $tmp[6] eq 'PASS' 
        ? "<input class=\"edit\" type=\"submit\" value=$tmp[6]>" 
        : "<input class=\"delete\" type=\"submit\" value=$tmp[6]>";

        my $link = "https://gnomad.broadinstitute.org/region/$tmp[0]-$tmp[1]-$End?dataset=gnomad_sv_r2_1";

        my ($AN) = grep($_=~/^AN=/, @info);
        if ($AN) {
            $AN =~s/AN=//;
        }
        else {
            $AN = '.';
        }

        my ($AC) = grep($_=~/^AC=/, @info);
        if ($AC) {
            $AC =~s/AC=//;
        }
        else {
            $AC = '.';
        }

        my ($AF) = grep($_=~/^AF/, @info);
        if ($AF) {
            $AF =~s/AF=//;
        }
        else {
            $AF = '.';
        }

        my ($AFR_AC) = grep($_=~/^AFR_AC/, @info);
        if ($AFR_AC) {
            $AFR_AC =~s/AFR_AC=//;
        }
        else {
            $AFR_AC = '.';
        }

        my ($AFR_AF) = grep($_=~/^AFR_AF/, @info);
        if ($AFR_AF) {
            $AFR_AF =~s/AFR_AF=//;
        }
        else {
            $AFR_AF = '.';
        }

        my ($AMR_AC) = grep($_=~/^AMR_AC/, @info);
        if ($AMR_AC) {
            $AMR_AC =~s/AMR_AC=//;
        }
        else {
            $AMR_AC = '.';
        }

        my ($AMR_AF) =  grep($_=~/^AMR_AF/, @info);
        if ($AMR_AF) {
            $AMR_AF =~s/AMR_AF=//;
        }
        else {
            $AMR_AF = '.';
        }

        my ($EAS_AC) = grep($_=~/^EAS_AC/, @info);
        if ($EAS_AC) {
            $EAS_AC =~s/EAS_AC=//;
        }
        else {
            $EAS_AC = '.';
        }

        my ($EAS_AF) = grep($_=~/^EAS_AF/, @info);
        if ($EAS_AF) {
            $EAS_AF =~s/EAS_AF=//;
        }
        else {
            $EAS_AF = '.';
        }

        my ($EUR_AC )= grep($_=~/^EUR_AC/, @info);
        if ($EUR_AC) {
            $EUR_AC =~s/EUR_AC=//;
        }
        else {
            $EUR_AC = '.';
        }

        #my $EUR_AN =
        my ($EUR_AF) = grep($_=~/^EUR_AF/, @info);
        if ($EUR_AF) {
            $EUR_AF =~s/EUR_AF=//;
        }
        else {
            $EUR_AF = '.';
        }

        my ($OTH_AC) = grep($_=~/^OTH_AC/, @info);
        if ($OTH_AC) {
            $OTH_AC =~s/OTH_AC=//;
        }
        else {
            $OTH_AC = '.';
        }
        my ($OTH_AF) = grep($_=~/^OTH_AF/, @info);
        if ($OTH_AF) {      
            $OTH_AF=~s/OTH_AF=//;
        }
        else {
            $OTH_AF = '.';
        }

        my $chromo = $tmp[0];
        $chromo =~s/chr//;

        my ($suportTag) = grep($_=~/^EV/, @info);
        $suportTag =~s/EV=//;

        #Check if variant has gzipped file for plotting 
        my $singleExonPlot = "";
        my $multipleExonPlot = "";
        my $plotDiv = "";

        my $tag = "$sample.$tmp[0].$tmp[1].$End";
        my ($savedPlot, $class) = getAvailablePlot("$tag.ratios.bed.gz");
        my $ideogram = returnIdeogram($tmp[0], $tmp[1], $End, $tag, $sample);
        #my $geneTrack = renderGeneTrck($tmp[0], $tmp[1], $End, $tag);
        my $geneTrack = renderGeneTrack2($tmp[0], $tmp[1], $End, $tag);

        my $snapshot = getUCSCsnapshot($tmp[0], $tmp[1], $End);
        my %HoG = fetchGenes($tmp[0], $tmp[1], $End);
        my $rectangle = createRectangle($tmp[1], $End, $tag, \%HoG);

        #print "$tag\n";
        #<div id =\"$tag-rectangle\"></div>
        #$rectangle
        if ($ROIs ne '.') {
        if ( $ROIs == 1) {
            $singleExonPlot = createSingleExonPlots($savedPlot, $sample, $tag) if $tag ne '.';
            if ($singleExonPlot) {
              
               $plotDiv = "
                <!-- Modal HTML embedded directly into document -->
                    <div id=\"$tmp[0]$tmp[1]$End\" class=\"modal\" style=\"float:center;text-align:center;max-width:600px;\">
                        <div style=\"float:center;width:600px;\">
                            <h3>$tmp[0]:$tmp[1]-$End $svtype <i>$GENE</i></h3>                                    
                            <div id=\"$tag\" style=\"float:center;width:500px;height:500px;margin: 0 auto;\"></div>
                            $singleExonPlot
                        </div>
                    </div>
                <!-- Link to open the modal -->
                 <p><a href=\"#$tmp[0]$tmp[1]$End\" rel=\"modal:open\"><i class=\"fas fa-eye fa-xs\"></i></a></p>\n";                           

            }
        }
        elsif ($ROIs > 1 && $ROIs < 15) {
            my $onlyOffDiv = "";
            if (defined $class eq 'on-target') {
                $multipleExonPlot = createMultipleExonPlots1($savedPlot, $sample, $tag, $tmp[1], $End) if $tag ne '.';
            }
            else {
                $multipleExonPlot = createMultipleExonPlots2($savedPlot, $sample, $tag, $tmp[0], $End) if $tag ne '.';

                $onlyOffDiv =
                "<div id=\"$tmp[0]$tmp[1]$End-ideogram\" style=\"display:block;float:center;width:700px;height:80px;margin: 0 auto;\"></div>
                $ideogram
                <div id =\"$tag-rectangle\" style=\"float:center;width:1200px;height:100px;margin: 0 auto;\"></div>
                $rectangle\n";
            }
            if ($multipleExonPlot) {
                $plotDiv = "
                <!-- Modal HTML embedded directly into document -->
                    <div id=\"$tmp[0]$tmp[1]$End\" class=\"modal\" style=\"float:center;text-align:center;max-width:1250px;\">
                        <div style=\"float:center;width:1200px;\">
                            <h3>$tmp[0]:$tmp[1]-$End $svtype <i>$GENE</i></h3>
                            $onlyOffDiv                               
                            <div id=\"$tag\" style=\"float:center;width:1200px;height:500px;margin: 0 auto;\"></div>
                            $multipleExonPlot
                        </div>
                    </div>
                <!-- Link to open the modal -->
                 <p><a href=\"#$tmp[0]$tmp[1]$End\" rel=\"modal:open\"><i class=\"fas fa-eye fa-xs\"></i></a></p>\n";
            }
        }                  

        elsif ($ROIs >=15) {
            $multipleExonPlot = createMultipleExonPlots2($savedPlot, $sample, $tag, $tmp[0], $End) if $tag ne '.';
 
            if ($multipleExonPlot) {
                $plotDiv = "

                <!-- Modal HTML embedded directly into document -->
                    <div id=\"$tmp[0]$tmp[1]$End\" class=\"modal\" style=\"background-color:white;float:center;text-align:center;max-width:1250px;\">
                        <div style=\"background-color:white;float:center;width:1200px;\">
                            <h3>$tmp[0]:$tmp[1]-$End $svtype <i>$GENE</i></h3>
                            <div id=\"$tmp[0]$tmp[1]$End-ideogram\" style=\"display:block;float:center;width:700px;height:80px;margin: 0 auto;\"></div>
                            $ideogram
                            <div id =\"$tag-rectangle\" style=\"float:center;width:1200px;height:100px;margin: 0 auto;\"></div>
                            $rectangle                                           
                            <div id=\"$tag\" style=\"float:center;width:1200px;height:400px;margin: 0 auto;\"></div>
                            $multipleExonPlot
                        </div>
                    </div>
                <!-- Link to open the modal -->
                 <p><a href=\"#$tmp[0]$tmp[1]$End\" rel=\"modal:open\"><i class=\"fas fa-eye fa-xs\"></i></a></p>\n";
            }
        }
        }
        if ($singleExonPlot eq ""  && $multipleExonPlot eq ""){
            $plotDiv = "<a href=\"#\" style=\"color:black\"><i class=\"fas fa-eye-slash fa-xs\"></i></a>";
        }
        #<div id=\"tnt-$tmp[1]$End\" style=\"float:center;width:1200px;height:200px;margin: 0 auto;\"></div>
        $accumulateModal .= "
             <script>
                // Get the modal
                var modal = document.getElementById(\"$sample-$tag\");

                // Get the button that opens the modal
                var btn = document.getElementById(\"$sample-$tag-button\");

                // Get the <span> element that closes the modal
                var span = document.getElementsByClassName(\"close\")[0];

                // When the user clicks the button, open the modal 
                btn.onclick = function() {
                modal.style.display = \"block\";
                }

                // When the user clicks on <span> (x), close the modal
                span.onclick = function() {
                    modal.style.display = \"none\";
                }

                // When the user clicks anywhere outside of the modal, close it
                window.onclick = function(event) {
                    if (event.target == modal) {
                        modal.style.display = \"none\";
                    }
                }
            </script> ";

        $vcfInfo{$sample}.= 
        "<tr>
            <td>   
                $plotDiv
            </td>
            <td><a href=\"$link\"> $tmp[0]:$tmp[1]-$End</a></td>
            <td>$svtype</td>
            <td>$size</td>
            <td>$GENE</td>
            <td>$ROIs</td>
            <td>$suportTag</td>
            <td>$info[0]</td>
            <td>$PE</td>
            <td>$BR</td>
            <td>$ASBR</td>
            <td>$RDRATIO</td>
            <td>$ZSCORE</td>
            <td>$SNR</td>
            <td>$GC</td>
            <td>$MAPQ</td>
            <td>
                $filterTag
            </td>
            <td>$AN</td>
            <td>$AC</td>
            <td>$AF</td>
            <td>$AFR_AC</td>
            <td>$AFR_AF</td>
            <td>$AMR_AC</td>
            <td>$AMR_AF</td>
            <td>$EAS_AC</td>
            <td>$EAS_AF</td>
            <td>$EUR_AC</td>
            <td>$EUR_AF</td>
            <td>$OTH_AC</td>
            <td>$OTH_AF</td>
        </tr>";

        if ($line =~/PASS/) {
            $varPass{$sv}++;
        }
        else {
            $varFilter{$sv}++;
        }
    }
    close IN;
    if (!$vcfInfo{$sample}) {
        $vcfInfo{$sample}.= 
        "<tr>
            <td></td>
            <td></td>
            <td></td>
            <td></td>
            <td></td>
            <td></td>
            <td></td>
            <td></td>
            <td></td>
            <td></td>
            <td></td>
            <td></td>
            <td></td>
            <td></td>
            <td></td>
            <td></td>
            <td></td>
            <td></td>
            <td></td>
            <td></td>
            <td></td>
            <td></td>
            <td></td>
            <td></td>
            <td></td>
            <td></td>
            <td></td>
            <td></td>
            <td></td>
            <td></td>
        </tr>";
    }
}

################## 
sub renderHeatmap {

    my $corrMatrix = "$::outDir/ON_TARGET/$::outName.cor_corrected.txt";

    my %hashCorr = ();
    my @rows = ();
    my $i = 0;
    my $xLab = ();
    open (CORR, "<", $corrMatrix) || die " ERROR: Unable to open $corrMatrix\n";
    while (my $line=<CORR>) {
        chomp $line;
        $line =~s/\"//g;
        my @tmp = split (/\t/, $line);
        if ($i == 0) {
            my @tmpLab;
            foreach my $sample (@tmp) {
               push @tmpLab,  "\'$sample\'";
            }
            $xLab = join (",", @tmpLab);
        }
        else {
            my $values = join (",", @tmp[1..@tmp-1]);
            push @rows, "[$values]";
        }
        $i++;
    }
    close CORR;

    my $yLab = $xLab;
    my @arrZ;
    foreach my $row (@rows) {
        push @arrZ, $row;
    }
    my $z = join (",", @arrZ);
    $z = "[$z]";
    my $x = "[$xLab]";
    my $y = $x;

    my $heatMap = " 
    <script>

    const letterWidth = 8;
    var data = [{
        z: $z,
        x: $x,
        y: $y,
        type: \'heatmap\',
        hoverongaps: false }];    
    var maxLabelLength = d3.max(data, d => d3.max(d.y, label => label.length)); 

    var layout = {
        margin:{
            l: maxLabelLength * letterWidth,
            b: maxLabelLength * letterWidth

        },

        title: {
            text:\'<b>Read depth correlation matrix</b>\',
            font: {
                family: \'Helvetica\',
                size: 18,
                color: \'#000000\'
            },
    }};

    Plotly.newPlot(\'plotly-div\', data, layout, {responsive: true});
    </script>\n";

    return $heatMap;
}
##################
sub printFooter {

    my $localtime = localtime();

    my $footer = "
    <footer style=\"position:absolute;margin-bottom:30px;margin-top:30px;width:100%;background-color:#F8F8F8;color:grey;text-align:center;\">
        <hr style=\"border: 0.5px solid #F8F8F8;\">
        <p>VCF generated by: <i>GRAPES $::version</i> on $localtime</p>
    </footer>    
    </body>
 </html>";
    return $footer;
}

##################
sub printBody {

    #my $humanVars = number2human( scalar keys %totalVars);
    my $heatMap       = renderHeatmap();
    my $summaryTable  = createSummaryTable();
    my $filteredTable = getFilteredRois();

    my $faInfo = "<i class=\"fas fa-info-circle fa-xs\"></i>";

    my $runName   = basename($::outDir);
    my $refGenome = basename($::genome);
    my $roiFile   = basename($::bed);
    my $nSamples  = scalar keys %::sampleHash;
    my $nROI      = `$::cat $::bed | $::wc -l`; chomp $nROI;
    my $filteredROIs = `$::cat $::onTargetFilteredROIs | $::wc -l`; chomp $filteredROIs;

    my $samplePages = "";
    foreach my $sample (natsort keys %::sampleHash) {
        $samplePages.= "<button class=\"tablinks\" style=\"border-bottom:1px solid white;font-size:13px;\" onclick=\"openCity(event, \'$sample\')\"><b>$sample</b></button>\n"
    }

    my $sampleResults = "";
    foreach my $sample (natsort keys %::sampleHash) {

        my $infoTable = $vcfInfo{$sample};

        my $scatterPlot = renderGenomeWideScatterPlot($sample);

        $sampleResults .= "
        <div class = \"container\" id=\"$sample\" style=\"display: none;margin-bottom:35px;\">
             <h3 style=\"margin-left:auto;margin-right:auto;text-align:center\">$sample SV Report</h3>
             <div id=\"$sample-scatter\" style=\"margin: 0 auto;\"></div>
             $scatterPlot
                <table id=\"\" class=\"display hover\" style=\"overflow:auto;margin: 0 auto;\">
                  <thead style=\"background-color:white;text-align:left;padding:3px;\">
                    <tr>
                        <th></th>
                        <th>Coordinates</th>
                        <th>Class</th>
                        <th>Size (bp)</th>
                        <th>ROI</th>
                        <th>#ROIs</th>
                        <th>Support</th>
                        <th>Breakpoint</th>
                        <th>PE
                         <div class=\"popup\" onclick=\"myFunction('Paired-End (PE) discordant reads are abnormally mapped reads that indicate the presence of medium/large SVs')\"> $faInfo
                         </div>
                         </th>
                        <th>SR
                         <div class=\"popup\" onclick=\"myFunction(\'Split-Reads (SR) are reads that directly span SV breakpoints, creating multi-part alignments\')\"> $faInfo
                         </div>
                        </th>
                        <th>AS
                         <div class=\"popup\" onclick=\"myFunction(\'Proportion of reads that can be assembled into contigs that support he SV\')\"> $faInfo
                         </div></th>
                        <th>Ratio</th>
                        <th>Z-score</th>
                        <th>SNR</th>
                        <th>\%GC</th>
                        <th>\%Map</th>
                        <th>Filter</th>
                        <th>AN</th>
                        <th>AC</th>
                        <th>AF</th>
                        <th>AFR_AC</th>
                        <th>AFR_AF</th>
                        <th>AMR_AC</th>
                        <th>AMR_AF</th>
                        <th>EAS_AC</th>
                        <th>EAS_AF</th>
                        <th>EUR_AC</th>
                        <th>EUR_AF</th>
                        <th>OTH_AC</th>
                        <th>OTH_AF</th>
                    </tr>
                  </thead>
                  <tbody>
                    $infoTable
                  </tbody>
                </table>
        </div>\n"; 
    }

    #my $circos = renderCircos();
    my $body = "
    <body style=\"margin-bottom:0px;\">
        <div class=\"tab\" style=\"margin-top:20px;margin-bottom:5px;border-bottom:1px solid #ddd;\">
            <button class=\"tablinks\" style=\"font-size:14px;color:white;background-color:#ff6633;\" onclick=\"openCity(event, 'summary')\" id=\"defaultOpen\"><b>$runName summary</b></button>
            $samplePages
        </div>
        <div class = \"container\" id=\"summary\" style=\"border-top:none\";>
                <div class=\"box\" style=\"width:92.5%;margin-bottom:0px;\">
                </div>
                <div class=\"box\" style=\"width:90%;\" >                    
                    <p style=\"font-size:16px;margin-right:2px;text-align:center;\"><b>Run:</b> $runName, <b>#Samples:</b> $nSamples, <b>Reference genome:</b> $refGenome, <b>ROI file:</b> $roiFile, <b>#ROIs:</b> $nROI, <b>Filtered ROIs:</b> $filteredROIs</p>
                    <h3 style=\";margin-top:35px;margin-left:auto;margin-right:auto;text-align:center\">Sequencing metrics</h3>
                        $summaryTable
                </div>
                <div class=\"box\" style=\"margin-left:20px;width:45%;margin-top:20px;\" >
                    <h3 style=\"margin-left:auto;margin-right:auto;text-align:center\">Filtered ROIs <div class=\"popup\" onclick=\"myFunction(\'<b>Filtering ROIs:</b> Low mappability, abnormal GC content (<25% or >80%) or small sizes can inflate the number of false positive CNVs.\')\"> $faInfo
                         </div></h3>
                        $filteredTable
                </div>
                <div class=\"box\" style=\"width:45%;margin-left:0px;margin-bottom:25px;margin-top:13px;margin-right:20px;\">
                    <div id=\"plotly-div\" style=\"margin-right:0px;\" ></div>
                    $heatMap
                </div>
        </div>
        $sampleResults    
        $accumulateModal
\n";


    return $body;
}

##################
sub printNavBar {

    my $navBar = "
    <ul>
        <li><a class = \"active\" href=\"#\"><b>Structural Variation</b></a></li>
        <li style=\"float:left\"><a href=\"#\"><b><i>reportHTML::Hybrid capture</i></b></i></a></li>

        <li style=\"float:right\"><a href=\"https://github.com/bdolmo/GRAPES\">Details</a></li>
        <li style=\"float:right\"><a href=\"mailto:bdelolmo\@gencardio.com\">Contact</a></li>
    </ul>";
    return $navBar;
}        
#<link rel=\"stylesheet\" type=\"text/css\" href=\"css/style.css\"\">
#<link rel=\"stylesheet\" type=\"text/css\" href=\"https://drive.google.com/uc?export=view&id=1F9vrJUlmU3Y1yLEu6MIwuvyTsT_unl_G\">

##################

sub renderGeneTrack2 {

    my $chr   = shift;
    my $start = shift;
    my $end   = shift;
    my $tag   = shift;

    my $chromo = $chr;
    $chromo =~s/chr//;

    my $json = "genes_" . $chromo . "_" . "$start-$end" . ".json?";

    my $geneTrack = "
    <script>

        var genome = tnt.board.genome()
            .species(\"human\")
            .chr($chromo)
            .width(1000)
            .min_coord (new Promise (function (res) {
                res($start);
            }))
            .max_coord (new Promise (function (res) {
                res($end);
            }));
        current_height = 200;

        var gene_track = tnt.board.track()
            .height(200)
            .color(\"white\")
            .display(tnt.board.track.feature.genome.gene()
                .color(\"#0000FF\")
            )
            .data(tnt.board.track.data.genome.gene());

        gene_track.display().layout()
            .fixed_slot_type(\"expanded\")
            .keep_slots(false)
            .on_layout_run (function (types, current) {
                var needed_height = types.expanded.needed_slots * types.expanded.slot_height;
                if (needed_height !== current_height) {
                    current_height = needed_height;
                    gene_track.height(needed_height);
                    genome.tracks(genome.tracks());
                }
            });

        var transcript_track = tnt.board.track()
            .height(200)
            .color(\"white\")
            .display(tnt.board.track.feature.genome.transcript()
                .color(\"#fff\")
            )
            .data(tnt.board.track.feature.genome.transcript());
            
        genome
            .zoom_in(100)
            .add_track(gene_track);

        genome(document.getElementById(\"tnt-$start$end\"));
        genome.start();
    </script>\n";

 return $geneTrack;
}

##################
sub renderGeneTrck {

    my $chr   = shift;
    my $start = shift;
    my $end   = shift;
    my $tag   = shift;

    my $chromo = $chr;
    $chromo =~s/chr//;

    my $json = "genes_" . $chromo . "_" . "$start-$end" . ".json?";

    my $geneTrack = "
    <script>
    var apiBase, data_sources;
    apiBase = \"https://portaldev.sph.umich.edu/api/v1/\";
    data_sources = new LocusZoom.DataSources()
        .add(\"gene\", [\"GeneLZ\", { url: apiBase + \"annotation/genes/\", params: { build: \'GRCh37\' } }]);

    var layoutGT = {
        width: 1200,
        height: 200,
        panels: [
                LocusZoom.Layouts.get(\"panel\", \"genes\", {namespace: { \"gene\":\"gene\"} })
        ]
        }
    var plot = LocusZoom.populate(\"#lz-$start$end\", data_sources, layoutGT); 

    </script>\n";

 return $geneTrack;
}

##################
sub printHead {

    my $title = shift;
    my $jsDatatables = "";

    $jsDatatables .= "\$(document).ready(function() {
        \$(\'table.display\').DataTable( {
            fixedColumns: true,
            \"order\": [[ 0, \"desc\" ]],
           \"scrollX\": true
        } ).columns.adjust().draw();
    } );
    \n";
#https://drive.google.com/file/d/1yJIL462ONHzawn-af6_QQOkJ5vdwbXSR/view?usp=sharing
#https://drive.google.com/file/d/1C_BQNQGF5ZoluU9nqPZVbHVnLvKOTq3O/view?usp=sharing
#https://drive.google.com/file/d/1S7X5dxrykQR8Wy6xUc1I3XOCs-_0yGH6/view?usp=sharing
#https://drive.google.com/file/d/1wgosiWnRcF0duu1LaHBqi4HH8fTKddlJ/view?usp=sharing

#https://drive.google.com/file/d/1hKnFy8VZZ5L52UbNkLw4sNxku1y9btX0/view?usp=sharing
#https://drive.google.com/file/d/1FPZfRO894DSox23_YmjQCxiVTnYkbbyg/view?usp=sharing

    my $head = "
    <!DOCTYPE html>
    <html style=\"position:relative;min-height:100vh;\">    
        <head>
        <title>$title</title>
        <link rel=\"stylesheet\" type=\"text/css\" href=\"https://drive.google.com/uc?export=view&id=1leVVXa8LPMM_IgVApwf1nbTv03jcToBg\">
        <link rel=\"stylesheet\" href=\"https://cdnjs.cloudflare.com/ajax/libs/font-awesome/4.7.0/css/font-awesome.min.css\">
        <link rel=\"stylesheet\" type=\"text/css\" href=\"https://cdn.datatables.net/1.10.19/css/jquery.dataTables.min.css\">
        <link rel=\"stylesheet\" type=\"text/css\" href=\"https://cdn.jsdelivr.net/npm/locuszoom\@0.10.0/dist/locuszoom.css\">

        <link rel=\"stylesheet\" href=\"https://use.fontawesome.com/releases/v5.5.0/css/all.css\" integrity=\"sha384-B4dIYHKNBt8Bc12p+WXckhzcICo0wtJAoU8YZTY5qE0Id1GSseTk6S+L3BlXeVIU\" crossorigin=\"anonymous\">
        <script src=\"http://code.jquery.com/jquery-1.11.3.min.js\"></script>
        <script src=\"https://d3js.org/d3.v3.min.js\"></script>        
        <script src=\"https://cdnjs.cloudflare.com/ajax/libs/jszip/3.1.3/jszip.min.js\"></script>
        <link href=\"https://nightly.datatables.net/css/jquery.dataTables.css\" rel=\"stylesheet\" type=\"text/css\"/>
        <script src=\"https://nightly.datatables.net/js/jquery.dataTables.js\"></script>
        <script type=\"text/javascript\" src=\"https://cdnjs.cloudflare.com/ajax/libs/Chart.js/2.7.1/Chart.bundle.min.js\"></script>
        <script src=\"https://cdn.datatables.net/1.10.20/js/jquery.dataTables.min.js\"></script>
        <script src=\"https://cdn.datatables.net/fixedheader/3.1.6/js/dataTables.fixedHeader.min.js\"></script>
        <script src=\"https://cdn.plot.ly/plotly-latest.min.js\"></script>
        <script src=\"https://cdn.datatables.net/buttons/1.5.2/js/dataTables.buttons.min.js\"></script>
        <script src=\"https://cdnjs.cloudflare.com/ajax/libs/jszip/3.1.3/jszip.min.js\"></script>
        <script src=\"https://cdnjs.cloudflare.com/ajax/libs/pdfmake/0.1.36/pdfmake.min.js\"></script>
        <script src=\"https://cdnjs.cloudflare.com/ajax/libs/pdfmake/0.1.36/vfs_fonts.js\"></script>
        <script src=\"https://cdn.datatables.net/buttons/1.5.2/js/buttons.html5.min.js\"></script>
        <script src=\"https://cdn.jsdelivr.net/npm/igv\@2.3.5/dist/igv.min.js\"></script>
        <script src=\"https://cdnjs.cloudflare.com/ajax/libs/d3/4.5.0/d3.js\"></script>
        <script src=\"https://cdnjs.cloudflare.com/ajax/libs/d3-queue/3.0.3/d3-queue.js\"></script>
        <script src=\"https://cdn.jsdelivr.net/gh/nicgirault/circosJS\@v2/dist/circos.js\"></script>
        <script src=\"https://cdn.jsdelivr.net/npm/sweetalert2\@9\"></script>
        <script src=\"https://code.jquery.com/ui/1.12.1/jquery-ui.js\"></script>
        <script src=\"https://cdn.jsdelivr.net/npm/ideogram\@1.15.0/dist/js/ideogram.min.js\"></script>
    
        <script src=\"https://tntvis.github.io/tnt.genome/build/tnt.genome.min.js\"></script>
        <script src=\"https://tntvis.github.io/tnt.genome/build/tnt.genome.js\" charset=\"utf-8\"></script>

        <script src=\"https://cdn.jsdelivr.net/npm/locuszoom\@0.10.0/dist/locuszoom.vendor.min.js\"></script>
        <script src=\"https://cdn.jsdelivr.net/npm/locuszoom\@0.10.0/dist/locuszoom.app.min.js\"></script>
        <script src=\"https://cdnjs.cloudflare.com/ajax/libs/jquery-modal/0.9.1/jquery.modal.min.js\"></script>
        <link rel=\"stylesheet\" href=\"https://cdnjs.cloudflare.com/ajax/libs/jquery-modal/0.9.1/jquery.modal.min.css\">
        <script>
            $jsDatatables
        </script>
        <script>

  \$( function() {
    \$( document ).tooltip();
  } );

    \"use strict\";

    var igvBrowser;

    document.addEventListener(\"DOMContentLoaded\", function () {
        var igvDiv = document.getElementById(\"igv-div\");
        var options = {
            locus: \"chr8:128,747,267-128,754,546\",
            genome: \"hg19\"
        };

        igv.createBrowser(igvDiv, options)
            .then(function (b) {
                igvBrowser = b;
            })
    })

    function load() {
        var fileWidget = document.getElementById(\"fileWidget\");
        var files = fileWidget.files;

        var fileTxt = \"<ul>\";
        for (let file of files) {
            fileTxt += \"<li>\" + file.name + \"</li>\";
        }
        fileTxt += \"</ul>\";
        document.getElementById(\"fileNameDiv\").innerHTML = fileTxt;

        // Find BAM files and cache index files.  Note there are 2 index naming conventions, .bam.bai and .bai
        // This scheme catches both.
        var bamFiles = [];
        var indexFiles = {};

        for (let file of files) {
            if (file.name.endsWith(\".bam\")) {
                bamFiles.push(file);
            }
            else if (file.name.endsWith(\".bai\")) {
                var key = getKey(file.name);
                indexFiles[key] = file;
            }
            else {
                alert(\"Unsupported file type: \" + file.name);
            }
        }

        // Create track objects
        var trackConfigs = [];

        for (let file of bamFiles) {

            var key = getKey(file.name);
            var indexFile = indexFiles[key];
            if (indexFile) {
                trackConfigs.push({
                    name: file.name,
                    type: \"alignment\",
                    format: \"bam\",
                    url: file,
                    indexURL: indexFile
                })
            }
            else {
                alert(\"No index file for: \" + file.name);
            }
        }

        if (trackConfigs.length > 0) {
            igvBrowser.loadTrackList(trackConfigs);
        }

        function getKey(filename) {

            var idx = filename.indexOf(\".\");
            if (idx < 0) {
                console.error(\"File with no extension: \" + filename);
            }
            else {
                return filename.substring(0, idx);
            }
        }

        // igv.browser.loadSession(sessionFile)
        //     .catch(function (error) {
        //         alert(\"Error loading session file\");
        //     })

    }
     </script>
     <script>
        function showPlot(myDiv) {
            var x = document.getElementById(myDiv);
            if (x.style.display == \"none\") {
                x.style.display = \"block\";
            } else {
                x.style.display = \"none\";
            }
        }

        function myFunction(mypopup) {
            Swal.fire({
                icon: 'info',
                text: mypopup
            })
        }

        function openCity(evt, section) {
            var i, container, tablinks;
            container = document.getElementsByClassName(\"container\");
            for (i = 0; i < container.length; i++) {
                container[i].style.display = \"none\";
            }
            tablinks = document.getElementsByClassName(\"tablinks\");
            for (i = 0; i < tablinks.length; i++) {
                tablinks[i].className = tablinks[i].className.replace(\" active\", \"\");
                tablinks[i].style.color=\"black\";
                if (i == 0) {
                    tablinks[0].style.color=\"white\";
                    tablinks[0].style.background=\"#ff9f80\";
                }
                else {
                    tablinks[i].style.color=\"black\";
                    tablinks[i].style.background=\"white\";
                    tablinks[i].style.borderColor=\"white\";    
  
                }
            }
            if (section == 'summary') {
               document.getElementById(section).style.display =\"flex\";
               document.getElementById(section).style.color =\"black\";
               evt.currentTarget.style.background=\"#ff6633\";
            }
            else {
               document.getElementById(section).style.display =\"block\";
               document.getElementById(section).style.color =\"black\";
               evt.currentTarget.style.color=\"black\";
               evt.currentTarget.style.borderColor=\"#ddd\";
               evt.currentTarget.style.background=\"#ddd\";

            }
            evt.currentTarget.className += \" active\";
        }
            // Get the element with id=\"defaultOpen\" and click on it            
            document.getElementById(\"defaultOpen\").style.background = \"grey\";
            document.getElementById(\"defaultOpen\").click();
        </script>
        </head>";

    return $head;
}

sub createRectangle {

my $start = shift;
my $end   = shift;
my $div   = shift;
my $hash  = shift;

my %HoG = %$hash;

my $n = 0;
my $shapes = "";
my $traces = "";
my @xarr;
my @yarr;
my @geneText;

foreach my $coordinate ( natsort keys %HoG ) {
    
    my ($chr, $Start, $End) = split (/\t/, $coordinate);
    my $sumPos = int (($End-$Start)/2); 
    my $position = $Start+$sumPos;

    my $is_even = $n % 2 == 0;
    my $yvalue;
    my $y0;
    my $y1;

    push @xarr, $position;    
    if ($is_even) {
        $yvalue = "140";
        $y0 = 200;
        $y1 = 250;
    }
    else {
        $yvalue = "390";
        $y0 = 300;
        $y1 = 360;  
    }
    push @yarr, $yvalue;

    push @geneText, "\'$HoG{$coordinate}{NAME}\'";

    $shapes .= "
        {
            type: \'rect\',
            x0: $Start,
            y0: $y0,
            x1: $End,
            y1: $y1,
            line: {
                color: \'#0000FF\',
                width: 1.5
            },
            fillcolor: \'#0000FF\'
        },
    \n";
    $n++;
}
my $x = join (",", @xarr);
my $y = join (",", @yarr);
my $gt= join (",", @geneText);

my $rectangle = "
<script>
    var trace1 = {
        x: [$x],
        y: [$y],
        text: [$gt],
        mode: \'text\'
    };

    var layout$start$end = {
        margin: {
            t: 4,
            b: 4
        },
        xaxis: {
            range: [$start, $end],
            showgrid: false,
            showgrid: false,
            zeroline: true,
            showline: true,
            linecolor: \'black\',
            linewidth: 2,
            mirror: \'ticks\',
            showticklabels: false
        },
        yaxis: {
            range: [0, 500],
            showgrid: false,
            zeroline: true,
            showline: true,
            linecolor: \'black\',
            linewidth: 2,
            mirror: \'ticks\',
            showticklabels: false,
            title: {
                text: \'Genes\',
                font: {
                    family: \'Helvetica\',
                    size: 16,
                    color: \'black\'
                }
            },
        },
        width: 1200,
        height: 100,
        
        shapes: [
            $shapes
        ]
    };

var data = [trace1];

Plotly.newPlot(\'$div-rectangle\', data, layout$start$end, {responsive: true});
</script>
\n";
return $rectangle;

}

########################################
sub getAvailablePlot {
    my $dataFile = shift;
    my $plotFile;
    my $class;

    if (-e "$::outDir/ON_TARGET/PLOT_DATA/$dataFile") {
        $plotFile = "$::outDir/ON_TARGET/PLOT_DATA/$dataFile";
        $class = "on-target";
    }
    elsif (-e "$::outDir/OFF_TARGET/PLOT_DATA/$dataFile") {
        $plotFile = "$::outDir/OFF_TARGET/PLOT_DATA/$dataFile";
        $class = "off-target";
    }
    elsif (-e "$::outDir/PLOT_DATA/$dataFile") {
        $plotFile = "$::outDir/PLOT_DATA/$dataFile";
        $class = "mixed";
    }

    return $plotFile, $class;
}


return 1;


