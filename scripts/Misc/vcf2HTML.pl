#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;
use File::Basename;

my $inputVCF = $ARGV[0]; 
my $bamFile  = $ARGV[1];

my $sampleName = basename($bamFile);
$sampleName=~s/.bam//;

if ( !$inputVCF || !$bamFile) {
    print "Usage: $0 <file.vcf> <file.bam>\n";
    exit;
}

my $bamDir = dirname($bamFile);
my $bamName = basename($bamFile);
my %var = ();

my %varPass = ();
my %varFilter = ();
my $totalVars = 0;

my %evidence = (
    'RP' => 0,
    'SR/AS'=> 0,
    'RP+SR/AS' => 0,
    'RD'=> 0,
);

my %svSize = (
    '50-100'  => 0,
    '100-150' => 0,
    '150-200' => 0,
    '200-250' => 0,
    '250-300' => 0,
    '300-350' => 0,
    '350-400' => 0,
    '400-450' => 0,
    '450-500' => 0,
    '500-600' => 0,
    '600-700' => 0,
    '700-800' => 0,
    '800-900' => 0,
    '900-1000' => 0,
    '1000-2000' => 0,
    '2000-5000' => 0,
    '5000-10000' => 0,
    '10000-50000' => 0,
    '50000-100000' => 0,
    '100000' => 0
);

#Goal is to print an STATIC html report
my $vcfInfo = readVCF();

my $head = printHead("VCF");
print "$head\n";

my $navBar = printNavBar();
print "$navBar\n";

my $body = printBody();
print "$body\n";

my $footer = printFooter();
print "$footer\n";

#############
sub readVCF {

    my $vcfInfo = "";
    open (IN, "<", $inputVCF) || die " ERROR: Unable to open $inputVCF\n";
    while (my $line=<IN>) {
        chomp $line;
        next if $line =~/^#/;
        my @tmp = split (/\t/, $line);
        my $sv = $tmp[4];
        $sv =~s/[<|>]//g;
        $var{$sv}++;
        $totalVars++;

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

        my $size = $End-$tmp[1];

        if ($size > 50 && $size <= 100) {
            $svSize{"50-100"}++;
        }
        if ($size > 100 && $size <= 150) {
            $svSize{"100-150"}++;
        }
        if ($size > 150 && $size <= 200) {
            $svSize{"150-200"}++;
        }
        if ($size > 200 && $size <= 250) {
            $svSize{"200-250"}++;
        }
        if ($size > 250 && $size <= 300) {
            $svSize{"250-300"}++;
        }
        if ($size > 300 && $size <= 350) {
            $svSize{"300-350"}++;
        }
        if ($size > 350 && $size <= 400) {
            $svSize{"350-400"}++;
        }
        if ($size > 400 && $size <= 450) {
            $svSize{"400-450"}++;
        }
        if ($size > 450 && $size <= 500) {
            $svSize{"450-500"}++;
        }
        if ($size > 500 && $size <= 600) {
            $svSize{"500-600"}++;
        }
        if ($size > 600 && $size <= 700) {
            $svSize{"600-700"}++;
        }
        if ($size > 700 && $size <= 800) {
            $svSize{"700-800"}++;
        }
        if ($size > 800 && $size <= 900) {
            $svSize{"800-900"}++;
        }
        if ($size > 900 && $size <= 1000) {
            $svSize{"900-1000"}++;
        }
        if ($size > 1000 && $size <= 2000) {
            $svSize{"1000-2000"}++;
        }
        if ($size > 2000 && $size <= 5000) {
            $svSize{"2000-5000"}++;
        }
        if ($size > 5000 && $size <= 10000) {
            $svSize{"5000-10000"}++;
        }
        if ($size > 10000 && $size <= 50000) {
            $svSize{"10000-50000"}++;
        }
        if ($size > 50000 && $size <= 100000) {
            $svSize{"50000-100000"}++;
        }
        if ($size >= 100000) {
            $svSize{"100000"}++;
        }

        my $filterTag;
        if ($tmp[6] eq 'PASS') {
            $filterTag = "<input class=\"edit\" type=\"submit\" value=$tmp[6]>";         
        }
        else {
            $filterTag = "<input class=\"delete\" type=\"submit\" value=$tmp[6]>";
        }

        my $suportTag;
        if ( $line =~/BR=0;ASBR=0;PE=0;PPE=0/) {
            $suportTag = "RD";
        }
        else {
            if ( $line =~/IMPRECISE/) {
                $suportTag = "RP";
            }
             else {
                if ( $line =~/PE=0/) {
                    $suportTag = "SR";
                }
                else {
                    $suportTag = "RP+SR";
                }
             }
        }

        $vcfInfo.= 
        "<tr>
            <td>$tmp[0]:$tmp[1]-$End</td>
            <td>$svtype</td>
            <td>$suportTag</td>
            <td>$info[0]</td>
            <td>$BR</td>
            <td>$PE</td>
            <td>$ASBR</td>
            <td>
                $filterTag
            </td>
        </tr>";

        if ($line =~/PASS/) {
            $varPass{$sv}++;

            if ( $line =~/BR=0;ASBR=0;PE=0;PPE=0/) {
                $evidence{"RD"}++;
            }
            else {
                if ( $line =~/IMPRECISE/) {
                    $evidence{"RP"}++;
                }
                else {
                    if ( $line =~/PE=0/) {
                        $evidence{"SR/AS"}++;
                    }
                    else {
                        $evidence{"RP+SR/AS"}++;
                    }
                }
            }
        }
        else {
            $varFilter{$sv}++;
        }
    }
    close IN;
    return $vcfInfo;
}

##################
sub renderCircos {

    my $circos = "
    <script>
        var myCircos = new Circos({
            container: '#chart',
            width: 800,
            height: 800,
        });

        var configuration = {
            innerRadius: 250,
            outerRadius: 300,
            cornerRadius: 10,
            gap: 0.04, // in radian
            labels: {
                display: true,
                position: 'center',
                size: '14px',
                color: '#000000',
                radialOffset: 20,
            },
            ticks: {
                display: true,
                color: 'grey',
                spacing: 10000000,
                labels: true,
                labelSpacing: 10,
                labelSuffix: 'Mb',
                labelDenominator: 1000000,
                labelDisplay0: true,
                labelSize: '10px',
                labelColor: '#000000',
                labelFont: 'default',
                majorSpacing: 5,
                size: {
                minor: 2,
                major: 5,
                }
            },
            events: {}
        };
        var data = [
            { len: 31, color: \"#8dd3c7\", label: \"January\", id: \"january\" },
            { len: 28, color: \"#ffffb3\", label: \"February\", id: \"february\" },
            { len: 31, color: \"#bebada\", label: \"March\", id: \"march\" },
            { len: 30, color: \"#fb8072\", label: \"April\", id: \"april\" },
            { len: 31, color: \"#80b1d3\", label: \"May\", id: \"may\" },
            { len: 30, color: \"#fdb462\", label: \"June\", id: \"june\" },
            { len: 31, color: \"#b3de69\", label: \"July\", id: \"july\" },
            { len: 31, color: \"#fccde5\", label: \"August\", id: \"august\" },
            { len: 30, color: \"#d9d9d9\", label: \"September\", id: \"september\" },
            { len: 31, color: \"#bc80bd\", label: \"October\", id: \"october\" },
            { len: 30, color: \"#ccebc5\", label: \"November\", id: \"november\" },
            { len: 31, color: \"#ffed6f\", label: \"December\", id: \"december\" }
        ];
        myCircos.layout(data, configuration);
        myCircos.render();
    </script>
 \n";

}

##################
sub renderPie {

    my $n = shift;

    my @labelTmp;
    my @evidenceTmp;
    foreach my $ev (sort keys %evidence) {
        push @labelTmp, "\'$ev\'";
        push @evidenceTmp, "\'$evidence{$ev}\'";
    }
    my $labels   = join (",", @labelTmp);
    my $data     = join (",", @evidenceTmp);

    my $pie = "
<canvas id=\'pie-chart\' width=\'400\' height=\'400\'></canvas>
<script>
    new Chart(document.getElementById(\'pie-chart\'), {
        type: 'pie',
        data: {
        labels: [$labels],
        datasets: [{
            label: \'Population (millions)\',
            backgroundColor: [\'#123456\', \'#28819f\',\'#3eb5c1\',\'#abede0\'],
            data: [$data]
        }]
        },
        options: {
        }
    });
</script>";
    return $pie;
}

##################
sub renderHistogram {

    my $n = shift;

    my @labTmp;
    my @sizeTmp;

    foreach my $size ( sort { $a <=> $b } keys %svSize) {

        my $resize = $size;
        if ($size eq 100000) {
            push @labTmp, "\'>=100000\'";
            push @sizeTmp, "\'$svSize{$size}\'";
        }
        else {
            push @labTmp, "\'$size\'";
            push @sizeTmp, "\'$svSize{$size}\'";
        }
    }
    my $labels   = join (",", @labTmp);
    my $data     = join (",", @sizeTmp);


    my $chart = "
    <canvas id=\'myChart$n\' width=\'400\' height=\'400\'></canvas>
    <script>
    var ctx = document.getElementById(\'myChart$n\')
    var chart = new Chart(ctx, {
        // The type of chart we want to create
        type: 'bar',

        // The data for our dataset
        data: {
            labels: [$labels],
            datasets: [{
                backgroundColor: 'rgb(0,0,255)',
                borderColor: 'rgb(0,0,255)',
                data: [$data]
            }]
        },

    });

    </script>
    ";
    return $chart;

}

##################
sub renderChart {

    my $n = shift;

    my @labTmp;
    my @passTmp;
    my @filteredTmp;

    foreach my $sv (sort keys %var) {
        push @labTmp, "\'$sv\'";
        push @passTmp, "\'$varPass{$sv}\'";

        $varFilter{$sv} = 0 if !$varFilter{$sv};
        push @filteredTmp, "\'$varFilter{$sv}\'";
    }
    my $labels   = join (",", @labTmp);
    my $pass     = join (",", @passTmp);
    my $filtered = join (",", @filteredTmp);

    my $chart = "
    <canvas id=\'myChart$n\' width=\'400\' height=\'400\'></canvas>
    <script>
    var ctx = document.getElementById(\'myChart$n\')
    var myChart$n = new Chart(ctx, {
        type: \'bar\',
        data: { 
            labels: [$labels],
            datasets: [
            {
                label: \'PASS\',
                data: [$pass],
                backgroundColor: [
                    \'rgba(220, 20, 60, 1)\',
                    \'rgba(0,128,0, 1)\',
                    \'rgba(255,165,0, 1)\',
                ],
                borderColor: [
                    \'rgba(220, 20, 60, 1)\',
                    \'rgba(0,128,0, 1)\',
                    \'rgba(255,165,0, 1)\',

                ],
                borderWidth: 1
            },
            {
                label: \'Filtered\',
                data: [$filtered],
                backgroundColor: [
                    \'rgba(220, 20, 60, 0.3)\',
                    \'rgba(0,128,0, 0.3)\',
                    \'rgba(255,165,0, 0.3)\',
                ],
                borderColor: [
                    \'rgba(220, 20, 60, 0.3)\',
                    \'rgba(0,128,0, 0.3)\',
                    \'rgba(255,165,0, 9.3)\',
                ],
                borderWidth: 1   
            }
        ]
        },
    options: {
        scales: {
        yAxes: [{
            stacked: true,
            ticks: {
            beginAtZero: true
            }
        }],
        xAxes: [{
            stacked: true,
            ticks: {
            beginAtZero: true
            }
        }]

        }
    }
    });
    </script>
    ";
    return $chart;

}

##################
sub printFooter {

    my $footer = "
    <div class=\"footer\">
            <h3></h3>
        </div>
 </html>";
    return $footer;
}

##################
sub printBody {

    my $chart1 = renderChart("1");
    my $pie2   = renderPie("2");
    my $chart3 = renderHistogram("3");

               #  <div id=\"chart\" style=\"overflow: visible; display=block\">
               #   <h3 style=\"text-align:center\">Circos</h3>

               #  $circos 
               # </div>

    my $circos = renderCircos();
    my $body = "
    <body>
        <div class=\"tab\">
            <button class=\"tablinks\" onclick=\"openCity(event, 'summary')\" id=\"defaultOpen\"><b>Summary</b></button>
            <button class=\"tablinks\" onclick=\"openCity(event, 'VCF')\"><b>Variants</b></button>
        </div>
        <div class = \"container\" id=\"summary\">
                <div class=\"box\">
                    <h3 style=\"text-align:center\">SV size distribution</h3>
                    $chart3
                </div>
                <div class=\"box\">
                    <h3 style=\"text-align:center\">SV signal contribution</h3>
                    $pie2
                </div>
                <div class=\"box\">
                    <h3 style=\"text-align:center\">PASS / LowQual</h3>
                    $chart1
                </div>

        </div>
        <div class = \"container\" id=\"VCF\" style=\"display: none\">
            <div style=\"float: center;margin-left:auto;width:90%; margin-right:auto;\">BAM file and index 
                <input id=\"fileWidget\" type=\"file\" multiple=\"true\" accept=\".bam,.bai\" onchange=\"load()\"/>
            </div>

    <p>
    <div id=\"fileNameDiv\" style=\"display:none\"></div>
    </p>                
    <div id=\"igv-div\" style=\"padding-top: 10px;width:90%; margin-bottom:20px; margin-left:auto; margin-left:auto; margin-right:auto; padding-bottom: 10px; border:1px solid lightgray\"></div>
                <table id=\"example\" class=\"display hover\" style=\"margin-top:10px;\">
                  <thead>
                    <tr>
                        <th>Location</td>
                        <th>Class</td>
                        <th>Support</td>
                        <th>Breakpoint</td>
                        <th>Split reads</td>
                        <th>Discordant pairs</td>
                        <th>Assembled reads</td>
                        <th>Filter</td>
                    </tr>
                  </thead>
                  <tbody>
                    $vcfInfo
                  </tbody>
                </table>
        </div>
    </body>";
    return $body;
}

##################
sub printNavBar {

    my $navBar = "
    <ul>
        <li><a class = \"active\" href=\"home\"><b>Structural Variation</b></a></li>
        <li style=\"float:left\"><a href=\"#\"><b>display VCF</b></a></li>
        <li style=\"float:right\"><a href=\"#\">Details</a></li>
        <li style=\"float:right\"><a href=\"#\">Contact</a></li>
    </ul>";
    return $navBar;
}

##################
sub printHead {

    my $title = shift;
    my $head = "
    <!DOCTYPE html>
    <html>    
        <head>
        <title>$title</title> 

        <link rel=\"stylesheet\" type=\"text/css\" href=\"css/style.css\"\">
        <link rel=\"stylesheet\" type=\"text/css\" href=\"https://cdn.datatables.net/1.10.19/css/jquery.dataTables.min.css\">
        <link rel=\"stylesheet\" href=\"https://use.fontawesome.com/releases/v5.5.0/css/all.css\" integrity=\"sha384-B4dIYHKNBt8Bc12p+WXckhzcICo0wtJAoU8YZTY5qE0Id1GSseTk6S+L3BlXeVIU\" crossorigin=\"anonymous\">
        <script src=\"http://code.jquery.com/jquery-1.11.3.min.js\"></script>
        <script src=\"https://cdnjs.cloudflare.com/ajax/libs/jszip/3.1.3/jszip.min.js\"></script>
        <link href=\"https://nightly.datatables.net/css/jquery.dataTables.css\" rel=\"stylesheet\" type=\"text/css\"/>
        <script src=\"https://nightly.datatables.net/js/jquery.dataTables.js\"></script>
        <script type=\"text/javascript\" src=\"https://cdnjs.cloudflare.com/ajax/libs/Chart.js/2.7.1/Chart.bundle.min.js\"></script>
        <script src=\"https://cdn.datatables.net/1.10.19/js/jquery.dataTables.min.js\"></script>
        <script src=\"https://cdn.datatables.net/buttons/1.5.2/js/dataTables.buttons.min.js\"></script>
        <script src=\"https://cdnjs.cloudflare.com/ajax/libs/jszip/3.1.3/jszip.min.js\"></script>
        <script src=\"https://cdnjs.cloudflare.com/ajax/libs/pdfmake/0.1.36/pdfmake.min.js\"></script>
        <script src=\"https://cdnjs.cloudflare.com/ajax/libs/pdfmake/0.1.36/vfs_fonts.js\"></script>
        <script src=\"https://cdn.datatables.net/buttons/1.5.2/js/buttons.html5.min.js\"></script>
        <script src=\"https://cdn.jsdelivr.net/npm/igv\@2.3.5/dist/igv.min.js\"></script>
        <script src=\"https://cdnjs.cloudflare.com/ajax/libs/d3/4.5.0/d3.js\"></script>
        <script src=\"https://cdnjs.cloudflare.com/ajax/libs/d3-queue/3.0.3/d3-queue.js\"></script>
        <script src=\"https://cdn.jsdelivr.net/gh/nicgirault/circosJS\@v2/dist/circos.js\"></script>        
        <script>
         \$(document).ready(function() {
        \$(\'#example\').DataTable( {
            dom: \'Bfrtip\',
            buttons: [
                \'copyHtml5\',
                \'excelHtml5\',
                \'csvHtml5\',
                \'pdfHtml5\'
            ]
        } );
        } );
        </script>
        <script>

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
        function openCity(evt, section) {
            var i, container, tablinks;
            container = document.getElementsByClassName(\"container\");
            for (i = 0; i < container.length; i++) {
                container[i].style.display = \"none\";
            }
            tablinks = document.getElementsByClassName(\"tablinks\");
            for (i = 0; i < tablinks.length; i++) {
                tablinks[i].className = tablinks[i].className.replace(\" active\", \"\");
            }
            if (section == 'summary') {
               document.getElementById(section).style.display =\"flex\";
            }
            if (section == 'VCF') {
               document.getElementById(section).style.display =\"block\";
            }
            evt.currentTarget.className += \" active\";
            }
            // Get the element with id=\"defaultOpen\" and click on it
            document.getElementById(\"defaultOpen\").click();
        </script>
        </head>";

    return $head;
}




