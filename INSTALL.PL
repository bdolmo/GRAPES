#!/usr/bin/env perl


use strict;
use warnings;
use Cwd;
use Term::ANSIColor qw(:constants);
use Parallel::ForkManager;
use Statistics::Descriptive;
use Config;

my $make  = `which make`;
chomp $make;
if (!$make) {
	print "ERROR: make command not found\n";
	exit;
}
my $date = localtime;

print " ## $date ##\n";
print "\n ## Wellcome to GRAPES installer ##\n";

my $currentDir = cwd();

if (!-e "$currentDir/mappability") {
  mkdir "$currentDir/mappability";
}

if (!-e "$currentDir/db/gnomad") {
  mkdir "$currentDir/db/gnomad";
}

my $cut  = `which cut`;
chomp $cut;
if (!$cut) {
	print "ERROR: cut command not found\n";
	exit;
}

my $cat  = `which cat`;
chomp $cat;
if (!$cat) {
	print "ERROR: cat command not found\n";
		exit;
}

my $wc  = `which wc`;
chomp $wc;
if (!$wc) {
	print "ERROR: wc command not found\n";
	exit;
}

my $awk  = `which awk`;
chomp $awk;
if (!$awk) {
	print "ERROR: awk command not found\n";
	exit;	
}

my $sed  = `which sed`;
chomp $sed;
if (!$sed) {
	print "ERROR: sed command not found\n";
	exit;
}

my $paste  = `which paste`;
chomp $paste;
if (!$paste) {
	print "ERROR: paste command not found\n";
	exit;
}

my $head  = `which head`;
chomp $head;
if (!$head) {
	print "ERROR: head command not found\n";
	exit;
}

my $uniq  = `which uniq`;
chomp $uniq;
if (!$uniq) {
	print "ERROR: uniq command not found\n";
	exit;
}
 
my $mv  = `which mv`;
chomp $mv;
if (!$mv) {
	print "ERROR: mv command not found\n";
	exit;
}

my $tail  = `which tail`;
chomp $tail;
if (!$tail) {
	print "ERROR: mv command not found\n";
	exit;
}

my $grep;
my $sort;
if ($Config{osname} =~/darwin/) {
 $grep  = `which egrep`; chomp $grep;
 $sort = `which gsort`; chomp $sort;
}
else {
 $grep  = `which grep`; chomp $grep;
 $sort = `which sort`; chomp $sort;
}

if (!$grep) {
	print " ERROR: grep command was not found\n";
	exit;
}
if (!$sort) {
	print " ERROR: sort command was not found\n";
	exit;
}

my $Rscript  = `which Rscript`;
chomp $Rscript;

if (!$Rscript) {
    print " ERROR: Rscript was not found on PATH\n"; 
    exit;
}

my $tabix = `which tabix`; 
chomp $tabix;

if (!$tabix) {
	print " ERROR: tabix was not found on PATH\n"; exit;
}

my $bedtools = `which bedtools`;
chomp $bedtools;

if (!$bedtools) {
    print " ERROR: BEDtools was not found on PATH\n";
	  exit;
}
else {
    my $bedtools_version = `$bedtools --version | $cut -d \" \" -f 2 | $awk '{ split (\$0, a, \".\");  print a[2] }'`;
    chomp $bedtools_version;
	if ($bedtools_version < 19) { 
    	print " ERROR: BEDtools version $bedtools_version is too old. Please, consider installing the newest version\n"; exit;
    }
}

my $samtools = `which samtools`;
chomp $samtools;

if (!$samtools) {
	print " ERROR: samtools was not found on PATH\n"; exit;
}

my $macs2 = `which macs2`;
chomp $macs2;

if (!$macs2) {
	print " ERROR: macs2 was not found on PATH\n"; exit;
}

 my %CPP = (
	"1_SEQLIB"      => "$currentDir/bin/SeqLib",
	"2_GRAPES_SV"   => "$currentDir/bin/grapes_sv",
	"3_TARGETDEPTH" => "$currentDir/bin/TargetDepth",
	"4_OFFTARGET"   => "$currentDir/bin/offtarget",
	"5_AnnFeat"     => "$currentDir/bin/AnnFeat"
 );

 print " ## Compiling CPP ##\n";
 foreach my $program (sort keys %CPP) {

  print YELLOW " ## Compiling $program ##\n", RESET;
  chdir $CPP{$program};

  if ($program eq "1_SEQLIB") {

	my $cmd = "make clean";
	system $cmd;

	$cmd = "./configure";
	system $cmd;

	$cmd = "make CXXFLAGS=\"-std=c++11\"";
	system $cmd;

	$cmd = "make install";
	system $cmd;
  }
  else {
	  my $cmd = "make clean";
	  system ($cmd);
	  
	  $cmd = "make CXXFLAGS=\"-std=c++11\"";
	  system ($cmd);
  }
 }

 print " ## Setting SeqLib ##\n";
 
 print YELLOW " ## Checking executables ##\n", RESET;
 my $main = 
 ( -x "$currentDir/GRAPES" )
 ? " SUCCESS: GRAPES main can be executed\n"
 : die " ERROR: GRAPES cannot be executed\n";

 my $cmd = "perl $currentDir/GRAPES";
 if (`$cmd`) {
	print " SUCCESS: GRAPES main can be executed\n"
 }
 else {
	print " ERROR: GRAPES cannot be executed\n"
 }
 my $grapes_sv = 
 ( -x "$currentDir/bin/grapes_sv/GRAPES" )
 ? print " SUCCESS: grapes_sv can be executed\n"
 : die " ERROR: grapes_sv cannot be executed\n";

 my $targetDepth = 
 ( -x "$currentDir/bin/TargetDepth/TargetDepth" )
 ? print " SUCCESS: TargetDepth can be executed\n"
 : die " ERROR: TargetDepth cannot be executed\n";

 my $AnnFeat = 
 ( -x "$currentDir/bin/AnnFeat/AnnFeat" )
 ? print " SUCCESS: AnnFeat can be executed\n"
 : die " ERROR: AnnFeat cannot be executed\n";

 #Downloading mappability tracks

 if ( !-e "$currentDir/mappability/GRCh38.mappability.100mer.bedGraph.gz" ) {
	    my $name = "GRCh38.mappability.100mer.bedGraph.gz?dl=0";
	    my $new_name = "GRCh38.mappability.100mer.bedGraph.gz";
	    print " INFO: Downloading GRCh38.mappability.100mer.bedGraph.gz\n";
	    `wget https://www.dropbox.com/s/3hnvpczq9sbown6/GRCh38.mappability.100mer.bedGraph.gz?dl=0 -O $currentDir/mappability/$new_name`;
	   # rename "$currentDir/$name", "$currentDir/$new_name";
 }

 if ( !-e "$currentDir/mappability/GRCh38.mappability.100mer.bedGraph.tar.gz" ) {
	    my $name = "GRCh38.mappability.100mer.bedGraph.tar.gz?dl=0";
	    my $new_name = "GRCh38.mappability.100mer.bedGraph.tar.gz";
	    print " INFO: Downloading GRCh38.mappability.100mer.bedGraph.tar.gz\n";
	     `wget https://www.dropbox.com/s/14fn751kd7tzzue/GRCh38.mappability.100mer.bedGraph.tar.gz?dl=0 -O $currentDir/mappability/$new_name`;
	  #  rename "$currentDir/$name", "$currentDir/$new_name";
 }

 if ( !-e "$currentDir/mappability/wgEncodeCrgMapabilityAlign100mer.chr.bedgraph.gz" ) {
	    my $name = "wgEncodeCrgMapabilityAlign100mer.chr.bedgraph.gz?dl=0";
	    my $new_name = "wgEncodeCrgMapabilityAlign100mer.chr.bedgraph.gz";
	    print " INFO: Downloading wgEncodeCrgMapabilityAlign100mer.chr.bedgraph.gz\n";
	    `wget https://www.dropbox.com/s/hjf0a9ga8je6ndz/wgEncodeCrgMapabilityAlign100mer.chr.bedgraph.gz?dl=0 -O $currentDir/mappability/$new_name`;
	   # rename "$currentDir/$name", "$currentDir/$new_name";
 }

 if ( !-e "$currentDir/mappability/wgEncodeCrgMapabilityAlign100mer.nochr.bedgraph.gz" ) {
	    my $name = "wgEncodeCrgMapabilityAlign100mer.nochr.bedgraph.gz?dl=0";
	    my $new_name = "wgEncodeCrgMapabilityAlign100mer.nochr.bedgraph.gz";
	    print " INFO: Downloading wgEncodeCrgMapabilityAlign100mer.nochr.bedgraph.gz\n";
		`wget https://www.dropbox.com/s/e1mtadamo1jr7um/wgEncodeCrgMapabilityAlign100mer.nochr.bedgraph.gz?dl=0 -O $currentDir/mappability/$new_name`;
	    #rename "$currentDir/$name", "$currentDir/$new_name";
 }

if ( !-e "$currentDir/db/gnomad/gnomad_v2.1_sv.sites.vcf.gz" ) {
	    my $name = "gnomad_v2.1_sv.sites.vcf.gz";
	    print " INFO: Downloading $name\n";
		`wget https://storage.googleapis.com/gnomad-public/papers/2019-sv/gnomad_v2.1_sv.sites.vcf.gz -O $currentDir/db/gnomad/$name`;
		`wget https://storage.googleapis.com/gnomad-public/papers/2019-sv/gnomad_v2.1_sv.sites.vcf.gz.tbi -O $currentDir/db/gnomad/$name.tbi`;
	    #rename "$currentDir/$name", "$currentDir/$new_name";
 } 

if ( !-e "$currentDir/db/gencode.v27/gencode.v27.annotation.gff3.gz") {
	my $name = "gencode.v27.annotation.gff3.gz";
	print " INFO: Downloading $name\n";
	`wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.annotation.gff3.gz -O $currentDir/db/gencode.v27/$name`;
}

checkRpackages();

sub checkRpackages {
    print " INFO: Checking R packages\n";
    open R, ">", "Rinstaller.R" || die "ERROR: Unable to open Rinstaller\n";
    print R "is.installed <- function(mypkg){ is.element(mypkg, installed.packages()[,1])}\n";

    print R "if (!is.installed(\"ggplot2\")){
print(\"Installing ggplot2:\n\"); install.packages(\"ggplot2\", repos = \"http://cran.us.r-project.org\")
}";
    print R "else {
print(\"INFO: ggplot2 is already installed\n\")
}\n";

    print R "if (!is.installed(\"RColorBrewer\")){
print(\"Installing RColorBrewer:\n\"); install.packages(\"RColorBrewer\", repos = \"http://cran.us.r-project.org\")
}";
    print R "else {
print(\"INFO: RColorBrewer is already installed\n\")
}\n";

    print R "if (!is.installed(\"ggsignif\")){
print(\"Installing ggsignif:\n\"); install.packages(\"ggsignif\", repos = \"http://cran.us.r-project.org\")
}";
    print R "else {
print(\"INFO: ggsignif is already installed\n\")
}\n";

    print R "if (!is.installed(\"gtools\")){
print(\"Installing gtools:\n\"); install.packages(\"gtools\", repos = \"http://cran.us.r-project.org\")
}";
    print R "else {
print(\"INFO: gtools is already installed\n\")
}\n";

    print R "if (!is.installed(\"PSCBS\")){
print(\"Installing PSCBS:\n\"); install.packages(\"PSCBS\", repos = \"http://cran.us.r-project.org\")
}";
    print R "else {
print(\"INFO: PSCBS is already installed\n\")
}\n";

    print R "if (!is.installed(\"egg\")){
print(\"Installing egg:\n\"); install.packages(\"egg\", repos = \"http://cran.us.r-project.org\")
}";
    print R "else {
print(\"INFO: egg is already installed\n\")
}\n";
    print R "if (!is.installed(\"corrplot\")){
print(\"Installing corrplot:\n\"); install.packages(\"corrplot\", repos = \"http://cran.us.r-project.org\")
}";
    print R "else {
print(\"INFO: corrplot is already installed\n\")
}\n";

    print R "if (!is.installed(\"gplots\")){
print(\"Installing gplots:\n\"); install.packages(\"gplots\", repos = \"http://cran.us.r-project.org\")
}";
    print R "else {
print(\"INFO: gplots is already installed\n\")
}\n";

    print R "if (!is.installed(\"gridExtra\")){
print(\"Installing gridExtra:\n\"); install.packages(\"gridExtra\", repos = \"http://cran.us.r-project.org\")
}";
    print R "else {
print(\"INFO: gridExtra is already installed\n\")
}\n";

    close R;
    my $cmd = "Rscript Rinstaller.R";
    system($cmd);
}

 print " ## DONE ##\n";
