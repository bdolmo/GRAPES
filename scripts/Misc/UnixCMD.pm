package UnixCMD;
use strict;
use warnings;
use List::MoreUtils qw(uniq);
use Term::ANSIColor qw(:constants);
use Cwd qw(cwd abs_path);
use Config;

our $dirname = cwd();
=head1 NAME

  Grapes::IO 

=head1 SYNOPSIS

  use Grapes::UnixCMD

=head1 DESCRIPTION

=head1 AUTHOR

Written by Bernat del Olmo, PhD

=cut

=head1 METHODS

=cut

our $tabix = `which tabix`; 
chomp $tabix;

if (!$tabix) {
	print " ERROR: tabix was not found on PATH\n"; exit;
}

our $bgzip = `which bgzip`; 
chomp $bgzip;

if (!$bgzip) {
	print " ERROR: bgzip was not found on PATH\n"; exit;
}

our $cat  = `which cat`;
chomp $cat;
if (!$cat) {
	print "ERROR: cat command not found\n";
	exit;
}

our $zcat  = `which zcat`;
chomp $zcat;
if (!$zcat) {
	print "ERROR: zcat command not found\n";
	exit;	
}

our $cut  = `which cut`;
chomp $cut;
if (!$cut) {
	print "ERROR: cut command not found\n";
	exit;
}

our $awk  = `which awk`;
chomp $awk;
if (!$awk) {
	print "ERROR: awk command not found\n";
	exit;	
}
our $grep;
our $sort;
if ($Config{osname} =~/darwin/) {
	$grep = `which egrep`; chomp $grep;
	$sort = `which gsort`; chomp $sort;
}
else {
	$grep = `which grep`; chomp $grep;
	$sort = `which sort`; chomp $sort;
}

our $wc  = `which wc`;
chomp $wc;
if (!$wc) {
	print "ERROR: wc command not found\n";
	exit;
}

our $bedtools = `which bedtools`;
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

our $samtools  = `which samtools`;
chomp $samtools;
if (!$samtools) {
	print "ERROR: samtools command not found\n";
	exit;	
}

  our $macs2 = `which macs2`;
  chomp $macs2;

  if (!$macs2) {
    print " ERROR: macs2 was not found on PATH\n"; exit;
  }

our $sed  = `which sed`;
chomp $sed;
if (!$sed) {
	print "ERROR: sed command not found\n";
	exit;
}

our $paste  = `which paste`;
chomp $paste;
if (!$paste) {
	print "ERROR: paste command not found\n";
	exit;
}

our $head  = `which head`;
chomp $head;
if (!$head) {
	print "ERROR: head command not found\n";
	exit;
}

our $uniq  = `which uniq`;
chomp $uniq;
if (!$uniq) {
	print "ERROR: uniq command not found\n";
	exit;
}
 
our $mv  = `which mv`;
chomp $mv;
if (!$mv) {
	print "ERROR: mv command not found\n";
	exit;
}

our $tail  = `which tail`;
chomp $tail;
if (!$tail) {
	print "ERROR: mv command not found\n";
	exit;
}

our $Rscript  = `which Rscript`;
chomp $Rscript;

if (!$Rscript) {
    print " ERROR: Rscript was not found on PATH\n"; exit;
}

return 1;