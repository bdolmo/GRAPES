#!/usr/bin/env perl

package Annotate;

use strict;
use Getopt::Long;
use File::Basename;
use List::MoreUtils qw(uniq);
use Misc::Utils;

# Purpose of this module is to annotate VCF files with gnomAD database

#Hardcoded default fields
our @popDefaults = (
    "EVIDENCE",
    "AF",
    "AN",
    "AC",
    "AFR_AC",
    "AFR_AN",
    "AFR_AF",
    "AMR_AC",
    "AMR_AN",
    "AMR_AF",
    "EAS_AC",
    "EAS_AN",
    "EAS_AF",
    "EUR_AC",
    "EUR_AN",
    "EUR_AF",
    "OTH_AC",
    "OTH_AN",
    "OTH_AF"
);

################################
sub doAnnotation {

    my $inputVCF = shift;
    my $outDir   = shift;
    my $minOlap  = shift;

    if (!-e $inputVCF) {
      print " ERROR: missing vcf $inputVCF\n";
      return;
    }

    # Get a hash with vcf fields from gnomad
    my %headerVCF = returnVcfFields();

    # Annotate Gene names (HUGO)
    my $annotVCF = annotateGenes($inputVCF, $outDir);

     # Annotate with gnomAD SV dataset (v2.1)
    $annotVCF = annotateGnomad($annotVCF, $outDir, $minOlap, \%headerVCF);

}

################################
sub getGenomicCode {
    my $chr    = shift;
    my $start  = shift;
    my $end    = shift;
    my $svtype = shift;
}

################################
sub rewriteHeader  {

  my $inputVCF     = shift;
  my $newAnnoField = shift;

  my $headerStr;
  if ( isGzipped($inputVCF) ) {
    $headerStr = `$::zcat $inputVCF | $::head -5000 | $::grep -e '^#'`;
  }
  else {
    $headerStr = `$::head -5000 $inputVCF | $::grep -e '^#'`;
  }
  chomp $headerStr;
  my @tmpHeader = split("\n", $headerStr);

  # Get an index of the latest INFO field
  my $idx = 0;
  foreach my $field (@tmpHeader)  {
    if ($field=~/^##INFO/) {
      last;
    }
    $idx++;
  }
  # Insert new field into last INFO position
  splice @tmpHeader, $idx, 0, $newAnnoField;
  my $newHeader = join("\n", @tmpHeader);

  return $newHeader;

}


################################
sub checkChrConvention  {

  my $inputVCF = shift;
  my $isChr = 0;
  my $fieldFound = 0;

  my $zcat = isGzipped($inputVCF) ? $::zcat : $::cat;
  open (IN, "$zcat $inputVCF | $::head -100 |") || die "ERROR: Unable to open $inputVCF\n";

  while (my $line=<IN>) {
    chomp $line;
    if ($line=~/^##contig/)  {
      $fieldFound = 1;
      if ($line=~/chr/) {
        $isChr = 1;
        last;
      }
      else {
        last;
      }
    }
  }
  close IN;

  if (!$fieldFound) {
    print " ERROR: Missing contig header fields on $inputVCF to extract chr naming convention (##contig)\n";
    exit;
  }
  return $isChr;
}

################################
sub annotateGenes {

  my $inputVCF = shift;
  my $outDir   = shift;

  my $basenameInputVcf = basename($inputVCF);
  my $tmpVCF   = "$outDir/$basenameInputVcf";
  my $annotVCF = "$outDir/$basenameInputVcf";

  if ($basenameInputVcf =~/.gz/) {
    $tmpVCF   =~s/.vcf.gz/.tmp.vcf/;
    $annotVCF =~s/.vcf.gz/.annotated.vcf/ if $annotVCF !~/.annotated/;
  }
  else {
    $tmpVCF   =~s/.vcf/.tmp.vcf/;
    $annotVCF =~s/.vcf/.annotated.vcf/ if $annotVCF !~/.annotated/;
  }

  if ( !checkChrConvention($inputVCF) ) {
    open(IN, "$::zcat $::geneList | $::awk \'{gsub(/chr/,\"\", \$1); print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$4\"\t\"\$5\"\t\"\$6\"\t\"\$7}\' |".
      " $::bedtools intersect -a $inputVCF -b stdin -wao |");
  }
  else {
    open(IN, "$::bedtools intersect -a $inputVCF -b $::geneList -wao |");
  }
  print "$::bedtools intersect -a $inputVCF -b $::geneList -wao\n" if $::verbose;
  my %HashAnno = ();
  while (my $line=<IN>) {
      chomp $line;
      next if $line =~/#/;
      my @tmp = split("\t", $line);
      my $coordinates = join("\t", @tmp[0..4]);
      my $gene = $tmp[-3];
      if ($gene ne ".") {
        push @{$HashAnno{$coordinates}}, $gene;
        @{$HashAnno{$coordinates}}= uniq(@{$HashAnno{$coordinates}});
      }
  }
  close IN;

  my $GeneHeader = "##INFO=<ID=GENES,Number=1,Type=String,Description=\"Gene name in HUGO notation\">";
  my $newHeader = rewriteHeader($inputVCF, $GeneHeader);
  my @newHeaderArr = split("\n", $newHeader);

  # Write new header
  open (OUT, ">", $tmpVCF) || die " ERROR: Unable to find $tmpVCF\n";
  foreach my $entry (@newHeaderArr) {
    print OUT "$entry\n";
  }

  # Now loop from initial unnanotated variants and append new fields
  if ( isGzipped($inputVCF) ) {
    open(IN, "$::zcat $inputVCF |");
  }
  else {
    open(IN, "<", $inputVCF) || die " ERROR: Unable to open $inputVCF\n";
  }
  while (my $line=<IN>) {
    chomp $line;
    # Skip header lines
    next if $line =~/#/;
    my @tmp = split("\t", $line);
    my $coordinates = join("\t", @tmp[0..4]);

    my $genes = '.';
    if (exists $HashAnno{$coordinates}) {
      my @geneArr = ();
      foreach my $val (@{$HashAnno{$coordinates}}) {
        push @geneArr, $val;
      }
      $genes = join(",", @geneArr);
    }

    # Fetch genes
    my @info = split(";", $tmp[7]);

    # Append new field and rejoin INFO
    push @info, "GENES=$genes";
    $tmp[7] = join(";", @info);

    # Dump new annotated variant
    my $newLine = join("\t", @tmp);
    print OUT "$newLine\n";
  }
  close IN;
  close OUT;

  # Rename new annotated vcf same as initial vcf
  rename $tmpVCF, $annotVCF;

  # Compress VCFs
  my $compressedVCF = $annotVCF . ".gz";
  Utils::compressFileBgzip($annotVCF);

  # Index compressed VCF with tabix
  Utils::tabixIndex($compressedVCF, "vcf");

  return $compressedVCF;
}

################################
sub annotateGnomad {

    my $inputVCF = shift;
    my $outDir   = shift;
    my $minOlap  = shift;
    my $hashRef  = shift;

    my %headerVCF = %$hashRef;

    my $basenameInputVcf = basename($inputVCF);
    my $tmpVCF   = "$outDir/$basenameInputVcf";
    my $annotVCF = "$outDir/$basenameInputVcf";

    if ($basenameInputVcf =~/.gz/) {
      $tmpVCF   =~s/.vcf.gz/.tmp.vcf/;
      $annotVCF =~s/.vcf.gz/.annotated.vcf/ if $annotVCF !~/.annotated./;
      $annotVCF =~s/.gz//;
    }
    else {
      $tmpVCF   =~s/.vcf/.tmp.vcf/;
      $annotVCF =~s/.vcf/.annotated.vcf/ if $annotVCF !~/.annotated./;
    }
    my $fh;
    if ( isGzipped($inputVCF) ) {
        open $fh, "$::zcat $inputVCF |";
    }
    else {
        open ($fh, "<", $inputVCF)
            || die "ERROR: Unable to open $inputVCF\n";
    }
    my $hasChrom = checkChrConvention($inputVCF);

    my %seenVar = ();
    my $flag = 0;

    open (VCF, ">", $tmpVCF) || die " ERROR: Unable to open $tmpVCF\n";
    while (my $line=<$fh>) {
        chomp $line;
        # Header section
        if ($line=~/^#/) {
            if ($line=~/##ALT/ && !$flag) {
                $flag = 1;
                foreach my $field (@popDefaults) {

                    my $outField = $headerVCF{$field}{DESCRIPTION};
                    $outField =~s/EVIDENCE/gnomAD_SV_method/;

                    print VCF "$outField\n";
                }
                print VCF "$line\n";
            }
            else {
                print VCF "$line\n";
            }
            next;
        }

        # Skipping repeated lines
        $seenVar{$line}++;
        next if $seenVar{$line} > 1;

        my @tmp = split (/\t/, $line);
        my @info= split (/;/, $tmp[7]);

        my $chr   = $tmp[0];
        my $start = $tmp[1];

        my ($end) = grep ($_=~/^END=/, @info);
        $end=~s/END=// if $end;

        my ($svtype) = grep ($_=~/^SVTYPE=/, @info);
        ($svtype) =~s/SVTYPE=// if $svtype;

        my ($precision) = grep ($_=~/PRECISE/, @info);
        #print "$line\n";
        if (!$precision) {
          $precision = 'PRECISE';
        }
        if (!$hasChrom) {
          #$chr = "chr" . $chr;
          #print "$chr\n";
        }
        my @gnomad = tabixQuery($chr, $start, $end, $svtype, $precision, $minOlap);

        if (@gnomad) {
            foreach my $entry (@gnomad) {
                my @info = split(";", $entry);
                my @gathered = ();

                foreach my $desired (@popDefaults) {
                    my ($x) = grep ($_=~/^$desired=/, @info);
                    push @gathered, $x;
                }

                my $annot = join(";", @gathered);
                $annot =~s/EVIDENCE/gnomAD_SV_method/;

                my $extendedVar = join("\t", @tmp[0..7]) . ";" . $annot . join ("\t", @tmp[8..@tmp-1]);

                $seenVar{$extendedVar}++;
                next if $seenVar{$extendedVar} > 1;

                print VCF join("\t", @tmp[0..7]) . ";" . $annot . "\t". $tmp[8] . "\t" . $tmp[9] . "\n";
            }
        }
        else {
            my @tmp = split (/\t/, $line);
            my %popHash = map { $popDefaults[$_] => '.' } 0..$#popDefaults;
            my @popData = map { $_ . '=' . $popHash{$_} } keys %popHash;
            print VCF join("\t", @tmp[0..7]) . ";" . join(";", @popData) ."\t". $tmp[8] . "\t" . $tmp[9] . "\n";
        }
    }
    close $fh;
    close VCF;

    my ($previousAnnVCF) = basename(glob ("$outDir/$annotVCF"));
    if ($previousAnnVCF eq basename($annotVCF) ) {
      unlink $annotVCF;
    }
    rename $tmpVCF, $annotVCF;

    # Compress VCFs
    my $compressedVCF = $annotVCF . ".gz";
    Utils::compressFileBgzip($annotVCF);

    # Index compressed VCF with tabix
    Utils::tabixIndex($compressedVCF, "vcf");

    return $compressedVCF;
}

################################
sub tabixQuery {

    my $chrQuery   = shift;
    my $startQuery = shift;
    my $endQuery   = shift;
    my $typeQuery  = shift;
    my $precision  = shift;
    my $minOlap    = shift;

    my @records = ();

    if (!$endQuery || !$typeQuery || !$precision) {
        return @records;
    }

    # Remove chr if present
    if ($chrQuery =~/^chr/) {
        $chrQuery =~s/chr//;
    }

    my $str = `$::tabix $::gnomADvcf $chrQuery:$startQuery-$endQuery`;
    chomp $str;
    print "$::tabix $::gnomADvcf $chrQuery:$startQuery-$endQuery\n" if $::verbose;

    my @tmpStr = split(/\n/, $str);

    foreach my $entry (@tmpStr) {

        my @tmp = split (/\t/, $entry);
        my @info= split (/;/, $tmp[7]);

        my $chrDB   = $tmp[0];
        my $startDB = $tmp[1];

        my ($endDB) = grep ($_=~/^END=/, @info);
        $endDB=~s/END=//;

        my ($typeDB) = grep ($_=~/^SVTYPE=/, @info);
        ($typeDB)=~s/SVTYPE=//;

        # Calculate overlap between query and database
        my ($overlapQuery, $overlapDB)= calculateOverlap($chrQuery,
        $startQuery, $endQuery, $chrDB, $startDB, $endDB);

        my $joinedInfo = join(";", @info);

        # Variant class must be the same
        if ($typeDB eq $typeQuery) {

            # Do not apply reciprocal overlap for IMPRECISE calls
            if ($precision eq 'IMPRECISE') {

                if ($overlapQuery == 1) {
                    push @records, $joinedInfo;
                }
            }
            # But do it for PRECISE calls
            else {
                if ($overlapQuery >= $minOlap && $overlapDB >= $minOlap ) {
                    push @records, $joinedInfo;
                }
            }
        }
    }
    return @records;
}

################################
sub isGzipped {
    my $inputFile = shift;
    if ($inputFile =~/.gz/) {
        return 1;
    }
    else {
        return 0;
    }
}

################################
sub returnVcfFields {

    my %headerVCF = ();

    my @infoArray = ();
    open (IN, "$::zcat $::gnomADvcf | ") || die " ERROR: Unable to open $::gnomADvcf\n";
    while (my $line=<IN>) {
        chomp $line;
        if ($line =~/^##INFO/) {
            push @infoArray, $line;
        }
        last if $line !~/^##/;
    }
    close IN;

    my $i = 0;
    foreach my $field (@infoArray) {
        my @info = split(/,/, $field);
        my $tag = $info[0];
        $tag =~s/INFO=<ID=//;
        $tag =~s/#//g;
        $headerVCF{$tag}{DESCRIPTION} = $field;
        $headerVCF{$tag}{IDX} = $i;
        $i++;
    }
    return %headerVCF;
}

################################
sub calculateOverlap {

    my $chr   = shift;
    my $start = shift;
    my $end   = shift;
    my $chrom = shift;
    my $st    = shift;
    my $ed    = shift;

    my $length_A = $end-$start;
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
    return ($olap_A, $olap_B);
}

return 1;
