#!/usr/bin/env perl

package Reference;

use strict;
use warnings;
use File::Basename;
use Sort::Key::Natural qw(natsort);
use List::MoreUtils qw(uniq);
use DBI;

# Package for managing sample references. Currently supporting:
# -Creation
# -Update
# -Importing

######################
sub main  {
    if (!-e $::database) {
        print " INFO: Creating $::roiName database\n";
        createReference();
    }
}

######################
sub createDB {

    my $normCounts = shift;
    my $function   = shift;

    if ($normCounts =~/.gz$/) {
        Utils::decompressFile($normCounts);
        $normCounts =~s/.gz//;
    }

    my $header = `$::head -1 $normCounts`;
    chomp $header;

    my @tmpHeader = split("\t", $header);
    my $count = 0;
    my @arrHeader;
    foreach my $value (@tmpHeader) {
        if ($count == 0) {
            push @arrHeader, "$value TEXT NOT NULL";
        }
        elsif ($count == 1) {
            push @arrHeader, " $value INTEGER NOT NULL";
        }
        elsif ($count == 2) {
            push @arrHeader, " $value INTEGER NOT NULL";
        }
        elsif ($count == 3) {
            push @arrHeader, " $value TEXT NOT NULL";
        }
        else {
            push @arrHeader, " \"$value INTEGER NOT NULL\"";
        }
        $count++;
    }

    # Create table if not existent
    my $column_list = join ', ', map $::dbh->quote_identifier($_), @arrHeader;
    my $headerSql = join (",", @arrHeader);

    my $sql = "DROP TABLE IF EXISTS ROI";
    my $rv  = $::dbh->do($sql);

    $sql = "CREATE TABLE ROI ($headerSql)";
    $rv = $::dbh->do($sql);

    # Insert values at each roi
    my @arrSearch =();
    foreach my $f (@arrHeader) {
        push @arrSearch, "?";
    }
    my $nFields = join (", ", @arrSearch);

    $sql = "INSERT INTO ROI VALUES($nFields)";    
    $rv = $::dbh->prepare($sql);

    open (NC, "<", $normCounts) || die " ERROR: Unable to open $normCounts\n";
    while (my $line=<NC>) {
        chomp $line;
        next if $line =~/^chr\tstart/; # skipping header line
        my @tmp = split(/\t/, $line);
        $rv->execute(split (/\t/, $line) );
    }
    close NC;
    $::dbh->commit() or die $::dbh->errstr;

    print " INFO: Successfully $function sqlite database\n";
    Utils::compressFile($normCounts);
}
######################
sub exportDB {

    open (OUTDB, ">", $::HoF{EXPORTED_DB}) 
    || die " ERROR: Unable to open $::HoF{EXPORTED_DB}\n";

    # Selecting header and exporting it
    my $sth = $::dbh->prepare("SELECT * FROM ROI WHERE 1=0");
    $sth->execute();
        $::dbh->commit() or die $::dbh->errstr;

    my $fields = $sth->{NAME};
    foreach my $field ( @$fields ) {
        my @tmp = split (/\s/, $field);
        $field = $tmp[0] if @tmp > 1;
        print OUTDB "$field\t";
    }
    print OUTDB "\n";

    # Exporting all body content 
    $sth = $::dbh->prepare("SELECT * FROM ROI");
    $sth->execute();
    my $result = $sth->fetchall_arrayref;

    foreach my $row ( @$result ) {
        my $tmp = join("\t", @$row);
        #my ($chr, $start, $end, $info) = split (/\t/, $tmp);
        print OUTDB "$tmp\n";
    }
    $sth->finish();
    #$::dbh->disconnect();
    close OUTDB;

    if (-s $::HoF{EXPORTED_DB}) {
        print " INFO: Database succesfully exported\n";
    }

    Utils::compressFile( $::HoF{EXPORTED_DB} );
}

######################
sub updateDB {

    print " INFO: Updating database\n";

    # Get all created references
    my @references = glob ("$::ontargetDir/*.ref");

    my %HoR = (); # Hash of references
    my @refs= ();
    foreach my $reference (@references ) {

        my $refName = ( split "\t", `$::head -1 $reference` )[-1];
        chomp $refName;
        push @refs, $refName;

        open (IN, "<", $reference) || die " ERROR: Unable to open $reference\n";
        my $n = 0;
        while (my $line=<IN>) {
            chomp $line;
            my @tmp = split (/\t/, $line);
            next if $line =~/^chr\tstart/;
            my $coordinate = join ("\t", @tmp[0..5]);

            # Now filling hash of references
            $HoR{$coordinate}{$refName} = $tmp[6];
        }
        close IN;
    }

    # Getting coordinate information
    my $header = `$::zcat $::HoF{NORM_COUNTS_ON}`;
    chomp $header;
    my @tmpHeader = split (/\t/, $header);

    open (OUT, ">", "$::ontargetDir/$::outName.merged.refs.txt")
     or die " ERROR: cannot open $::ontargetDir/$::outName.merged.refs.txt\n";

    # Printing header
    print OUT join("\t", @tmpHeader[0..5]) . "\t" . join("\t", @refs) . "\n";

    # Printing body
    foreach my $region ( natsort keys %HoR ) {
        print OUT "$region";

        # Dumping normalized coverage results
        foreach my $ref (natsort @refs) {
            print OUT "\t$HoR{$region}{$ref}";
        }
        print OUT "\n";
    }
    close OUT;

    # Now, creating a new database
    createDB("$::ontargetDir/$::outName.merged.refs.txt", "updated");

    # remove temporal reference files
    unlink (@references);
}

######################
sub mergeRunWithDB {

    # Goal is to merge *NormalizedCounts.bed with References.sqlite.txt

    my %HoS = (); # Hash of samples

    my ($hrefNC, $sampRefNC) = allocateFile($::HoF{NORM_COUNTS_ON});
    my ($hrefDB, $sampRefDB) = allocateFile($::HoF{EXPORTED_DB});

    # derreferencing sample arrays
    my @sampNC = @$sampRefNC;
    my @sampDB = @$sampRefDB;

    # Get uniq samples from both files
    my @samples = uniq (@sampNC, @sampDB);

    my %HoNC = %$hrefNC; # Hash of Normalized Coverage
    my %HoDB = %$hrefDB; # Hash of Database

    my $header = `$::zcat $::HoF{NORM_COUNTS_ON}`;
    chomp $header;
    my @tmpHeader = split (/\t/, $header);

    open (OUT, ">", $::HoF{MERGED_NORM_COV}) || die " ERROR: Unable to open $::HoF{MERGED_NORM_COV}\n";

    # Printing header
    print OUT join("\t", @tmpHeader[0..5]) . "\t" . join("\t", @samples) . "\n";

    # Now printing full body
    foreach my $coordinate (natsort keys %HoNC ){
        print OUT "$coordinate";

        # Dumping normalized coverage results
        foreach my $sample (natsort @samples) {
            if (exists $HoNC{$coordinate}{$sample}) {
                print OUT "\t$HoNC{$coordinate}{$sample}";
                next;
            }
            elsif (exists $HoDB{$coordinate}{$sample}) {
                print OUT "\t$HoDB{$coordinate}{$sample}";
                next;
            }
        }
        print OUT "\n";
    }
    close OUT;
}

###################### 
sub allocateFile {

    # Goal is to read HoF{EXPORTED_DB} and HoF{NORM_COUNTS_ON} 
    # and return a hash reference to merge later

    my $inputFile = shift;

    my @fileArray = split ("\n", `$::zcat $inputFile`);
    my %hash = ();
    my $i = 0;
    my @samples = ();
    foreach my $line (@fileArray) {
        my @tmp = split ("\t", $line);
        my $coordinate = join ("\t", @tmp[0..5]);
        if ($i == 0) { # header
            for (my $i = 6; $i < @tmp; $i++) {
                push @samples, $tmp[$i];
            }
        }
        else {
            my $j = 0;
            for (my $i= 6; $i < @tmp; $i++) {
                $hash{$coordinate}{$samples[$j]} = $tmp[$i];
                $j++;
            }
        }
        $i++;
    }

    return \%hash, \@samples;
}

########################
 sub clusterBatches {

  my %nodes = ();
  my %seen  = ();
  my $count = 0;
  my @samples = ();

  open (IN, "<", $::HoF{CORRELATIONS_ON}) || die " ERROR: Cannot open $::HoF{CORRELATIONS_ON}\n";
  while (my $line =<IN>) {

	chomp $line;
	$line=~s/"//g;
	my @tmp = split (/\t/, $line);
    my $sample = $tmp[0];

    # Getting sample names from the header
	if ($count == 0) {
		@samples = @tmp;
		$count++;
		next;
	}

    if (!exists $::sampleHash{$sample}) {
        $::referenceHash{$sample} = 1;
    }

	my $j = 0;
	for (my $i = 1; $i < @tmp; $i++) {

		# skipping self correlation
		if ($sample eq $samples[$j]) {
			$j++;
			next;
		}
		# Creating nodes
	 	$nodes{$sample}{$samples[$j]} = $tmp[$i];
	 	$nodes{$samples[$j]}{$sample} = $tmp[$i];
		$j++;
	}
	$seen{$sample}++;
	$count++;
 }
 close IN;

 open (LOG, ">", $::HoF{REFERENCES}) || die " ERROR: Cannot open $::HoF{REFERENCES}\n";
 print LOG "N\tSAMPLE\tREFERENCE_SET\tMEAN_CORRELATION\n";
 my @cluster = ();
 $count = 0;
 my $sum = 0;
 my $mean_corr;
 my $sampleCount = 0;

 foreach my $sample1 (@samples) {
    #print"$sample1\n";
    next if exists $::referenceHash{$sample1};

    if ($::doCaseControl) {
        next if $::sampleHash{$sample1}{CONTROL};
    }

	$count++;
	my %sampleCorr = ();
	foreach my $sample2 (@samples) {
        if ($::doCaseControl) {
            next if $::sampleHash{$sample2}{CASE};
        }
		next if $sample1 eq $sample2;
		if ($nodes{$sample1}{$sample2} >= $::minCorrelation) {
			my $corr = $nodes{$sample1}{$sample2};
			$sampleCorr{$corr} = $sample2;
		}
	}

	# Selecting the 15th most well correlated samples
	foreach my $corr (reverse sort keys %sampleCorr) {

		$sampleCount++;

		if ($::doPooled) {
			last if $sampleCount > $::maxSampleSizeCluster;
		}

		push @{$::sampleHash{$sample1}{REFERENCE}}, $sampleCorr{$corr};
		$sum+=$corr;
		my $corrdecimal = sprintf "%.3f", ($corr); 
		push @cluster, "$sampleCorr{$corr}($corrdecimal)";
	}
	$sampleCount = 0;
	my $reference;
	if (@cluster == 0) {
		$mean_corr = 1;	
		$reference = "none";
		push @{$::sampleHash{$sample1}{REFERENCE}}, $reference;
	}
	else {
		my $n = scalar @cluster > 0 ? scalar@cluster : 1;
		$mean_corr = sprintf "%.3f", ($sum/$n);
		$reference = join (",", @cluster);
	}

	$::sampleHash{$sample1}{CORRELATION} = $mean_corr;

	print LOG "$count\t$sample1\t" . join ("\t", @cluster) . "\t" . $mean_corr . "\n";
	@cluster = ();
	$sum = 0;
 }

}


########################
 sub getCorrAndPlot {

	my $outputDir = $::ontargetDir;

	my $libraries = "library(ggplot2)\nlibrary(RColorBrewer)\nlibrary(corrplot)\nlibrary(gplots)
    \nlibrary(grid)\nlibrary(gridExtra)\nlibrary(reshape2)\n";

	my $myFile = qq("$::HoF{MERGED_NORM_COV}\");
	open (R, ">", "$outputDir/plotCorrelationMatrix.R") || die " ERROR: Cannot open plotCorrelationMatrix.R\n";
    print R "$libraries\n";
	print R "mydata<-read.table(file=$myFile ,sep ='\t', check.names = FALSE, header=TRUE)\n";
	print R "attach(mydata)\n";
	print R "corrected_seq <- as.matrix(mydata[seq (7, ncol(mydata))])\n";

	# Calculating Pearson correlation
	print R "cor_corrected<-cor(corrected_seq, method=\"pearson\")\n";
	my $greyStr;
	for (my $i = 8; $i <=87; $i++) {
		$greyStr .= qq(\"grey$i\",);
	}

    if (!-e "$outputDir/$::outName.heatmap.png") {
        print R "png(\"$outputDir/$::outName.heatmap.png\", res = 300, height=2800, width=2800)\n";
        
        # Adjusting margin size based on the longest label
        #https://stackoverflow.com/questions/38249567/r-figure-being-cropped
        print R "margin_width = max(strwidth(rownames(cor_corrected), units=\'inches\')) * par(\'fin\')[1]\n";

        # Defining a color pallete
        print R "my_palette <- colorRampPalette(c(\"darkblue\", \"#196aff\",\"#6f6fff\", \"#327aff\", 
        \"white\", \"#ffdfc0\", \"#ffccb3\", \"#CC483A\", \"#790000\"))(n = 150)\n";

        # Creating heatmap
        print R "heatmap.2(x = cor_corrected,colsep=1:ncol(cor_corrected),
        rowsep=1:nrow(cor_corrected), trace=\"none\", sepwidth=c(0.01, 0.01), sepcolor=\"white\", scale=\"none\", 
        col = my_palette, symm = TRUE, margins=c(margin_width, margin_width))\n";
        print R "dev.off()\n";
    }

    # Writing a correlation table
	print R "write.table(cor_corrected, \"$outputDir/$::outName.cor_corrected.txt\", sep=\"\\t\")\n";
	close R;
    if ($::verbose) {
	    `$::Rscript $outputDir/plotCorrelationMatrix.R`;
    }
    else {
        `$::Rscript $outputDir/plotCorrelationMatrix.R $::devNull`;
    }
	#unlink("$outputDir/bias_info.tmp.txt");
 }



return 1;