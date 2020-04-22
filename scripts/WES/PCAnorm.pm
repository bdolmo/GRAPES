#!/usr/bin/env perl

package PCAnorm;

use strict;
use Getopt::Long;
use File::Basename; 
use Data::Dumper;
use Parallel::ForkManager;
use Sort::Key::Natural qw(natsort);

sub doPCA {

    my $inputFile    = shift;
	my $outputDir    = shift;
    my $hashOfRegion = shift;

    my $isCov = 0;
    my $normalizedRD = "$outputDir/$::outName.NormalizedCounts.bed";
    if ($inputFile =~/Coverage/) {
        $isCov = 1;
        $normalizedRD = "$outputDir/$::outName.NormalizedCoverage.bed";
    }

	my $tmpPCAnorm = "$outputDir/$::outName.normalized_PCA.tmp.bed";
    my %HashReg = fillExonHash($inputFile . ".gz");

    fillMappability();
    Utils::decompressFile($inputFile . ".gz");

    #my %HoS = getSampleIdxFromRD($inputFile);
    open (R, ">", "$::outDir/PCA.R");
    print R "library(ggplot2)\n";
    print R "library(gridExtra)\n";
    print R "mydata<-read.table(file=\"$inputFile\", sep =\"\t\", header = TRUE)\n";
    print R "attach(mydata)\n";
    print R "mydata.subset<-mydata[,6:ncol(mydata)]\n";
    print R "newdata<-t(mydata.subset)\n";
    print R "pc.use <- nrow(newdata)\n";

    print R "newdata<-newdata[ , apply(newdata, 2, var) != 0]\n";
    print R "res.pca <- prcomp(newdata, scale=T, center=T)\n";
    print R "VE <- res.pca\$sdev^2\n";
    print R "PVE <- VE / sum(VE)\n";
    print R "PVE <- round(PVE, 2)\n";
    print R "PVE<- as.data.frame(PVE)\n";
    print R "newPVE <-PVE[PVE>0]\n";
    print R "nComp <-length(newPVE)\n";
    print R "newPVE<- as.data.frame(newPVE)\n";

    print R "std_dev <- res.pca\$sdev\n";
    print R "pr_var <- std_dev^2\n";
    print R "prop_varex <- pr_var/sum(pr_var)\n";
    print R "prop_varex[1:nComp]\n";
    print R "cumulativeVariance <- 0\n";
    print R "count <- 0\n";
    print R "numPC <- 2\n";
    print R "for (val in prop_varex) {
                cumulativeVariance<-cumulativeVariance+(val)
                count<-count+1
                if (cumulativeVariance > $::PCAvariance) {
                    numPC <- count-1
                    break
                }
            }\n";

    print R "png(\"$outputDir/ScreePlot.png\", width=1000, height=500, res=190)\n";
    print R "screePlot<-ggplot(newPVE) + geom_col(aes(x=seq_along(newPVE), y=newPVE*100)) + xlab(\"Principal Component\") + ylab(\"Variation explained(%)\") +ggtitle(\"Scree Plot\")+ theme_classic() + scale_fill_manual(values=(\"red\")) + scale_x_continuous(breaks = seq(0, nComp, by = 1)) \n";
    print R "Cumulative <- ggplot(newPVE) + geom_line(aes(x=seq_along(newPVE), y=cumsum(newPVE*100))) + geom_point(aes(x=seq_along(newPVE), y=cumsum(newPVE*100)), colour = \"black\", fill = \"white\") + xlab(\"Principal Component\") + ylab(\"Cumulative var (%)\") +ggtitle(\" \")  + theme_classic()+ scale_x_continuous(breaks = seq(0, nComp, by = 1)) \n";
    print R "grid.arrange(screePlot, Cumulative, nrow = 1)\n";
    print R "dev.off()\n";
    print R "trunc <- res.pca\$x[,numPC:pc.use] %*% t(res.pca\$rotation[,numPC:pc.use])\n";
    print R "trunc <- scale(trunc, center = FALSE , scale=1/res.pca\$scale)\n";
    print R "trunc <- scale(trunc, center = -1 * res.pca\$center, scale=FALSE)\n";
    print R "trunc <- t(trunc)\n";
    print R "write.table(trunc, file = \"$tmpPCAnorm\", sep=\"\t\", row.names = TRUE, col.names = TRUE)\n";
    close R;

    `$::Rscript $::outDir/PCA.R`;

    # Add header and exon columns
    open (IN, "<", $tmpPCAnorm);
    open (OUT, ">", $normalizedRD ) || die " ERROR: Unable to open $normalizedRD\n";
    my $n = 0;
	foreach my $coordinate( natsort keys %HashReg ) {
        if ($n == 0) {
            my $line = <IN>;
            chomp $line;
            my @tmp    = split (/\t/, $line);
            my @subset = join("\t", @tmp[0..@tmp-1]);
            my $newLine= join ("\t", @subset);
            $newLine =~s/\"//g;
            if ($isCov) {
                print OUT "chr\tstart\tend\texon\t\%GC\t$newLine\n";
            }
            else {
                print OUT "chr\tstart\tend\texon\t\%GC\t\%MAP\t$newLine\n";
            }
            $n++;
            next;
        }
        my $line = <IN>;
        chomp $line;
        my @tmp = split (/\t/, $line);
        my @subset = @tmp[1..@tmp-1];
        @subset = map { int($_) } @subset;        
        my $newLine = join ("\t", @subset);

        if ($isCov) {
            print OUT "$coordinate\t$HashReg{$coordinate}{GC}\t$newLine\n";
        }
        else {
            print OUT "$coordinate\t$HashReg{$coordinate}{GC}\t$::ExonFeatures{$coordinate}{MAP}\t$newLine\n";
        }
        $n++;
    }
    close IN;
    close OUT;

    unlink($tmpPCAnorm);

    # Compress normalized counts
    Utils::compressFile($normalizedRD);
    Utils::compressFile($inputFile);

}

##########################
sub fillExonHash {

    my $inputFile = shift;
    my %hash = ();
    open (IN, "$::zcat $inputFile |") || die " ERROR: Unable to open $inputFile\n";
    while (my $line=<IN>) {
        chomp $line;
        next if $line =~/^chr\tstart/;
        my @tmp = split (/\t/, $line);
        my $coordinate = join ("\t", @tmp[0..3]);
        my $gc = $tmp[4];
        $hash{$coordinate}{GC} = int($gc);
    }
    close IN;
    return %hash;
}

##########################
sub fillMappability {

    #my %hash = ();
    open (IN, "<", "$::ontargetDir/mappability_ontarget.bed") || die " ERROR: Unable to open $::ontargetDir/mappability_ontarget.bed\n";
    while (my $line=<IN>) {
        chomp $line;
        my @tmp = split (/\t/, $line);
        my $coordinate = join ("\t", @tmp[0..3]);
        my $map = $tmp[4];
        $::ExonFeatures{$coordinate}{MAP} = $map;
    }
    close IN;
}

##########################
 sub getSampleIdxFromRD {

	my $ratioFile = shift;
	my $strHeaderRatios = `$::head -1 $ratioFile`;
	my @tmpHeaderRatios = split (/\t/, $strHeaderRatios );

	my %HashPosAll = ();
	for (my $i = 5; $i < @tmpHeaderRatios; $i++) {
		$HashPosAll{$i} = $tmpHeaderRatios[$i];
	} 
	return %HashPosAll;
 }
 return 1;