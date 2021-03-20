#!/usr/bin/env Rscript  
library(tidyr)
library(dplyr)
library(argparse)
library(GenomicRanges)
library(reshape2)
library(ggplot2)
library(data.table)
parser<-ArgumentParser(description="Analyze single-exon CNV")

parser$add_argument('--raw_calls', dest='raw_calls', type="character", help="raw calls")
parser$add_argument('--coverage_data', dest='coverage_data', type="character", help="Per base coverage data")
parser$add_argument('--references', dest='references', type="character", help="Sample references")


parser$add_argument('--min_zscore', dest='min_zscore', type="double", help="Minimum z-score")
parser$add_argument('--min_s2n', dest='min_s2n', type="double", help="Minimum signal to noise")
parser$add_argument('--lower_del_cutoff', dest='lower_del_cutoff', type="double", help="Lower deletion ratio cutoff")
parser$add_argument('--upper_del_cutoff', dest='upper_del_cutoff', type="double", help="Upper deletion ratio cutoff")
parser$add_argument('--lower_dup_cutoff', dest='lower_dup_cutoff', type="double", help="Lower duplication ratio cutoff")
parser$add_argument('--output_dir', dest='output_dir', type="character", help="Output directory")

args <- parser$parse_args()

# Getting arguments
raw_calls <- args$raw_calls
coverage_data <- args$coverage_data
references <- args$references
min_zscore <- args$min_zscore
min_s2n    <- args$min_s2n
lower_del_cutoff <- args$lower_del_cutoff
upper_del_cutoff <- args$upper_del_cutoff
lower_dup_cutoff <- args$lower_dup_cutoff

output_dir <- args$output_dir
output_dir <- normalizePath(output_dir)

raw_calls_df <- read.table(raw_calls, sep="\t", header=F, fill=T, check.names=FALSE)
coverage_data_df <- read.table(coverage_data, sep="\t", header=T, check.names=FALSE)
references_df <- read.table(references, sep="\t", header=F, fill=T, check.names=FALSE)
references_df <- references_df[-1,]

header = names(coverage_data_df)
samples = header[6:length(header)] 

na.omit.list <- function(y) { 
    return(y[!sapply(y, function(x) all(is.na(x)))]) 
}

# Removing unwanted strings
references_df[] <- lapply(references_df, function(x) gsub("\\(.*\\)", "",x)) 
reference_dict = list()
for (i in 1:nrow(references_df)) {
    sample <- unlist(references_df[i,2])
    others <- unname(unlist(references_df[i, 3:ncol(references_df)]))
    others<- others[!grepl("^ref*", others)]
    others<- others[!grepl("^0\\.[0-9]+", others)]
    others<- others[!grepl("^[[:space:]]*$", others)] 
    others <- na.omit.list(others)
    reference_dict[[sample]] = others
}

gr_obj <- makeGRangesFromDataFrame(coverage_data_df, keep.extra.columns=TRUE)
 
calculateRatio <- function(Data, case, controls) {
    cov_case     <- as.numeric(Data[case])
    cov_controls <- Data[controls] 
    median_cov   <- as.numeric(median(na.omit(cov_controls)))
   
    ratio_case <- 0
    if (median_cov > 0) {
        ratio_case <- cov_case/median_cov
    }
    else {
        ratio_case <- 0
    }
    return (ratio_case) 
} 

zscore <- function(data, x) {
    mean_data <- as.numeric(mean(data))
    sd_data <- as.numeric(sd(data))
    if (sd_data == 0) {
        sd_data = 0.1
    } 
    z_score <- (x-mean_data)/sd_data
    return(z_score)
} 

calculateRatioControls <- function(Data, controls, reference_dict ) {

    control_ratios <- c()
    sample_ratios <- list()
    for(control in controls) {
        references    <- reference_dict[[control]]  
        cov_control   <- as.numeric(Data[control])
        cov_reference <- as.numeric(Data[unlist(references)])
        median_reference <- as.numeric(median(na.omit(cov_reference)))
        ratio <- 0
        if (median_reference > 0) {
            ratio <- cov_control/median_reference
        }
        else {
            ratio <- 0
        }
        control_ratios<-append(control_ratios, ratio)
        sample_ratios[[control]] <- ratio
    }
    return(sample_ratios)
} 


signal2noise<-function(ratio_data) 
{
    median_ratio <- as.numeric(median(ratio_data))
    sd_ratio     <- as.numeric(sd(ratio_data))

    if (sd_ratio == 0) {
        sd_ratio = 0.1
    }
    signal2noise <- as.numeric(median_ratio/sd_ratio)
    return(signal2noise)
} 

plotSingleExon<-function(df, old_sample_name, sample_name, exon, output_dir) {

    output_png <- paste(output_dir, "/", old_sample_name, ".", exon, ".png", sep="")
    png(output_png, res =180, width = 1200, height=600)
    max_ratio <- max(df$Ratio)

    #class_factor <- factor(df$Classe, levels=c("Controls", sample_name))
    #group.colors <- c("Controls" = "darkgrey", sample_name = "#ffb2b2")
    myplot<-ggplot(df, aes(x=start, y=Ratio, group=Samples )) + 
        geom_line(aes(color=Classe, size=Classe))+ ylim(0, max_ratio+0.5) +
        theme(plot.title = element_text(hjust = 0.5))+ xlab("Genomic position") +
        ylab("Ratio")+ theme_bw()+ scale_colour_manual("Class", labels=c("Controls", old_sample_name), values=c("darkgrey", "red")) +
        ggtitle(exon)+ scale_size_manual(values=c(0.5, 1)) + guides(size=FALSE) +
        geom_hline(yintercept = 0.5, colour="red", linetype="dashed", size=0.2) + 
        geom_hline(yintercept = 1.5, colour="blue", linetype="dashed", size=0.2) + 
        theme(plot.title = element_text(hjust = 0.5)) 
    myplot
    print(myplot)
    dev.off()
} 

for (row in 1:nrow(raw_calls_df)){

    chromosome <- raw_calls_df[row,1]
    start      <- raw_calls_df[row,2]
    end        <- raw_calls_df[row,3]
    exon       <- raw_calls_df[row,4]
    sample_name<- raw_calls_df[row,7] 

    q = GRanges(seqnames=chromosome, ranges=IRanges(start=start, end=end))
    results<-subsetByOverlaps(gr_obj , q)
    gr_df <- as.data.frame(results, check.names=FALSE, header=TRUE)
    #attach(gr_df)
    #print(sample_name)
    names(gr_df)<- c("chr", "start", "end" , "width", "strand", "exon", "X.GC", samples)

    # Get the case sample
    case_sample <- gr_df$sample_name
 
    vec_case_sample <- c(sample_name)

    # Get the control samples for baseline calculation 
    control_samples <- setdiff(samples, vec_case_sample)

    output_list <- list()

    ratio_data_case <- c()
    ratio_data_controls <- c()
    refs_case <- reference_dict[[sample_name]]

    sample_header <- c(vec_case_sample, refs_case)

    for (i in 1:nrow(gr_df)){

        # Case ratio
        cov_case <- as.numeric(gr_df[i, sample_name]) 
        cov_controls <- as.numeric(gr_df[i, refs_case]) 
        median_controls <- as.numeric(median(na.omit(cov_controls)))
        ratio_case <- 0
        if (median_controls > 0){
            ratio_case <- cov_case/median_controls
        } 
        ratio_data_case <- append(ratio_data_case, ratio_case)

        # Now for controls 
        control_ratios_row <- c()
        for (control in refs_case) {
            refs <- reference_dict[[control]]
            cov_control  <- as.numeric(gr_df[i, control])
            cov_ref      <- as.numeric(gr_df[i, refs])
            median_ref   <- as.numeric(median(na.omit(cov_ref)))
            ratio_control<- 0
            if (median_ref > 0) {
                ratio_control <- cov_control/median_ref
            }
            control_ratios_row <- append(control_ratios_row, ratio_control)
            ratio_data_controls <- append(ratio_data_controls, ratio_control)
        } 
        all_ratios <- c(ratio_case, control_ratios_row)
        output_list[[i]]<-all_ratios 
    }

    # For the case sample 
    median_ratio_case <- round(median(ratio_data_case), 3)
    signal_to_noise_case <- round(signal2noise(ratio_data_case), 3)
    # For the control samples 
    median_ratio_controls <- round(median(ratio_data_controls), 3)
    signal_to_noise_controls <- round(signal2noise(ratio_data_controls), 3)

    zscore_case <- round(zscore(ratio_data_controls, median_ratio_case), 3)

    cat(" INFO: ", exon, " z-score=", zscore_case, " case_ratio=", median_ratio_case, " control_ratio=", median_ratio_controls, "\n")
    if ( (signal_to_noise_controls > min_s2n) & (signal_to_noise_case > min_s2n) & (abs(zscore_case) > min_zscore) ){

        if(((median_ratio_case < upper_del_cutoff) & (median_ratio_case > lower_del_cutoff))|((median_ratio_case > lower_dup_cutoff))) {


            out_file <- paste(output_dir, "/", sample_name, ".single.exon.cnv" , ".bed", sep="")
            out_df <- do.call("rbind",output_list)
            colnames(out_df) <- sample_header

            out_df <- cbind(chr=gr_df$seqnames, start=gr_df$start, end=gr_df$end, out_df)
            out_df <- as.data.frame(out_df)

            # Transform to long format
            out_df <- gather_(data=out_df, key="Samples", value="Ratio", sample_header)


            # Add a new columns assigning class category 
            out_df <- out_df %>%         
                mutate(Classe = if_else(out_df$Samples == sample_name, paste("2_", sample_name, sep=""), '1_Controls'))

            new_sample_name <- paste("2_", sample_name, sep="")

            # Plot single exon CNV 
            plotSingleExon(out_df, sample_name, new_sample_name, exon, output_dir)
            svtype <- ""
            if (median_ratio_case < upper_del_cutoff) {
                svtype <- "DEL"
            }
            else{
                svtype <- "DUP"
            }  
            variant_df<-data.frame(chromosome, start, end, svtype, exon, median_ratio_case, signal_to_noise_case, signal_to_noise_controls, zscore_case)
            variant <-c(chromosome, start, end, svtype, exon, median_ratio_case, signal_to_noise_case, signal_to_noise_controls, zscore_case)
            write.table(variant_df, file=out_file, col.names=FALSE, row.names=FALSE, append=TRUE, quote=FALSE, sep="\t")
        } 
    }   


} 