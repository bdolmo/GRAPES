# What is GRAPES?
 GRAPES stands for "Genomic Rearrangement Analysis of Panel Enrichment Sequencing".
 It is yet another SV and CNV caller, but with new interesting features (Keep reading!)

## Features
* SV and CNV calling for targeted sequencing: using data from on-target, off-target and breakpoint-informative reads.
* Clustering of samples with highly correlated read depth for improved CNV calling specificity.
* Dynamic update of sample references using SQLite.

#### What's new in GRAPES 0.9.5 ? March-2021
* New scoring metric to evaluate the confidence of each CNV.
* Improved single-exon CNV detection.
* Improved speed performance.
* BAF extraction on the fly.
* Several bug fixes.

![Slide](img/Figure1.png)

## Installation

### Recommended Installation (Docker)
 There is a docker repository available at:
```
docker pull bdolmo/grapes
```
And run commands with:
```
docker run -it bdolmo/grapes:latest GRAPES
```

### Source installation
 GRAPES will only work on Unix-based systems and supports only human genome.
 Before building and installing GRAPES, make sure you have available on path:
* Perl
* R (>= 3.3). In addition Rscript must be accessible
* g++ compiler (>= 4.7)
* Boost C++ library
* openMP C++ API
* [BEDtools](https://github.com/arq5x/bedtools2)
* [SAMtools](http://www.htslib.org/)
* [Tabix](https://github.com/samtools/tabix)
* [macs2](https://github.com/taoliu/MACS)
* SQLite3
* GNU core utils:  wget, awk, sort, cat, grep, head, tail, sed, cut, paste, uniq, wc and mv.
  macOS users: need to install Homebrew (https://brew.sh/) and the latest XCode.

#### Perl modules
* Parallel::ForkManager
* Sort::Key::Natural
* Statistics::Descriptive
* JSON::MaybeXS
* Excel::Writer::XLSX

You can install them through CPAN:
```
cpan Parallel::ForkManager Sort::Key::Natural Statistics::Descriptive JSON::MaybeXS Excel::Writer::XLSX
```
Then you can download and install the latest release:
```
 git clone --recursive https://github.com/bdolmo/GRAPES.git
 cd GRAPES
 ./INSTALL.PL
```
INSTALL.PL will compile all C++ code along with the required R packages for segmentation and plotting.
In addition it will download both 75-mer and 100-mer mappability tracks for GRCh37 and GRCh38 and gnomaAD annotation files.


## Commands
### Targeted sequencing analysis: ```GRAPES wes```

#### Workflow 1: Pooled Analysis (creates a reference using all available samples):
 ```
 GRAPES wes --pooled <bam_dir> --outdir <output_dir> --bed <roi> --genome <genome_fa> --all
 ```

(docker)
```
docker run -t -i \
-v $HOME/BAM_FOLDER:/bam_folder \
-v $HOME/BED_FOLDER:/bed_folder \
-v $HOME/GENOME_FOLDER:/genome_folder \
-v $HOME/OUTPUT_FOLDER:/output_folder \
-it bdolmo/grapes:latest GRAPES wes \
-all -pooled /bam_folder/ -b /bed_folder/targets.bed -g /genome_folder/genome.fa -o /out_dir/ -t 4
```

#### Workflow 2: Case-Control analysis:
```
GRAPES wes --cases <bam_dir> --control <bam_dir> --outdir <output_dir> --bed <roi> --genome <genome_fa> --all
```

(docker)
```
docker run -t -i \
-v $HOME/CASE_FOLDER:/case_folder \
-v $HOME/CONTROL_FOLDER:/control_folder \
-v $HOME/BED_FOLDER:/bed_folder \
-v $HOME/GENOME_FOLDER:/genome_folder \
-v $HOME/OUTPUT_FOLDER:/output_folder \
-it bdolmo/grapes:latest GRAPES wes \
-all -cases /case_folder/ -controls /control_folder/ -b /bed_folder/targets.bed -g /genome_folder/genome.fa -o /out_dir/ -t 4
```

##### I/O:
```
-o,--outdir	        STRING	  Output directory
-g,--genome_fasta	  STRING    Genome reference in FASTA format
-r,--genome_version	STRING    Genome version. Choose: hg19, hg38 (default = hg19)
-b,--bed	          STRING	  Regions file in BED format
-t,--threads	        INT	    Number of CPUs (default = 1)
```

##### Options:
```
--all	                     Perform all steps below
--refdir	                     Input directory where DB references will be stored.
--breakpoint	                     Perform Breakpoint analysis
  --nobreakpoint	                     Turn off breakpoint analysis
--extract	                     Extract Depth, GC and Mappability
--offtarget	                     Perform Off-target analysis
  --noofftarget	                     Turn off offtarget analysis
--buildref	                     Build a reference from a pool of samples
--callcnv	                 Segment and call CNVs
  --nocallcnv	             Turn off CNV calling
--normalize	               Normalize read depth. Choose from 'median', 'PCA'. (default='median')
--plotsingleexon	         Plot single exon CNVs
  --noplotcnv	             Turn off CNV plotting
--plotlargecnv	           Plot segmented CNVs
  --noplotcnv	             Turn off CNV plotting
--plotscatter	             Plot genome-wide CNV scatter plot
--vaf	                     Include Variant-Allele Frequency (VAF) analysis
  --novaf	                 Turn off VAF analysis
--samtools	               Default if --vaf set. Perform variant call with samtools
--freebayes	               if --vaf set, perform variant call with freebayes
--annotate	               Annotate VCF
  --noannotate	           Turn off VCF annotation
--filtervcf	               Filter low qual VCF entries
  --nofiltervcf	           Turn off VCF filtering
--reporthtml	             Write results to an HTML file (Only for gene panels)
  --noreporthtml	         Turn off report HTML creation
--filterdiscordantonly     Filter discordant-only SV predictions (default = true)
  --nofilterdiscordantonly Turn off discordant-pair only filtering
--verbose	                 Print sub-command messages

 ```
##### Tuning parameters:
  ```
   --mincorr            FLOAT	  Minimum pairwise-correlation to build a reference set (default = 0.91)
   --minrefsize         INT 	  Default minimum number of samples to build a single baseline (default = 2)
   --maxrefsize         INT   	Default maximum number of samples to build a single baseline (default = 15)
   --minzscore          FLOAT   Minimum Z-score required to output a CNV prediction (default = 2.58)
   --pcavariance	      FLOAT 	Variance to remove when normalizing by PCA (default = 0.7)
   --lowerdelcutoff     FLOAT	  Lower-bound deletion cutoff ratio (default = 0.35)
   --upperdelcutoff     FLOAT	  Upper-bound deletion cutoff ratio (default = 0.71)
   --lowerdupcutoff     FLOAT	  Lower-bound duplication cutoff ratio (default = 1.24)
   --minofftargetreads	INT	    Default minimum number of offtarget reads required to trigger off-target analysis (default = 1e6)
   --minofftargetsd	    FLOAT	  Default minimum number of std.dev from off-target rartios trigger off-target analysis	(default = 0.2)
   --minsvsize          INT     Minimum SV size to report a breakpoint call (default = 15)
   --maxsvsize          INT     Maximum SV size to report a breakpoint call (default = 5000000)
   --mindiscordants     INT	    Minimum number of discordant read pairs (default = 5)
   --mindiscordantssd   INT     Minimum std.deviations from the mean insert size to consider discordant pairs (default = 10)
   --breakreads         INT	    Minimum number of break reads (default = 5)
  ```

### SV annotation: ```GRAPES annotate```
 ```
 GRAPES annotate --input_vcf <vcf> --output_dir <output_dir> --r_overlap <reciprocal_overlap>
 ```
##### Tuning parameters:
 ```
 -i,--input_vcf  STRING  Input VCF file
 -o,--output_dir STRING  Output directory (default=.)
 -l,--r_overlap  FLOAT   Reciprocal overlap fraction (default=0.5)
 ```
