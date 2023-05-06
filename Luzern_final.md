# Characterization of genital microbiota in infertile couples

```
require(knitr)
opts_knit$set(root.dir = normalizePath('../'))
opts_chunk$set(tidy.opts=list(width.cutoff=80),tidy=TRUE)
knitr::opts_chunk$set(message=FALSE,echo=TRUE,eval=FALSE)
```

## Description

In this study, we want the characterize genital microbiota of infertile couples.

Samples include:

\- vaginal swab

\- follicular liquid

\- semen

\- penis swab

The analysis was performed by sequencing the variable regions V1-V2 of the 16S rRNA gene amplified using custom barcoded primers (F-27/R-338) containing Illumina sequencing adapters.

Illumina sequencing was performed at the Lausanne Genomic Technologies Facility (GTF) of the Lausanne University using an Illumina MiSeq instrument in paired-end mode 2 x 250 nt.

Demultiplexing of the raw sequencing data was performed with illumina-utils package. Reads were processed with the DADA2 pipeline. Decontam package was used to remove potential contaminations based on negative controls.Phyloseq and ampvis2 R packages were used for analysis and graphical visualization of the microbiota results. R package maaslin2 was used to infer differentially abundant bacteria.

## Data

\- Illumina adapter file

\- Raw data location contains:

\- fastq.gz

\- metadata

# 1 - Reads processing

## 1.1 - Demultiplex reads

#Install and activate illumina utils - create a Python 3 virtual environment:
```
mkdir -p ~/virtual-envs/
virtualenv ~/virtual-envs/illumina-utils-v2.7
source ~/virtual-envs/illumina-utils-v2.7/bin/activate
python --version
##Make sure the output starts with a '3'. Then continue with installation:

conda install illumina-utils
#pip install illumina-utils
```
#Activate illumina-utils in virtualenv:
```
echo 'alias illumina-utils-activate-v2.7="source ~/virtual-envs/illumina-utils-v2.7/bin/activate"' >> ~/.bash_profile

cd luzern2021
```
#Check the correct filename and decompress files
```
gunzip Luzern1_L1_R1_001.fastq.gz
##tar -zxvf Unaligned-BIRTH01.tar
```

#create output folder based on the run number
```
mkdir demultiplexed_2

gunzip *R1*.fastq.gz *R2*.fastq.gz *I1*.fastq.gz

iu-demultiplex -s ~/luzern2021/barcode_to_sample_4.txt --r1 *R1*.fastq.gz --r2 *R2*.fastq.gz -i *I1*.fastq.gz -x -o ~/luzern2021/demultiplexed_4
```
#Repeat for all the runs

#Move files to the correct location

```
mkdir ~/luzern2021/01_raw_sequences

cp ~/luzern2021/demultiplexed*/*fastq ~/luzern2021/01_raw_sequences

```

## 1.2 - Reads count

#Count the number of reads in each files and store them in ReadsCount.csv. This is performed in bash.

```
cd ~/luzern2021/01_raw_sequences

ls *.fastq| while read LINE; do

        echo -n $LINE ;

        echo $(cat < $LINE | wc -l)/4|bc ;

done |sed 's/\.fastq/,/g' > ~/luzern2021/ReadsCount_rawsequences.csv
```

## 1.3 - Quality control of the reads

Assessment of the quality of the raw reads with fastqc.

Install and activate multiqc in conda

```{bash}
mkdir ~/luzern2021/02_FastQC_RAW
for files in *.fastq; do fastqc $files -o ~/luzern2021/02_FastQC_RAW ; done
multiqc ~/luzern2021/02_FastQC_RAW
```

Comments:

\- Overall bases quality is good

\- presence of sequence duplication levels ; over represented sequences

\- Per base sequence content and per sequence GC content failed QC

## 1.4 - DADA2 pipeline

Load required packages

```
# Load the libraries
library(dada2)
library(ShortRead)
library(Biostrings)
library(ggplot2)
library(biomformat)
library(microbiome)
library(DESeq2)
```

Save

```
getwd()
save.image("~/luzern2021/luzern2021_DADA2_.RData")
load("~/luzern2021/luzern2021_DADA2.RData")
```

## 1.5 - Get files path and sample names

```
library(dada2); packageVersion("dada2")
path<-"~/luzern2021/01_raw_sequences/"
# List all the files in trimmed directory
list.files(path)
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.trim.fastq and SAMPLENAME_R2_001.fastq
FWDfiles <- sort(list.files(path, pattern="-R1.fastq", full.names = TRUE))
REVfiles <- sort(list.files(path, pattern="-R2.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(FWDfiles), "-R"), `[`, 1)
str(sample.names)
write.csv(sample.names,"samplenames.csv")
```

## 1.6 - Quality scores

Median quality score is the green line. Quartile quality scores are the orange lines.

The red line (bottom) is the proportion of reads that reach the position (length).

The overall quality of the reads is good, median and quartiles quality scores are above 30 (phred score).

Generally, reverse reads have 'lower' scores at the 5' and 3' ends, but for the rest of the sequence scores were good.

```
# Quality scores of R1 reads (partial)
plotQualityProfile(FWDfiles[1:3]) 
# Quality scores of R2 reads (partial)
plotQualityProfile(REVfiles[24:26])
```

## 1.7 - Trim the data

Comments:

\- Trimming should be adapted to the data-type and quality of the reads.

\- \`truncLen\` must be large enough to maintain an overlap between forward and reverse reads of at least \`20 + biological.length.variation\` nucleotides.

\- We will try two trimming settings: truncLen=c(180,140) and truncLen=c(200,160)

```
# Place filtered files in filtered/subdirectory
path<-"~/luzern2021"
filtFWD <- file.path(path,"05_filtered", paste0(sample.names, "_F_filt.fastq"))
filtREV<- file.path(path,"05_filtered", paste0(sample.names, "_R_filt.fastq"))
names(filtFWD) <- sample.names
names(filtREV) <- sample.names
length(filtFWD)
out1 <- filterAndTrim(FWDfiles, filtFWD, REVfiles, filtREV, truncLen=c(180,140), 
                     maxN=0, maxEE=c(2,4), truncQ=2, rm.phix=TRUE,      
                     compress=TRUE, multithread=TRUE) 

out2 <- filterAndTrim(FWDfiles, filtFWD, REVfiles, filtREV,truncLen=c(200,160), 
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE)
# On Windows set multithread=FALSE
# all parameters but truncLen are default DADA2 params

#any(duplicated(c(FWDfiles, filtFWD)))
head(out1)
head(out2)

save.image("luzern2021_DADA2.RData")
```

Move filtered reads to a new folder

```{bash}
mkdir ~/luzern2021/05_filtered_sequences

mv ~/luzern2021/04_trimmed_sequences/filtered/* ~/luzern2021/05_filtered_sequences
```

## 1.8 - \`derepFastq\` Dereplication step

All identical sequences are combined in "unique sequences" that are associated with "abundance" (number of reads that have this unique sequence)

```
setwd("~/luzern2021/05_filtered")
derepFWD <- derepFastq(filtFWD)
derepREV <- derepFastq(filtREV)
sam.names <- sapply(strsplit(basename(filtFWD),"_"),`[`,1)
names(derepFWD) <- sam.names
names(derepREV) <- sam.names
save.image("~/luzern2021/luzern2021_DADA2.RData")

```

## 1.9 - Learn the error rates

Data is used to model the probability of transitions and transversions (errors)

in function of the read quality.

Each run has its specific error rates (cannot combine data from two different runs)

\- black dots : observed error rates for each consensus quality score.

\- black lines : estimates error rate after convergence of the algorithm

\- red line : error rates expected under the nominal definition of the Q-score

! Parameter learning is computationally intensive, so by default the learnErrors function uses only a subset of the data (the first 100M bases = 1e8). If you are working with a large dataset and the plotted error model does not look like a good fit, you can try increasing the nbases parameter to see if the fit improves !

```
sys_str <- Sys.time()
errF <- learnErrors(derepFWD, randomize=TRUE,nbases = 5e+08 ,multithread=TRUE)  
errR <- learnErrors(derepREV, randomize=TRUE,nbases = 5e+08, multithread=TRUE) 
plotErrors(errR, nominalQ=TRUE)
# In the plots, the black line is the error model, the dots are the actual errors
sys_str[2] <- Sys.time()
sys_str
#save.image("luzern2021_DADA2.RData")
rm(sys_str)
```

## Comments: Overall very good reads quality

## 1.10 - Sample inference

The DADA2 algorithm divides the reads in ASVs

```
sys_str <- Sys.time()
dadaFs <- dada(filtFWD, err=errF, multithread=TRUE) # we need to incorporate "selfconsist" and "pool=TRUE"
dadaRs <- dada(filtREV, err=errR, multithread=TRUE) # we need to incorporate "selfconsist" and "pool=TRUE"
sys_str[2] <- Sys.time()
sys_str
rm(sys_str)
#save.image("luzern2021_DADA2.RData")
dadaFs[[1]]
dadaRs[[1]]
```

## 1.11 - Merging paired reads

Merging reads to obtain full denoised sequences. Merged sequences are output if the overlap is at least of 12 \*identical\* nucleotides.

! Most of your reads should successfully merge. If that is not the case upstream parameters may need to be revisited: Did you trim away the overlap between your reads? !

merger contains a list of data.frames. Each data.frame contains the merged \`\$sequence\`, \`\$abundance\`, the indices of FWD and REV sequences variant that were merged. Paired-reads that did not exactly match were removed by the \`mergePairs\` function.

```
mergers <- mergePairs(dadaFs, derepFWD, dadaRs, derepREV, verbose=TRUE, trimOverhang=TRUE)
#save.image("luzern2021_DADA2.RData")
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
```

## 1.12 - Construct sequence table

! some sequences may be shorter or longer than what is expected - here \~250 bp.

```
seqtab <- makeSequenceTable(mergers)

dim(seqtab)

# Inspect distribution of sequence lengths
plot(table(nchar(getSequences(seqtab))))
```

## 1.13 - Remove chimeras

Chimeric sequences are identified if they can be exactly reconstructed by combining a left-segment and a right-segment from two more abundant "parent" sequences.

```
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

dim(seqtab.nochim)

rownames(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

# save sequences

sequences <- data.frame(colnames(seqtab.nochim))

colnames(sequences) <- 'sequences'
write.csv(sequences, 'SB_sequences.csv')
save.image("luzern2021_DADA2.RData")
```

## 1.14 - Track the number of reads after each filtering steps (Figure 2)

Create one file containing read counts from raw data and post-trimmomatic

```{bash}
#awk 'BEGIN{FS=","; OFS=","} FNR==NR{a[FNR]=$2;next};{print $0, a[FNR]}' 01_trimmed/readsCounts2.csv  SOURCE DIRECTORY/ReadsCount.csv | grep 'R1' > #ReadsCounts.csv
```

```
library(ggplot2)
library(reshape)

setwd("~/luzern2021")

preDADA2 <-read.table('ReadsCount_rawsequences.csv', sep=',', header = FALSE)

colnames(preDADA2) <- c('sample', 'raw')

getN <- function(x) sum(getUniques(x))

reads_counts <- cbind(out1, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

colnames(reads_counts) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")

rownames(reads_counts) <- sample.names

head(reads_counts)

# check dimensions
 dim(reads_counts)
 colSums(reads_counts)
 a<-colSums(reads_counts)
 cat("Percentage of sequences removed after filtering:", 1-(a[6]/a[1]))

write.csv(reads_counts,"~/luzern2021/reads_tracking.csv", row.names = TRUE)

reads_counts.melt<-melt(reads_counts)

colnames(reads_counts.melt)<-c("sample","step","reads")

medians<-as.numeric(apply(reads_counts,2,median))

retained_reads<-100/medians[1]*medians[6]

ggplot(reads_counts.melt, aes(x=step, y=reads)) + 
    geom_boxplot() + scale_y_continuous(trans='log10') 
    ggtitle("Read counts")

reads_counts$input

class(reads_counts)
```

## The filtering steps (trimmings, denoising, removal of artifact and chimeras) removed \~15% of the initial reads.

## 1.15 - Track reads through the pipeline

```
ReadsTracking <- read.csv("~/luzern2021/reads_tracking.csv")

div <- function(x,y) (x/y)*100

lostReads <- (1-(ReadsTracking[,-c(1)]/ReadsTracking$input))*100

averageLost <- mean(lostReads$nonchim)

lostperSpecies <- lostReads$nonchim

names(lostperSpecies) <- ReadsTracking$sample.names

averageLost

lostperSpecies
```

# 2 - Assign taxonomy

Fasta release files from the UNITE ITS database can be used as is. To follow along, download the silva_nr_v132_train_set.fa.gz

```
#assignTaxonomy using DADA2/Silva
setwd("~/luzern2021")

taxa <- assignTaxonomy(seqtab.nochim, "~/SILVA/silva_nr_v132_train_set.fa", multithread=TRUE)

taxa <- addSpecies(taxa,"~/SILVA/silva_species_assignment_v132.fa")

taxa.print <- taxa # Removing sequence rownames for display only

rownames(taxa.print) <- NULL

head(taxa.print)

path_trim<-paste(path, "06_Taxonomy", sep="/")

dir.create(path_trim)

write.csv2(file=paste(path_trim, "Taxtable_dada2.csv", sep="/"),taxa)

write.csv2(file=paste(path_trim, "ASV_sequences.csv", sep="/"),seqtab.nochim)

save.image("~/luzern2021/luzern2021_DADA2.RData")
```

## 2.1 - Convert to fasta and run on SILVA

Fasta files have sequence as header (so that header is unique)

This was done locally, as the limit for online upload was reached

```{bash eval=FALSE}
awk -F';'  'NR>1{ print  ">ASV"++i "\n" $1 }' ~/luzern2021/06_Taxonomy/Taxtable_dada2.csv  > ASV_sequences.fasta

sed 's/\"//g' ASV_sequences.fasta | sed 's/NA//g' > ASV_sequences2.fasta

#run locally
sina -i ASV_sequences2.fasta -o ASV_sequences2_aligned.csv --meta-fmt csv --db SILVA_138.1_SSURef_NR99_12_06_20_opt.arb --search --search-db SILVA_138.1_SSURef_NR99_12_06_20_opt.arb --lca-fields tax_slv

# Download result file - trim file for SINA taxonomy

# arb-silva.de_align_resultlist_867879.csv

cd ~/luzern2021/06_Taxonomy

echo '"","Kingdom","Phylum","Class","Order","Family","Genus","Species"' > sina_taxonomy.csv

cut -f 1,8 -d , ASV_sequences2_aligned.csv | sed 's/;/,/g' | sed '1d' >> sina_taxonomy.csv
```

# 3 - Create Phyloseq object

```
library(ggplot2)
library(vegan) 
library(dplyr)
library(scales) 
library(grid)
library(reshape2) 
library(cowplot)
library(phyloseq)
library(tidyverse)
library(Biostrings)
library(remotes)
library(stringr)
library(gridExtra)
library(microbiome)

setwd("~/luzern2021/07_output")

# Set plotting theme
theme_set(theme_bw())

#Data frame containing sample information
samdf = read.table(file="~/luzern2021/00_metadata/luzern_metadata_final_290322.csv", 
                   sep=",",header = T, fill=TRUE) 
                  # fill=TRUE allows to read a table with missing entries

rownames(samdf) = samdf$sample
head(samdf)
samdf$qpcr = as.numeric(samdf$qpcr)

# Import SINA taxonomy

taxa.sina <- read.csv(file="~/luzern2021/06_Taxonomy/sina_taxonomy.csv",sep=",",header = TRUE, row.names = 1)

#check dimensions
dim(taxa.sina)
dim(seqtab.nochim)

namestochange<-row.names(seqtab.nochim)
namestochange<-str_remove(namestochange, "_F_filt.fastq")
row.names(seqtab.nochim)<-namestochange

#Create a phyloseq object
ps_raw <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=F), 
               sample_data(samdf), 
               tax_table(taxa))

#check phyloseq object
otu_table(ps_raw)
sample_data(ps_raw)
sample_names(ps_raw)
tax_table(ps_raw)

# save sequences as refseq and give new names to ASV's
dna <- Biostrings::DNAStringSet(taxa_names(ps_raw))

names(dna) <- taxa_names(ps_raw)

ps_raw <- merge_phyloseq(ps_raw, dna)

taxa_names(ps_raw) <- paste0("ASV", seq(ntaxa(ps_raw)))
```

## 3.1 - Add tree to the phyloseq object

## 3.1.1 - Select sequences and perform Sina alignment

```{bash}
conda activate sina

sina -i ~/luzern2021/06_Taxonomy/ASV_sequences2.fasta -o ~/luzern2021/06_Taxonomy/ASV_aligned_sequences.fna –meta-fmt csv –db ~/luzern2021/06_Taxonomy/SILVA_138.1_SSURef_NR99_12_06_20_opt.arb –search –search-db ~/luzern2021/06_Taxonomy/SILVA_138.1_SSURef_NR99_12_06_20_opt.arb –lca-fields tax_slv

#The alignment file in FastTree using aligned sequences

FastTree -nt ~/luzern2021/06_Taxonomy/ASV_aligned_sequences.fna > ~/luzern2021/06_Taxonomy/FastTree_phyloseq

```

## 3.1.2 - Import tree file in the Phyloseq object

```

library("ape")

#import FastTree file
tree_phyloseq <- read.tree("~/luzern2021/06_Taxonomy/FastTree_phyloseq")

#import the tree file in the phyloseq object

ps = phyloseq::merge_phyloseq(ps_raw, sampledata, tree_phyloseq)
#remove community sample controls, mitochondria and chloroplast

ps <- phyloseq::subset_samples(ps, sample_type!="community")
ps <- ps %>% subset_taxa( Family!= "mitochondria" | is.na(Family) & Class!="Chloroplast" | is.na(Class) )

rank_names(ps)
```

## 3.2. - Export tax_table

```
# Export ASV table

table = merge(tax_table(ps),t(otu_table(ps)), by="row.names")

write.csv(table, "~/luzern2021/06_Taxonomy/ASVtable.csv")

# Export to FASTA with Biostrings

writeXStringSet(refseq(ps), "phyloseq_ASVs.fasta",append=FALSE, format="fasta")

# Then align it with SINA/SILVA, and edit the taxonomy table to get it

# To save phyloseq objects:

# With Biostrings

writeXStringSet(refseq(ps), "outfile.fasta",append=FALSE, format="fasta")

# With seqRFLP

#install.packages("seqRFLP")

#names <- table$Genus

#sequences <- table$Row.names

#dfasta <- data.frame(names,sequences)

#df.fasta = dataframe2fas(dfasta, file="df.fasta")

# Export Taxtable

write.csv2(data.frame(tax_table(ps)),file="taxtable.csv")
```

## 3.3. - Plot gel_band vs qPCR result

```
table.band<-samdf[,c(1,3,7,8,14)]

plot.gelband.qpcr <- ggplot(table.band, aes(x=gel_band, y=qpcr))+
      geom_boxplot()+
      geom_point(aes(fill = sample_type), size = 2, shape = 21, na.rm=TRUE,position="dodge")+
      scale_y_log10() + expand_limits(y = c(0, 1e10))+
      scale_fill_brewer(palette="Set1")+
      ggtitle("Correlation gel_band-qpcr")

plot.gelband.qpcr
```

## 3.4. - Plot gel_band vs qPCR result (Figure 1)

```
table.band<-samdf[,c(1,3,7,8,14)]

plot.gelband.qpcr <- ggplot(table.band, aes(x=gel_band, y=qpcr)) +
  geom_boxplot()+
  geom_point(aes(fill = sample_type), size = 2, shape = 21, na.rm=TRUE,position="dodge")+
  scale_y_log10() + expand_limits(y = c(0, 1e10))+
  scale_fill_brewer(palette="Set1")+
  ggtitle("Correlation gel_band-qpcr")

plot.gelband.qpcr
table.band<-samdf[,c(1,3,7,8,14)]
plot.gelband.qpcr <- ggplot(table.band, aes(x=gel_band, y=log10(qpcr))) +
  geom_boxplot()+
  geom_point(aes(fill = sample_type), size = 2, shape = 21, na.rm=TRUE,position="dodge")+
  #scale_y_log10() + expand_limits(y = c(0, 1e10))+
  scale_fill_brewer(palette="Set1")+
  ggtitle("Correlation gel_band-qpcr")

plot.gelband.qpcr
```

## 3.5. - Plot qPCR for each sample_type (Figure 4)

```
qpcr_type<-data.frame(samdf$sample_type,samdf$qpcr,samdf$is_sample)
colnames(qpcr_type)<-c("sample_type","qpcr","is_sample")
qpcr_type<-subset(qpcr_type,is_sample=="yes")
qpcr_type<-subset(qpcr_type,qpcr>200)

qpcr_type_plot <- ggplot(qpcr_type, aes(x=sample_type, y=qpcr,fill=sample_type)) +
  geom_boxplot() + scale_y_log10() + expand_limits(y = c(0, 1e10))+
  scale_fill_brewer(palette="Set1") + ggtitle("qpcr - sample_type") +
  geom_hline(yintercept = 300,colour = "red",linetype = "dashed")

qpcr_type_plot
```

# 4 - Check control community samples (Figure 3)

```
#subset
ps_comm<-subset_samples(ps_raw, sample_type=="community")
ps_comm<-tax_glom(ps_comm, taxrank="Genus")

comm_data<-as.data.frame(t(ps_comm@otu_table@.Data))
comm_data$ASV<-row.names(comm_data)
ASVtable<-select(table,Row.names,Genus)
ASVtable$ASV<-ASVtable$Row.names

comm_genus<-full_join(comm_data,ASVtable,by="ASV")
comm_genus$ASV<-NULL
comm_genus$Row.names<-NULL

comm.genera<-c("Lactobacillus","Listeria","Pseudomonas","Staphylococcus",
"Escherichia/Shigella","Salmonella", "Enterococcus","Bacillus")

comm.samples<-c("run1","run2","run3","run4")

select.genera<-filter(comm_genus,Genus %in% comm.genera)
select.genera<-select.genera[1:8,]
tot.all<-(comm_genus[,1:4])
sums.tot<-colSums(tot.all,na.rm = TRUE)
tot.comm<-(select.genera[,1:4])
sums.comm<-colSums(tot.comm,na.rm = TRUE)
new.row<-c((sums.tot-sums.comm),"z_other")
final.table<-rbind(select.genera,new.row)

final.table.temp<-final.table
final.table.temp$Genus<-NULL
sapply(final.table.temp, class)
final.table.temp<-as.data.frame(sapply(final.table.temp, as.numeric))
all.sum<-colSums(final.table.temp,na.rm = TRUE)

#relative abundance table
final.table.temp<-as.data.frame(t(apply(final.table.temp, 1, function(x) x/all.sum))*100)
final.table.temp$expected<-c(15.5,9.9,17.4,14.1,10.1,18.4,10.4,4.1,0)
row.names(final.table.temp)<-final.table$Genus
final.table.temp$Genus<-final.table$Genus
melted.final.table <- melt(final.table.temp, id = "Genus")

runs<-vector()
runs[1:36]<-"runs"

expected<-vector()
expected[1:9]<-"expected"

type<-c(runs,expected)
melted.final.table$type<-type

comm.plot<-ggplot(melted.final.table, aes(x = variable, y = value, fill = Genus)) +
  geom_bar(position = "fill",stat = "identity") +
  scale_y_continuous(labels = scales::percent_format()) +
  facet_grid(. ~ type,scales = "free", space='free') +
  scale_fill_brewer(palette="Set1")

pdf("community_samples.pdf",width = 8,height = 4)
comm.plot
dev.off()

#Check proportion proteobacteria vs firmicutes
ps_phylum<-tax_glom(ps_genus, taxrank=rank_names(ps_genus)[2], NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))
phylum_reads<-ps_phylum@otu_table@.Data
phylum_reads<-colSums(phylum_reads)
phylum_sums<-sum(phylum_reads)
phylum_percent<-phylum_reads/phylum_sums*100
```

# 5 - Run decontam package (Figure 5)

```
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(decontam); packageVersion("decontam")
library(ggpubr)
library(grid)
library(dplyr)
library(reshape)
library(RColorBrewer)

setwd("~/luzern2021/07_output")

#Prepare phyloseq object - remove samples with no reads
phylobj <- prune_samples(sample_sums(ps) >= 1, ps)
phylobj <- subset_samples(phylobj,qpcr!="na")
phylobj <- subset_samples(phylobj,sample_type!="community")


#Inspect Library Sizes
df <- as.data.frame(sample_data(phylobj)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(phylobj)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))

ggplot(data=df, aes(x=Index, y=LibrarySize, color=sample_type)) + geom_point() +   
    scale_y_continuous(trans='log10')

ggplot(data=df, aes(x=Index, y=LibrarySize, color=is_control)) + geom_point() + 
    scale_y_continuous(trans='log10')

sample_data(phylobj)$is.neg <- sample_data(phylobj)$is_control == "yes"

#prevalence method
contamdf.prev <- isContaminant(phylobj, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)
head(which(contamdf.prev$contaminant))
sample_data(phylobj)$is.neg <- sample_data(phylobj)$is_control == "yes"

#combined methods
contamdf.comb2<-isContaminant(
  phylobj,
  conc = "qpcr",
  neg = "is.neg",
  method = "either", #c("auto", "frequency", "prevalence","combined","minimum",both"),
  batch = sample_data(phylobj)$run,
  batch.combine = "product", #"minimum", "product", "fisher"),
  threshold = 0.05,
  normalize = TRUE,
  detailed = TRUE)

ps.pa <- transform_sample_counts(phylobj, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$is_control == "yes", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$is_control == "no", ps.pa)

# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
contaminant=contamdf.comb2$contaminant)
df.pa$asv<-row.names(df.pa)

pdf("decontam_graph.pdf",width = 6,height = 3)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
dev.off()

contaminant_toremove<-subset(contamdf.comb2,contaminant==TRUE)
contaminant_toremove<-row.names(contaminant_toremove)

#check samples vs controls in the most abundant 5 ASV
ASV1.df<-otu_table(phylobj)[,1]
ASV1.df<-as.data.frame(ASV1.df)
ASV1.df$sample<-rownames(ASV1.df)
ASV1.df$type<-"sample"
ASV1.df$type[259:278]<-"ctrl"

plot1<-ggplot(data=ASV1.df, aes(x=sample,y=ASV1)) +
  geom_bar(stat="identity",aes(fill = type)) +
  scale_y_continuous(trans='log10') +
  theme(text = element_text(size=5),axis.text.x = element_text(angle = 90))+
  ggtitle("ASV1")

ASV2.df<-otu_table(phylobj)[,2]
ASV2.df<-as.data.frame(ASV2.df)
ASV2.df$sample<-rownames(ASV2.df)
ASV2.df$type<-"sample"
ASV2.df$type[259:278]<-"ctrl"

plot2<-ggplot(data=ASV2.df, aes(x=sample,y=ASV2)) +
geom_bar(stat="identity",aes(fill = type)) +
scale_y_continuous(trans='log10') +
theme(text = element_text(size=5),axis.text.x = element_text(angle = 90))+
ggtitle("ASV2")

ASV3.df<-otu_table(phylobj)[,3]
ASV3.df<-as.data.frame(ASV3.df)
ASV3.df$sample<-rownames(ASV3.df)
ASV3.df$type<-"sample"
ASV3.df$type[259:278]<-"ctrl"

plot3<-ggplot(data=ASV3.df, aes(x=sample,y=ASV3)) +
geom_bar(stat="identity",aes(fill = type)) +
scale_y_continuous(trans='log10') +
theme(text = element_text(size=5),axis.text.x = element_text(angle = 90))+
ggtitle("ASV3")

ASV4.df<-otu_table(phylobj)[,4]
ASV4.df<-as.data.frame(ASV4.df)
ASV4.df$sample<-rownames(ASV4.df)
ASV4.df$type<-"sample"
ASV4.df$type[259:278]<-"ctrl"

plot4<-ggplot(data=ASV4.df, aes(x=sample,y=ASV4)) +
geom_bar(stat="identity",aes(fill = type)) +
scale_y_continuous(trans='log10') +
theme(text = element_text(size=5),axis.text.x = element_text(angle = 90))+
ggtitle("ASV4")

ASV5.df<-otu_table(phylobj)[,5]
ASV5.df<-as.data.frame(ASV2.df)
ASV5.df$sample<-rownames(ASV2.df)
ASV5.df$type<-"sample"
ASV5.df$type[259:278]<-"ctrl"

plot5<-ggplot(data=ASV5.df, aes(x=sample,y=ASV5)) +
geom_bar(stat="identity",aes(fill = type)) +
scale_y_continuous(trans='log10') +
theme(text = element_text(size=5),axis.text.x = element_text(angle = 90))+
ggtitle("ASV5")
```

## 5.1 - Remove taxa identified by decontam

```
head(which(contamdf.comb2$contaminant))
badTaxa = contaminant_toremove
allTaxa = taxa_names(ps)
allTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
ps = phyloseq::prune_taxa(allTaxa, ps)
```

# 6 - Remove samples with less then 1000 reads

```
ps_1000 = phyloseq::prune_samples(sample_sums(ps)>=1000, ps)
samdf_1000<-meta(ps_1000)
```

## 6.1. - Create plot kept vs removed samples (Figure 6)

```
OTUtable.filt<-ps@otu_table@.Data
OTUtable.filt<-as.data.frame(t(OTUtable.filt))

OTUtable.filt.sums<-colSums(OTUtable.filt)
OTUtable.filt.sums<-sort(OTUtable.filt.sums, decreasing = TRUE)
OTUtable.filt.sums<-as.data.frame(OTUtable.filt.sums)
OTUtable.filt.sums$sample<-rownames(OTUtable.filt.sums)

OTUtable.1000<-ps_1000@otu_table@.Data
OTUtable.1000<-as.data.frame(t(OTUtable.1000))

OTUtable.1000.sums<-colSums(OTUtable.1000)
OTUtable.1000.sums<-sort(OTUtable.1000.sums, decreasing = TRUE)
OTUtable.1000.sums<-as.data.frame(OTUtable.1000.sums)
OTUtable.1000.sums$sample<-rownames(OTUtable.1000.sums)

sample_table<-data.frame(samdf$sample,samdf$sample_type)
colnames(sample_table)<-c("sample","sample_type")

OTUtable.all<-left_join(OTUtable.filt.sums,sample_table,by="sample")
OTUtable.1000<-left_join(OTUtable.1000.sums,sample_table,by="sample")
OTUtable.removed<-anti_join(OTUtable.all, OTUtable.1000, by = "sample")

OTUtable.1000$result<-"kept"
colnames(OTUtable.1000)[1] <-"value"
OTUtable.1000$result2<-paste(OTUtable.1000$sample_type, OTUtable.1000$result, sep="_")
OTUtable.removed$result<-"rem"
colnames(OTUtable.removed)[1] <-"value"
OTUtable.removed$result2<-paste(OTUtable.removed$sample_type, OTUtable.removed$result, sep="_")

OTUtable.1000<-melt(OTUtable.1000)
OTUtable.removed<-melt(OTUtable.removed)

OTUtable.final<-bind_rows(OTUtable.1000, OTUtable.removed)
OTUtable.final$result[OTUtable.final$sample_type=="community"]<-"rem"

level_order <- c("vag","fol","spe","pen","ctrl","fol_ctrl","tech_ctrl","water")

graph.filt <- ggplot(OTUtable.final, aes(x=factor(sample_type, level = level_order), y=value)) + #, group=dose)) + 
  geom_boxplot(aes(fill=sample_type)) + scale_y_log10() +
  geom_hline(yintercept=3000, linetype="dashed", color = "red") +
  geom_point() +
  scale_x_discrete("sample_type", breaks=factor(1:8), drop=FALSE) +
  facet_grid(result ~ .,scales="free_y")

level_order <- c("vag_kept","vag_rem","fol_kept","fol_rem",
                 "spe_kept","spe_rem","pen_kept","pen_rem", 
                 "ctrl_kept","ctrl_rem","fol_ctrl_kept","fol_ctrl_rem",
                 "tech_ctrl_kept","tech_ctrl_rem","water_rem")

#kept vs removed for all the samples
ggplot(OTUtable.final, aes(x=factor(result2, level = level_order), y=value)) + 
  geom_boxplot(aes(fill=sample_type)) + scale_y_log10() +
  geom_hline(yintercept=500, linetype="dashed", color = "red") +
  geom_point() +
  xlab("Samples") + ylab("Reads") +
  ggpubr::rotate_x_text()
```

## 6.2. - Keep only samples in ps object

```
ps_1000<-subset_samples(ps_1000,is_sample=="yes")
```

## 6.3. - Transform data to relative abundance

```
phylo.ab <- transform_sample_counts(ps_1000, function(x){x/sum(x)})

phylo.abs <- phylo.ab

        for(n in 1:nsamples(phylo.ab))
        {otu_table(phylo.abs)[,n] <- otu_table(phylo.ab)[,n]*sample_data(phylo.ab)$qpcr[n]}

ps2<-phylo.abs
```

# 7 - Most prevalent ASVs (Figure 7)

Prevalence and abundance of top20 ASVs with maximum-likelihood phylogenetic tree

```
library(tidyverse)
library(ggplot2)
library(ggtree)
library(treeio)
library(phyloseq)
library(gridExtra)
library(grid)
library(ggpubr)
library(dplyr)
library(cowplot)
library(gplots)
library(RColorBrewer)
library(reshape2)
library(fantaxtic)

#extract legend for combined plots
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

#small legend function
addSmallLegend <- function(myPlot, pointSize = 1, textSize = 3, spaceLegend = 0.1) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
  }

ps_1000<-subset_samples(ps_1000,is_sample=="yes")

ps_1000_f<-subset_samples(ps_1000,sample_type!="pen")
ps_1000_f<-subset_samples(ps_1000_f,sample_type!="spe")

ps_1000_m<-subset_samples(ps_1000,sample_type!="vag")
ps_1000_m<-subset_samples(ps_1000_m,sample_type!="fol")

ps_20<-get_top_taxa(ps_1000,20)

#SWITCH
phylo.obj<-ps_1000_f
phylo.obj<-get_top_taxa(phylo.obj,19)

#Extraxt OTU table from all the samples
OTU1 = as(otu_table(phylo.obj), "matrix")
if(taxa_are_rows(phylo.obj)){OTU1 <- t(OTU1)}
OTUdf = as.data.frame(OTU1)

#Prevalence
prev<-colSums(OTUdf!=0)
prev<-(100/nrow(OTUdf)*prev)

#Abundance
abund<-colSums(OTUdf)
c<-sum(abund)
abund<-(abund/c*100)

#Create table
table<-data.frame(prev,abund)

#Select ASV that have prevalence above 25%
table_top<-subset(table, prev > 1)

table_top$ASV<-row.names(table_top)

#Select top taxa
ps.top = prune_taxa(table_top$ASV, ps_1000)

#subset phyloseq object:
ps.top.noneg<-subset_samples(ps.top, is_sample=="yes")

ps.top.noneg_vag<-subset_samples(ps.top.noneg, sample_type=="vag")
ps.top.noneg_fol<-subset_samples(ps.top.noneg, sample_type=="fol")
ps.top.noneg_spe<-subset_samples(ps.top.noneg, sample_type=="spe")
ps.top.noneg_pen<-subset_samples(ps.top.noneg, sample_type=="pen")

#Create OTU tables
##vag
OTU.vag = as(otu_table(ps.top.noneg_vag), "matrix")
if(taxa_are_rows(ps.top.noneg_vag)){OTU.vag <- t(OTU.vag)}
OTUdf.vag = as.data.frame(t(OTU.vag))

##fol
OTU.fol = as(otu_table(ps.top.noneg_fol), "matrix")
if(taxa_are_rows(ps.top.noneg_fol)){OTU.fol <- t(OTU.fol)}
OTUdf.fol = as.data.frame(t(OTU.fol))

##spe
OTU.spe = as(otu_table(ps.top.noneg_spe), "matrix")
if(taxa_are_rows(ps.top.noneg_spe)){OTU.spe <- t(OTU.spe)}
OTUdf.spe = as.data.frame(t(OTU.spe))

##pen
OTU.pen = as(otu_table(ps.top.noneg_pen), "matrix")
if(taxa_are_rows(ps.top.noneg_pen)){OTU.pen <- t(OTU.pen)}
OTUdf.pen = as.data.frame(t(OTU.pen))

#selection of most prevalent ASV in home and hospital data
selection<-rownames(table_top)

write(selection, "~/luzern2021/06_Taxonomy/prevalentASV.txt")

#calculate prevalences and abundance
#Prevalence
prev.vag<-rowSums(OTUdf.vag!=0)
prev.vag<-100/nsamples(ps.top.noneg_vag)*prev.vag

prev.fol<-rowSums(OTUdf.fol!=0)
prev.fol<-100/nsamples(ps.top.noneg_fol)*prev.fol

prev.spe<-rowSums(OTUdf.spe!=0)
prev.spe<-100/nsamples(ps.top.noneg_spe)*prev.spe

prev.pen<-rowSums(OTUdf.pen!=0)
prev.pen<-100/nsamples(ps.top.noneg_pen)*prev.pen

#Abundance
##calculate number of reads for each subgroup
totabund.vag <- sum(colSums(OTUdf[grep("vag", rownames(OTUdf)), ]))
totabund.fol <- sum(colSums(OTUdf[grep("fol", rownames(OTUdf)), ]))
totabund.spe <- sum(colSums(OTUdf[grep("spe", rownames(OTUdf)), ]))
totabund.pen <- sum(colSums(OTUdf[grep("pen", rownames(OTUdf)), ]))

##abundance
abund.vag<-rowSums(OTUdf.vag)/totabund.vag*100
abund.fol<-rowSums(OTUdf.fol)/totabund.fol*100
abund.spe<-rowSums(OTUdf.spe)/totabund.spe*100
abund.pen<-rowSums(OTUdf.pen)/totabund.pen*100

#create a dataframe with prevalences and abundances
table.prev<-data.frame(prev.vag,prev.fol,prev.spe,prev.pen)
table.prev$ASV<-rownames(table.prev)

table.abund<-data.frame(abund.vag,abund.fol,abund.spe,abund.pen)
table.abund$ASV<-rownames(table.abund)

#retreive phylum
taxa.table_top<- tax_table(ps)
taxa.table_top<-as.matrix(taxa.table_top@.Data)
taxa.table_top<-as.data.frame(taxa.table_top[,2])
colnames(taxa.table_top)<-"Phylum"
taxa.table_top$ASV<-rownames(taxa.table_top)

#add Phylum column
metadata<-dplyr::inner_join(table_top, taxa.table_top, by = "ASV")

selection_vag<-rownames(as.data.frame(abund.vag))
selection_fol<-rownames(as.data.frame(abund.fol))
selection_spe<-rownames(as.data.frame(abund.spe))
selection_pen<-rownames(as.data.frame(abund.pen))

selection_m<-intersect(selection_spe,selection_pen)
```

```{bash}
cd /Users/mstojano/luzern2021/06_Taxonomy

seqtk subseq ASV_sequences2.fasta prevalentASV.txt > prevalentASV.fasta

conda activate sina

sina -i ~/luzern2021/06_Taxonomy/prevalentASV.fasta -o ~/luzern2021/06_Taxonomy/aligned_prevalentASV.fna --meta-fmt csv --db ~/luzern2021/06_Taxonomy/SILVA_138.1_SSURef_NR99_12_06_20_opt.arb --search --search-db ~/luzern2021/06_Taxonomy/SILVA_138.1_SSURef_NR99_12_06_20_opt.arb --lca-fields tax_slv

#The alignment file in FastTree using aligned sequences

FastTree -nt ~/luzern2021/06_Taxonomy/aligned_prevalentASV.fna > ~/luzern2021/06_Taxonomy/FastTree_prevalent

#Blast selected ASV to find closest matches
cd /Users/mstojano/blast/16SMicrobial_v4/

blastn -db /Users/mstojano/blast/16SMicrobial_v4/16SMicrobial -query ~/luzern2021/06_Taxonomy/prevalentASV.fasta -task blastn -dust no -outfmt "7 delim=, scomname qacc sacc evalue bitscore qcovus pident" -max_target_seqs 2 > ~/luzern2021/06_Taxonomy/blast_preval_result.csv

#clean the blast output file
cd /Users/mstojano/luzern2021/06_Taxonomy/

grep -A 1 "found" blast_preval_result.csv > blast_preval_result_clean.csv

sed -n '/# 2 hits found/!p' blast_preval_result_clean.csv > blast_preval_result_cleann.csv
sed -n '/--/!p' blast_preval_result_cleann.csv > blast_preval_result_cleannn.csv

#bacteria names were manually cleaned
```

```
#import blast results
blast<-read.csv("~/luzern2021/06_Taxonomy/blast_preval_result_cleannn.csv",header = F,sep="\t")
colnames(blast)<-c("hit", "ASV", "accession" ,"evalue" ,"bitscore", "coverage", "identity")

#round the numbers of identity
blast<-blast %>% mutate_at(vars(coverage, identity), funs(round(., 0)))

#import FastTree file
tree <- read.tree("~/luzern2021/06_Taxonomy/FastTree_prevalent")
tree1<- ggtree(tree)
plot_tree<-tree1 + geom_tiplab(align=TRUE, linetype='dashed', linesize=.4)

#extract tree1 data information - ASV order
tree1_data<-tree1$data
tree1_data<-tree1_data[, c(2,4,7)]   
tree1_data<- tree1_data[grep("ASV", tree1_data$label),]
tree1_data<-tree1_data %>% arrange(desc(y))
#tree1_data$y<-NULL
colnames(tree1_data) <- c("node", "ASV", "position")

#Rearrange tables based on ASV positions in the tree
table.prev.ord<-table.prev[tree1_data$ASV,]
table.prev.ord<-table.prev.ord %>% mutate_at(vars(prev.vag, prev.fol, prev.spe, prev.pen), funs(round(., 1)))

table.abund.ord<-table.abund[tree1_data$ASV,]
table.abund.ord<-table.abund.ord %>% mutate_at(vars(abund.vag, abund.fol, abund.spe, abund.pen), funs(round(., 1)))

#Separate sample types
table.prev.ord.vag<-table.prev.ord[,c(1,5)]
table.prev.ord.fol<-table.prev.ord[,c(2,5)]
table.prev.ord.spe<-table.prev.ord[,c(3,5)]
table.prev.ord.pen<-table.prev.ord[,c(4,5)]

table.abund.ord.vag<-table.abund.ord[,c(1,5)]
table.abund.ord.fol<-table.abund.ord[,c(2,5)]
table.abund.ord.spe<-table.abund.ord[,c(3,5)]
table.abund.ord.pen<-table.abund.ord[,c(4,5)]

#Rearrange blast table based on ASV positions in the tree
rownames(blast)<-blast$ASV
blast.ord<-blast[tree1_data$ASV,]

tree1_data$species<-paste(blast.ord$hit,"_",blast.ord$ASV, sep = "")

for(i in 1:nrow(tree1_data)){
  
  for(j in 1:nrow(tree1_data)){
    if(tree$tip.label[i]==tree1_data$ASV[j]){
      tree$tip.label[i]=tree1_data$species[j]
    }
  }
}  

tree1<- ggtree(tree)
plot_tree<-tree1 + geom_tiplab(align=TRUE, linetype='dashed', linesize=.3,size=3) + xlim(NA, 1)

pdf("tree_most_preval.pdf",width = 8,height = 8)
plot_tree
dev.off()

#print only ASVs
mytheme <- gridExtra::ttheme_minimal(col.just="left",base_size = 10,
                                     core = list(fg_params=list(hjust=0, x=0.1)),
                                     rowhead=list(fg_params=list(hjust=0, x=0)),
                                     padding=unit(c(10, 5), "mm"))


blast.ord$OTU<-gsub("ASV", "OTU", blast.ord$ASV)
pASVs<-tableGrob(blast.ord$ASV,theme=mytheme)
grid.draw(pASVs)

#print only hits
phits<-tableGrob(blast.ord$hit,theme=mytheme)
grid.draw(phits)

pdf("tree_most_preval.pdf",width = 8,height = 8)
plot_tree
dev.off()
```

```
heatmap.prev<-table.prev.ord %>% 
  pivot_longer(!ASV, names_to = "type", values_to = "value")

heatmap.abund<-table.abund.ord %>% 
  pivot_longer(!ASV, names_to = "type", values_to = "value")

prev.vag<- melt(table.prev.ord.vag)
prev.fol<- melt(table.prev.ord.fol)
prev.spe<- melt(table.prev.ord.spe)
prev.pen<- melt(table.prev.ord.pen)

abund.vag<- melt(table.abund.ord.vag)
abund.fol<- melt(table.abund.ord.fol)
abund.spe<- melt(table.abund.ord.spe)
abund.pen<- melt(table.abund.ord.pen)

#PLOTS PREVALENCE
#Prepare the legend
px<-ggplot(prev.vag, aes(x = variable, y = factor(ASV), fill = value)) +
  geom_raster(aes(fill=value))+
  scale_fill_distiller(palette = "OrRd", direction=2,breaks=c(0, 25,50,75,100),limits = c(0,100)) +
  geom_text(aes(label = value),size=3) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_y_discrete(drop = F,limits=rev(heatmap.prev$ASV))

mylegend<-g_legend(px + theme(legend.position="bottom") + theme(legend.title=element_blank()))

#Plots
p1<-ggplot(prev.vag, aes(x = variable, y = factor(ASV), fill = value)) +
  geom_raster(aes(fill=value))+
  scale_fill_distiller(palette = "OrRd", direction=2,breaks=c(0, 25,50,75,100),limits = c(0,100)) +
  geom_text(aes(label = value),size=3) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_y_discrete(drop = F,limits=rev(heatmap.prev$ASV))+
  theme_nothing()+
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

p2<-ggplot(prev.fol, aes(x = variable, y = factor(ASV), fill = value)) +
  geom_raster(aes(fill=value))+
  scale_fill_distiller(palette = "OrRd", direction=2,breaks=c(0, 25,50,75,100),limits = c(0,100)) +
  geom_text(aes(label = value),size=3) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_y_discrete(drop = F,limits=rev(heatmap.prev$ASV))+
  theme_nothing()+
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

p3<-ggplot(prev.spe, aes(x = variable, y = factor(ASV), fill = value)) +
  geom_raster(aes(fill=value))+
  scale_fill_distiller(palette = "OrRd", direction=2,breaks=c(0, 25,50,75,100),limits = c(0,100)) +
  geom_text(aes(label = value),size=3) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_y_discrete(drop = F,limits=rev(heatmap.prev$ASV))+
  theme_nothing()+
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

p4<-ggplot(prev.pen, aes(x = variable, y = factor(ASV), fill = value)) +
  geom_raster(aes(fill=value))+
  scale_fill_distiller(palette = "OrRd", direction=2,breaks=c(0, 25,50,75,100),limits = c(0,100)) +
  geom_text(aes(label = value),size=3) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_y_discrete(drop = F,limits=rev(heatmap.prev$ASV))+
  theme_nothing()+
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

#Combine the plots
p_preval <- grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                                    p2 + theme(legend.position="none"),
                                    p3 + theme(legend.position="none"),
                                    p4 + theme(legend.position="none"),
                                    nrow=1),
                                    mylegend, nrow=2,heights=c(20, 10))

#PLOTS ABUNDANCE TOT
#Prepare the legend
py<-ggplot(abund.vag, aes(x = variable, y = factor(ASV), fill = value)) +
  geom_raster(aes(fill=value))+
  scale_fill_distiller(palette = "Blues", direction=2,breaks=c(0, 10,20,30,40),limits = c(0,40)) +
  geom_text(aes(label = round(value, digits=1)),size=3,show.legend = F) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_y_discrete(drop = F,limits=rev(heatmap.abund$ASV))

#run p8 without the last two paragraphs to memorize the legend
addSmallLegend(py)

mylegend<-g_legend(py + theme(legend.position="bottom") + theme(legend.title=element_blank())+ guides(color = guide_legend(override.aes = list(size = 0.01))))

#Plots
p5<-ggplot(abund.vag, aes(x = variable, y = factor(ASV), fill = value)) +
  geom_raster(aes(fill=value))+
  scale_fill_distiller(palette = "Blues", direction=2,breaks=c(0, 10,20,30,40),limits = c(0,40))+
  geom_text(aes(label = round(value, digits=1)),size=3) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_y_discrete(drop = F,limits=rev(heatmap.abund$ASV))+
  theme_nothing()+
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

p6<-ggplot(abund.fol, aes(x = variable, y = factor(ASV), fill = value)) +
  geom_raster(aes(fill=value))+
  scale_fill_distiller(palette = "Blues", direction=2,breaks=c(0, 10,20,30,40),limits = c(0,40)) +
  geom_text(aes(label = round(value, digits=1)),size=3) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_y_discrete(drop = F,limits=rev(heatmap.abund$ASV))+
  theme_nothing()+
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

p7<-ggplot(abund.spe, aes(x = variable, y = factor(ASV), fill = value)) +
  geom_raster(aes(fill=value))+
  scale_fill_distiller(palette = "Blues", direction=2,breaks=c(0, 10,20,30,40),limits = c(0,40))+
  geom_text(aes(label = round(value, digits=1)),size=3) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_y_discrete(drop = F,limits=rev(heatmap.abund$ASV))+
  theme_nothing()+
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

p8<-ggplot(abund.pen, aes(x = variable, y = factor(ASV), fill = value)) +
  geom_raster(aes(fill=value))+
  scale_fill_distiller(palette = "Blues", direction=2,breaks=c(0, 10,20,30,40),limits = c(0,40)) +
  geom_text(aes(label = round(value, digits=1)),size=3,show.legend = F) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_y_discrete(drop = F,limits=rev(heatmap.abund$ASV))+
  theme_nothing()+
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

#Combine the plots
p_abund <- grid.arrange(arrangeGrob(p5 + theme(legend.position="none"),
                                     p6 + theme(legend.position="none"),
                                     p7 + theme(legend.position="none"),
                                     p8 + theme(legend.position="none"),
                                     nrow=1),
                                    mylegend, nrow=2,heights=c(20, 10))

# Combine the plots
g_prevalabund = cbind(p_preval,p_abund, size = "last")

# Draw it
grid.newpage()
grid.draw(g_prevalabund)

#Save figure
pdf("figure_tree.pdf",width = 5,height = 14)
grid.newpage()
grid.draw(g_prevalabund)
dev.off()
```

# 8 - Shared ASVs between sample types. (Figure 8)

## Venn diagram of shared ASVs across the four sample types

```
library(tidyverse)
library(caret)
library(leaps)
library(dplyr)
library(fantaxtic)
library(ggVennDiagram)

ps_spe<-phyloseq::subset_samples(ps_1000, sample_type=="spe")
ps_pen<-phyloseq::subset_samples(ps_1000, sample_type=="pen")
ps_fol<-phyloseq::subset_samples(ps_1000, sample_type=="fol")
ps_vag<-phyloseq::subset_samples(ps_1000, sample_type=="vag")

#Extract ASVs present in all sample types
#spe samples
ps_spe_100<-get_top_taxa(ps_spe, 99, relative = TRUE, discard_other = FALSE,other_label = "Other")
spe_asv100<-colSums(ps_spe_100@otu_table@.Data)
spe_asv100<-spe_asv100[spe_as100v>0]
spe_asv100<-names(spe_asv100)

spe_asv<-colSums(ps_spe_100@otu_table@.Data)
spe_asv<-spe_asv[spe_asv>0]
spe_asv<-names(spe_asv)

#pen samples
ps_pen_100<-get_top_taxa(ps_pen, 99, relative = TRUE, discard_other = FALSE,other_label = "Other")
pen_asv100<-colSums(ps_pen_100@otu_table@.Data)
pen_asv100<-pen_asv100[pen_asv100>0]
pen_asv100<-names(pen_asv100)

pen_asv<-colSums(ps_pen@otu_table@.Data)
pen_asv<-pen_asv[pen_asv>0]
pen_asv<-names(pen_asv)

#vag samples
ps_vag_100<-get_top_taxa(ps_vag, 99, relative = TRUE, discard_other = FALSE,other_label = "Other")
vag_asv100<-colSums(ps_vag_100@otu_table@.Data)
vag_asv100<-vag_asv100[vag_asv100>0]
vag_asv100<-names(vag_asv100)

vag_asv<-colSums(ps_vag@otu_table@.Data)
vag_asv<-vag_asv[vag_asv>0]
vag_asv<-names(vag_asv)

#fol samples
ps_fol_100<-get_top_taxa(ps_fol, 99, relative = TRUE, discard_other = FALSE,other_label = "Other")
fol_asv100<-colSums(ps_fol_100@otu_table@.Data)
fol_asv100<-fol_asv100[fol_asv100>0]
fol_asv100<-names(fol_asv100)

fol_asv<-colSums(ps_fol@otu_table@.Data)
fol_asv<-fol_asv[fol_asv>0]
fol_asv<-names(fol_asv)

all_asv_100<-list("vag_asv"=vag_asv100,"fol_asv"=fol_asv100,"pen_asv"=pen_asv100,"spe_asv"=spe_asv100)
all_asv<-list("vag_asv"=vag_asv,"fol_asv"=fol_asv,"pen_asv"=pen_asv,"spe_asv"=spe_asv)

#find ASVs present in all samples (top100)
intersect_f<-intersect(vag_asv,fol_asv)
intersect_m<-intersect(pen_asv,spe_asv)
intersect_all<-as.data.frame(intersect(intersect_f,intersect_m))
colnames(intersect_all)<-"ASV"

intersect_all<-left_join(intersect_all,ASVtable)

intersect_all_melt<-melt(intersect_all)

bp<- ggplot(intersect_all_melt, aes(x="", y=Genus, fill=Genus))+
  geom_bar(width = 1, stat = "identity") 

pie <- bp + coord_polar("y", start=0)
pie

intersect_100_asv<-left_join(intersect_all,asv_table)

venn_all<-ggVennDiagram(all_asv) +
  ggplot2::scale_fill_gradient(low="cornsilk",high = "orange")

pdf("~/luzern2021/07_output/venn_all_asv.pdf",width = 10,height = 10)
ggVennDiagram(all_asv_100) +
  ggplot2::scale_fill_gradient(low="cornsilk",high = "orange")
dev.off()
```

# 9 - Alpha diversity (Figures 9 & 10)

```
ps_1000<-subset_samples(ps_1000,is_sample=="yes")

ps_m<-subset_samples(ps_1000, sample_type!="fol")
ps_m<-subset_samples(ps_m, sample_type!="vag")

ps_f<-subset_samples(ps_1000, sample_type!="pen")
ps_f<-subset_samples(ps_f, sample_type!="spe")

ps_vag<-subset_samples(ps_1000,sample_type=="vag")
ps_fol<-subset_samples(ps_1000,sample_type=="fol")
ps_spe<-subset_samples(ps_1000,sample_type=="spe")
ps_pen<-subset_samples(ps_1000,sample_type=="pen")

ps_vag_pen<-subset_samples(ps_1000,sample_type!="fol")
ps_vag_pen<-subset_samples(ps_vag_pen,sample_type!="spe")

ps_fol_spe<-subset_samples(ps_1000,sample_type!="vag")
ps_fol_spe<-subset_samples(ps_fol_spe,sample_type!="pen")

#plot alpha diversity
alpha<- plot_richness(ps_1000,title="Alpha diversity",x="sample_type",measures=c("Observed","Shannon","InvSimpson")) + 
  geom_boxplot()

#find correlations between alpha diversities
alpha.<-data.frame(sample=alpha$data$sample,sample_type=alpha$data$sample_type,value=alpha$data$value,couple=alpha$data$couple)

alpha.vag <- subset(alpha., sample_type=="vag")
alpha.fol <- subset(alpha., sample_type=="fol")
alpha.pen <- subset(alpha., sample_type=="pen")
alpha.spe <- subset(alpha., sample_type=="spe")

alpha.pen.vag<-left_join(alpha.pen,alpha.vag,by="couple")
alpha.vag.fol<-left_join(alpha.fol,alpha.vag,by="couple")
alpha.spe.pen<-left_join(alpha.spe,alpha.pen,by="couple")
alpha.spe.vag<-left_join(alpha.spe,alpha.vag,by="couple")
alpha.spe.fol<-left_join(alpha.spe,alpha.fol,by="couple")
alpha.pen.fol<-left_join(alpha.pen,alpha.fol,by="couple")

ggscatter(alpha.pen.fol, x = "value.y", y = "value.x",color = "black", shape = 21, size = 2, # Points color, shape and size
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE, # Add confidence interval
   cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
   cor.coeff.args = list(method = "pearson", label.x = 0, label.sep = "\n"),
   font.label = c(8, "plain"),
   xlim = c(0, 3.5),
   xlab="fol",
   ylab="pen"
   #label = "couple"
   )

#function to display some useful information on boxplots
stat_box_data <- function(y) {
  return( 
    data.frame(
      y = 0.95 * upper_limit,
      label = paste('count =', length(y), '\n',
                    'mean =', round(mean(y), 1), '\n')
    )
  )
}

#infertility cause
plot_richness(ps_1000,title="Alpha diversity",x="infretility",measures=c("Chao1", "Shannon","Observed")) + geom_boxplot()

#sperm vs normosperm
plot_richness(ps_spe,title="Alpha diversity",x="normosperm",measures=c("Chao1", "Shannon","Observed")) + geom_boxplot()

#pen vs normosperm
plot_richness(ps_pen,title="Alpha diversity",x="normosperm",measures=c("Chao1", "Shannon","Observed")) + geom_boxplot()

#sperm vs oligosperm
plot_richness(ps_spe,title="Alpha diversity",x="oligosperm",measures=c("Chao1", "Shannon","Observed")) + geom_boxplot()

#sperm vs asthenozoosperm
plot_richness(ps_spe,title="Alpha diversity",x="asthenozoosperm",measures=c("Chao1", "Shannon","Observed")) + geom_boxplot()

#sperm vs teratozoosperm
plot_richness(ps_spe,title="Alpha diversity",x="teratozoosperm",measures=c("Chao1", "Shannon","Observed")) + geom_boxplot()
```

## 9.1. - Perform ANOVA

Find significant results for metadata variables

```
variables<-colnames(samdf_1000)

#vagina
alpha.div.anova.vag<-data.frame()

for (z in 1:length(variables)) {
rich_vag<-estimate_richness(ps_vag, measures=c("Observed", "InvSimpson", "Shannon", "Chao1"))
rownames(rich_vag)<-sample_names(ps_vag)
rich_vag[[variables[z]]]<-sample_data(ps_vag)[[variables[z]]] 
try(anova_sample_type <- aov(Shannon ~ rich_vag[[6]], rich_vag))
new.table<-as.data.frame(do.call(rbind, summary(anova_sample_type)))
rownames(new.table)[1]<-variables[z]
new.table<-new.table[1,]
if(nrow(new.table) == 0) next
alpha.div.anova.vag <- bind_rows(alpha.div.anova.vag,new.table)
}

#Plot significant results
plot_richness(ps_vag,title="Shannon index",x="endometriosis",measures="Shannon") + 
  geom_boxplot() + 
  stat_summary(fun.data = stat_box_data,geom = "text", hjust = 0.5, vjust = 5.0)

    

#fol
alpha.div.anova.fol<-data.frame()
for (z in 1:length(variables)) {
rich_fol<-estimate_richness(ps_fol, measures=c("Observed", "InvSimpson", "Shannon", "Chao1"))
rownames(rich_fol)<-sample_names(ps_fol)
rich_fol[[variables[z]]]<-sample_data(ps_fol)[[variables[z]]] 
try(anova_sample_type <- aov(Shannon ~ rich_fol[[6]], rich_fol))
new.table<-as.data.frame(do.call(rbind, summary(anova_sample_type)))
rownames(new.table)[1]<-variables[z]
new.table<-new.table[1,]
if(nrow(new.table) == 0) next
alpha.div.anova.fol <- bind_rows(alpha.div.anova.fol,new.table)
}

#Plot significant results
plot_richness(ps_fol,title="Shannon index",x="ICSI",measures="Shannon") +   
    geom_boxplot() +
    stat_summary(fun.data = stat_box_data, geom = "text", hjust = 0.5, vjust = 5.0)

#pen
alpha.div.anova.pen<-data.frame()
for (z in 1:length(variables)) {
rich_pen<-estimate_richness(ps_pen, measures=c("Observed", "InvSimpson", "Shannon", "Chao1"))
rownames(rich_pen)<-sample_names(ps_pen)
rich_pen[[variables[z]]]<-sample_data(ps_pen)[[variables[z]]] 
try(anova_sample_type <- aov(Shannon ~ rich_pen[[6]], rich_pen))
new.table<-as.data.frame(do.call(rbind, summary(anova_sample_type)))
rownames(new.table)[1]<-variables[z]
new.table<-new.table[1,]
if(nrow(new.table) == 0) next
alpha.div.anova.pen <- bind_rows(alpha.div.anova.pen,new.table)
}

#Plot significant results - pen
plot_richness(ps_pen,title="Shannon index",x="number_blastocystes",measures="Shannon") + 
  geom_boxplot() +
  stat_summary(fun.data = stat_box_data, geom = "text", hjust = 0.5, vjust = 5.0)

#spe
alpha.div.anova.spe<-data.frame()
for (z in 1:length(variables)) {
rich_spe<-estimate_richness(ps_spe, measures=c("Observed", "InvSimpson", "Shannon", "Chao1"))
rownames(rich_spe)<-sample_names(ps_spe)
rich_spe[[variables[z]]]<-sample_data(ps_spe)[[variables[z]]] 
try(anova_sample_type <- aov(Shannon ~ rich_spe[[6]], rich_spe))
new.table<-as.data.frame(do.call(rbind, summary(anova_sample_type)))
rownames(new.table)[1]<-variables[z]
new.table<-new.table[1,]
if(nrow(new.table) == 0) next
alpha.div.anova.spe <- bind_rows(alpha.div.anova.spe,new.table)
}

#Plot significant results - spe
plot_richness(ps_spe,x="azoosperm",measures="Shannon") + 
  geom_boxplot() +
  stat_summary(fun.data = stat_box_data, geom = "text", hjust = 0.5, vjust = 5.0)

plot_richness(ps_spe,title="myometritis",x="myometritis",measures="Observed") + 
  geom_boxplot() +
  stat_summary(fun.data = stat_box_data, geom = "text", hjust = 0.5, vjust = 200.0) +
  annotate("text", x = 1.5, y=300, label = "p = 0.02")

plot_richness(ps_spe,title="age",x="age",measures="Shannon") + 
  geom_boxplot() +
  stat_summary(fun.data = stat_box_data, geom = "text", hjust = 0.5, vjust = 5.0) +
  annotate("text", x = 1.5, y=7.5, label = "p = 0.0")

plot_richness(ps_spe,title="azoosperm",x="azoosperm",measures="Shannon") + 
  geom_boxplot() +
  stat_summary(fun.data = stat_box_data, geom = "text", hjust = 0.5, vjust = 5.0) +
  annotate("text", x = 1.5, y=7.5, label = "p = 0.0")

##Male samples
alpha.div.anova.m<-data.frame()
for (z in 1:length(variables)) {
rich_m<-estimate_richness(ps_m, measures=c("Observed", "InvSimpson", "Shannon", "Chao1"))
rownames(rich_m)<-sample_names(ps_m)
rich_m[[variables[z]]]<-sample_data(ps_m)[[variables[z]]] 
try(anova_sample_type <- aov(Shannon ~ rich_m[[6]], rich_m))
new.table<-as.data.frame(do.call(rbind, summary(anova_sample_type)))
rownames(new.table)[1]<-variables[z]
new.table<-new.table[1,]
if(nrow(new.table) == 0) next
alpha.div.anova.m <- bind_rows(alpha.div.anova.m,new.table)
}

#Plot significant results

plot_richness(ps_m,title="Shannon index",x="oligosperm",measures="Shannon") + 
  geom_boxplot() +
  stat_summary(fun.data = stat_box_data, geom = "text", hjust = 0.5, vjust = 5.0)

##Female samples
alpha.div.anova.f<-data.frame()
for (z in 1:length(variables)) {
rich_f<-estimate_richness(ps_f, measures=c("Observed", "InvSimpson", "Shannon", "Chao1"))
rownames(rich_f)<-sample_names(ps_f)
rich_f[[variables[z]]]<-sample_data(ps_f)[[variables[z]]] 
try(anova_sample_type <- aov(Shannon ~ rich_f[[6]], rich_f))
new.table<-as.data.frame(do.call(rbind, summary(anova_sample_type)))
rownames(new.table)[1]<-variables[z]
new.table<-new.table[1,]
if(nrow(new.table) == 0) next
alpha.div.anova.m <- bind_rows(alpha.div.anova.m,new.table)
}

#Plot significant results
plot_richness(ps_m,title="Shannon index",x="oligosperm",measures="Shannon") + 
  geom_boxplot() +
  stat_summary(fun.data = stat_box_data, geom = "text", hjust = 0.5, vjust = 5.0)

##All samples
alpha.div.anova<-data.frame()
for (z in 1:length(variables)) {
rich_1000<-estimate_richness(ps_1000, measures=c("Observed", "InvSimpson", "Shannon", "Chao1"))
rownames(rich_1000)<-sample_names(ps_1000)
rich_1000[[variables[z]]]<-sample_data(ps_1000)[[variables[z]]] 
try(anova_sample_type <- aov(Shannon ~ rich_1000[[6]], rich_1000))
new.table<-as.data.frame(do.call(rbind, summary(anova_sample_type)))
rownames(new.table)[1]<-variables[z]
new.table<-new.table[1,]
if(nrow(new.table) == 0) next
alpha.div.anova <- bind_rows(alpha.div.anova,new.table)
}

#Plot significant results
plot_richness(ps_1000,title="Shannon index",x="previous_life_birth",measures="Shannon") + 
  geom_boxplot() +
  stat_summary(fun.data = stat_box_data, geom = "text", hjust = 0.5, vjust = 5.0)

#Perform Wicoxon test
pairwise.wilcox.test(rich_spe$Observed,rich_spe$sample_type)
```

# 10 - Beta diversity (Figure 11)

```
library("phyloseq")
library("ggplot2")
library("vegan")
library("data.table"); packageVersion("data.table")

ps_noneg<-ps_1000

dist.bin<-phyloseq::distance(ps_noneg,"bray",binary=T)
hn.dist<-phyloseq::distance(ps_noneg,"horn")
uni.dist<-phyloseq::distance(ps_noneg,"unifrac")
wuni.dist<-phyloseq::distance(ps_noneg,"wunifrac")

ordination = ordinate(ps_noneg, method="PCoA", distance=uni.dist)

#ordination = ordinate(ps_noneg, method="NMDS", distance=uni.dist)
#stressplot(ordination)

plot_ordination(ps_noneg, ordination, color="sample_type") + 
  theme(aspect.ratio=1) + stat_ellipse() +
  theme_classic()
```

# 11 - Pairwise distance intra- and inter-group (Figure 12)

```
#subset phyloseq object
ps_vag_fol<-subset_samples(ps_1000, sample_type!="pen")
ps_vag_fol<-subset_samples(ps_vag_fol, sample_type!="spe")

ps_vag_spe<-subset_samples(ps_1000, sample_type!="pen")
ps_vag_spe<-subset_samples(ps_vag_spe, sample_type!="fol")

ps_vag_pen<-subset_samples(ps_1000, sample_type!="spe")
ps_vag_pen<-subset_samples(ps_vag_pen, sample_type!="fol")

ps_spe_pen<-subset_samples(ps_1000, sample_type!="vag")
ps_spe_pen<-subset_samples(ps_spe_pen, sample_type!="fol")

ps_spe_fol<-subset_samples(ps_1000, sample_type!="vag")
ps_spe_fol<-subset_samples(ps_spe_fol, sample_type!="pen")

ps_pen_fol<-subset_samples(ps_1000, sample_type!="spe")
ps_pen_fol<-subset_samples(ps_pen_fol, sample_type!="vag")

###################

#VAG_FOL
#get couples number
couples<-as.vector(c(unique(metadata$couple)))
unique<-unique(row.names(metadata))

sample.names<-sample_names(ps_vag_fol)

unique_vag <- grep("vag", sample.names)
unique_vag <-sample.names[unique_vag]
unique_vag<-str_replace(unique_vag, "_vag", "")

unique_fol <- grep("fol", sample.names)
unique_fol <-sample.names[unique_fol]
unique_fol<-str_replace(unique_fol, "_fol", "")

couples_vag_fol<-as.integer(intersect(unique_vag,unique_fol))
names(couples_vag_fol)<-couples_vag_fol

if("3112" %in%couples_vag_fol){
  couples_vag_fol["3112"]<-312
}

#create an empty t vector
t<-vector()

#loop: 
for(i in 1:length(couples_vag_fol)){
  ps.subs <- subset_samples(ps_vag_fol, couple %in% couples_vag_fol[i])
  metadata_sub <- data.frame(sample_data(ps.subs))
  
  distance_pairwise<-phyloseq::distance(ps.subs, method="bray")
  
  
  t<-c(t,distance_pairwise)
}

t<-as.vector(t)
table_vag_fol<-data.frame(couples_vag_fol,t)
colnames(table_vag_fol)<-c("couple","vag_fol")

############

#VAG_SPE
#get couples number
couples<-as.vector(c(unique(metadata$couple)))
unique<-unique(row.names(metadata))

sample.names<-sample_names(ps_vag_spe)

unique_vag <- grep("vag", sample.names)
unique_vag <-sample.names[unique_vag]
unique_vag<-str_replace(unique_vag, "_vag", "")

unique_spe <- grep("spe", sample.names)
unique_spe <-sample.names[unique_spe]
unique_spe<-str_replace(unique_spe, "_spe", "")

couples_vag_spe<-as.integer(intersect(unique_vag,unique_spe))
names(couples_vag_spe)<-couples_vag_spe

if("3112" %in%couples_vag_spe){
  couples_vag_spe["3112"]<-312
}

#create an empty t vector
t<-vector()

#loop: 
for(i in 1:length(couples_vag_spe)){
  ps.subs <- subset_samples(ps_vag_spe, couple %in% couples_vag_spe[i])
  metadata_sub <- data.frame(sample_data(ps.subs))
  
  distance_pairwise<-phyloseq::distance(ps.subs, method="bray")
  
  
  t<-c(t,distance_pairwise)
}

t<-as.vector(t)
table_vag_spe<-data.frame(couples_vag_spe,t)
colnames(table_vag_spe)<-c("couple","vag_spe")

############

#VAG_PEN
#get couples number
couples<-as.vector(c(unique(metadata$couple)))
unique<-unique(row.names(metadata))

sample.names<-sample_names(ps_vag_pen)

unique_vag <- grep("vag", sample.names)
unique_vag <-sample.names[unique_vag]
unique_vag<-str_replace(unique_vag, "_vag", "")

unique_pen <- grep("pen", sample.names)
unique_pen <-sample.names[unique_pen]
unique_pen<-str_replace(unique_pen, "_pen", "")

couples_vag_pen<-as.integer(intersect(unique_vag,unique_pen))
names(couples_vag_pen)<-couples_vag_pen

if("3112" %in%couples_vag_pen){
  couples_vag_pen["3112"]<-312
}

#create an empty t vector
t<-vector()

#loop: 
for(i in 1:length(couples_vag_pen)){
  ps.subs <- subset_samples(ps_vag_pen, couple %in% couples_vag_pen[i])
  metadata_sub <- data.frame(sample_data(ps.subs))
  
  distance_pairwise<-phyloseq::distance(ps.subs, method="bray")
  
  
  t<-c(t,distance_pairwise)
}

t<-as.vector(t)
table_vag_pen<-data.frame(couples_vag_pen,t)
colnames(table_vag_pen)<-c("couple","vag_pen")

############

#SPE_FOL
#get couples number
couples<-as.vector(c(unique(metadata$couple)))
unique<-unique(row.names(metadata))

sample.names<-sample_names(ps_spe_fol)

unique_spe <- grep("spe", sample.names)
unique_spe <-sample.names[unique_spe]
unique_spe<-str_replace(unique_spe, "_spe", "")

unique_fol <- grep("fol", sample.names)
unique_fol <-sample.names[unique_fol]
unique_fol<-str_replace(unique_fol, "_fol", "")

couples_spe_fol<-as.integer(intersect(unique_spe,unique_fol))
names(couples_spe_fol)<-couples_spe_fol

if("3112" %in%couples_spe_fol){
  couples_spe_fol["3112"]<-312
}

#create an empty t vector
t<-vector()

#loop: 
for(i in 1:length(couples_spe_fol)){
  ps.subs <- subset_samples(ps_spe_fol, couple %in% couples_spe_fol[i])
  metadata_sub <- data.frame(sample_data(ps.subs))

  distance_pairwise<-phyloseq::distance(ps.subs, method="bray")
  
  
  t<-c(t,distance_pairwise)
}

t<-as.vector(t)
table_spe_fol<-data.frame(couples_spe_fol,t)
colnames(table_spe_fol)<-c("couple","spe_fol")

############

#SPE_PEN
#get couples number
sample.names<-sample_names(ps_spe_pen)

unique_spe <- grep("spe", sample.names)
unique_spe <-sample.names[unique_spe]
unique_spe<-str_replace(unique_spe, "_spe", "")

unique_pen <- grep("pen", sample.names)
unique_pen <-sample.names[unique_pen]
unique_pen<-str_replace(unique_pen, "_pen", "")

couples_spe_pen<-as.integer(intersect(unique_spe,unique_pen))
names(couples_spe_pen)<-couples_spe_pen

if("3112" %in%couples_spe_pen){
couples_spe_pen["3112"]<-312
}

#create an empty t vector
t<-vector()
count=0

#loop: 
for(i in 1:length(couples_spe_pen)){
  ps.subs <- subset_samples(ps_spe_pen, couple %in% couples_spe_pen[i])
  metadata_sub <- data.frame(sample_data(ps.subs))
  
  distance_pairwise<-phyloseq::distance(ps.subs, method="bray")
  
  t<-c(t,distance_pairwise)
  
  count=count+1
  
  print(count)
}

t<-as.vector(t)
table_spe_pen<-data.frame(couples_spe_pen,t)
colnames(table_spe_pen)<-c("couple","spe_pen")

#PEN_FOL
#get couples number
couples<-as.vector(c(unique(metadata$couple)))
unique<-unique(row.names(metadata))

sample.names<-sample_names(ps_pen_fol)

unique_pen <- grep("pen", sample.names)
unique_pen <-sample.names[unique_pen]
unique_pen<-str_replace(unique_pen, "_pen", "")

unique_fol <- grep("fol", sample.names)
unique_fol <-sample.names[unique_fol]
unique_fol<-str_replace(unique_fol, "_fol", "")

couples_pen_fol<-as.integer(intersect(unique_pen,unique_fol))
names(couples_pen_fol)<-couples_pen_fol

if("3112" %in%couples_pen_fol){
  couples_pen_fol["3112"]<-312
}

#create an empty t vector
t<-vector()

#loop: 
for(i in 1:length(couples_pen_fol)){
  ps.subs <- subset_samples(ps_pen_fol, couple %in% couples_pen_fol[i])
  metadata_sub <- data.frame(sample_data(ps.subs))
  
  distance_pairwise<-phyloseq::distance(ps.subs, method="bray")
  
  
  t<-c(t,distance_pairwise)
}

t<-as.vector(t)
table_pen_fol<-data.frame(couples_pen_fol,t)
colnames(table_pen_fol)<-c("couple","pen_fol")

############

##pairwise_distance between samples of the same type

#subset for sample_type
ps_vag<-subset_samples(ps_1000, sample_type=="vag")
ps_fol<-subset_samples(ps_1000, sample_type=="fol")
ps_spe<-subset_samples(ps_1000, sample_type=="spe")
ps_pen<-subset_samples(ps_1000, sample_type=="pen")

#get couples combinations
#cbn_couples<-combn(x=unique(metadata$couple),m=2)

#SPE-SPE
metadata <- data.frame(sample_data(ps_spe))
cbn_couples<-combn(x=unique(metadata$couple),m=2)
only_spe<-vector()

#loop: 
for(i in 1:ncol(cbn_couples)){
  ps.subs <- subset_samples(ps_spe, couple %in% cbn_couples[,i])
  metadata_sub <- data.frame(sample_data(ps.subs))
  
  distance_pairwise<- phyloseq::distance(ps.subs, method="bray")
  
  only_spe<-c(only_spe,distance_pairwise)
}

table_spe_spe<-data.frame(only_spe)
colnames(table_spe_spe)<-"spe_spe"

#PEN-PEN
metadata <- data.frame(sample_data(ps_pen))
cbn_couples<-combn(x=unique(metadata$couple),m=2)
only_pen<-vector()

#loop: 
for(i in 1:ncol(cbn_couples)){
  ps.subs <- subset_samples(ps_pen, couple %in% cbn_couples[,i])
  metadata_sub <- data.frame(sample_data(ps.subs))
  
  distance_pairwise<- phyloseq::distance(ps.subs, method="bray")
  
  only_pen<-c(only_pen,distance_pairwise)
}

table_pen_pen<-data.frame(only_pen)
colnames(table_pen_pen)<-"pen_pen"

#VAG-VAG
metadata <- data.frame(sample_data(ps_vag))
cbn_couples<-combn(x=unique(metadata$couple),m=2)
only_vag<-vector()

#loop: 
for(i in 1:ncol(cbn_couples)){
  ps.subs <- subset_samples(ps_vag, couple %in% cbn_couples[,i])
  metadata_sub <- data.frame(sample_data(ps.subs))
  
  distance_pairwise<- phyloseq::distance(ps.subs, method="bray")
  
  only_vag<-c(only_vag,distance_pairwise)
}

table_vag_vag<-data.frame(only_vag)
colnames(table_vag_vag)<-"vag_vag"

#create an empty vector
metadata <- data.frame(sample_data(ps_fol))
cbn_couples<-combn(x=unique(metadata$couple),m=2)
only_fol<-vector()

#loop: 
for(i in 1:ncol(cbn_couples)){
  ps.subs <- subset_samples(ps_fol, couple %in% cbn_couples[,i])
  metadata_sub <- data.frame(sample_data(ps.subs))
  
  distance_pairwise<- phyloseq::distance(ps.subs, method="bray")
  
  only_fol<-c(only_fol,distance_pairwise)
}

table_fol_fol<-data.frame(only_fol)
colnames(table_fol_fol)<-"fol_fol"

#intra_groups<-data.frame(only_vag,only_)

#Create table
beta_table<-as.data.frame(rbind(melt(table_fol_fol),melt(table_pen_pen),melt(table_spe_spe),melt(table_vag_vag),
                  melt(table_pen_fol[2]),melt(table_spe_fol[2]),melt(table_spe_pen[2]),melt(table_vag_fol[2]),
                  melt(table_vag_pen[2]),melt(table_vag_spe[2]),melt(table_vag_pen[2])))
  
beta_table_fol<-as.data.frame(rbind(melt(table_fol_fol),melt(table_vag_fol[2]),melt(table_pen_fol[2]),melt(table_spe_fol[2])))

beta_table_vag<-as.data.frame(rbind(melt(table_vag_vag),melt(table_vag_fol[2]),melt(table_vag_pen[2]),melt(table_vag_spe[2])))

beta_table_spe<-as.data.frame(rbind(melt(table_spe_spe),melt(table_spe_pen[2]),melt(table_spe_fol[2]),melt(table_vag_spe[2])))

beta_table_pen<-as.data.frame(rbind(melt(table_pen_pen),melt(table_spe_pen[2]),melt(table_pen_fol[2]),melt(table_vag_pen[2])))

test1<-beta_table$variable
test2<-as.numeric(factor(beta_table$value))
test<-data.frame(test1,test2)
                  
ggplot(beta_table_vag, aes(x=variable, y=value)) +
  geom_boxplot(alpha=0.6) +
  theme(legend.position="none") +
  theme_bw() +
  geom_boxplot(fill=c("#56b1e9", "#3da822", "#e9568e","#e9568e"))
```

# 12 - Cluster samples based on their community composition (Bray-Curtis dissimilarities) (Figures 13 and 14)

## 

12.1 - Functions

```
# This R script is an extension of vegan library's bioenv()
# function and uses the bio.env() and bio.step() of
#	http://menugget.blogspot.co.uk/2011/06/clarke-and-ainsworths-bioenv-and-bvstep.html
#	The original author suggested these functions to overcome
#	the inflexibility of the bioenv() function which uses
#	a similarity matrix based on normalized "euclidean" distance.
# The new functions are given below and implement the following algorithms: 
# Clarke, K. R & Ainsworth, M. 1993. A method of linking multivariate community structure to environmental variables. Marine Ecology Progress Series, 92, 205-219.
# Clarke, K. R., Gorley, R. N., 2001. PRIMER v5: User Manual/Tutorial. PRIMER-E, Plymouth, UK.
# Clarke, K. R., Warwick, R. M., 2001. Changes in Marine Communities: An Approach to Statistical Analysis and Interpretation, 2nd edition. PRIMER-E Ltd, Plymouth, UK.
# Clarke, K. R., Warwick, R. M., 1998. Quantifying structural redundancy in ecological communities. Oecologia, 113:278-289. 

bv.step <- function(fix.mat, var.mat, 
                    fix.dist.method="bray", var.dist.method="euclidean", correlation.method="spearman",
                    scale.fix=FALSE, scale.var=TRUE,
                    max.rho=0.95,
                    min.delta.rho=0.001,
                    random.selection=TRUE,
                    prop.selected.var=0.2,
                    num.restarts=10,
                    var.always.include=NULL,
                    var.exclude=NULL,
                    output.best=10
){
  
  if(dim(fix.mat)[1] != dim(var.mat)[1]){stop("fixed and variable matrices must have the same number of rows")}
  if(sum(var.always.include %in% var.exclude) > 0){stop("var.always.include and var.exclude share a variable")}
  require(vegan)
  
  if(scale.fix){fix.mat<-scale(fix.mat)}else{fix.mat<-fix.mat}
  if(scale.var){var.mat<-scale(var.mat)}else{var.mat<-var.mat}
  
  fix.dist <- vegdist(as.matrix(fix.mat), method=fix.dist.method)
  
  #an initial removal phase
  var.dist.full <- vegdist(as.matrix(var.mat), method=var.dist.method)
  full.cor <- suppressWarnings(cor.test(fix.dist, var.dist.full, method=correlation.method))$estimate
  var.comb <- combn(1:ncol(var.mat), ncol(var.mat)-1)
  RES <- data.frame(var.excl=rep(NA,ncol(var.comb)), n.var=ncol(var.mat)-1, rho=NA)
  for(i in 1:dim(var.comb)[2]){
    var.dist <- vegdist(as.matrix(var.mat[,var.comb[,i]]), method=var.dist.method)
    temp <- suppressWarnings(cor.test(fix.dist, var.dist, method=correlation.method))
    RES$var.excl[i] <- c(1:ncol(var.mat))[-var.comb[,i]]
    RES$rho[i] <- temp$estimate
  }
  delta.rho <- RES$rho - full.cor
  exclude <- sort(unique(c(RES$var.excl[which(abs(delta.rho) < min.delta.rho)], var.exclude)))
  
  if(random.selection){
    num.restarts=num.restarts
    prop.selected.var=prop.selected.var
    prob<-rep(1,ncol(var.mat))
    if(prop.selected.var< 1){
      prob[exclude]<-0
    }
    n.selected.var <- min(sum(prob),prop.selected.var*dim(var.mat)[2])
  } else {
    num.restarts=1
    prop.selected.var=1  
    prob<-rep(1,ncol(var.mat))
    n.selected.var <- min(sum(prob),prop.selected.var*dim(var.mat)[2])
  }
  
  RES_TOT <- c()
  for(i in 1:num.restarts){
    step=1
    RES <- data.frame(step=step, step.dir="F", var.incl=NA, n.var=0, rho=0)
    attr(RES$step.dir, "levels") <- c("F","B")
    best.comb <- which.max(RES$rho)
    best.rho <- RES$rho[best.comb]
    delta.rho <- Inf
    selected.var <- sort(unique(c(sample(1:dim(var.mat)[2], n.selected.var, prob=prob), var.always.include)))
    while(best.rho < max.rho & delta.rho > min.delta.rho & RES$n.var[best.comb] < length(selected.var)){
      #forward step
      step.dir="F"
      step=step+1
      var.comb <- combn(selected.var, RES$n.var[best.comb]+1, simplify=FALSE)
      if(RES$n.var[best.comb] == 0){
        var.comb.incl<-1:length(var.comb)
      } else {
        var.keep <- as.numeric(unlist(strsplit(RES$var.incl[best.comb], ",")))
        temp <- NA*1:length(var.comb)
        for(j in 1:length(temp)){
          temp[j] <- all(var.keep %in% var.comb[[j]]) 
        }
        var.comb.incl <- which(temp==1)
      }
      
      RES.f <- data.frame(step=rep(step, length(var.comb.incl)), step.dir=step.dir, var.incl=NA, n.var=RES$n.var[best.comb]+1, rho=NA)
      for(f in 1:length(var.comb.incl)){
        var.incl <- var.comb[[var.comb.incl[f]]]
        var.incl <- var.incl[order(var.incl)]
        var.dist <- vegdist(as.matrix(var.mat[,var.incl]), method=var.dist.method)
        temp <- suppressWarnings(cor.test(fix.dist, var.dist, method=correlation.method))
        RES.f$var.incl[f] <- paste(var.incl, collapse=",")
        RES.f$rho[f] <- temp$estimate
      }
      
      last.F <- max(which(RES$step.dir=="F"))
      RES <- rbind(RES, RES.f[which.max(RES.f$rho),])
      best.comb <- which.max(RES$rho)
      delta.rho <- RES$rho[best.comb] - best.rho 
      best.rho <- RES$rho[best.comb]
      
      if(best.comb == step){
        while(best.comb == step & RES$n.var[best.comb] > 1){
          #backward step
          step.dir="B"
          step <- step+1
          var.keep <- as.numeric(unlist(strsplit(RES$var.incl[best.comb], ",")))
          var.comb <- combn(var.keep, RES$n.var[best.comb]-1, simplify=FALSE)
          RES.b <- data.frame(step=rep(step, length(var.comb)), step.dir=step.dir, var.incl=NA, n.var=RES$n.var[best.comb]-1, rho=NA)
          for(b in 1:length(var.comb)){
            var.incl <- var.comb[[b]]
            var.incl <- var.incl[order(var.incl)]
            var.dist <- vegdist(as.matrix(var.mat[,var.incl]), method=var.dist.method)
            temp <- suppressWarnings(cor.test(fix.dist, var.dist, method=correlation.method))
            RES.b$var.incl[b] <- paste(var.incl, collapse=",")
            RES.b$rho[b] <- temp$estimate
          }
          RES <- rbind(RES, RES.b[which.max(RES.b$rho),])
          best.comb <- which.max(RES$rho)
          best.rho<- RES$rho[best.comb]
        }
      } else {
        break()
      }
      
    }
    
    RES_TOT <- rbind(RES_TOT, RES[2:dim(RES)[1],])
    print(paste(round((i/num.restarts)*100,3), "% finished"))
  }
  
  RES_TOT <- unique(RES_TOT[,3:5])
  
  
  if(dim(RES_TOT)[1] > output.best){
    order.by.best <- RES_TOT[order(RES_TOT$rho, decreasing=TRUE)[1:output.best],]
  } else {
    order.by.best <-  RES_TOT[order(RES_TOT$rho, decreasing=TRUE), ]
  }
  rownames(order.by.best)<-NULL
  
  order.by.i.comb <- c()
  for(i in 1:length(selected.var)){
    f1 <- which(RES_TOT$n.var==i)
    f2 <- which.max(RES_TOT$rho[f1])
    order.by.i.comb <- rbind(order.by.i.comb, RES_TOT[f1[f2],])
  }
  rownames(order.by.i.comb)<-NULL
  
  if(length(exclude)<1){var.exclude=NULL} else {var.exclude=exclude}
  out <- list(
    order.by.best=order.by.best,
    order.by.i.comb=order.by.i.comb,
    best.model.vars=paste(colnames(var.mat)[as.numeric(unlist(strsplit(order.by.best$var.incl[1], ",")))], collapse=","),
    best.model.rho=order.by.best$rho[1],
    var.always.include=var.always.include,
    var.exclude=var.exclude
  )
  out
  
}

bio.env <- function(fix.mat, var.mat, 
                    fix.dist.method="bray", var.dist.method="euclidean", correlation.method="spearman",
                    scale.fix=FALSE, scale.var=TRUE,
                    output.best=10,
                    var.max=ncol(var.mat)
){
  if(dim(fix.mat)[1] != dim(var.mat)[1]){stop("fixed and variable matrices must have the same number of rows")}
  if(var.max > dim(var.mat)[2]){stop("var.max cannot be larger than the number of variables (columns) in var.mat")}
  
  require(vegan)
  
  combn.sum <- sum(factorial(ncol(var.mat))/(factorial(1:var.max)*factorial(ncol(var.mat)-1:var.max)))
  
  if(scale.fix){fix.mat<-scale(fix.mat)}else{fix.mat<-fix.mat}
  if(scale.var){var.mat<-scale(var.mat)}else{var.mat<-var.mat}
  fix.dist <- vegdist(fix.mat, method=fix.dist.method)
  RES_TOT <- c()
  best.i.comb <- c()
  iter <- 0
  for(i in 1:var.max){
    var.comb <- combn(1:ncol(var.mat), i, simplify=FALSE)
    RES <- data.frame(var.incl=rep(NA, length(var.comb)), n.var=i, rho=0)
    for(f in 1:length(var.comb)){
      iter <- iter+1
      var.dist <- vegdist(as.matrix(var.mat[,var.comb[[f]]]), method=var.dist.method)
      temp <- suppressWarnings(cor.test(fix.dist, var.dist, method=correlation.method))
      RES$var.incl[f] <- paste(var.comb[[f]], collapse=",")
      RES$rho[f] <- temp$estimate
      if(iter %% 100 == 0){print(paste(round(iter/combn.sum*100, 3), "% finished"))}
    }
    
    order.rho <- order(RES$rho, decreasing=TRUE)
    best.i.comb <- c(best.i.comb, RES$var.incl[order.rho[1]])
    if(length(order.rho) > output.best){
      RES_TOT <- rbind(RES_TOT, RES[order.rho[1:output.best],])
    } else {
      RES_TOT <- rbind(RES_TOT, RES)
    }
  }
  rownames(RES_TOT)<-NULL
  
  if(dim(RES_TOT)[1] > output.best){
    order.by.best <- order(RES_TOT$rho, decreasing=TRUE)[1:output.best]
  } else {
    order.by.best <- order(RES_TOT$rho, decreasing=TRUE)
  }
  OBB <- RES_TOT[order.by.best,]
  rownames(OBB) <- NULL
  
  order.by.i.comb <- match(best.i.comb, RES_TOT$var.incl)
  OBC <- RES_TOT[order.by.i.comb,]
  rownames(OBC) <- NULL
  
  out <- list(
    order.by.best=OBB,
    order.by.i.comb=OBC,
    best.model.vars=paste(colnames(var.mat)[as.numeric(unlist(strsplit(OBB$var.incl[1], ",")))], collapse=",") ,
    best.model.rho=OBB$rho[1]
  )
  out
}
```

## 12.2 - Subset samples (female, male or selected couples)

```
library(ggplot2)
library(ggdendro)
library(ggtext)
library(gridExtra)
library(grid)

ps_1000<-subset_samples(ps_1000, is_sample=="yes")
#ps_100<-get_top_taxa(ps_1000,100)

#SWITCH
#female samples
ps_switch<-subset_samples(ps_1000,sample_pairs=="fv")

#male samples
ps_switch<-subset_samples(ps_1000,sample_pairs=="ps")

#couples with possible microbiota interaction
ps_327<-subset_samples(ps_1000,couple=="327")
ps_415<-subset_samples(ps_1000,couple=="415")
ps_31<-subset_samples(ps_1000,couple=="31")
ps_340<-subset_samples(ps_1000,couple=="340")
ps_416<-subset_samples(ps_1000,couple=="416")
ps_317<-subset_samples(ps_1000,couple=="317")
ps_339<-subset_samples(ps_1000,couple=="339")
ps_311<-subset_samples(ps_1000,couple=="311")
ps_406<-subset_samples(ps_1000,couple=="406")

ps_switch<-merge_phyloseq(ps_327,ps_416,ps_317)

ps_switch<-ps_1000
```

## 12.3 - Make clusters

```
abund_table <- as.data.frame(otu_table(ps_switch))
abund_table <-abund_table[, colSums(abund_table != 0) > 0]
abund_table<-abund_table[ , colSums(is.na(abund_table)) == 0]


meta_table<-samdf_1000
meta_table <- meta_table[rownames(abund_table),]

#select specific metadata
#keep only numeric columns in meta_table
meta_table <- meta_table[,c(14,15,37,49,54), drop=FALSE]

meta_table <- sapply(meta_table, as.numeric)
meta_table <- as.data.frame(meta_table[ , apply(meta_table, 2, function(x) !any(is.na(x)))])
meta_table<-as.data.frame(meta_table[ , colSums(is.na(meta_table)) == 0])
meta_table <-meta_table[, colSums(meta_table != 0) > 0]

#Get grouping information
grouping_info<-data.frame(row.names=rownames(abund_table),t(as.data.frame(strsplit(rownames(abund_table),"_"))))
head(grouping_info)

#Parameters
cmethod<-"pearson" #Correlation method to use: pearson, pearman, kendall
fmethod<-"bray" #Fixed distance method: euclidean, manhattan, gower, altGower, canberra, bray, kulczynski, morisita,horn, binomial, and cao
vmethod<-"bray" #Variable distance method: euclidean, manhattan, gower, altGower, canberra, bray, kulczynski, morisita,horn, binomial, and cao
nmethod<-"bray" #NMDS distance method:  euclidean, manhattan, gower, altGower, canberra, bray, kulczynski, morisita,horn, binomial, and cao


res <- bio.env(wisconsin(abund_table), meta_table,fix.dist.method=fmethod, 
               var.dist.method=vmethod, correlation.method=cmethod,
               scale.fix=FALSE, scale.var=TRUE) 


envNames<-colnames(meta_table)
bestEnvFit<-""
for(i in (1:length(res$order.by.best$var.incl)))
{
  bestEnvFit[i]<-paste(paste(envNames[as.numeric(unlist(strsplit(res$order.by.best$var.incl[i], split=",")))],collapse=' + '), " = ",res$order.by.best$rho[i],sep="")
}
bestEnvFit<-data.frame(bestEnvFit)
colnames(bestEnvFit)<-"Best combination of environmental variables with similarity score"

betad<-vegdist(abund_table,method="bray")
res_adonis <- adonis2(betad ~ X2, grouping_info) 

hc <- hclust(betad)

#We will color the labels according to countries(group_info[,1])
hc_d <- dendro_data(as.dendrogram(hc))
hc_d$labels$Type<-grouping_info[as.character(hc_d$labels$label),2]

#Coloring function
gg_color_hue<-function(n){
  hues=seq(15,375,length=n+1)
  hcl(h=hues,l=65,c=100)[1:n]
}

# cols=gg_color_hue(length(unique(hc_d$labels$Type)))
cols=gg_color_hue(4)

colors<-c()

for(i in 1:length(hc_d$labels$Type)){
  
  if(hc_d$labels$Type[i]=="fol"){
    colors<-c(colors,cols[1])
    print(colors)
  }
  
  #else if
  else if(hc_d$labels$Type[i]=="pen"){
    colors<-c(colors,cols[2])
    print(colors)
  }
  else if(hc_d$labels$Type[i]=="spe"){
    colors<-c(colors,cols[3])
    print(colors)
  }
  else if(hc_d$labels$Type[i]=="vag"){
    colors<-c(colors,cols[4])
    print(colors)
  }
}


hc_d$labels$color=colors

hc_d$leaf_labels<-hc_d$labels$label

element_text <- function(family = NULL, face = NULL, colour = NULL,
                         size = NULL, hjust = NULL, vjust = NULL, angle = NULL, lineheight = NULL,
                         color = NULL, margin = NULL, debug = NULL, inherit.blank = FALSE) {
  
  if (!is.null(color))  colour <- color
  
  n <- max(
    length(family), length(face), length(colour), length(size),
    length(hjust), length(vjust), length(angle), length(lineheight)
  )
  if (n > 1) {
    cli::cli_warn(c(
      "Vectorized input to {.fn element_text} is not officially supported.",
      "i" = "Results may be unexpected or may change in future versions of ggplot2."
    ))
  }
  
  
  structure(
    list(family = family, face = face, colour = colour, size = size,
         hjust = hjust, vjust = vjust, angle = angle, lineheight = lineheight,
         margin = margin, debug = debug, inherit.blank = inherit.blank),
    class = c("element_text", "element")
  )
}

tip_names<-hc_d$labels$label
tip_names<-rev(tip_names)


## ######Plot cluster beta diversity ########
p1 <- ggplot(data = segment(hc_d)) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
  scale_x_discrete(limits=label(hc_d)$label)+
  coord_flip() +
  ylab("Distance (beta diversity = bray)") + theme_void() +
  theme(axis.text.y = element_text(color = hc_d$labels$color),
        axis.title.y = element_blank(),text = element_text(size = 6))
        
p1 <- p1 + geom_point(data=hc_d$label, aes(x = x, y = y, color = Type), inherit.aes =F, alpha = 0)
p1 <- p1 + guides(colour = guide_legend(override.aes = list(size=3, alpha = 1)))+
  scale_color_manual(values = cols)

p1 + theme(axis.text.y=element_text(size=rel(0.5))) + theme(
  panel.background = element_rect(fill='transparent'), #transparent panel bg
  plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
  panel.grid.major = element_blank(), #remove major gridlines
  panel.grid.minor = element_blank(), #remove minor gridlines
  legend.background = element_rect(fill='transparent'), #transparent legend bg
  legend.box.background = element_rect(fill='transparent') #transparent legend panel
)

pdf("Cluster_samples.pdf",height=10)
print(p1)
dev.off()
```

## 12.4 - Find the dominant ASV

```
otu <- otu_table(t(ps_switch))
otu<-otu@.Data

df_all = data.frame()

for(i in 1:ncol(otu)){
  result<-t(sort(otu[,i],decreasing = TRUE)[1])
  result<-result/(sum(otu[,i]))*100
  rownames(result)<-colnames(otu)[i]
  result<-melt(result)
  df_all <- rbind(df_all,result)
}

colnames(df_all)<-c("sample","asv","prevalence")

df_all<-left_join(df_all,asv_table,by="asv")
df_all$type<- str_sub(df_all$sample,-3)
df_all$dominant_species<-df_all$name

#melt and reorder
df_all_melt<-melt(df_all)
df_all_melt<-df_all_melt[match(label(hc_d)$label, df_all_melt$sample),]

df_all_melt$sample <- as.character(df_all_melt$sample)

#Then turn it back into a factor with the levels in the correct order
df_all_melt$sample <- factor(df_all_melt$sample, levels=unique(df_all_melt$sample))

#plot
prevplot<-ggplot(df_all_melt, aes(fill=sample, y=value, x=sample))+  
  geom_point() + theme(legend.position="none") + coord_flip()

###########dominant species - create table#####################
dominant_species<-data.frame(rev(df_all_melt$sample), rev(df_all_melt$dominant_species))
colnames(dominant_species)<-c("sample","dominant")
dominant_species$dominant<-gsub("NA", "sp.", dominant_species$dominant)
```

## 12.5 - Plot colors of each sample

```
###################plot colors - samples ########
plot(seq_len(length(colors)), rep_len(1, length(colors)),
     col = colors, pch = 15, cex = 1, xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
```

## 12.6 - Plot Shannon index for each sample

```
mytheme <- gridExtra::ttheme_minimal(col.just="left",base_size = 10,
                                     core = list(fg_params=list(hjust=0, x=0.1)),
                                     rowhead=list(fg_params=list(hjust=0, x=0)),
                                     padding=unit(c(4, 2), "mm"))

#p2<-tableGrob(dominant_species$dominant,theme=ttheme_minimal(base_size = 4))

p2<-tableGrob(dominant_species$dominant,theme=mytheme)

pdf("dominant_species.pdf", height=30, width=8)
grid.draw(p2)
dev.off()

###############alpha diversity####################
p<-plot_richness(ps_switch, measures=("Shannon")) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_x_discrete(limits = label(hc_d)$label) +
  coord_flip() +
  theme(
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent') #transparent legend panel
  )

pdf("shannon_samples.pdf",width= 2,height=18)
print(p)
dev.off()
```

# 13 - LEfSeR analysis (Figure 15)

```
## Genus level
library(lefser)
library(readr)
library(tibble)
library(dplyr)

setwd("~/luzern2021")

#agglomerate phyloseq object at the genus level
ps_genus<-subset_samples(ps_1000,is_sample=="yes")
ps_genus<-tax_glom(ps_genus, taxrank=rank_names(ps_genus)[6], NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))

#extract otu table (be sure to not include controls)
otu_table_genus<-ps_genus@otu_table@.Data
otu_table_genus<-as.data.frame(t(otu_table_genus))

#create table with asv
change_asv<-as.data.frame(colnames(t(otu_table_genus)))
change_asv$asv<-change_asv$`colnames(t(otu_table_genus)`

#retreive genera corresponding to each asv
taxa_table<-as.data.frame((as.matrix(ps_genus@tax_table@.Data)))
taxa_table$asv<-row.names(taxa_table)
keeps <- c("asv", "Genus")
change_genus<-taxa_table[keeps]

#join tables asv and genus table
genus_names<-left_join(change_asv, change_genus, by = "asv")

#change taxa names
rownames(otu_table_genus)<-genus_names$Genus
otu_table_genus<-as.data.frame(otu_table_genus)

keeps<-colnames(otu_table_genus)

#import metadata file
metadata<-read.csv(
  file = "~/luzern2021/00_metadata/luzern_metadata_final_290322.csv",header = T)
row.names(metadata)<-metadata$sample

#keep only samples (no controls)
metadata<-metadata[keeps,]

#create SummarizedExperiment
se_genus <- SummarizedExperiment(otu_table_genus, colData = metadata)
se_genus<-se_genus[, se_genus$age_opu != "na"]

#subset for each sample
se_genus_vag<-se_genus[, se_genus$sample_type == "vag"]
se_genus_fol<-se_genus[, se_genus$sample_type == "fol"]
se_genus_pen<-se_genus[, se_genus$sample_type == "pen"]
se_genus_spe<-se_genus[, se_genus$sample_type == "spe"]

se_genus_f<-se_genus[, se_genus$sample_pairs == "fv"]
se_genus_m<-se_genus[, se_genus$sample_pairs == "ps"]

#perform lefser

###comparison between female samples###
res_lefser_f <- lefser(se_genus_f, groupCol = "sample_type", blockCol = NULL)
lefserPlot(res_lefser_f)

#plot relative abundance
ps_genus_f<-subset_samples(ps_genus,sample_pairs=="fv")
ps_genus_rel= transform_sample_counts(ps_genus_f, function(x) x/sum(x))

test<-as.data.frame(t(ps_genus_rel@otu_table@.Data))
test$ASV<-rownames(test)
test<-left_join(test,ASVtable,by="ASV")
rownames(test)<-test$Genus
test<-test[c(1:111)]

keeps<-res_lefser_f$Names
test<-as.data.frame(t(test[keeps,]))

test$sample_type<-(str_split_fixed(rownames(test), "_", 2)[,2])
test<-melt(test) 

test %>% ggplot(aes(x=variable, y=value)) +
  geom_boxplot()+
  coord_flip() +
  theme_light()

###male samples###
res_lefser_m <- lefser(se_genus_m, groupCol = "sample_type", blockCol = NULL)
lefserPlot(res_lefser_m)

#Relative abundance plot
ps_genus_m<-subset_samples(ps_genus,sample_pairs=="ps")
ps_genus_rel= transform_sample_counts(ps_genus_m, function(x) x/sum(x))

test<-as.data.frame(t(ps_genus_rel@otu_table@.Data))
test$ASV<-rownames(test)
test<-left_join(test,ASVtable,by="ASV")
rownames(test)<-test$Genus
test<-test[c(1:92)]

keeps<-res_lefser_m$Names
test<-as.data.frame(t(test[keeps,]))

test$sample_type<-(str_split_fixed(rownames(test), "_", 2)[,2])
test<-melt(test) 

test %>% ggplot(aes(x=variable, y=value)) +
  geom_boxplot()+
  coord_flip() +
  theme_light()


######comparison female vs male#####
res_lefser_FvsM <- lefser(se_genus, groupCol = "sample_pairs", blockCol = NULL)
lefserPlot(res_lefser_f)
```

# 14 - Supplementary figures

## 14.1 - Percent stacked barplot of the most abundant genera (Supplementary figure 1)

```
library(phyloseq)
library(microbiome)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(visibly)
library(grid)
library(ggpubr)
library(randomcoloR)
library(V8)

#Agglomerate to genus level and subset
ps_genus<-tax_glom(ps_1000, taxrank=rank_names(ps_1000)[6], 
                   NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))

taxa_table<- ps_genus@tax_table
taxa_table<-taxa_table@.Data
taxa_table<-as.matrix(taxa_table)
taxa_table<-as.data.frame(taxa_table)

genus_vag<-subset_samples(ps_genus,sample_type=="vag")
genus_fol<-subset_samples(ps_genus,sample_type=="fol")
genus_pen<-subset_samples(ps_genus,sample_type=="pen")
genus_spe<-subset_samples(ps_genus,sample_type=="spe")

l<-list("genus_vag"=genus_vag,"genus_fol"=genus_fol,"genus_pen"=genus_pen,"genus_spe"=genus_spe)

#Make graph using for loop
graphs<-list()
graph.names<-c("Vagina","Follicular fluid","Penis","Sperm")

for(i in 1:length(l)){
  #Extraxct otu table
  otu.table<-as(otu_table(l[[i]]), "matrix")
  otu.table<-as.data.frame(otu.table)

  #top.genera<-top_n(as.data.frame(otu.table), 10)
  top.genera<-names(sort(taxa_sums(l[[i]]), TRUE)[1:10])
  
  #Extract matching rows
  genera10 <- as.data.frame(t(otu.table[,top.genera]))
  genera10$ASV<-row.names(genera10)
  genera10<-left_join(genera10,ASVtable,by="ASV")
  row.names(genera10)<-genera10$Genus
  genera10<-dplyr::select(genera10, -c('ASV', 'Row.names',"Genus"))
  
  
  genera_plot<-c("Lactobacillus","Gardnerella","Prevotella","Bifidobacterium",
                "Atopobium","Alloscardovia","Ureaplasma","Dialister","Campylobacter",
                "Porphyromonas","Staphylococcus","Streptococcus","Sphingomonas",
                "Aquabacterium","Caulobacter","Afipia","Ezakiella",
                "Escherichia/Shigella","Peptoniphilus","Prevotella_6","Pelomonas") 
  
  #Rename unknown
  rownames(genera10)[rownames(genera10) == "Unknown"] <- "z_Other"
  genera10_df<-as.data.frame(genera10)
  
  for(j in 1:length(diff[[i]])){
    name<-diff[[i]][j]
    
    genera10_df <- genera10_df %>%
      add_column(name = c(0,0,0,0,0,0,0,0,0,0))
  }
    
  genera10_df$genus<-rownames(genera10_df)
  
  #Reshape table for ggplot2
  ncols<-(ncol(genera10_df)-1)
  genera10.ggplot<-pivot_longer(genera10_df,cols=1:all_of(ncols),
                                names_to = "sample",values_to = "reads")
  
  #Plot stacked
  bar.plot<-ggplot(genera10.ggplot, aes(fill=genus, y=reads, x=sample)) + 
    geom_bar(position="fill", stat="identity") +
    theme(text = element_text(size=8),axis.text.x = element_text(angle = 90)) 

  graphs[[i]] <-ggplot(genera10.ggplot, aes(fill=genus, y=reads, x=sample)) + 
    geom_bar(position="fill", stat="identity") +
    theme(text = element_text(size=8),axis.text.x = element_text(angle = 90)) +
    ggtitle(print(graph.names[i]))
}
  
#find all the genera 
uniques<-c(unique(graphs[[1]]$data$genus),unique(graphs[[2]]$data$genus),unique(graphs[[3]]$data$genus),unique(graphs[[4]]$data$genus))
uniques<-c(unique(uniques))

combinations<-vector()

for(i in 1:length(uniques)){
  combinations<-c(combinations,uniques[i],code[i])
}

combinations<-c("Lactobacillus"="#d45dda","Gardnerella"="#FF0000",
                "Prevotella"="#05B00B","Bifidobacterium"="#FFC3BD",
                "Atopobium"="#7C30FF","Alloscardovia"="#8bfc96",      
                "Ureaplasma"="#75B8FF","Dialister"="#F57F18",
                "Campylobacter"="#a53a0b","Porphyromonas"="#e8a13a",
                "Staphylococcus"="#dcf223","Streptococcus"="#844c73",
                "Sphingomonas"="#3f6960","Aquabacterium"="#c4b8d6",
                "Caulobacter"="#96ae24","Afipia"="#cf8f23",
                "Ezakiella"="#beea82","Escherichia/Shigella"="#5bf4f2",
                "Peptoniphilus"="#d1b598","Prevotella_6"="#709047",
                "Pelomonas"="#5cab93","z_Other"="#8c8c8c") 
  
#all sample names
all<-row.names(ps@otu_table@.Data)

#samples names after filtration
filtered<-row.names(ps_1000@otu_table@.Data)

#all missing samples
diff<-setdiff(all, filtered)

#create missing samples for each sample_type
diff_vag <- grep("vag", diff, value = TRUE)
diff_fol <- grep("fol", diff, value = TRUE)
diff_spe <- grep("spe", diff, value = TRUE)
diff_pen <- grep("pen", diff, value = TRUE)

#add pen missing samples (n=3)
diff_pen<-c(diff_pen,"0320_pen","0408_pen","3112_pen")

#create a list with all missing samples
diff<-list(diff_vag,diff_fol,diff_pen,diff_spe)

for(i in 1:length(l)){

  #Extraxct otu table
  otu.table<-as(otu_table(l[[i]]), "matrix")
  otu.table<-as.data.frame(otu.table)
  
  #top.genera<-top_n(as.data.frame(otu.table), 10)
  top.genera<-names(sort(taxa_sums(l[[i]]), TRUE)[1:10])
  
  #Extract matching rows
  genera10 <- as.data.frame(t(otu.table[,top.genera]))
  genera10$ASV<-row.names(genera10)
  genera10<-left_join(genera10,ASVtable,by="ASV")
  row.names(genera10)<-genera10$Genus
  genera10<-dplyr::select(genera10, -c('ASV', 'Row.names',"Genus"))
  
  #Rename unknown
  rownames(genera10)[rownames(genera10) == "Unknown"] <- "z_Other"
  genera10_df<-as.data.frame(genera10)
  
  #add genus column
  for(a in 1:length(rownames(genera10_df))){
        if (rownames(genera10)[a] %in% genera_plot) {
            genera10_df$genus[a]<-rownames(genera10)[a]
        } else {
          genera10_df$genus[a]<-"z_Other"
        }
    }

  #add missing samples  
  genera10_df[,diff[[i]]]<-c(0,0,0,0,0,0,0,0,0,0)
    
  #reorder columns
  genera10_df<-genera10_df[ , order(names(genera10_df))]
  
  #melt table
  genera10.ggplot<-melt(genera10_df)
  
  #Plot stacked
  bar.plot<-ggplot(genera10.ggplot, aes(fill=genus, y=value, x=variable)) + 
    geom_bar(position="fill", stat="identity") +
    theme(text = element_text(size=8),axis.text.x = element_text(angle = 90)) 
  
  graphs[[i]] <-ggplot(genera10.ggplot, aes(fill=genus, y=value, x=variable)) + 
  geom_bar(position="fill", stat="identity") +
  theme(text = element_text(size=8),axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values=combinations) +
  ggtitle(print(graph.names[i]))

}

pdf("luzern_stackedplot.pdf",width = 30,height = 16)
grid.newpage()
grid.draw(rbind(ggplotGrob(graphs[[1]]),ggplotGrob(graphs[[2]]), 
                ggplotGrob(graphs[[3]]), ggplotGrob(graphs[[4]]),
                size = "last"))
dev.off()  
```

## 14.2 - Dominant ASVs

```
top.taxa <- tax_glom(ps_1000, taxrank=rank_names(ps_1000)[7], 
                     NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))
otu <- otu_table(t(top.taxa))
otu<-otu@.Data

df_all = data.frame()

for(i in 1:ncol(otu)){
result<-t(sort(otu[,i],decreasing = TRUE)[1])
result<-result/(sum(otu[,i]))*100
rownames(result)<-colnames(otu)[i]
result<-melt(result)
df_all <- rbind(df_all,result)
}

colnames(df_all)<-c("sample","asv","relative_abundance")
df_all_sample<-df_all

df_all<-left_join(df_all,asv_table,by="asv")
df_all$type<- str_sub(df_all$sample,-3)
df_all$dominant_species<-df_all$name
df_all$count <- table(df_all$name)[df_all$name]
df_all$name2 <- ifelse(df_all$count>=3, df_all$name, "z_other")

#Plot dominant ASVs counts
p <- ggplot(df_all, aes(name2)) 
p + theme_bw() + theme(axis.text.x = element_text(angle = 75, vjust = 1, hjust=1)) +
  geom_bar(aes(fill=type))

#Plot the relative abundance of the dominant ASV in the samples
q <-ggplot(df_all, aes(name2)) + theme_bw() +
    geom_boxplot(aes(y=relative_abundance)) +
    theme(axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank()) 
```
