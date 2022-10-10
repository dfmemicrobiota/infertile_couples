# infertile_couples

In this study, we want the characterize genital microbiota of infertile couples.

Samples include:

- vaginal swab

- follicular liquid

- semen

- penis swab

The analysis was performed by sequencing the variable regions V1-V2 of the 16S rRNA gene amplified using custom barcoded primers (F-27/R-338) containing Illumina sequencing adapters.

Illumina sequencing was performed at the Lausanne Genomic Technologies Facility (GTF) of the Lausanne University using an Illumina MiSeq instrument in paired-end mode 2 x 250 nt.

Demultiplexing of the raw sequencing data was performed with illumina-utils package. Reads were processed with the DADA2 pipeline. Decontam package was used to remove potential contaminations based on negative controls.Phyloseq and ampvis2 R packages were used for analysis and graphical visualization of the microbiota results. R package maaslin2 was used to infer differentially abundant bacteria.
