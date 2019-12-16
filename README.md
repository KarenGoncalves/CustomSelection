# CustomSelection - package for selection of reference genes from RNAseq data

This package calculates the Transcripts per Million (TPM) data frame from the read count matrix, calculates the minimum expressin level for a gene to be considered as expressed in each sample and selects as reference genes those with lowest covariance.

It contains three main functions: __Counts_to_tpm__, __DAFS__ and __gene_selection__. The function __customReferences__ merges the three main functions and only returns the data frame with the selected genes, their average TPM and the covariance of their TPM.

The package only needs the __read count matrix__ and the __length of the genes__ in the read count matrix.
We offer for a test a data frame of read counts for three samples (4 replicates per sample) of *Arabidopsis* plants expressing GFP (control) or one of two candidate effectors from *Melampsora larici-populina* (Mlp37347-GFP and Mlp124499-GFP). We also offer the gene lengths from Arabidopsis thaliana TAIR10 obtained from [Ensembl](plants.ensembl.org) with __biomaRt__.

### Count matrix example
Deposited in NCBI GEO under accession [GSE136038](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136038)

### Counts_to_tpm
Modified from [Slowkow gist](https://gist.github.com/slowkow/c6ab0348747f86e2748b)

### DAFS
Obtained from [George and Chang (2014)](doi:10.1186/1471-2105-15-92)

## Workflow

Starting with raw reads run:
1. Quality filtering [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
2. Alignment ([HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml), [TopHat2](https://ccb.jhu.edu/software/tophat/index.shtml)) or do a De novo assembly ([Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki), [Trans-abyss](https://github.com/bcgsc/transabyss)), then align the reads to the assembly
3. Sort the alignment ([samtools sort](http://samtools.sourceforge.net/))
4. Get the read count for each gene in each sample (proposed workflow in R bellow)
5. Run Counts_to_tpm, DAFS and gene_selection OR run customReferences


### Getting read counts

1. Use the package "GenomicFeatures" to construct a transcript database from biomart. Example, for _Arabidopsis thaliana_:

> library("GenomicFeatures")

> ath <- makeTxDbFromBiomart(biomart="plants_mart",
                        dataset="athaliana_eg_gene",
                        transcript_ids=NULL,
                        circ_seqs=DEFAULT_CIRC_SEQS,
                        filter=NULL,
                        id_prefix="ensembl_",
                        host="plants.ensembl.org",
                        taxonomyId=3702,
                        miRBaseBuild=NA)

> tx <- transcriptsBy(ath)

2. Create a table with the names of the samples as they are in the bam and bai files. Then use the package "Rsamtools" to create the bam file list. Example, for three bam files of three replicates of the control sample (CNT-1, CNT-2, CNT-3):

> samples <- c(paste0("CNT-", 1:3))

> sampleTable <- data.frame(sampleName = samples, fileNameBam = paste0(samples, ".bam"), fileNameBai = paste0(samples, ".bai"))

> library("Rsamtools")

> bfl <- BamFileList(file.path(sampleTable$fileNameBam), file.path(sampleTable$fileNameBai), yieldSize = 20000)

3. Use the package "GenomicAlignments" to count the number of reads aligned to each gene in the transcript database. Then, create a data frame of read counts per sample. Example for paired reads:

> txMat=list()\nfor (i in names(bfl)){

> txMat[[i]] <- assays(summarizeOverlaps(tx, bfl[[i]], mode = Union, singleEnd=F, ignore.strand=F, fragments=T))$counts
  
> write.csv(as.data.frame(txMat[[i]]), file=paste(gsub(".bam","",i), ".csv", sep=""))}

> counts <- do.call(cbind, txMat)

> colnames(counts) <- samples
