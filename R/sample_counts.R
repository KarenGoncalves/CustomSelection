#' Example - Counts data frame
#'
#' Counts of 3 samples (4 replicates per sample) of Arabidopsis thaliana genes.
#'
#' Transgenic Arabidopsis thaliana Columbia-0 plants expressing GFP alone (Control) or fused to a candidate secreted effector protein of the fungus Melampsora larici-populina (Mlp37347 or Mlp124499) were used for the transcriptome analysis.
#'
#' RNA was extracted from pooled aerial tissue of 2-week-old soil-grown plants, doing four replicates per genotype. Libraries were generated using the TruSeq Stranded mRNA Library Prep kit (Illumina) and 100 ng of total RNA. The libraries were sequenced with Illumina HiSeq 4000 Sequencer paired-end reads of 100nt.
#'
#' Trimmomatic (LEADING:4 TRAILING:4 SLIDINGWINDOW:4:20 MINLEN:20) and then the surviving paired reads were aligned to the TAIR10 assembly of the genome of A. thaliana with TopHat v2.0.14 in Galaxy (default options, with average mate inner distance varying for each replicate and standard deviation of distance between pairs of 50 base pairs).
#'
#' Further analyses were done using R software v.3.2.5. Genomic ranges of Arabidopsis transcripts were obtained from Ensembl plants with GenomicFeatures and overlaps of sequencing reads with the transcripts were counted using GenomicAlignments, using options for paired-end reads and union mode.
#' @format A data frame with 32833 rows (genes) and 12 columns (samples). Each column shows the read count of the genes for a replicate of a sample.
#' @source NCBI GEO \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136038}
"sample_counts"
