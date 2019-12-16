# CustomSelection - package for selection of reference genes from RNAseq data

This package calculates the Transcripts per Million (TPM) data frame from the read count matrix, calculates the minimum expressin level for a gene to be considered as expressed in each sample and selects as reference genes those with lowest covariance.

It contains three main functions: __Counts_to_tpm__, __DAFS__ and __gene_section__. The function __customReferences__ merges the three main functions and only returns the data frame with the selected genes, their average TPM and the covariance of their TPM.

The package only needs the __read count matrix__ and the __length of the genes__ in the read count matrix.
We offer for a test a data frame of read counts for three samples (4 replicates per sample) of *Arabidopsis* plants expressing GFP (control) or one of two candidate effectors from *Melampsora larici-populina* (Mlp37347-GFP and Mlp124499-GFP). We also offer the gene lengths from Arabidopsis thaliana TAIR10 obtained from [Ensembl](plants.ensembl.org) with __biomaRt__.
