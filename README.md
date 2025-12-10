# PARTAGE: Parallel analysis of replication timing and gene expression
The human genome is partitioned into functional compartments that replicate at specific times during the S-phase. This temporal program, referred to as replication timing (RT), is co-regulated with the 3D genome organization, is cell type-specific, and changes during development in coordination with gene expression. Moreover, RT alterations are linked to abnormal gene expression, genome instability, and structural variation in multiple diseases, including cancer. However, mechanistic links between RT, large-scale 3D genome architecture, and transcriptional regulation remain poorly understood. A major limitation is that current approaches require the separate profiling of RT and transcriptomes from independent batches of samples, obscuring the complex co-regulation between the epigenome and transcriptome. Here, we developed PARTAGE, a multiomics approach that enables joint profiling of copy number variation (CNV), RT, and gene expression from the same sample, providing a more accurate integrative view of the complex relationships between RT and gene regulation.


<img width="3323" height="1606" alt="Figure1_v1" src="https://github.com/user-attachments/assets/816b104f-0672-4301-9d4d-4c32c12d1352" />



Scripts for the analysis of PARTAGE data
Consists of four jupyter notebook files

PARTAGE_RT.ipynb: Used for calculating and plotting partage data. 
Requires RPKM data for G1 samples and other timepoints. Will output 
correlation graph and PARTAGE heatmap representing the replication 
timing.

rna_enrichment_corr.ipynb: Used for calculating clusters and 
enrichment analysis of expression data, along with self-correlation and 
correlation between pseudo-RT and expression.

CNV_Analysis_plot.ipynb: Used for plotting copy number variation along 
with called deletions and duplications called by CNVpytor. Outputs calls
text file for use with ProOvErlap

convert_calls_to_bed.ipynb: convert text file output to bed for
ProOvErlap
