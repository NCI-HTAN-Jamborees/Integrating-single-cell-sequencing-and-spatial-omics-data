# Evaluating methods for integrating scRNA-seq and spatial transcriptomics data using publicly available cancer datasets

### Human Tumor Atlas Network (HTAN) Data Jamboree | Dec. 4-7, 2023 

## Background

scRNA-seq allows for profiling of the entire transcriptome at the single cell level, but fails to preserve spatial information within the respective tissue. The development of spatial transcriptomic assays solves this issue by directly profiling of RNA transcripts from a tissue, thus preserving the spatial information. However, [current image-based spatial transcriptomic (ST) assays](https://www.nature.com/articles/s41592-022-01409-2) are limited to the number of genes they can profile and sequencing depth. There is motivation to integrate scRNA-seq data and ST data to leverage the strengths of each technology. Many computational tools and algorithms have been developed to carry out the integration process, but use different approaches. For instance, [MaxFuse](https://www.nature.com/articles/s41587-023-01943-0) uses weakly linked features for it's integration to allow integration across distinct modalities. In contrast, Seurat uses an anchor-based integration method, which depends on strongly linked features to transfer label annotations from a reference to query dataset. Here, we compare the computational methods MaxFuse and Seurat for the integration of scRNA-seq and ST data using data collected for the Human Tumor Atlas Network HTAPP atlas. 

## About the data 

The data used in this project were accessed using the [HTAN Data Portal](https://humantumoratlas.org/explore), and were generated as part of the [Human Tumor Atlas Pilot Project](https://humantumoratlas.org/hta1) (HTAPP). A manuscript describing and analyzing some of the samples from this dataset is currently in [pre-print on BioRxiv](https://doi.org/10.1101/2023.03.21.533680).

Analyzed here are paired scRNA-seq and MERFISH datasets from 9 samples derived from patients with various subtypes of metastatic breast cancer. scRNA-seq datasets were processed using either a fresh tissue or frozen tissue protocol, which dictated the generation of scRNA-seq or snRNA-seq datasets, respectively. Experimental parameters such as isolation of cells or nuclei, fresh or frozen will affect both RNA quality and detectable transcriptome population.


![plot](./figures/HTAPP_MERFISH_sample_cellcount_barplot.png)

## Analysis goals 

- Integrate the data 
- Transfer cell type labels between modalities 
- Determine which method works best for this particular problem 
- Draw novel biological conclusions from the integrated dataset 

## Basic workflow 

- Data ingress and pre-processing 
- Integration: `Method 1`
- Integration: `Method 2`
- Evaluation
- Downstream analysis

## Results 

## Future directions 

## Team 

- **John Duc Nguyen**, Genentech 
- **Jon Akutagawa**, University of California, San Francisco
- **Xiyue Zhao**, Oregon Health and Science University
- **Mark Dane**, Oregon Health and Science University
- **Cameron Watson**, Oregon Health and Science University
- **Gautam Machiraju**, Stanford University
