# Evaluating methods for integrating scRNA-seq and MERFISH spatial transcriptomics data using publicly available cancer datasets

### Human Tumor Atlas Network (HTAN) Data Jamboree | Dec. 4-7, 2023 

---

Here, we explore integration methods to leverage the complementary advantages of scRNA-seq and MERFISH spatial transcriptomic technologies. scRNA-seq enables whole transcriptome profiling at the single cell level, but does not preserve the spatial context of the tissue. MERFISH provides transcript signatures for single cells while preserving cell spatial information within a tissue, but is limited in transcipt panel size and sequencing depth. Integrating the two assays that were run on paired samples from an HTAN metastatic breast cancer cohort allows us to 1. evaluate he best tools for integrating this type of data, and 2. elucidate new biological insights not previously established by analyzing the modalities separately. 

*Disclaimer* This work was done over a very short period for a hackathon, and represents a preliminary effort. With this in mind, code is not optimized, and there is ample room for further work. 

---

## Background

scRNA-seq allows for profiling of the entire transcriptome at the single cell level, but fails to preserve spatial information within the respective tissue. The development of spatial transcriptomic assays solves this issue by directly profiling of RNA transcripts from a tissue, thus preserving the spatial information. However, [current image-based spatial transcriptomic (ST) assays](https://www.nature.com/articles/s41592-022-01409-2) are limited to the number of genes they can profile and sequencing depth. There is motivation to integrate scRNA-seq data and ST data to leverage the strengths of each technology. Many computational tools and algorithms have been developed to carry out the integration process, but use different approaches. For instance, [MaxFuse](https://www.nature.com/articles/s41587-023-01943-0) was developed generally to integrate weakly-linked modalities through iterative joint embeddings. In contrast, [Tangram](https://www.nature.com/articles/s41592-021-01264-7) was more specifically developed for integrating scRNA-seq and spatial transcriptomics data. Here, we use both MaxFuse and Tangram to integrate scRNA-seq and MERFISH ST paired datasets, and begin to leverage the integrated data to draw novel biological conclusions. 

## About the data 

The data used in this project were accessed using the [HTAN Data Portal](https://humantumoratlas.org/explore), and were generated as part of the [Human Tumor Atlas Pilot Project](https://humantumoratlas.org/hta1) (HTAPP). A manuscript describing and analyzing some of the samples from this dataset is currently in [pre-print on BioRxiv](https://doi.org/10.1101/2023.03.21.533680).

Analyzed here are paired scRNA-seq and MERFISH datasets from 9 samples derived from patients with various subtypes of metastatic breast cancer. scRNA-seq datasets were processed using either a fresh tissue or frozen tissue protocol, which dictated the generation of scRNA-seq or snRNA-seq datasets, respectively. Experimental parameters such as isolation of cells or nuclei, fresh or frozen will affect both RNA quality and detectable transcriptome population.

Figure 1 shows total cell counts from the MERFISH samples. Each MERFISH dataset had a paired scRNA-seq dataset derived from the same biospecimen. For this project, we successfully integrated all datasets using both Tangram and MaxFuse; however, analysis was only performed on the sample **HTAPP-313** due to time contraints. Future work would involve a comprehensive comparison of the entire cohort. 

# Figure 1.
![plot](./figures/HTAPP_MERFISH_sample_cellcount_barplot.png)

## Analysis goals 

- Integrate the data using both MaxFuse and Tangram 
- Transfer cell type labels between modalities 
- Map scRNA-seq transcriptional data into the MERFISH spatial context 
- Determine pro's and con's of each integration method  
- Draw novel biological conclusions from the integrated dataset 

## Basic workflow 

- **Data ingress**
    - HTAN level 3 and 4 data ([Learn more about the HTAN data model](https://humantumoratlas.org/standards)) for scRNA-seq and MERFISH were
      accessed using syanpseclient. See notebook: `data_ingress/Data_extraction.ipynb`
    - Metadata was accessed and merged with datasets using the HTAN data portal and Google BigQuery 
- **Pre-processing**
    - Separate pre-processing pipelines were set up for MERFISH and scRNA-seq. In short these processes are described below:
    - MERFISH
        - Cell centroids were calculated for each cell from the middle slice of the 3D z-stack using the Shapely Python package. See
          notebook: `pre_processing/MERFISH_preprocessing.ipynb`
        - Control probes were removed from the expression matrix, all transcript counts were normalized, and cell types were annotated
          manually using Leiden clustering implemented in Scanpy. 
          See notebook for an example on one sample: `pre_processing/MERFISH_celltype_514.ipynb`
    - scRNA-seq
        - Using Scanpy, cells were filtered based on mitochondrial counts, genes were filtered based on expression, and transcript counts
          were normalized. Clustering and gene rankings, also implemented in Scanpy, were used to validate cell type labels that were available in the HTAN metadata. 
          See notebook: `pre_processing/preprocess_scrna.ipynb`
- **Integration**
- **Evaluation & Downstream analysis**

# Figure 2.
![overview](./figures/pipeline_overview.png)

## Results 

![results1](./figures/MaxFuse_integration_results.png)

![results2](./figures/tangram_313_scRNAfeatures_MERFISHlabels.png)

![results3](./figures/MERFISH_integrated_spatial.png)

## Future directions 

## Team 

- **John Duc Nguyen**, Genentech 
- **Jon Akutagawa**, University of California, San Francisco
- **Xiyue Zhao**, Oregon Health and Science University
- **Mark Dane**, Oregon Health and Science University
- **Cameron Watson**, Oregon Health and Science University
- **Gautam Machiraju**, Stanford University
