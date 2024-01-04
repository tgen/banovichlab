This directory contains custom scripts to process and analyze the data presented in

## Cell type-specific and disease-associated eQTL in the human lung

Heini M Natri<sup>1</sup>, Christina B Del Azodi<sup>2</sup>, Lance Peter<sup>1</sup>, Chase Taylor<sup>3</sup>, Sagrika Chugh<sup>2</sup>, Robert Kendle<sup>1</sup>, Jonathan Kropski<sup>3</sup>, Davis McCarthy<sup>2</sup> & Nicholas E Banovich<sup>1</sup>

<sup>1</sup>Translational Genomics Research Institute (TGen), Phoenix, AZ;
<sup>2</sup>St. Vincent’s Institute of Medical Research, Melbourne, Australia;
<sup>3</sup>Vanderbilt University Medical Center, Nashville, TN

Corresponding author information: Nicholas E Banovich, nbanovich@tgen.org

For detailed information, see the BiorXiv preprint: URL

### Single-cell sequence data processing and cell type annotation

10x data were processed and integrated using <i>Seurat</i> v4 (processing_integration.R) and cell types were annotated based on marker gene expression (annotation.R).

### eQTL calling with <i>LIMIX</i> and <i>mashr</i>
A snakemake pipeline for reproducing the eQTL results is located in [a separate repository](https://gitlab.svi.edu.au/biocellgen-public/musj_2021_multi-omics-lung-CBA). Plotting and multi-cell type eQTL analysis was carried out using custom software (url pending).

### Colocalization analysis using <i>coloc</i>

Colocalization with GTEx bulk-eQTL and lung trait GWAS was carried out using a bayesian method implemented in <i>coloc</i> (coloc.R).

### Enrichment testing

To test for the enrichment of cell type eQTL among IPF GWAS implicated variants, a null set of non-significant eQTL was selected using <i>nullranges</i> and significance was calculated using Fisher's test (enrichment.R).
