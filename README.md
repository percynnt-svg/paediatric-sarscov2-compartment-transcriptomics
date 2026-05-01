\# Paediatric SARS-CoV-2 compartment-specific transcriptomic analysis



This repository contains R scripts used for a secondary transcriptomic analysis of paediatric acute SARS-CoV-2 infection using BRAVE cohort RNA-sequencing data from GSE231409.



The analysis compares upper respiratory tract and peripheral blood responses in children and adolescents under 21 years old. The workflow includes differential gene expression analysis, Module 61 immune pathway enrichment, symptom-associated immune module analysis, and paired URT–blood correlation analysis.



\## Data availability



Raw RNA-sequencing data are publicly available from the Gene Expression Omnibus under accession number GSE231409. De-identified clinical metadata are available from the associated COVID\_RNASeq\_Age GitHub repository.



Raw data and metadata are not included in this repository.



\## Scripts



Scripts are provided in the `scripts/` folder:



1\. `01\_volcano\_differential\_expression.R`  

&#x20;  Performs differential gene expression analysis and generates volcano plots.



2\. `02\_module61\_fgsea\_enrichment.R`  

&#x20;  Performs Module 61 immune pathway enrichment using ranked limma results.



3\. `03\_symptom\_module\_heatmap.R`  

&#x20;  Performs symptom-associated immune module analysis in acute SARS-CoV-2-positive participants.



4\. `04\_paired\_urt\_blood\_correlation.R`  

&#x20;  Performs paired URT–blood ssGSEA module scoring and Pearson correlation analysis.



\## Software



Analyses were performed in R version 4.5.2 using packages including edgeR, limma, DESeq2, fgsea, GSVA, ggplot2, dplyr, readxl and writexl.



\## Notes



These scripts document the analytical workflow used for the dissertation figures. File paths may need to be adjusted depending on the local directory structure.

