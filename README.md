# covid_paper
 
This repository contains codes for our protein-COVID MR and coloc analysis from Open Targets Genetics.

The GSMR_coloc_COVID_analysis_codes.Rmd contains R and shell codes that were used to format the protein (Sun et al) and COVID outcome datasets (Oct 2020, release 4) for GSMR analysis.

The base file_for_cisMR_coloc_analysis.Rmd contains R codes that were used to generate a file (sun_proteins_ensid_coloc.csv) that contained eligible proteins from Sun et al (excluding proteins in X/Y chromosomes and proteins that have no chromosome number) with their cis-region parameters.

The make_coloc_dataset_covid.R, make_coloc_dataset_sun.R, and make_cis_dataset_sun.R used the sun_proteins_ensid_coloc.csv file to generate datasets for coloc and cis-MR analysis

The coloc_pqtl_covid.R creates a single file containing pair-wise coloc and single-SNP MR analysis results
