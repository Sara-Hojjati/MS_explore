# System requirements
Operating system: Ubuntu 20.04.4 LTS

Python 3.9.12
pandas 1.5.2
numpy 1.23.5
scikit-learn 1.1.2
statsmodels 0.13.2
scipy 1.9.1

R 4.2.1
ChAMP 2.26.0
ggplot2 3.4.0
limma 3.52.4
clusterProfiler 4.4.4
org.Hs.eg.db 3.15.0
stats 3.6.2
verification 1.42
cutpointr 1.1.2

# Installation guide
Typical install time on "normal" desktop computer: 2h

Python installation: https://www.python.org/downloads/
R installation: https://www.r-project.org/

## Install python packages from command line using pip
pip install pandas numpy scikit-learn statsmodels scipy

## Install R packages from R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ChAMP")
BiocManager::install("limma")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
install.packages("ggplot2")
install.packages("verification")
install.packages("cutpointr")

# Demo
## Input files:
Files provided in the folder data:
- Proteins below limit of detection (LOD): data/proteins_below_LOD.csv
- Example file for MS enrichment analysis: data/known_ms_biomarkers.txt
- Proteins in the MS network: data/ms_network.txt
- Proteins in PPI-network from STRINGdb (ENTREZ IDs): data/string_proteins_background.txt
- R data with CSF DEPs in discovery cohort: data/proteinnames.RData

Files not provided due to data access restrictions:
- CSF protein data: data/CSF_proteindata.csv. Contains protein levels measured using the Olink Explore platform which uses Proximity Extension Assay (PEA) technology. Comma separated file with the headers SampleID, OlinkID, UniProt, Assay, LOD, NPX.
- Plasma protein data: data/plasma_proteindata.csv. Contains protein levels measured using the Olink Explore platform which uses Proximity Extension Assay (PEA) technology. Comma separated file with the headers SampleID, OlinkID, UniProt, Assay, LOD, NPX.
- Clinical data: data/clinical_data.csv. Contains clinical data associated with each sample. Tab separated file with the headers sample_ID, sex, age, NEDA_EDA_2years, nARMSS, treatment_duration_index.

## Run pipeline
The scripts should be run in the outlined order to guarantee that all necessary files are generated.

Open a terminal and set the working directory to the folder containing the scripts.

1. Preprocessing
python preprocessing_control_correction.py
Rscript preprocessing_combat.R
python preprocessing_add_clinical_info.py

2. Differential expression analysis and MS enrichment analysis
Rscript dea.R
python MS_enrichment_analysis.py

3. Diagnosis model
Rscript logistic_regression_singlevariable_models_diagnosis.R
Rscript logistic_regression_stepwise_models_diagnosis.R

4. NEDA model
Rscript logistic_regression_singlevariable_models_disease_activity.R
Rscript logistic_regression_disease_activity_model.R

5. nARMSS model
python linearregression_model_narmss.py

6. GO enrichment analysis of MS network
Rscript GO_enrichment_analysis.R

## Expected output:
- Data with NPX values stored in data/CSF_data.tsv and data/Plasma_data.tsv
- Batch corrected data stored in data/CSF_corrected_data.tsv and Plasma_corrected_data.tsv
- Differentially expressed proteins in CSF dicovery cohort stored in results/Discovery_DEP.txt
- Results from MS enrichment analysis on example file with known MS biomarkers.
- Diagnosis: Performance of individual proteins for diagnosis are stored in results/Disease_activity_model_AUCs.tsv.
- Diagnosis: Selected stepwise model and performance score
- Disease activity: Performance of individual proteins for diagnosis are stored in results/Diagnosis_model_AUCs.tsv.
- Disease activity: Selected stepwise model and performance score
- Disability worsening: nARMSS model proteins and perforance score
- Significant GO terms from enrichment analysis of MS network proteins stored in data/ms_network_GO_enrichment.tsv.

# Instructions for use
The scripts can be used as outlined in Demo. To run on other input files, the file names need to be replaced within the scripts.
The scripts provided are not intended for use on other datasets then those provided.
