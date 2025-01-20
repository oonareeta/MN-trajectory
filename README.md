# MN-Trajectory

## Overview
This repository contains the code and supplementary materials associated with our publication: 

**"Early Detection of Myeloid Neoplasm Using Blood Count Trajectories"**

### Abstract
Despite the discovery of precursor states to myeloid neoplasm (MNs), their large-scale screening is not yet possible. Here, we compared longitudinal blood count data from individuals who subsequently developed *de novo* acute myeloid leukemia (AML, N=361), myelodysplastic syndrome (MDS, N=749), and primary myelofibrosis (MF, N=93) to controls without later cancer diagnosis (N=75 714).

By characterizing patient age, gender and pre-diagnostic blood count trajectories, we collected over 67M data points between 2002-2024 and developed machine learning models that accurately predicted the onset of MNs up to five years before clinical diagnosis. Notably, these models achieved high predictive accuracy without genomic profiling in contrast to other available pre-MN risk scores. We validated the models in both a hold-out and an independent test set(N=183 AML, N=468 MDS, N=49 MF, N=10 058 controls).

By exploring the connection between pre-MN blood counts and diagnostic genomic and chromosomal alterations, we revealed the connection between elevated erythroid RDW, bone marrow erythroid dysplasia, *TP53*mut, *SRSF2*mut and *NPM1*wt as well as poor prognosis. These findings demonstrate the feasibility of earlier and accessible detection of pre-MN, which is a prerequisite to study and apply personalized surveillance plans and early therapy.

## Repository Contents
This repository includes:
1. **Data Preparation**:
   - Scripts for processing blood count trajectory data.
   - Scripts to handle missing data and ensure data consistency.
   - Feature engineering scripts to derive trajectories from raw blood counts.

2. **Machine Learning Models**:
   - Code to train, validate, and test predictive models.

3. **Results**:
   - Performance metrics of the models (e.g., accuracy, sensitivity, and specificity).
   - Scripts to visualize predictions.

# Application
We developed an application for easy use of our models. The app calculates the risk of developing acute myeloid leukemia (AML), myelodysplastic syndrome (MDS), and primary myelofibrosis (MF) based on laboratory data trajectories, age, and sex.

### Try the App:
[**MN-Trajectory Risk Prediction Application**](<https://hematoscopelab.shinyapps.io/mn-trajectory/>)

## How to Use the Repository
### Requirements
- Python 3.8 or higher, Python dependencies are describes at `./requirements.txt` (for machine learning models)
- R 4.4.2 and R libraries are described at `./sessioninfo.txt` (for data preparation, statistical analyses and some plots)
  
