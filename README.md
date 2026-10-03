# MN-Trajectory

## Overview
This repository contains the code and supplementary materials associated with our publication: 

**"Early Detection of Myeloid Neoplasm Using Blood Count Trajectories"**

### Abstract
Myeloid neoplasm (MNs) precursor states cannot be effectively screened at large scale. Here, we compared longitudinal blood count data from individuals who subsequently developed de novo acute myeloid leukemia (AML, N=365), myelodysplastic syndrome (MDS, N=629), and primary myelofibrosis (MF, N=69) to controls (N=121,739). We collected over 814M data points covering patient age, gender and blood count trajectories between 2002-2024 and developed machine learning models predicting the onset of MNs up to five years before clinical diagnosis. We validated the models in both an internal and an external test set demonstrating high predictive accuracy without reliance on genomic profiling. By exploring the connection between blood counts and genomic hallmarks, we describe persistent granulocytosis as a promising precursor state of JAK2-mutated MF. These findings demonstrate the feasibility of early detection of pre-MN, facilitating personalized surveillance and early intervention.


## Repository Contents
This repository includes:
1. **Data Preparation**:
   - Scripts for processing blood count trajectory data.
   - Scripts to handle missing data and ensure data consistency.
   - Feature engineering scripts to derive trajectories from raw blood counts.

2. **Machine Learning Models**:
   - Code to train, validate, and test predictive models.


# Application
We developed an application for easy use of our models. The app calculates the risk of developing acute myeloid leukemia (AML), myelodysplastic syndrome (MDS), and primary myelofibrosis (MF) based on laboratory data trajectories, age, and sex.

### Try the App:
[**MN-Trajectory Risk Prediction Application**](<https://hematoscopelab.shinyapps.io/mn-trajectory/>)

## How to Use the Repository
### Requirements
- Python 3.8 or higher, Python dependencies are described at `./requirements.txt` (for machine learning models)
- R 4.4.2 and R libraries are described at `./sessioninfo.txt` (for data preparation, statistical analyses and some plots)
  
