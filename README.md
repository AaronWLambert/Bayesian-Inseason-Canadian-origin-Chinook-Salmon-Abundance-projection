---
output:
  pdf_document: default
  html_document: default
---
# Bayesian Inseason Canadian-origin Chinook Salmon Abundance Projection
## Purpose
A management tool created to generate an inseason projection for each day during the summer Chinook salmon run on the Yukon River to estimate the end-of-season escapement of Canadian-origin Chinook salmon into Canda. This method utilizes a Bayesian updating approach using Rstan in R. Users can generate projections retrospectively for model testing and selection purposes or to generate inseason projections for the current year.
## Folders
|File|Description|
|:---:|---|
|Data|Data files that are used in the projection model. These data are generally preprocessed in an R script and fed into the Stan model. The two subfolders labeled as ADFG Eagle/PSS Daily Reports is used to download passage observations inseason to run the r-markdown scripts.|
|R|Scripts used to preprocess data and to call the Stan model using Rstan. |
|Stan| Stan scripts that are called in the Rstan model. These include different versions of the model utilizing methods that are currently being assessed for their predictive performance.|
|Mathematical Model Descriptions| Descriptions of the model versions with the accompanying math and equations.|
|How to run projection models| Brief directions for running the model inseason and retrospectively.|
|ADFG Normal Daily Projection| R markdown script for generating daily projection reports using model *PSSnormal_ESprop*. Instructions in the How To Run Projection Model PDF.|
|ADFG Proportion Model Daily Projection| R markdown script for generating daily projection reports using model *PSSprop_ESprop*. Instructions in the How To Run Projection Model PDF.|
## Note
This project is in progress and all results are preliminary.
