 # A joint model for the estimation of species distributions and environmental characteristics from point-referenced data
 
This repository contains the data, codes, and results for our paper. 

Abstract:
> Background: Predicting and explaining species occurrence using environmental characteristics is essential for nature conservation and management. Species distribution models consider species occurrence as the dependent variable and environmental conditions as the independent variables. Suitable conditions are estimated based on a sample of species observations, where one assumes that the underlying environmental conditions are known. This is not always the case, as environmental variables at broad spatial scales are regularly extrapolated from point-referenced data. However, treating the predicted environmental conditions as accurate surveys of independent variables at a specific point does not take into account their uncertainty. 
> Methods: We present a joint hierarchical Bayesian model where models for the environmental variables, rather than a set of predicted values, are input to the species distribution model. All models are fitted together based only on point-referenced observations, which results in a correct propagation of uncertainty. We use 50 plant species representative of the Dutch flora in natural areas with 8 soil condition predictors taken during field visits in the Netherlands as a case study. We compare the proposed model to the standard approach by studying the difference in associations, predicted maps, and cross-validated accuracy. 
> Findings: We find that there are differences between the two approaches in the estimated association between soil conditions and species occurrence (correlation 0.64-0.84), but the predicted maps are quite similar (correlation 0.82-1.00). The differences are more pronounced in the rarer species. The cross-validated accuracy is substantially better for 5 species out of the 50, and the species can also help to predict the soil characteristics. The estimated associations tend to have a smaller magnitude with more certainty. 
> Conclusion: These findings suggests that the standard model is often sufficient for prediction, but effort should be taken to develop models which take the uncertainty in the independent variables into account for interpretation.

The codes compare the two-stage model and the joint model in identical interpretation and prediction problems.

 ## Data curation

The file **paper_data.Rmd** describes the curation and combination steps of the two primary data sets:
 - The National Flora Monitoring Network - Environment and Nature Quality (LMF-M\&N) [link](https://www.rivm.nl/publicaties/ontwerp-landelijk-meetnet-flora-milieu-natuurkwaliteit-lmf-mn)
 - Wageningen University \& Research abiotic factors [link](https://library.wur.nl/WebQuery/wurpubs/reports/367477)

These data sets are not directly available to download by public but may be obtained for research from their respective authors. 

These data sets are linked into additional spatial open data sets in the Netherlands: 
- Provinces (Bestuurlijke Gebieden) [link](https://www.pdok.nl/geo-services/-/article/bestuurlijke-gebieden)
- FGR regions (Fysisch Geografische Regio’s) [link](https://www.pdok.nl/introductie/-/article/fysisch-geografische-regio-s)
- Landuse (Bestand Bodemgebruik) [link](https://www.pdok.nl/-/bestand-bodemgebruik-2015-van-cbs-nu-bij-pdok)
- Soiltype (Grondsoorten) [link](http://www.geodesk.nl/Grondsoorten.htm)

This results in the following files, which are required to run the experiments:
- **data_paper/data_plots.csv** publishes a subset of the data used in our paper with the author's permission. 
- **data_paper/data_grid.csv** contains the Netherlands grid for which predictions were made.
- **data_paper/NL.shp** contains a shapefile for the land boundary of the Netherlands.

## Codes

The file **paper_codes.Rmd** contains the data set statistics, experiments, and visualization of the results.

The experiments take a very long time to run, so we provide an alternative [LSF cluster](https://www.ibm.com/support/pages/what-lsf-cluster) implementation in the folder **lsf/**.
See the file **lsf/run.sh** for submission of the jobs and combination of the results, where each job fits a species specific SDM to a given data set.

The results from running the experiments are saved in the following files:
- **data_paper/results.csv** contains the model parameters and predictions for the 50 species in the entire data set.
- **data_paper/predictions_[province/fgr].csv** contains the model predictions for the 50 species in the validation data.
- **data_paper/prevalences_[province/fgr].csv** contains posterior prevalences for the 50 species in the validation data.

These files are sufficient to reproduce all of the tables and visualizations in the paper.

The file **simulation.R** contains the experiment to generate simulated data, fit the models, calculate accuracy and visualize the simulation. 

