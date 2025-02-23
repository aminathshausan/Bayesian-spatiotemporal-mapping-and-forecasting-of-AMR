# Bayesian-spatiotemporal-mapping-and-forecasting-of-AMR

This repository makes available the source code for the manuscript: 
"Bayesian spatiotemporal mapping and forecasting of antimicrobial resistant infection with application to northern Australia".
The computation concerns model fits and forecast validations for 4 types of Generalised Linear Mixed Effects Models (GLMM): 
Poisson, negative binomial, and their zero-inflated versions. 
The computation is performed using Integrated Laplace Approximation (INLA) method. 
All codes are implemented using the R software.

For ethical reasons, the actual data has been provided here with noise. 
The actual data can be visualised and requested from the HOTspots website: https://amr-hotspots.net 

# Prerequisites

To run the codes you will require R(>=3.30) and latest version of RStudio.

# Directory structure
```
├── project
│   ├── data
│   ├── results
│   ├── src
│   │   ├── fit_models_noCovt.R
│   │   ├── plots.R
│   │   ├── posterior_linr_combs.R
```
- The `data` folder contains all processed data required to run the codes in the `src` folder. 
- The `results` folder contains all model fits from 5-fold-temporal cross validations (as .RData format) and summary of these fits.
- The `src` folder contains all source codes required to produce the results and images in the manuscript.

# How to run scripts

1. Run the script `fit_models_noCovt.R` to  fit the models.
2. Run the script `plots.R` to generate all figures and summary of results.  
 