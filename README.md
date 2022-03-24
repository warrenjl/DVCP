# DVCP: Directionally-Varying Change Points Model

This package implements a Directionally-Varying Change Points (DVCP) model that aims to estimate the magnitude of the impact of a point source as well as its range of influence across the spatial domain.  DVCP includes a Gaussian process with directionally-defined correlation structure nested within a change point framework to introduce unique change point parameters in every direction extending from the point source.  The Gaussian predictive process approximation is used to facilitate model fitting for large datasets, and DVCP is fit using Markov chain Monte Carlo sampling techniques. Please see the "DVCP_Model_Details" and "DVCP_Example" folders for more specific information regarding the statistical model and package use details, respectively.

# Reference
* Song J and Warren JL (2022). A directionally-varying change points model for quantifying the impact of a point source. Journal of Agricultural, Biological and Environmental Statistics, 27(1):46-62.
* https://link.springer.com/article/10.1007/s13253-021-00466-y

