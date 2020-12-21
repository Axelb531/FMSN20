
# 1. Classification of fMRI data using MRFs and Gibbs (MCMC):
The task is to classify the pixels in an image, assuming an MRF structure for the pixel classes, and conditionally independent Gaussian observations,
􏰉y(si)|z(si, k) = 1􏰊 ∈ N 􏰇μk, Σk􏰈
Study the fMRI data and the initial analysis presented in fmri class init.m. The idea is to first perform data reduction using either PCA directly on the fMRI data or by first regressing the data onto the temporal regressors in X and then performing PCA on the regression coefficients. Here we see the image as having either 160 (the time points) or 11 (the regression coefficients) “colour layers”.
For the classification use PCA to reduce the data dimension, and then clas- sify the images using the increasingly complex models provided by kmeans.m, normmix gibbs.m and the MRF model.

Sampling procedures for the MRF-model can be constructed using the func- tions for simulation of the field (using Gibbs sampling) and parameters (using MCMC) in MRFs provided in the course files (see list below). A suitable start- ing point for the MRF-estimation are the parameter estimates and classification obtained from normmix gibbs.m.
Study the difference between the classifications and investigate the effects of: 1) Different initial data reductions; 2) Different neighbourhood sizes (4 or 8); and 3) Having the same or different β:s for the classes.
Use the special functions for MRF:s
mrf gaussian post.m (posterior alpha-parameters for Gaussian data)
mrf sim.m (simulate samples from an MRF model)
gibbs alpha beta.m (pseudo-likelihood MH for MRF-parameters)
gibbs mu Sigma.m (Sampling from μ and Σ given a set of observations from a Gaussian-distribution.)
           
Spatial statistics with image analysis, Home assignment 3 2
 Note that the posterior for α, β is tricky, and it is typically not possible to fully explore the distribution. However, it should converge to a parameter set for which the results of classification works.
For a similar model used in computer tomography see

  https://arxiv.org/abs/1607.02188
