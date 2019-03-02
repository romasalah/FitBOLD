# FitBOLD

This is a Matlab toolbox used to predict brain regions correlated to timed behaviour.

fitBOLD extracts BOLD signal from subjects' fMRI scans, preprocesses it, clusters brain regions, then uses BOLD signal to perform various multivariate linear and quadratic regression and outputs 3D maps of correlated brain regions with their estimated beta coefficients in addition to different regression statistics e.g. RMSE, AIC, BIC. 

What it assumes:
- Behavioural experiment done under fMRI scanning. 
- First level GLM is already estimated by SPM: http://www.fil.ion.ucl.ac.uk/spm/

What can this toolbox do?


- Parallel regression algorithms built from scratch.
- Backward stepwise regression, optimized for improving adjusted R-square, decreasing AIC and increasing F-statistic of the model.
- Allows for different regression methods, backwards stepwise regression, OLS multilinear regression, lasso, ridge regression, fused lasso. 
- Allows for bootstrapping the data, subject-based cross validation of data and bootstrapping of the cross validated training data.
- K-means clustering of brain regions. Distance measure, min and max k are adjustable, Average silhouette of clusters determines the optimal K.
- Calculation and plotting of Confidence Intervals of regression statistics and parameter estimates.
-Added option for FIR filtering of the data during preprocessing.
- removed low signal drift in preprocessing.
- Code cleanse, organization, and bug fixes.
