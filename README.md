# FitBOLD

This is a Matlab toolbox used to predict brain regions correlated to timed behaviour.

fitBOLD extracts BOLD signal from subjects' fMRI scans, preprocesses it, clusters brain regions, then uses BOLD signal to perform various multivariate linear and quadratic regression and outputs 3D maps of correlated brain regions with their estimated beta coefficients in addition to different regression statistics e.g. RMSE, AIC, BIC.

Please see the example workflow in the root folder. xMind is need to view the file. you can download it  https://www.xmind.net/download/

What it assumes:
- Behavioural experiment done under fMRI scanning. 
- First level GLM is already estimated by SPM: http://www.fil.ion.ucl.ac.uk/spm/ (only of if dynamic BOLD signal shift is used).

What can this toolbox do?

Preprocessing:

- Z-scoring.
- Pooling subjects' data.
- Combining bilateral masks of the same brain region into one.
- K-means clustering of brain regions. Distance measure, min and max k are adjustable, Average silhouette of clusters determines the optimal K.
- Low signal drift BOLD filtering and bandpass FIR filtering.
- Shifting BOLD signal in each sub-region by either a standard (2-scans) or using a dynamic shift estimated by SPM's first level model.

BOLD Modelling:
- Building desing matrix of BOLD signals, behavioural data, and categorical variables.
- Parallel regression algorithms built from scratch.
- Backward stepwise regression, optimized for improving adjusted R-square, decreasing AIC and increasing F-statistic of the model.
- Allows for different regression methods, backwards stepwise regression, OLS multilinear regression, lasso, ridge regression, fused lasso. 
- Allows for bootstrapping the data, subject-based cross validation of data and bootstrapping of the cross validated training data.

Output:
- Regression statistics: correlated brain regions and their beta estimates, ordinary and adjusted R-squares, RMSE, AIC, BIC , and 95% confidence intervals.
- Plotting of Confidence Intervals of regression statistics and parameter estimates.
- plotting of tested cross validation, R-square, and RMSE.


