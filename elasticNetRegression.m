% this calls Matlab's lasso.m function with optional input "Lambda" set to
% 0.5, so that lasso runs an elastic net regression. See the help section
% in lasso.m for implementation explanation, and below for motivation:

% From Katha's paper (August 2013 version):

% In case of multicollinearity, the prediction accuracy of ordinary least
% squares regression can be reduced and the results are difficult to
% interpret (e.g., see Hoerl & Kennard, 1970). Compared to ordinary least
% squares regression, regularized regression obtains higher prediction
% accuracy (in particular if multicollinearity exists) and provides a
% simpler model by selecting the most informative predictors. Regularized
% regression methods include an additional regularization term in the cost
% function (e.g., an l1-norm regularizer as in the ?Lasso method?; see
% Tibshirani, 1996; or an l2-norm as in ?ridge regression?; see Hoerl &
% Kennard, 1970) for which a regularization parameter ? defines the degree
% to which small coefficients are penalized. Recently, a regularized
% regression method was proposed which combines Lasso and ridge regression
% (called ?Elastic net?; see Zou & Hastie, 2005). Elastic net selects the
% most important predictors under consideration of multicollinearity (Zou &
% Hastie, 2005) where an additional parameter ? (0 < ? <= 1) defines the
% weight of lasso (? = 1) versus ridge regression (? = 0).

% With increasing values of ?, elastic net retains optic flow and facial
% action activation as nonzero coefficients while the latter is set to zero
% last. Note that for ? = 0, the coefficients are equivalent to ordinary
% least squares regression. For a ? value of 0.06 the mean prediction
% squared error as calculated by cross-validation was minimal (MSE = 0.03;
% see black dashed vertical line in Figure 8). The best fitting model with
% ? = 0.06 explained 61% of the variance in the behavioral choices (R =
% 0.78; F(13) = 33.77; p < 0.0001). The fitted coefficients for this model
% were Beta = 0.88 for facial action activation, Beta = 0.46 for optic flow
% and Beta = 0 for Gabor similarity.

%% Specify data
[data, text] = xlsread('ClareData.xlsx')
idx = [2:6 8:12]; % we only want the parts of the questionnaires, not the total scores
X = data(:,idx);
Y = data(:,end-2);
predictorNames = text(idx);

% % Example data:
% load imports-85
% X = X(~any(isnan(X(:,1:16)),2),:);
% Y = X(:,16);
% Y = log(Y);
% X = X(:,3:15);
% predictorNames = {'wheel-base' 'length' 'width' 'height' ...
%   'curb-weight' 'engine-size' 'bore' 'stroke' 'compression-ratio' ...
%   'horsepower' 'peak-rpm' 'city-mpg' 'highway-mpg'};

%% Run elastic net regression, create best models
% Compute the default sequence of lasso fits.
[B,S] = lasso(X,Y,'CV',10,'PredictorNames',predictorNames,'Alpha',.5);
% lassoPlot(B,S,'PlotType','Lambda');
lassoPlot(B,S);
ylabel('Predictor weights')
set(gcf,'Name','Trace plot of coefficients')
printfig

% What variables are in the model corresponding to minimum
% cross-validated MSE, and in the sparsest model within one
% standard error of that minimum.
minMSEModel = S.PredictorNames(B(:,S.IndexMinMSE)~=0)
sparseModel = S.PredictorNames(B(:,S.Index1SE)~=0)

figure('name','PredictorWeightsInBestModels');
bar(B(:,[S.IndexMinMSE S.Index1SE])); xlabel_oblique(predictorNames)
set(gca,'position',[0.1300    0.1500    0.7750    0.7750],'xlim',[0 11])
ylabel('Predictor weights')
legend({'model with min MSE';'sparsest model within 1SE of min'})
title('Predictor weights for different best-fitting models')
printfig
% Fit the best model, compare to original data and examine residuals.
Xplus = [ones(size(X,1),1) X];
fitSparse = Xplus * [S.Intercept(S.Index1SE); B(:,S.Index1SE)];
fitMinMSE = Xplus * [S.Intercept(S.IndexMinMSE); B(:,S.IndexMinMSE)];
r2 = corr(Y,fitMinMSE)^2;
r_res = corr(fitMinMSE,Y-fitMinMSE);

% Same for sparsest model:
% r2 = corr(Y,fitSparse)^2;
% r_res = corr(fitSparse,Y-fitSparse);

%% display results
figure
subplot(1,2,1)
plot(Y,fitMinMSE,'o'); xlabel('Data'); ylabel('Fit')
title(sprintf('Min-MSE model: lambda = %.3f, r^2 = %.2f',S.LambdaMinMSE,r2))

subplot(1,2,2)
plot(fitMinMSE,Y-fitMinMSE,'o')
title(sprintf('Residuals of MinMSE model fit (r = %.2f)',r_res))
