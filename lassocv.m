function b=lassocv(XTRAIN,ytrain,XTEST)
[lso,fitstats]=lasso(XTRAIN,ytrain,'CV',10);
b=XTEST*(lso(:,fitstats.IndexMinMSE));
end