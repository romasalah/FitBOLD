function [resulttbl,stats,namesquad]=cv_parstepwiselm(Opts,fit)
%% Create Polynomials
[termsquad,namesquad]=Polynom_fitBOLD(fit.design,2);

S=(1:size(fit.sub,2))';
kfolds=Opts.kfolds;
CVO=cvpartition(S,'k',kfolds);
for sample_i= 1:CVO.NumTestSets
    trsub = CVO.training(sample_i);trIdx=[fit.S_i{find(trsub==1),1}]';
    tesub = CVO.test(sample_i);teIdx=[fit.S_i{find(tesub==1),1}]';
    design_train=fit.design(trIdx,:);
    termstest=termsquad(teIdx ,:); tgttest=fit.design{teIdx,end};
    [fitstats,~,nloops,niter,namesquad{sample_i},idxfrominput]=parstepwiselm(design_train,Opts.power,Opts.Incriterion,Opts.Outcriterion,Opts.upper,Opts.fit.cpu_parallel,Opts.fit);
%% Collect training statistics
stats(sample_i).fitstats=fitstats;
stats(sample_i).nloops=nloops;
stats(sample_i).niter=niter;
stats(sample_i).idxfrominput=idxfrominput;
AV=anova(fitstats,'summary');
fstat=AV.F(2);
if fstat==Inf; fstat=0; end
pVal=AV.pValue(2);
if isnan(pVal); pVal=1; end
[acfc,~,e2]=autocorr_c(fitstats.Residuals.Raw,10);
for h= 1: size(acfc,1)
    ar_ar{sample_i}(h)=(acfc(h)>e2(1)) + (acfc(h)<e2(2));
end
max_ar=find([ar_ar{sample_i}]==1,1,'last');
if isempty(max_ar) max_ar=0; end
for findquad= 1:size(fitstats.Coefficients.Properties.RowNames,1)
    quads_i=find(~cellfun(@isempty,strfind(fitstats.Coefficients.Properties.RowNames,'^2')));
    quads_names={fitstats.Coefficients.Properties.RowNames(quads_i)};
    quads_num=length(quads_i);
    lin_num= fitstats.NumCoefficients-quads_num;
end
result(sample_i,:)=[fitstats.ModelCriterion.AIC fitstats.ModelCriterion.BIC fitstats.Rsquared.Ordinary...
    fitstats.Rsquared.Adjusted fitstats.RMSE fitstats.NumCoefficients quads_num lin_num fstat pVal];

stats(sample_i).yhat_i_sample=trIdx;
stats(sample_i).yhat_i_test=teIdx;

%% test prediction performance
betas=zeros(1,size(termstest,2));
betas(1,idxfrominput)=fitstats.Coefficients.Estimate(1:end-1);
betas(1,end)=fitstats.Coefficients.Estimate(end);
stats(sample_i).testpred= betas*termstest';
stats(sample_i).predrsqr= corr(stats(sample_i).testpred',tgttest).^2;

end
resulttbl=array2table(result); resulttbl.Properties.VariableNames={'AIC','BIC','rsqr','adjrsqr','RMSE','coefnum'...
    'quadsnum','linnum','fstat','pval'};
end