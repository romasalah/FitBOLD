function [resulttbl,stats]=strp_parstepwiselm(strpopts,fit,idxatinitialmodel,trainsubjects_i)
allsubs=1:size(fit.sub,2);
S=allsubs(trainsubjects_i);
S_i_train=fit.S_i(trainsubjects_i);
nboots=strpopts.nboots;
[~,b]=bootstrp(nboots,[],S);
ar_ar=cell(1,nboots);
hbar = parfor_progressbar(nboots,['Bootstrap sampling']);  %create the progress bar
tic
yhat_i_sample2=cell(1,nboots);
yhat_i_sample=cell(1,nboots);
if fit.regress.feature_selection==1 % Don't redo feature selection  
    parfor i = 1:nboots
        stats(i).niter=nboots;
        stats(i).nloops=nboots;
        stats(i).idxfrominput=idxatinitialmodel;
        S_low=b(:,i);
        for unpackterms_i=1:length(S_low)
            yhat_i_sample{i}=vertcat(yhat_i_sample{i}, S_i_train{S_low(unpackterms_i),1}');
        end
        sample=array2table([fit.designquad(yhat_i_sample{i},idxatinitialmodel) fit.design{yhat_i_sample{i},end}]);
        fitstats{i}=fitlm(sample,'intercept',false);
        hbar.iterate(1);   % update progress by one iteration
        stats(i).fitstats=fitstats{i};
        stats(i).yhat_i_sample=yhat_i_sample{i};
    end
else
    parfor i2 = 1:nboots
        S_low=b(:,i2);
        for unpackterms_i=1:length(S_low)
            yhat_i_sample2{i2}=vertcat(yhat_i_sample2{i2}, fit.S_i{S_low(unpackterms_i),1}');
        end
        sample=fit.design(yhat_i_sample2{i2},:);
        [fitstats{i2},~,nloops{i2},niter{i2},~,idxfrominput{i2}]=parstepwiselm(sample,strpopts.power,strpopts.Incriterion,strpopts.Outcriterion,strpopts.upper,strpopts.fit.cpu_parallel,strpopts.fit);
        stats(i2).idxfrominput=idxfrominput{i2};
        stats(i2).nloops=nloops{i2};
        stats(i2).niter=niter{i2};
        stats(i2).fitstats=fitstats{i2};
        stats(i2).yhat_i_sample=yhat_i_sample{i2};
    end
end
close(hbar);   %close progress bar
hbar = parfor_progressbar(nboots,['Getting statistics']);  %create the progress bar
parfor i3=1:nboots
hbar.iterate(1);   % update progress by one iteration
AV=anova(fitstats{i3},'summary');
fstat=AV.F(2);
if fstat==Inf; fstat=0; end
pVal=AV.pValue(2);
if isnan(pVal); pVal=1; end
[acfc,~,e2]=autocorr_c(fitstats{i3}.Residuals.Raw,10);
for h= 1: size(acfc,1)
    ar_ar{i3}(h)=(acfc(h)>e2(1)) + (acfc(h)<e2(2));
end
max_ar=find([ar_ar{i3}]==1,1,'last');
if isempty(max_ar) max_ar=0; end
rnames=fitstats{i3}.Coefficients.Properties.RowNames;
coefnum=fitstats{i3}.NumCoefficients;
for findquad= 1:size(fitstats{i3}.Coefficients.Properties.RowNames,1)
    quads_i=find(~cellfun(@isempty,strfind(rnames,'^2')));
    quads_names={rnames(quads_i)};
    quads_num=length(quads_i);
    lin_num= coefnum-quads_num;
end
result(i3,:)=[fitstats{i3}.ModelCriterion.AIC fitstats{i3}.ModelCriterion.BIC fitstats{i3}.Rsquared.Ordinary...
    fitstats{i3}.Rsquared.Adjusted fitstats{i3}.RMSE fitstats{i3}.NumCoefficients quads_num lin_num fstat pVal];

end
resulttbl=array2table(result); resulttbl.Properties.VariableNames={'AIC','BIC','rsqr','adjrsqr','RMSE','coefnum'...
    'quadsnum','linnum','fstat','pval'};
close(hbar);
elapsed=toc;
  fprintf(['\n Elapsed time:' num2str(elapsed)])
end
