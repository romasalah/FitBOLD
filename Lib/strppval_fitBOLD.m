function pval=strppval_fitBOLD(fitstatsall,fitstats0,idxfrominput0,trIdx,fit)
nboots=1000;

for boot_i=1:nboots
fitstatssample=fitstatsall(boot_i).fitstats;
tinitstrp=zeros(1,fit.quadsize);
idxfrominputsample=fitstatsall(boot_i).idxfrominput;
tstatsample=fitstatssample.Coefficients.tStat;
for get_tnull=1:size(idxfrominputsample,2)
    tinitstrp(idxfrominputsample(get_tnull))=tstatsample(get_tnull);
end
if fitstatssample.Coefficients{end,1}~=0
    tinitstrp(end)=tstatsample(end);
else
    tinitstrp(end)=0;
end
tstatall(boot_i,:)=tinitstrp;
end

tinit=zeros(1,fit.quadsize);
tstat=fitstats0.Coefficients.tStat;
for get_tnull=1:size(idxfrominput0,2)
    tinit(idxfrominput0(get_tnull))=tstat(get_tnull);
end
if fitstats0.Coefficients{end,1}~=0
    tinit(end)=tstat{end,1};
else
    tinit(end)=0;
end
if fit.regress.useCImu
    tstatmu=mean(tstatall,1);
else
tstatmu=tinit;
end

%% get tnull distribution
data=fit.tgt(trIdx);n_var=fit.quadsize;
hbar=parfor_progressbar(nboots,'Approximating a null distribution');
    %Use random permutations
    tnull=zeros(nboots,n_var);
   parfor perm=1:nboots
        hbar.iterate(1)
        tgt_perm=shuffle(data);
        fitstatsnull=fitlm([fit.designquad(trIdx,:) fit.design{trIdx,end-1}],tgt_perm,'Intercept',false);
        tstattmp=fitstatsnull.Coefficients.tStat;
        for i=1:n_var
        tnull(perm,i)=tstattmp(i);
        end
    end
    close(hbar)
    nans= isnan(tnull);
    tnull(nans)=0;
        for i=1:n_var    
        pval(i)=(sum(tnull(:,i)>abs(tstatmu(i)))/nboots + sum(tnull(:,i)< -1*abs(tstatmu(i)))/nboots);
        end
    
end