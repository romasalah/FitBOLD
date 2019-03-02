function [p,tbl]=anova_fitBOLD(folder,models)
if folder==0
    folder=pwd;
end
conds={'active','passive','n-soc'};
mats={'*results.mat','*Best_features_stats.mat'};
crit1={'IS','RL'};
crit2=conds;
cd(folder);
for whichres=2:size(mats,2)
    clear results
    mtt=mats{whichres};
    query= ['*' models mtt];
    files= dir(query);
    for load_i=1: size(files,1)
        try results{1,load_i}=load(files(load_i).name);
            results{2,load_i}=files(load_i).name;
        catch failed_loads; 
        end
    end
    loaded_i=find(~cellfun(@(x) isempty(fields(x)),{results{1,1:end}}));
    results=results(:,loaded_i);
    files=files(loaded_i);
    names=erase_fitBOLD(files,models,mtt(2:end),conds);
    %names and order
lbl={names.name};
for conds_i=1:size(conds,2)
    conds{2,conds_i}=find(~cellfun(@isempty,strfind(lbl,conds{1,conds_i})));
end
lbl_ord=[conds{2,1:end}];
    if strcmp(mtt,'*results.mat')
        for m_i=1:size(results,2)
            res{m_i}=   results{1,m_i}.results;
            summary(m_i,1)=res{m_i}.summaryci.adjrsqr(1);
            summary(m_i,2)=res{m_i}.summaryci.RMSE(1);
            summary(m_i,3)=res{m_i}.summaryci.AIC(1);
            summary(m_i,4)=res{m_i}.summaryci.rsqr(1);
            summary(m_i,5)=res{m_i}.summaryci.BIC(1);
            summary(m_i,6)=res{m_i}.summaryci.coefnum(1);
            fitstats{m_i,1}=res{m_i}.summaryci.fstat;
            fitstats{m_i,2}=res{m_i}.summaryci.pval(1);
        end
        allcrit={crit1,crit2};
        
        for critall_i=1:size(allcrit,2)
            crit=allcrit{critall_i};
        for fcrit1=1:size(crit,2)
            crit1_i{fcrit1}=find(~cellfun(@isempty,strfind({results{2,1:end}},crit{fcrit1}))==1);
            for bld_gp=1:size(crit1_i{fcrit1},2)
                ctgp{critall_i}{crit1_i{fcrit1}(bld_gp)}=crit{fcrit1};
            end
        end
        end
        
       for crit_i=1:size(summary,2)
           [p{crit_i},tbl{crit_i},stats{crit_i},terms{crit_i}]=anovan(summary(:,crit_i),{ctgp{1},ctgp{2}});
        end
        
    elseif strcmp(mtt,'*Best_features_stats.mat')
        for fcrit1=1:size(crit1,2)
        crit1_i{fcrit1}=find(~cellfun(@isempty,strfind({results{2,1:end}},crit1{fcrit1}))==1);
        end
        for fcrit2=1:size(crit2,2)
        crit2_i{fcrit2}=find(~cellfun(@isempty,strfind({results{2,1:end}},crit2{fcrit2}))==1);
        end
        for m_i=1:size(results,2)
            alldist=(results{1,m_i}(1:end).fitstats);
            for iboot=1:size(alldist,2)
                dist{m_i,1}(iboot)=[alldist(iboot).fitstats.Rsquared.adjusted];
                dist{m_i,2}(iboot)=alldist(iboot).fitstats.RMSE;
                dist{m_i,3}(iboot)=alldist(iboot).fitstats.ModelCriterion.AIC;
                dist{m_i,4}(iboot)=alldist(iboot).fitstats.Rsquared.Ordinary;
                dist{m_i,5}(iboot)=alldist(iboot).fitstats.LogLikelihood;
                dist{m_i,6}(iboot)=alldist(iboot).fitstats.ModelCriterion.BIC;
                dist{m_i,7}(iboot)=alldist(iboot).fitstats.NumPredictors;
                dist{m_i,8}(iboot)=alldist(iboot).fitstats.MSE;
            end
            
            for crit_i=1:size(dist,2)
                med(m_i,crit_i)=median(dist{m_i,crit_i});
                mu(m_i,crit_i)=mean(dist{m_i,crit_i});
                distfit{m_i,crit_i}=allfitdist(dist{m_i,1}','BIC');
                distname{m_i,crit_i}=distfit{m_i,crit_i}(1).DistName;
            end
        end
        clear ctgp
        allcrit={crit1,crit2};
        for critall_i=1:size(allcrit,2)
            crit=allcrit{critall_i};
        for fcrit1=1:size(crit,2)
            crit0_i{fcrit1}=find(~cellfun(@isempty,strfind({results{2,1:end}},crit{fcrit1}))==1);
            for bld_gp=1:size(crit0_i{fcrit1},2)
                ctgp{critall_i}{crit0_i{fcrit1}(bld_gp)}=crit{fcrit1};
            end
        end
        end
      
        factors={'Theory','Condition'};
        clear anovastats factsmu factsmed
         for crit_i=1:3
             for fact1_i=1:size(crit1_i,2)
                 for subf_i=1:size(crit1_i{fact1_i},2)
                  allfact(:,subf_i)=dist{(crit1_i{fact1_i}(subf_i)),crit_i};
                 end
                 factsall(:,:,fact1_i)=allfact;
                 mufact1=mu((crit1_i{fact1_i}),crit_i);
                 medfact1=med((crit1_i{fact1_i}),crit_i);
                 factsmu{crit_i}(:,fact1_i)=mufact1;
                 factsmed{crit_i}(:,fact1_i)=medfact1;
             end
             [out,str]=anova2_repmeas_autodesign(factsall,{'Theory','Condition'},1);
             [p.mu{crit_i},tbl.mu{crit_i}]=anovan_js(factsmu{crit_i},'residuals','factornames',{'Theory','Condition'},'display','on');
             [p.med{crit_i},tbl.med{crit_i}]=anovan_js(factsmed{crit_i},'residuals','factornames',{'Theory','Condition'},'display','on');
             methods={tbl.mu{crit_i},tbl.med{crit_i}};
           
             for fact_i=1:size(tbl.med{crit_i},2)
                 for method_i=1:2
                     meth=methods{method_i};
                     indci=strfind(meth{1,fact_i},'=');
                     sepi=strfind(meth{1,fact_i},', ');
                     ends=[sepi-1,size(meth{1,fact_i},2)];
                     for indcii=1:size(indci,2)
                     val=str2num(meth{1,fact_i}(indci(indcii)+1:ends(indcii)));
                      anovastats{crit_i}{indcii}(method_i,fact_i)=val;
                     end
                 end
             end
         end
         
         indclbl={'F-statitic', 'p-Value'};
         lbls={'Adjusted R-Squared','RMSE','AIC'};
         
          numcrit=size(tbl.med,2);
         for indcii2=1:size(indclbl,2)
             figure;
         for crit_i2=1:numcrit
         subplot(numcrit,1,crit_i2)
         bar(anovastats{crit_i2}{1,indcii2});
         xlabel_oblique(factors)
         title([lbls{crit_i2} ' ' indclbl{indcii2}])
         if indcii2==2
             hold on;
             xpval=0:size(factors,2)+1;
             plot(xpval,zeros(1,size(xpval,2))+0.05,'r')
         end
         end
         end
        
end


end