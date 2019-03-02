function [FinalDesign,fitcrit,nloops,niter,namesquad,ROIs_i,termsquad]=parstepwiselm(design,power,Incriterion,outcriterion,upper,workers,fit)
%% get terms and response
ROIs_i=0;
warning off
catname={'Gender'}; tgtname={'Happiness ratings'};

%% Create polynomials
[termsquad,namesquad,contterms,names,tgt]=Polynom_fitBOLD(design,power);
catterms=design{:,end-1};
%% do first fit
initfit=fitlm(design,'linear','intercept',false);
switch Incriterion
    case 'aic'
        critold=initfit.ModelCriterion.AIC;
    case 'bic'
        critold=initfit.ModelCriterion.BIC;
    case 'adjrsquared'
        critold=initfit.Rsquared.Adjusted;
    case 'rsquared'
        critold=initfit.Rsquared.Ordinary;
    case 'fstat'
        initANOVA=anova(initfit,'summary');
        critold=initANOVA.F(2);
end

comb_num=size(contterms,2);
%% Claculate number of all combinattions
c=0;
for n=1:comb_num
    k=n-1;
    c=c+nchoosek(n,k);
end

%% start parallel fitting
nloops=0;
niter=0;
lw_vars=1:comb_num;
oldsamplsize=comb_num+1;
lastvalid_idx=comb_num;
Rednames={namesquad{1,1:size(names,2)}};
gcf; title(['Improvement in ' Incriterion])
xlabel('Model Size iterations'); ylabel(Incriterion); hold on;
for samplesize= comb_num-1:-1:1 %loop through all possible k(sample sizes);
    npos=samplesize+1;
    niter=niter+nchoosek(npos,samplesize);
    plot(nloops,critold,'ro'); hold on; xlabel(['Model Size iterations: Current size= ' num2str(size(Rednames,2)) '   Iteration No. ' num2str(niter) '.']);
    pause(0.05);
    %% choose the best power
    %here the construction of polyterms is done by the model itself.
    clear newvarssinglesidx Powerdlwvarsidx choosepowerfit ...
        tstat powernames quadnamesidx doubleidx ...
        newvarsidx2 newpowernames Rednamessingles Comb  newvarsidx ...
        newvarsnames  idxnewquadraticnames newvarsquad newvarsquadhalf
    if nloops==0
        Powerdlwvarsidx=1:comb_num;
    else
        %Rednamessingles=cellfun(@(x) replace(x,'^2',''),Rednames,'un',0);
        for getidx=1:size(Rednames,2)
            idxtmp=find(strcmp(namesquad, Rednames{1,getidx})==1);
            if ~isempty(idxtmp)
                newvarsidx(getidx)=idxtmp;
            end
        end
        switch power
            case 1
                Powerdlwvarsidx=newvarsidx;
            case 2
                if nloops==1
                    newvarsnames=namesquad(:,newvarsidx);
                    % create names to do the quadratic regression
                    for getquad = 1: size(newvarsnames,2)
                        if ~isempty(strfind(newvarsnames{getquad},'^2'))
                            newvarsquadhalf{getquad}= newvarsnames{getquad}(1:end-2);
                        else
                            newvarsquadhalf{getquad}= [newvarsnames{getquad} '^2'];
                        end
                    end
                    newvarsquad=horzcat(newvarsnames,newvarsquadhalf);
                    %find indices of those names in the original data
                    for getidx3=1:size(newvarsquad,2)
                        idxtmp=find(strcmp(namesquad, newvarsquad{1,getidx3})==1);
                        if ~isempty(idxtmp)
                            idxnewquadraticnames(getidx3)=idxtmp; clear idxtmp;
                        end
                    end
                    idxnewquadraticnames=sort(idxnewquadraticnames);
                    %% make the quadratic fit
                    choosepowerfit=fitlm(termsquad(:,idxnewquadraticnames),tgt,'linear','VarNames',...
                        horzcat(namesquad{:,idxnewquadraticnames},tgtname),'intercept',false); %we don't add categorical vars here
                    tstat=choosepowerfit.Coefficients.tStat;tstat(isnan(tstat))=0;
                    powernames=choosepowerfit.Coefficients.Properties.RowNames;
                    % quadnamesidx=find(~cellfun(@isempty,strfind(powernames,'^2'))==1);
                    doubleidx=size(powernames,1)/2;
                    for getpower=1:doubleidx
                        if abs(tstat(getpower))>=abs(tstat(getpower+doubleidx))
                            newpowernames{getpower,1}=powernames{getpower,1};
                        else
                            newpowernames{getpower,1}=powernames{getpower+doubleidx,1};
                        end
                    end
                    for getidx2=1:size(newpowernames,1)
                        idxtmp2=find(strcmp(namesquad, newpowernames{getidx2,1})==1);
                        if ~isempty(idxtmp2)
                            newvarsidx2(getidx2)=idxtmp2;
                        end
                    end
                    Powerdlwvarsidx=newvarsidx2;
                else
                    Powerdlwvarsidx=newvarsidx;
                end
                %             else
                %             Powerdlwvarsidx=newvarssinglesidx;
                %             end
        end
    end
    
    %% Calculate Combinations
    V= Powerdlwvarsidx; %all possible variables
    nvars=length(V);
    Comb=nchoosek(V,length(V)-1); %create the logic matrix for combinations
    fitcrit=zeros(nvars,1);
    clear fitlist fitcrit
    %% do Regression
    %hbar = parfor_progressbar(nvars,['Trying model size of ' num2str(samplesize) ' Variables: Progress...']); %create the progress bar
    switch Incriterion
        case 'aic'
            parfor (cr_i=1:nvars,workers*2)
                %for cr_i=1:nvars
                warning off
                fitlist{cr_i}=fitlm(termsquad(:,Comb(cr_i,:)),tgt,'linear','intercept',false); %construct and fit the temporary designs for each samplesize
                fitcrit(cr_i)=fitlist{cr_i}.ModelCriterion.AIC;
                fitlist{cr_i}={};
                %hbar.iterate(1);
            end
        case 'bic'
            parfor (cr_i=1:nvars,workers*2)
                %for cr_i=1:nvars
                warning off
                fitlist{cr_i}=fitlm(termsquad(:,Comb(cr_i,:)),tgt,'linear','intercept',false); %construct and fit the temporary designs for each samplesize
                fitcrit(cr_i)=fitlist{cr_i}.ModelCriterion.BIC;
                fitlist{cr_i}={};
                %hbar.iterate(1);
            end
            
        case 'fstat'
            parfor (cr_i=1:nvars,workers*2)
                %for cr_i=1:nvars
                warning off
                fitlist{cr_i}=fitlm(termsquad(:,Comb(cr_i,:)),tgt,'linear','intercept',false); %construct and fit the temporary designs for each samplesize
                ANOVAlist{cr_i}=anova(fitlist{cr_i},'summary');
                fitcrit(cr_i)=ANOVAlist{cr_i}.F(2);
                fitlist{cr_i}={};ANOVAlist{cr_i}={};
                %hbar.iterate(1);
            end
        case 'adjrsquared'
            parfor (cr_i=1:nvars,workers*2)
                %for cr_i=1:nvars
                warning off
                fitlist{cr_i}=fitlm(termsquad(:,Comb(cr_i,:)),tgt,'linear','intercept',false); %construct and fit the temporary designs for each samplesize
                fitcrit(cr_i)=fitlist{cr_i}.Rsquared.Adjusted;
                fitlist{cr_i}={};
                %hbar.iterate(1);
            end
        case 'rsquared'
            parfor (cr_i=1:nvars,workers*2)
                %for cr_i=1:nvars
                warning off
                fitlist{cr_i}=fitlm(termsquad(:,Comb(cr_i,:)),tgt,'linear','intercept',false); %construct and fit the temporary designs for each samplesize
                fitcrit(cr_i)=fitlist{cr_i}.Rsquared.Ordinary;
                fitlist{cr_i}={};
                %hbar.iterate(1);
            end
    end
    %close(hbar);
    %% find best model
    if fit.regress.failed>=3; fit.scaleoldcrit=0.02+fit.scaleoldcrit; else fit.scaleoldcrit=0;  end
    
    if strcmp(Incriterion,'adjrsquared') || strcmp(Incriterion,'rsquared') || strcmp(Incriterion,'fstat')
        crit=max(fitcrit);
        if crit>critold-abs(fit.scaleoldcrit*critold)
            improvement=1;
        else
            improvement=0;
        end
        crit_best_i=find(fitcrit==max(fitcrit));
    else
        crit=min(fitcrit);
        if crit<critold+abs(fit.scaleoldcrit*critold)
            improvement=1;
        else
            improvement=0;
        end
        crit_best_i=find(fitcrit==min(fitcrit));
    end
        lw_vars=Comb(crit_best_i(1,end),:); %get what variables in the localwinner
        critold=crit;
        nloops=nloops+1;
        Rednames={namesquad{1,lw_vars}};
        Rednames=horzcat(Rednames,horzcat(catname,tgtname));
        Redterms=[termsquad(:,lw_vars) catterms];
        pval=getpval(Redterms,Rednames,tgt,upper,outcriterion);

    if improvement==0
                if pval< fit.regress.threshold/size(fit.design,2)
        [FinalDesign,ROIs_i,~]=designfinal(Redterms,Rednames,tgt,upper,outcriterion,namesquad);
                return
        end
        if  mod(fit.regress.failed,3)==0
            figure;
            Incriterion='adjrsquared';
            fit.regress.failed=fit.regress.failed+1;
           [FinalDesign,fitcrit,nloops,niter,namesquad,ROIs_i,termsquad]=parstepwiselm(design,power,Incriterion,outcriterion,upper,workers,fit);
        elseif fit.regress.failed==1 || mod(fit.regress.failed,3)==1
            figure;
            fit.regress.failed=fit.regress.failed+1;
            newcrit=fit.regress.incriterion2;
            [FinalDesign,fitcrit,nloops,niter,namesquad,ROIs_i,~]=change_criterion(design,power,newcrit,outcriterion,upper,workers,fit);
         elseif fit.regress.failed==2 || mod(fit.regress.failed,3)==2
             figure;
             fit.regress.failed=fit.regress.failed+1;
             newcrit='fstat';
             [FinalDesign,fitcrit,nloops,niter,namesquad,ROIs_i,~]=change_criterion(design,power,newcrit,outcriterion,upper,workers,fit);
        end
    elseif improvement==1 && size(lw_vars,2)<5
        [FinalDesign,ROIs_i,~]=designfinal(Redterms,Rednames,tgt,upper,outcriterion,namesquad);
        return
    end
end
end
function [FinalDesign,fitcrit,nloops,niter,namesquad,ROIs_i,pval]=change_criterion(design,power,newcrit,outcriterion,upper,workers,fit)
Incriterion=newcrit;
[FinalDesign,fitcrit,nloops,niter,namesquad,ROIs_i]=parstepwiselm(design,power,Incriterion,outcriterion,upper,workers,fit);
finalanova=anova(FinalDesign,'summary');
pval=finalanova.pValue(2);
end
function pval=getpval(Red2terms,Red2names,tgt,upper,outcriterion)

if strcmp(upper,'interactions') || isa(outcriterion,'char')
    Red2Design=stepwiselm(Red2terms,tgt,'linear','intercept',false,'VarNames',Red2names,'criterion',outcriterion,'upper','linear','verbose',0);
else
    Red2Design=fitlm(Red2terms,tgt,'linear','intercept',false,'VarNames',Red2names); %take the last design as the winner
end
finalanova=anova(Red2Design,'summary');
pval=finalanova.pValue(2);
end
function [FinalDesign,ROIs_i,pval]=designfinal(Redterms,Rednames,tgt,upper,outcriterion,namesquad)
        if strcmp(upper,'interactions') || isa(outcriterion,'char')
            RedDesign=stepwiselm(Redterms,tgt,'linear','intercept',false,'VarNames',Rednames,'criterion',outcriterion,'upper','linear','verbose',0);
            final_lw_vars=find(RedDesign.Formula.InModel==1);
            Finnames={Rednames{1,final_lw_vars}};Finnames=horzcat(Finnames,tgtname);
            FinalDesign=fitlm(Redterms(:,final_lw_vars),tgt,upper,'intercept',false,'VarNames',Finnames); %take the last design as the winner
        else
            Finnames=Rednames;
            FinalDesign=fitlm(Redterms,tgt,'linear','intercept',false,'VarNames',Finnames); %take the last design as the winner
            finalanova=anova(FinalDesign,'summary');
            pval=finalanova.pValue(2);
            Finalnames=FinalDesign.Coefficients.Properties.RowNames;
            for des_loop=1:size(Finalnames,1)
                for power_loop=1:size(namesquad,2)
                    if strcmp(Finalnames{des_loop,1},namesquad{1,power_loop})
                        ROIs_i(des_loop)=power_loop;
                    end
                end
            end
        end
end