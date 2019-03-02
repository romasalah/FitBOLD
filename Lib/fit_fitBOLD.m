function [result,tgt,fitstats] = fit_fitBOLD(fit)
terms=fit.terms; tgt=fit.tgt; sub_coord=fit.sub_coord;
fit.designquad=Polynom_fitBOLD(fit.design,2);
%% *****Do Stepwise fit******
if strcmp(fit.dependence,'stepwise') || strcmp(fit.dependence,'OLS')
    fprintf(['\n Dependent explanatory variables: Using ' fit.dependence ' regression \n'])
    clear fitstats betas inmodel pval se
    switch fit.dependence
        case 'stepwise'
            if fit.regress.validate
            S=(1:size(fit.sub,2))';
            CVO=cvpartition(S,'k',fit.regress.validationfolds);
            fitnum=CVO.NumTestSets;
            else
            fitnum=1;
            end
            for sample_i= 1:fitnum
                if fit.regress.validate
                trsub = CVO.training(sample_i);trIdx=[fit.S_i{find(trsub==1),1}]';
                tesub = CVO.test(sample_i);teIdx=[fit.S_i{find(tesub==1),1}]';
                termstest=fit.design(teIdx ,:); tgttest{sample_i}=fit.design{teIdx,end};
                else
                    trsub=1:size(fit.sub,2);
                trIdx=1:size(fit.design,1);
                end
                termstrain=fit.design(trIdx,:);
                
                if fit.pooled==1 modelspec='linear'; Incriterion='adjrsquared'; upper='linear'; Outcriterion=0; power=2; else modelspec='linear'; criterion='sse'; end
                strpOpts.power=power;
                strpOpts.nboots=1000;
                strpOpts.Incriterion=Incriterion;
                strpOpts.Outcriterion=Outcriterion;
                strpOpts.upper=upper; strpOpts.fit=fit;
                if fit.cpu_parallel>0
                    stream = RandStream('mlfg6331_64');  % Random number stream
                    options = statset('UseParallel',1,'UseSubstreams',1,...
                        'Streams',stream);
                else options='default';
                end
                tic
                figure;
                %% Select features;
                [fitstats0,~,~,~,coefnames,idxattrainquadterms,termstrainquad]=parstepwiselm(termstrain,strpOpts.power,strpOpts.Incriterion,strpOpts.Outcriterion,strpOpts.upper,strpOpts.fit.cpu_parallel,strpOpts.fit);
                %termsquad=termsquad(trIdx,:);
                fit.regress.feature_selection=1;
                figure;
                fit.featurestgttbl=array2table([termstrainquad(:,idxattrainquadterms) termstrainquad(:,end-1) fit.design{trIdx,end}]);
                %% estimate CI
                switch fit.resample.method
                    case 'cv'
                        strpOpts.kfolds=10;
                        [summary,fullstats]=cv_parstepwiselm(strpOpts,fit);
                    case 'bootstrap'
                        [summary,fullstats]=strp_parstepwiselm(strpOpts,fit,idxattrainquadterms,trsub);
                        yhat=zeros(size(fullstats,2),size(termstrain,1));
                        pvalnull=strppval_fitBOLD(fullstats,fitstats0,idxattrainquadterms,trIdx,fit);
                end
                elapsed=toc;
                fprintf(['\n Elapsed time:' num2str(elapsed) '\n Number of iterations: ' num2str(sum([fullstats.niter],2)) '\n regressions, through ' num2str(sum([fullstats(:).nloops],2)) ' design sizes \n']);
                %% set reference structure for indices
                if fit.regress.feature_selection==1; fitstats_ref.fitstats=fitstats0; else; fitstats_ref=fullstats;end
                %% remove duplicate subjects
                catbetas=zeros(strpOpts.nboots,1);
                betastrap=zeros(1,fit.quadsize);yhat=zeros(strpOpts.nboots,size(termstrain,1));
                signames={};allbetasstrap=zeros(1,fit.quadsize);
                for get_ci=1:size(fullstats,2)
                    [yhat_un_idx{get_ci},yhat_un_idxinfitstats{get_ci}]=unique(fullstats(get_ci).yhat_i_sample);
                    %% get all and significant betas and their names
                    estimates=fullstats(get_ci).fitstats.Coefficients.Estimate;
                    coefnum=size(fullstats(get_ci).fitstats.Coefficients,1)-1;
                    for fill_betas=1:coefnum %separate the Categorical Variable;
                        allbetasstrap(get_ci,fullstats(get_ci).idxfrominput(fill_betas))=estimates(fill_betas);
                        if  strcmpi(fit.resample.method,'bootstrap') && pvalnull(fill_betas)< (fit.regress.threshold/(coefnum+1)) && estimates(fill_betas)~=0
                            betastrap(get_ci,fullstats(get_ci).idxfrominput(fill_betas))=estimates(fill_betas);
                        end
                    end
                    thismodelnames=fullstats(get_ci).fitstats.CoefficientNames;
                    if ~isempty(strfind(thismodelnames,'Gender')) && pvalnull(end)<=(0.05/(coefnum+1))
                        catbetas(get_ci,1)=fullstats(get_ci).fitstats.Coefficients{end,1};
                    end
                    yhat(get_ci,yhat_un_idx{get_ci})= (fullstats(get_ci).fitstats.Fitted(yhat_un_idxinfitstats{get_ci},1))';
                    if strcmpi(fit.resample.method,'cv')
                        predrsqr(get_ci,1)=fullstats(get_ci).predrsqr;
                    end
                end
                
                figure; subplot(3,1,1); imagesc(allbetasstrap);colorbar;
                subplot(3,1,2); imagesc(betastrap);colorbar;
                subplot(3,1,3); imagesc(mean(betastrap,1));colorbar;
                %% Calculate means and CI
                catdist=fitdist(catbetas,'Normal');
                catci=ci(catbetas);
                if strcmpi(fit.resample.method,'cv')
                    predrsqrdist=fitdist(predrsqr,'Normal');
                    predrsqrci=paramci(predrsqrdist,'Parameter','mu');
                    result.predrsqrci=predrsqrci;
                end
                catcitbl=array2table([mean(catbetas);catci]);
                yhatmean=mean(yhat,1);
                cirownames={'mean','lower','upper'};
                for yhat_ci=1:size(yhat,2)
                    yhatci(:,yhat_ci)=ci(yhat(:,yhat_ci));
                end
                parfor i=1:size(betastrap,2)
                    coefci(:,i)=ci(betastrap(:,i));
                end
                parfor i=1:size(summary,2)
                    sumci(:,i)=ci(summary{:,i});
                end
                
                sumcitbl=array2table([mean(table2array(summary),1); sumci],'VariableNames',summary.Properties.VariableNames,...
                    'RowNames',cirownames);
                coefcitbl=array2table([mean(betastrap,1);coefci],'RowNames',cirownames);
                
                result.allbetasstrap=allbetasstrap;
                result.BOLDpred=mean(yhat,1)';
                result.summaryci=sumcitbl;
                result.catci=catcitbl;
                result.yhatci=yhatci;
                result.model_tstat=mean([summary.fstat(:)],2);
                result.model_pval=mean([summary.pval(:)],2);
                result.ROIs=coefnames;
                result.predictorROIs=coefnames(sum(betastrap(:,1:end-1),1)~=0);
                result.coefci=coefcitbl(:,(sum(betastrap,1)~=0));
                linearityidx=zeros(2,size(result.predictorROIs,2));
                for linidx=1:size(result.predictorROIs,2)
                    if ~isempty(strfind(result.predictorROIs{1,linidx},'^2'))
                        linearityidx(2,linidx)=1;
                    else
                        linearityidx(1,linidx)=1;
                    end
                end
                result.model_in=1;
                result.linearityidx=linearityidx;
                if fit.use_subROIs==1
                    coords={};
                    for get_coords=1:size(sub_coord,2)
                        coords_rois=horzcat(coords,sub_coord{get_coords}{1,1:end});
                    end
                    result.coords_rois=coords_rois;
                    result.sub_coord=sub_coord;
                end
                if fitnum>1
                % test predictors
                meanbetas=mean(betastrap,1)';
                testterms=[fit.designquad(teIdx,:) fit.tgt(teIdx)];
               predictedtgt_CV{1,sample_i}=testterms*meanbetas;
               validation_rsquared(sample_i)=corr(tgttest{1,sample_i},predictedtgt_CV{1,sample_i}).^2;
               validation_RMSE(sample_i)=sqrt(mean((predictedtgt_CV{1,sample_i}-tgttest{1,sample_i}).^2));
                end
            end
        case 'OLS'
            fitstats=fitglm(terms,tgt,'Intercept',false);
            result.BOLDpred=fitstats.Fitted.Response;
            result.ROIs=(fitstats.Coefficients.Properties.RowNames)';
            result.model_in=1;
            
        end  
            fitstats=fullstats;
            result.rsqr=sumcitbl.rsqr(1);
            if fitnum>1
            result.validationR2=validation_rsquared;
            result.validationRMSE=validation_RMSE;
            result.testedresponse=tgttest;
            result.predictedresponse=predictedtgt_CV;
            end
            result.betas=coefcitbl{1,:};
            result.bic=sumcitbl.BIC(1);
            result.aic=sumcitbl.AIC(1);
            result.nlag=0;
            result.sub=fit.sub;
            result.y=tgt;
            
            
    
    
    %%    %********Do Ridge regression*******
elseif strcmp(fit.dependence,'ridge')
    fprintf('\n Dependent explanatory variables: Using Ridge regression \n')
    terms2=x2fx(terms,'interaction');
    terms2(:,1)=[];
    k = 0:10e-7:1;
    zsctgt=zscore(tgt);
    betahat = ridge(zsctgt,terms2,k);
    %Plot Ridge regression
    figure
    plot(k,betahat,'LineWidth',2)
    ylim([-100 100])
    grid on
    xlabel('Ridge Parameter')
    ylabel('Standardized Coefficient')
    title('{\bf Ridge Trace}')
    legend('x1','x2','x3','x1x2','x1x3','x2x3')
    %Calculate R-squared
    betas=betahat(:,end);
    %betas=betas';
    
    BOLDpred_ridge=nansum(terms2*betas(1:end,1),2)/size(terms2,2);
    lse_ridge = nansum((tgt-BOLDpred_ridge).^2); %sum least-squares error (SSR)
    BOLDr2_ridge = 1-(lse_ridge/re);
    r2=r2_ridge;
    b=betas;
    
    %%    %********DO Fused lasso followed by OLSAR(n)*********
elseif (strfind(fit.dependence,'fused'))
    fprintf('\n Dependent explanatory variables: \n Using Fused lasso followed by Autoregressive OLS regression \n')
    cd '/Applications/spm12/toolbox/spams-matlab-v2.6'
    start_spams
    addpath(genpath('/Applications/spm12/toolbox/jplv7'))
    fprintf('\n Dependent explanatory variables & autocorrelated observations \n Using Fused lasso regression \n')
    inx=inx';
    if fit.zscored==0
        max_lambda=10000;
        nvar=size(terms,2);
        nobs=size(terms,1);
        lambda_1=nobs:max_lambda/29:max_lambda;
        lambda_2=nobs:max_lambda/29:max_lambda;
        lambda_3=nobs:max_lambda/29:max_lambda;
    else
        
        nvar=size(terms,2);
        nobs=size(terms,1);
        max_lambda= nobs/nvar;
        lambda_1=(max_lambda*sqrt((nvar/nobs)))/29:(max_lambda*sqrt((nvar/nobs)))/29:(max_lambda*sqrt((nvar/nobs)));
        lambda_2=(max_lambda*sqrt((nvar/nobs)))/29:(max_lambda*sqrt((nvar/nobs)))/29:(max_lambda*sqrt((nvar/nobs)));
        lambda_3=(max_lambda*sqrt((nvar/nobs)))/29:(max_lambda*sqrt((nvar/nobs)))/29:(max_lambda*sqrt((nvar/nobs)));
        
    end
    lm=zeros(0,3);
    for lambda_i1=1:size(lambda_1,2)
        % fprintf(['\n Trying Fused lasso Lambda = ' num2str(lambda_1(lambda_i1))])
        for lambda_i2=1:size(lambda_2,2)
            for lambda_i3=1:size(lambda_3,2)
                lm(end+1,1)=lambda_1(lambda_i1);
                lm(end,2)=lambda_2(lambda_i2);
                lm(end,3)=lambda_3(lambda_i3);
            end
        end
    end
    fprintf('\nFISTA + Regression Fused-Lasso\n');
    %% %**********Iterate over lambdas to get best BIC**************
    %     %start with amore agressive regularization
    all_const_idx=size(all_corct_reg,2)+1:size(terms,2);
    S=(1:size(fit.sub,2))';
    if fit.pooled==1 folds=round(size(fit.sub,2)/5)+1; else S=[S;1]; folds=2; end
    CVS = cvpartition(S,'k',folds);
    for i = 1:CVS.NumTestSets
        if fit.pooled==1
            trsub = CVS.training(i);trIdx=[S_i{find(trsub==1),1}]';
            tesub = CVS.test(i);teIdx=[S_i{find(tesub==1),1}]';
            terms_train=terms(trIdx,:);tgt_train=tgt(trIdx,1);
            terms_test=terms(teIdx ,:);tgt_test=tgt(teIdx,1);
        else
            terms_train=terms;trIdx=(1:size(tgt,1))';
            tgt_train=tgt;trsub=1;
        end
        j=size(lm,1);
        param=cell(1,j);out=cell(1,j);fitstats_lm=cell(1,j);lambda_BIC=cell(j,1);
        lso=cell(1,j);terms_train_const=cell(1,j); highest_likelihood_i=cell(1,j);
        terms_train_in =cell(1,j);out2=cell(1,j);sig_lag=cell(1,j);nlag=cell(1,j); nonempt3=cell(1,j);
        g=cell(1,j);
        if fit.pooled==0 && i>1
            break
        else
            parfor (parlm=1:size(lm,1),fit.cpu_parallel)
                %fprintf(['\n Trying Fused lasso Lambda = ' num2str(lambda_1(lambda_i1))])
                param{parlm}.loss='square';
                param{parlm}.numThreads=-1; % all cores (-1 by default)
                param{parlm}.verbose=true;   % verbosity, false by default
                param{parlm}.it0=10;      % frequency for duality gap computations
                param{parlm}.max_it=200; % maximum number of iterations
                param{parlm}.L0=0.1;
                param{parlm}.tol=1e-3;
                param{parlm}.intercept=false;
                param{parlm}.pos=false;
                param{parlm}.ista=false;
                param{parlm}.regul='fused-lasso';
                lm1{parlm}=lm(parlm,1);
                lm2{parlm}=lm(parlm,2);
                lm3{parlm}=lm(parlm,3);
                param{parlm}.lambda=lm1{parlm};
                param{parlm}.lambda2=lm2{parlm};
                param{parlm}.lambda3=lm3{parlm};
                %********get features****
                [lso{parlm},~]=mexFistaFlat(tgt_train,terms_train(:,1:size(all_corct_reg,2)),inx(1:size(all_corct_reg,2),1),param{parlm});
                terms_train_in{parlm} = terms_train(:,lso{parlm}(:,1)~=0);
                %get**** the optimal time lag only in non-pooled modelling****
                minlag=1;maxlag=10;min_dof=5; nlag{parlm}=0;
                if fit.AR==0
                    if fit.zscored==0
                        terms_train_const{parlm}=horzcat(terms_train_in{parlm},terms(trIdx,all_const_idx(find(trsub==1))));
                        fitstats_lm{parlm}=fitglm(terms_train_const{parlm},tgt_train);
                    else
                        fitstats_lm{parlm}=fitglm(terms_train_in{parlm},tgt_train);
                    end
                    g{parlm}.rsqr=fitstats_lm{parlm}.Rsquared.Ordinary;
                else
                    while ((size(tgt,1)-maxlag)-(maxlag+size(terms_train_in{parlm},2)+1))<min_dof %have at least dof of 2 when testing (number of observations must be at least 2 more than number of explanatory variables)
                        maxlag=maxlag-1;
                    end
                    if maxlag<=minlag
                        maxlag=minlag;
                    end
                    out{parlm}=lrratio(tgt_train,maxlag,minlag,1,terms_train_in{parlm});
                    %get best lag index
                    nonempt3{parlm}=find(~isnan(out{parlm}(:,4)));
                    out2{parlm}=out{parlm}(nonempt3{parlm},:);
                    sig_lag{parlm}=out2{parlm}(find(real(out2{parlm}(:,4))<0.001),:);
                    if isempty(sig_lag{parlm})
                        sig_lag{parlm}=out2{parlm}(find(real(out2{parlm}(:,4))<0.05),:);
                        if isempty(sig_lag{parlm})
                            sig_lag{parlm}=out2{parlm}(find((real(out2{parlm}(:,4))==min(real(out2{parlm}(:,4))))),:);
                        end
                    end
                    highest_likelihood_i{parlm}=find(abs(sig_lag{parlm}(:,3))==max(abs(sig_lag{parlm}(:,3))));
                    nlag{parlm}=sig_lag{parlm}(highest_likelihood_i{parlm},1);
                    
                    
                    if size(tgt_train,1)~=size(terms_train_in{parlm},2)
                        try
                            %get Pramater estimates with the optimal time lag
                            if fit.zscored==0
                                terms_train_const{parlm}=horzcat(terms_train_in{parlm},terms(trIdx,all_const_idx(find(trsub==1))));
                                fitstats_lm{parlm}=vare(tgt_train,terms_train_const{parlm});
                            else
                                fitstats_lm{parlm}=vare(tgt_train,nlag{parlm},terms_train_in{parlm});
                            end
                            g{parlm}.rsqr=fitstats_lm{parlm}.rsqr;
                            %some specific lags cause an error in calculating f-statistic and fprob. we will ignore these nlags
                        end
                    end
                end
                lambda_BIC{parlm,1}=[param{parlm}.lambda,param{parlm}.lambda2,param{parlm}.lambda3, g{parlm}.rsqr,size(terms_train_in{parlm},2)]
            end
            
            lmbc=[lambda_BIC{:,1}];
            lmbc_features=lmbc(1,5:5:size(lmbc,2));
            empt_features=[find(lmbc_features==0) find(lmbc_features==nvar)];
            lmbc_r2=lmbc(1,4:5:size(lmbc,2));
            if ~isempty(empt_features)
                lmbc_r2(empt_features)=0;
            end
            max_r2=find(lmbc_r2==max(lmbc_r2));
            max_r2_i=max_r2(1,end);
            best_lambda_1= lambda_BIC{max_r2_i,1}(1,1);
            best_lambda_2= lambda_BIC{max_r2_i,1}(1,2);
            best_lambda_3= lambda_BIC{max_r2_i,1}(1,3);
            param2.loss='square'; param2.numThreads=-1; % all cores (-1 by default)
            param2.verbose=true; param2.max_it=200; param2.L0=0.1;
            param2.tol=1e-3; param2.intercept=false;
            param2.pos=false; param2.ista=false; param2.regul='fused-lasso';
            param2.lambda=best_lambda_1;  param2.it0=10;
            param2.lambda2=best_lambda_2;
            param2.lambda3=best_lambda_3;
            [lasso_train,~]=mexFistaFlat(tgt_train,terms_train(:,1:size(all_corct_reg,2)),inx(1:size(all_corct_reg,2),1),param2);
            terms_train_in = terms_train(:,lasso_train(:,1)~=0);
            if fit.pooled==1
                terms_test_in = terms_test(:,lasso_train(:,1)~=0);
            end
            param_in = {all_param{1,lasso_train(1:size(all_param,2),1)~=0}};
            timings_in={all_timings{1,lasso_train(1:size(all_param,2),1)~=0}};
            minlag=1;maxlag=10;min_dof=5;
            betas=zeros(0,0);
            nlag_tr=0;
            try
                if fit.AR==0
                    fitstats_train=fitglm(terms_train_in,tgt_train);
                    if isfield(fitstats_train,'Coefficients')
                        betas=fitstats_train.Coefficients.Estimate;
                    else
                        betas=zeros(size(fitstats_train.Variables,2)-1,1);
                    end
                    nlag_tr=0;
                else
                    %get**optimal lag
                    while ((size(tgt_train,1)-maxlag)-(maxlag+size(terms_train_in,2)+1))<min_dof %have at least dof of 2 when testing (number of observations must be at least 2 more than number of explanatory variables)
                        maxlag=maxlag-1;
                    end
                    if maxlag<=minlag
                        maxlag=minlag;
                    end
                    out=lrratio(tgt_train,maxlag,minlag,1,terms_train_in);
                    %get best lag index
                    nonempt3=find(~isnan(out(:,4)));
                    out2=out(nonempt3,:);
                    sig_lag=out2(find(real(out(:,4))<0.001),:);
                    if isempty(sig_lag)
                        sig_lag=out2(find(real(out2(:,4))<0.05),:);
                        if isempty(sig_lag)
                            sig_lag=out2(find((real(out2(:,4))==min(real(out2(:,4))))),:);
                        end
                    end
                    highest_likelihood_i=find(abs(sig_lag(:,3))==max(abs(sig_lag(:,3))));
                    nlag_tr=sig_lag(highest_likelihood_i,1);
                    if size(tgt_train,1)~=size(terms_train_in,2)
                        
                        %get Pramater estimates with the optimal time lag
                        if fit.zscored==0
                            terms_train_const=horzcat(terms_train_in,terms(trIdx,all_const_idx(find(trsub==1))));
                            fitstats_train=ols(tgt_train,terms_train_const);
                            beta_size=size(all_corct_reg,2)+1;
                        else
                            fitstats_train=vare(tgt_train,nlag_tr,terms_train_in);
                            beta_size=size(all_corct_reg,2);
                        end
                        %some specific lags cause an error in calculating f-statistic and fprob. we will ignore these nlags
                        betas=fitstats_train.beta;
                        fitstats_train.lse = sum((fitstats_train.y-fitstats_train.yhat).^2);
                        fitstats_train.aic = length(fitstats_train.yhat)*log(fitstats_train.lse/length(fitstats_train.yhat)) + 2*beta_size;
                        fitstats_train.bic = length(fitstats_train.yhat)*log(fitstats_train.lse/length(fitstats_train.yhat)) + beta_size*log(length(fitstats_train.yhat));
                    end
                end
                
            end
            b=betas;
            if fit.AR~=0
                dof1=size(tgt_train,1)-size(betas,1);
                if isempty(param_in)
                    pval_in=1;
                else
                    for get_pval=1:size(param_in,2)
                        tstat(1,get_pval)=fitstats_train.tstat(get_pval,1);
                        p_value(1,get_pval)=tcdf(tstat(1,get_pval),dof1,'upper');
                    end
                    pval_in=p_value;
                end
                fitstats_train.pval=pval_in;
            end
            
            
            %***test on test set***
            if fit.pooled==1
                test_pred=terms_test_in*betas(nlag_tr+1:size(terms_test_in,2)+nlag_tr,1);
                g2.test_yhat=test_pred;
                g2.test_y=tgt_test;
                result(i).test_yhat=test_pred;
                result(i).test_y=tgt_test;
                
                if fit.AR~=0
                    %get the optimal lag
                    out_test=lrratio(tgt_test,maxlag,minlag,1,test_pred);
                    sig_lag_test=out_test(find(real(out_test(:,4))<0.001),:);
                    if isempty(sig_lag_test)
                        sig_lag_test=out_test(find(real(out_test(:,4))<0.05),:);
                        if isempty(sig_lag_test)
                            sig_lag_test=out_test(find(min(real(out_test(:,4)))),:);
                        end
                    end
                    highest_likelihood_i_test=find(abs(sig_lag_test(:,3))==max(abs(sig_lag_test(:,3))));
                    nlag_tr_test=sig_lag_test(highest_likelihood_i_test,1);
                    tgt_test_lag=mlag(tgt_test,nlag_tr_test);
                    
                    try
                        fitstats_test=vare(tgt_test,nlag_tr_test,test_pred);
                        lag_pred=tgt_test_lag*fitstats_test.beta(1:nlag_tr_test,1);
                        test_pred=test_pred*fitstats_test.beta(nlag_tr_test+1:end-1,1);
                        result(i).fitstats_test=fitstats_test;
                    end
                end
                
                if fit.zscored==0
                    [test_pred_r,test_pred_pval]=corr(zscore(test_pred),zscore(tgt_test));
                    if fit.AR~=0
                        [test_pred_r_lagcorct,test_pred_pval_lagcorct]=corr(zscore(test_pred),zscore(tgt_test-lag_pred));
                    end
                else
                    [test_pred_r,test_pred_pval]=corr(test_pred,tgt_test);
                    if fit.AR~=0
                        [test_pred_r_lagcorct,test_pred_pval_lagcorct]=corr(test_pred,tgt_test-lag_pred);
                    end
                end
                g2.prediction_r2=(test_pred_r)^2;
                g2.predictionpval=test_pred_pval;
                
                if fit.AR~=0
                    fitstats_train.prediction_r2_lagcorct=(test_pred_r_lagcorct)^2;
                    fitstats_train.prediction_pval_lagcorct=test_pred_pval_lagcorct;
                    fitstats_train.test_lag_pred=lag_pred;
                    result(i).prediction_r2_lagcorct=fitstats_train.prediction_r2_lagcorct;
                    result(i).prediction_pval_lagcorct=test_pred_pval_lagcorct;
                    result(i).lag_stats=out;
                end
            end
            
            
            result(i).lambda1=param2.lambda;
            result(i).lambda2=param2.lambda2;
            result(i).lambda3=param2.lambda3;
            result(i).betas = betas;
            result(i).ROIs = param_in;
            result(i).timings=timings_in;
            result(i).nlag=nlag_tr;
            result(i).lambda_BIC=lambda_BIC;
            result(i).model_in=0; %dummy
            result(i).fitstats_train=fitstats_train;
            
            if fit.AR~=0
                if fit.pooled==1
                    result(i).prediction_r2=(test_pred_r)^2;
                    result(i).prediction_pval=(test_pred_pval);
                end
                result(i).rsqr=fitstats_train.rsqr;
                result(i).pval=fitstats_train.pval;
                result(i).train_yhat=fitstats_train.yhat;
                result(i).train_y=fitstats_train.y;
                result(i).bic=fitstats_train.bic;
                result(i).aic=fitstats_train.aic;
            else
                if fit.pooled==1
                    result(i).prediction_r2= g2.prediction_r2;
                    result(i).prediction_pval=g2.predictionpval;
                end
                result(i).rsqr=fitstats_train.Rsquared.Ordinary;
                if isfield(fitstats_train,'Coefficients')
                    try result(i).pval=fitstats_train.Coefficients.pValue; end
                else
                    result(i).pval=ones(size(fitstats_train.Variables,2)-1,1);
                end
                result(i).train_yhat=fitstats_train.Fitted.Response;
                result(i).train_y=fitstats_train.Variables.y;
                result(i).bic=fitstats_train.ModelCriterion.BIC;
                result(i).aic=fitstats_train.ModelCriterion.AIC;
                result(i).sub=fit.sub;
            end
            clear fitstats_train
        end
        
    end
    
    
    
    
    %%
    %%******Do the actual fitting with the best lambda********
    if fit.pooled==1
        if fit.AR~=0
            model_b=find([result(:).prediction_r2_lagcorct]==max([result(:).prediction_r2_lagcorct]));
        else
            model_b=find([result(:).prediction_r2]==max([result(:).prediction_r2]));
        end
        param{1,1}.lambda=result(model_b).lambda1;
        param{1,1}.lambda2=result(model_b).lambda2;
        param{1,1}.lambda3=result(model_b).lambda3;
        %  nlagB=result(model_b).nlag;
        clear fitstats
        all_constB_idx=size(all_corct_reg,2)+1:size(terms,2);
        S=(1:size(fit.sub,2))';
        CVB = cvpartition(S,'k',folds);
        for iB = 1:CVB.NumTestSets
            trsubB = CVB.training(iB);trIdxB=[S_i{find(trsubB==1),1}]';
            tesubB = CVB.test(iB);teIdxB=[S_i{find(tesubB==1),1}]';
            terms_trainB=terms(trIdxB,:);tgt_trainB=tgt(trIdxB,1);
            terms_testB=terms(teIdxB ,:);tgt_testB=tgt(teIdxB,1);
            [lasso_trainB,~]=mexFistaFlat(tgt_trainB,terms_trainB,inx,param{1,1});
            terms_trainB_in = terms_trainB(:,lasso_trainB(:,1)~=0);
            terms_testB_in= terms_testB(:,lasso_trainB(:,1)~=0);
            paramB_in = {all_param{1,lasso_trainB(1:size(all_param,2),1)~=0}};
            timingsB_in={all_timings{1,lasso_trainB(1:size(all_param,2),1)~=0}};
            if fit.AR==0
                fitstats_trainB=fitglm(terms_train_in,tgt_train);
                % beta_size=size(all_corct_reg,2);
                if isfield(fitstats_train,'Coefficients')
                    betasB=fitstats_trainB.Coefficients.Estimate;
                else
                    betasB=zeros(size(fitstats_trainB.Variables,2)-1,1);
                end
                %              fitstats_trainB.lse= fitstats_trainB.SSE;
                %               fitstats_trainB.aic= fitstats_trainB.Modelcriterion.AIC;
                %               fitstats_trainB.bic= fitstats_trainB.Modelcriterion.BIC;
                %             fitstats_trainB.lse = sum((fitstats_trainB.y-fitstats_trainB.yhat).^2);
                %             fitstats_trainB.aic = length(fitstats_trainB.yhat)*log(fitstats_trainB.lse/length(fitstats_trainB.yhat)) + 2*fitstats_trainB.beta_size;
                %             fitstats_trainB.bic = length(fitstats_trainB.yhat)*log(fitstats_trainB.lse/length(fitstats_trainB.yhat)) + fitstats_trainB.beta_size*log(length(fitstats_trainB.yhat));
                nlag_trB=0;
            else
                % % %     %estimate the best time lag in the time series
                minlag=1;maxlag=10;min_dof=5;
                while ((size(tgt_trainB,1)-maxlag)-(maxlag+size(terms_trainB_in,2)+1))<min_dof %have at least dof of 2 when testing (number of observations must be at least 2 more than number of explanatory variables)
                    maxlag=maxlag-1;
                end
                if maxlag<=minlag
                    maxlag=minlag;
                end
                evalc('outB=lrratio(tgt_trainB,maxlag,minlag,1,terms_trainB_in)');
                %     %get best lag index
                sig_lagB=outB(find(real(outB(:,4))<0.001),:);
                if isempty(sig_lagB)
                    sig_lagB=outB(find(real(outB(:,4))<0.05),:);
                    if isempty(sig_lagB)
                        sig_lagB=outB(find(min(real(outB(:,4)))),:);
                    end
                end
                highest_likelihood_iB=find(abs(sig_lagB(:,3))==max(abs(sig_lagB(:,3))));
                nlag_trB=sig_lagB(highest_likelihood_iB,1);
                
                % % %     %get Pramater extimates with the optimal time lag
                
                if size(tgt_trainB,1)~=size(terms_trainB_in,2)
                    try
                        %get Pramater estimates with the optimal time lag
                        if fit.zscored==0
                            terms_trainB_const=horzcat(terms_trainB_in,terms(trIdxB,all_constB_idx(find(trsubB==1))));
                            fitstats_trainB=ols(tgt_trainB,terms_trainB_const);
                            beta_size=size(all_corct_reg,2)+1;
                        else
                            fitstats_trainB=vare(tgt_trainB,nlag_trB,terms_trainB_in);
                            beta_size=size(all_corct_reg,2);
                        end
                        %some specific lags cause an error in calculating f-statistic and fprob. we will ignore these nlags
                        betasB=fitstats_trainB.beta;
                        fitstats_trainB.lse = sum((fitstats_trainB.y-fitstats_trainB.yhat).^2);
                        fitstats_trainB.aic = length(fitstats_trainB.yhat)*log(fitstats_trainB.lse/length(fitstats_trainB.yhat)) + 2*beta_size;
                        fitstats_trainB.bic = length(fitstats_trainB.yhat)*log(fitstats_trainB.lse/length(fitstats_trainB.yhat)) + beta_size*log(length(fitstats_trainB.yhat));
                    end
                end
            end
            
            bB=betasB;
            dofB=size(tgt_trainB,1)-size(betasB,1);
            if fit.AR~=0
                if isempty(paramB_in)
                    pval_in=1;
                else
                    for get_pval=1:size(paramB_in,2)
                        tstat(1,get_pval)=fitstats_trainB.tstat(get_pval,1);
                        p_value(1,get_pval)=tcdf(tstat(1,get_pval),dofB,'upper');
                    end
                    pval_in=p_value;
                end
                
                fitstats_trainB.pval=pval_in;
            end
            %Estimate the nlag in the test data
            
            %***test on test set***
            test_predB=terms_testB_in*betasB(nlag_trB+1:size(terms_testB_in,2)+nlag_trB,1);
            g3.testB_yhat=test_predB;
            g3.testB_y=tgt_testB;
            
            if fit.AR~=0
                outB_test=lrratio(tgt_testB,maxlag,minlag,1,test_predB);
                sig_lagB_test=outB_test(find(real(outB_test(:,4))<0.001),:);
                if isempty(sig_lagB_test)
                    sig_lagB_test=outB_test(find(real(outB_test(:,4))<0.05),:);
                    if isempty(sig_lagB_test)
                        sig_lagB_test=outB_test(find(min(real(outB_test(:,4)))),:);
                    end
                end
                highest_likelihood_iB_test=find(abs(sig_lagB_test(:,3))==max(abs(sig_lagB_test(:,3))));
                nlag_trB_test=sig_lagB_test(highest_likelihood_iB_test,1);
                tgt_testB_lag=mlag(tgt_testB,nlag_trB_test);
                try
                    %             tgt_testB_lagcorct=tgt_testB_lagcorct(nlag_trB_test+1:end,end);
                    fitstats_testB=vare(tgt_testB,nlag_trB_test,test_predB);
                    lag_predB=tgt_testB_lag*fitstats_testB.beta(1:nlag_trB_test,1);
                    %           test_predB=test_predB*fitstats_testB.beta(nlag_trB_test+1:end-1,1);
                end
            end
            if fit.zscored==0
                [test_predB_r,test_predB_pval]=corr(zscore(test_predB),zscore(tgt_testB));
                if fit.AR~=0
                    [test_predB_r_lagcorct,test_predB_pval_lagcorct]=corr(zscore(test_predB),zscore(tgt_testB-lag_predB));
                end
            else
                [test_predB_r,test_predB_pval]=corr(test_predB,tgt_testB);
                if fit.AR~=0
                    [test_predB_r_lagcorct,test_predB_pval_lagcorct]=corr(test_predB,tgt_testB-lag_predB);
                end
            end
            g3.prediction_r2=(test_predB_r)^2;
            g3.prediction_pval=test_predB_pval;
            if fit.AR~=0
                fitstats_trainB.prediction_r2_lagcorct=(test_predB_r_lagcorct)^2;
                fitstats_trainB.prediction_pval_lagcorct=test_predB_pval_lagcorct;
                fitstats_trainB.lag_predB=lag_predB;
            end
            fitstats(iB).testB_yhat=g3.testB_yhat;
            fitstats(iB).testB_y=g3.testB_y;
            fitstats(iB).prediction_r2=g3.prediction_r2;
            fitstats(iB).prediction_pval=g3.prediction_pval;
            fitstats(iB).trainB=fitstats_trainB;
            if fit.AR==0
                fitstats(iB).rsqr=fitstats_trainB.Rsquared.Ordinary;
                fitstats(iB).bic=fitstats_trainB.ModelCriterion.BIC;
                fitstats(iB).aic=fitstats_trainB.ModelCriterion.AIC;
                
                if isfield(fitstats_trainB,'Coefficients')
                    try fitstats(iB).pval=fitstats_trainB.Coefficients.pValue; end
                else
                    fitstats(iB).pval=ones(size(fitstats_trainB.Variables,2)-1,1);
                end
                fitstats(iB).beta=fitstats_trainB.Coefficients.Estimate;
                fitstats(iB).y=fitstats_trainB.Variables.y;
                fitstats(iB).yhat=fitstats_trainB.Fitted.Response;
            end
            clear fitstats_trainB;
        end
    else
        
        % % %
        %Pred,lse, R2, tstat, p_val of uncorrected Predictors
        fitstats=result;
        r2=fitstats.rsqr;
        BOLDpred=fitstats.train_yhat;
        fitstats.lse=sum((fitstats.train_y-fitstats.train_yhat).^2);
        lse=fitstats.lse;
        dof1=size(fitstats.train_y,1)-2;
        fitstats.model_tstat= (r2.*(sqrt(dof1)))./sqrt(1-r2.^2);
        try fitstats.model_pval=tcdf(fitstats.model_tstat,dof1,'upper'); end
        
        if  fit.AR==1
            lagged_BOLD_fit=ols(fitstats.y,fitstats.xmat(:,nlag+1:end));
            lagged_BOLDpred_fused=lagged_BOLD_fit.yhat;
            fitstats.lagged_yhat=lagged_BOLD_fit.yhat;
            meanresult = mean(lagged_BOLD_fit.y,1);
            re_lagged=sum((lagged_BOLD_fit.y-meanresult).^2);
            lagged_lse_fused = sum((lagged_BOLD_fit.y-lagged_BOLDpred_fused).^2); %sum least-squares error (SSR)
            fitstats.lagged_BOLDr2_fused = 1-(lagged_lse_fused/re_lagged);
            rt1=fitstats.lagged_BOLDr2_fused;
            dof1=size(lagged_BOLD_fit.y,1)-2;
            r_stat1= (rt1.*(sqrt(dof1)))./sqrt(1-rt1.^2);
            fitstats.lagged_BOLD_pval=tcdf(r_stat1,dof1,'upper');
            fitstats.lagged_BOLD_lse=lagged_lse_fused;
            fitstats.lagged_BOLD_tstat=r_stat1;
            fitstats.lagged_BOLD_bic = length(fitstats.lagged_yhat)*log(fitstats.lagged_BOLD_lse/length(fitstats.lagged_yhat)) + size(lagged_BOLD_fit.beta,1)*log(length(fitstats.lagged_yhat));
            
            %Pred,lse, R2,tstat, p_val of lag correction terms
            lagterms_BOLD_fit=ols(fitstats.y,fitstats.xmat(:,1:nlag));
            %lag_betas=fitstats.beta(1:nlag,1);
            fitstats.lagterms_pred=lagterms_BOLD_fit.yhat;
            meanresult2 = mean(lagterms_BOLD_fit.y,1);
            re_lagterms=sum((lagterms_BOLD_fit.y-meanresult2).^2);
            lse_lagterms=sum((lagterms_BOLD_fit.y-lagterms_BOLD_fit.yhat).^2);
            lagterms_r2=1-(lse_lagterms/re_lagterms);
            fitstats.lagterms_r2=lagterms_r2;
            rt2=fitstats.lagterms_r2;
            if rt2>1 || rt2<-1
                rt2=0;
            end
            r_stat2= (rt2.*(sqrt(dof1)))./sqrt(1-rt2.^2);
            fitstats.lagterms_tstat=r_stat2;
            fitstats.lagterms_pval=tcdf(r_stat2,dof1,'upper');
            fitstats.lagterms_lse=lse_lagterms;
            fitstats.lagterms_bic = length(lagterms_BOLD_fit.yhat)*log(fitstats.lagterms_lse/length(lagterms_BOLD_fit.yhat)) + size(lagterms_BOLD_fit.beta,1)*log(length(lagterms_BOLD_fit.yhat));
            
            %save results in the output stucture
            
            
            if isempty(param_in)
                pval_in=1;
            else
                for get_pval=1+nlag:size(param_in,2)+nlag
                    p_value(get_pval-nlag,1)=fitstats.tprob(get_pval);
                end
                pval_in=p_value;
            end
        end
    end
    
    %*******Do Lasso Regression*********
    cd '/Applications/spm12/toolbox/spams-matlab-v2.6'
    start_spams
    addpath(genpath('/Applications/spm12/toolbox/jplv7'))
    fprintf('\n Dependent explanatory variables & autocorrelated observations \n Using Fused lasso regression \n')
elseif ~isempty(strfind(fit.dependence,'NCA'))
    NCAmdl=fsrnca(terms,tgt,'Standardize',1);
    ft_in= NCAmdl.FeatureWeights > 0.001;
    result.betas= NCAmdl.FeatureWeights(ft_in);
    result.BOLDpred=terms(:,ft_in')*result.betas;
    result.rsqr=corr(result.BOLDpred,tgt)^2;
    result.train_yhat=result.BOLDpred;
    result.lse=sum(result.train_yhat-tgt)^2;
    result.bic = length(result.BOLDpred)*log(result.lse/length(result.BOLDpred)) + length(result.betas)*log(length(result.BOLDpred));
    result.aic = length(result.BOLDpred)*log(result.lse/length(result.BOLDpred)) + 2*length(result.betas);
    result.nlag=0;
    
    %     if isempty(ft_in)
    %         pval_in=1;
    %     else
    %         for get_pval=1+result.nlag:size(result.ROIs,2)+result.nlag
    %             p_value(get_pval-result.nlag,1)=fitstats.tprob(get_pval);
    %         end
    %         pval_in=p_value;
    %     end
    %     result.sub=sub;
    %     if isempty(result.ROIs)
    %     result.model_in=-1;
    %     else
    %     result.model_in=1;
    %     end
    
    
elseif ~isempty(strfind(fit.dependence,'lasso')) && isempty(strfind(fit.dependence,'fused'))
    if ~isempty(strfind(fit.dependence,'ridge'))
        alpha=0.2;
    elseif ~isempty(strfind(fit.dependence,'-el'))
        alpha=0.5;
    else
        alpha=1;
    end
    if fit.pooled==1
        terms_pool=terms;
        terms_pool(:,end+1)=tgt;
        terms_pool(:,end+1)=rand(nobs,1);
        terms_poolsort=sortrows(terms_pool,size(terms_pool,2));
        terms=terms_poolsort(:,1:end-2);
        tgt=terms_poolsort(:,end-1);
    end
    tic
    alpha_v=0.1:0.1:1;
    model_CVM=zeros(size(alpha_v,2),2);
    for alpha_i=1:size(alpha_v,2)
        alpha=alpha_v(1,alpha_i);
        CVO = cvpartition(CV,'k',10);
        CVM=cell(CVO.NumTestSets,1);
        lso=cell(1,CVO.NumTestSets);
        terms_train=tcell(1,CVO.NumTestSets);
        fitstat=cell(1,CVO.NumTestSets);
        betas=cell(1,CVO.NumTestSets);
        dof=cell(1,CVO.NumTestSets);
        x=CVO.NumTestSets;
        parfor i = 1:CVO.NumTestSets
            x=30;
            if i>x
                fprintf('wrong');
            end
            trIdx = CVO.training(i);
            teIdx = CVO.test(i);
            terms_train{i}=terms(trIdx,:);tgt_train=tgt(trIdx,1);
            terms_test=terms(teIdx ,:);tgt_test=tgt(teIdx,1);
            %getting features
            [lso{i},fitstat{i}]=lasso(terms_train{i}(:,1:size(all_corct_reg,2)),tgt_train,'CV',10,'PredictorNames',all_param(1,1:size(all_corct_reg,2)),'Alpha',alpha);
            %         param_in = fitstats.PredictorNames(lso(:,fitstats.IndexMinMSE)~=0);
            %         timings_in={all_timings{1,lso(:,fitstats.IndexMinMSE)~=0}};
            terms_train_in=terms_train{i}(:,lso{i}(:,fitstat{i}.IndexMinMSE)~=0);
            terms_test=horzcat(terms_test(:,lso{i}(:,fitstat{i}.IndexMinMSE)~=0),terms_test(:,size(all_corct_reg,2)+1:end));
            %         sparsemodel = fitstats.PredictorNames(lso(:,fitstats.Index1SE)~=0);
            %         lassoPlot(lso,fitstats);
            %         lassoPlot(lso,fitstats,'PlotType','CV');
            %         legend('show')
            %getting parameters
            
            terms_train_in_const=horzcat(terms_train_in,terms_train{i}(:,size(all_corct_reg,2)+1:end));
            fitstat{i}=ols(tgt_train,terms_train_in_const);
            %Calculate R-squared
            betas{i}=fitstat{i}.beta; %getting the betas at minimum MSE
            BOLDpred_training=fitstat{i}.yhat;
            lse_training = nansum((fitstat{i}.y-BOLDpred_training).^2); %sum least-squares error (SSR)
            BOLDr2_training = fitstat{i}.rsqr;
            dof{i}=nobs-size(terms_train_in_const,2);
            %Apply on the test set
            BOLDpred_test=terms_test*betas{i};
            [training_test_r,training_test_pval]=corr(tgt_test,BOLDpred_test);
            training_test_r2=training_test_r^2;
            %         BOLDr2_test=fitstats_test.rsqr;
            dof1=size(BOLDpred_test,1)-size(BOLDpred_test,2);
            test_tstat= (BOLDr2_test.*(sqrt(dof1)))./sqrt(1-BOLDr2_test.^2);
            test_pval=tcdf(test_tstat,dof1,'upper');
            %         fitstats_test.lse=sum((fitstats_test.y-BOLDpred_test).^2);
            %         if fit.zscored==0
            %             fitstats_test.beta_size=size(all_corct_reg,2)+1;
            %         end
            CVM{i,1}=[i,training_test_r2,training_test_pval];
        end
        CVM_v=[CVM{:,1}];
        model_pred_r2= nanmean(CVM_v(2:3:size(CVM_v,2)));
        model_CVM(alpha_i,1)=alpha;
        model_CVM(alpha_i,2)=model_pred_r2;
        clear CVM_v
    end
    toc
    best_alpha= CVM(find(CVM(:,2)==max(CVM(:,2))),1);
    
    %DO a lasso regression based on CV results.
    alpha=best_alpha;
    [lso,fitstats]=lasso(terms_train(:,1:size(all_corct_reg,2)),tgt,'CV',10,'PredictorNames',all_param(1,1:size(all_corct_reg,2)),'Alpha',alpha);
    param_in = fitstats.PredictorNames(lso(:,fitstats.IndexMinMSE)~=0);
    timings_in={all_timings{1,lso(:,fitstats.IndexMinMSE)~=0}};
    terms_in=terms(:,lso(:,fitstats.IndexMinMSE)~=0);
    lassoPlot(lso,fitstats);
    lassoPlot(lso,fitstats,'PlotType','CV');
    legend('show')
    
    terms_in_const=horzcat(terms_in,terms(:,size(all_corct_reg,2)+1:end));
    fitstats=ols(tgt,terms_in_const);
    %Calculate R-squared
    betas=fitstats.beta; %getting the betas at minimum MSE
    BOLDpred=fitstats.yhat;
    lse = nansum((fitstats.y-BOLDpred).^2); %sum least-squares error (SSR)
    r2 = fitstat.rsqr;
    if isempty(param_in)
        pval_in=1;
    else
        for get_pval=1:size(terms_in_const,2)
            p_value(get_pval,1)=tcdf(abs(fitstats.tstat(get_pval)),dof,'upper');
        end
        pval_in=p_value;
    end
    b=betas(1:size(all_corct_reg,2));
    if fit.zscored==0 %adds one constant term to the explanatory varaiables (irrespective of the number of subjects)
        %to calculate BIC and AIC properly
        b=betas(1:size(all_corct_reg,2)+1);
    end
    
    %Plot Lasso predictor weights*****
    if fit.pooled==0
        figure('name','PredictorWeightsInBestModels');
        bar(lso(:,[fitstats.IndexMinMSE fitstats.Index1SE])); try xlabel_oblique(ROIs_names); end
        set(gca,'position',[0.1300    0.1500    0.7750    0.7750],'xlim',[0 11])
        ylabel('Predictor weights')
        legend({'model with min MSE';'sparsest model within 1SE of min'})
        title('Predictor weights for different best-fitting models')
        printfig
        
        % Fit the best model, compare to original data and examine residuals.
        % % %     Xplus = [ones(size(terms,1),1) terms];
        % % %     fitSparse = Xplus * [fitstats.Intercept(fitstats.Index1SE); lso(:,fitstats.Index1SE)];
        % % %     fitMinMSE = Xplus * [fitstats.Intercept(fitstats.IndexMinMSE); lso(:,fitstats.IndexMinMSE)];
        % % %     r2 = corr(tgt,fitMinMSE)^2;
        % % %     r_res = corr(fitMinMSE,tgt-fitMinMSE);
        %     %Display results
        figure('name','fittings of data and residuals');
        subplot(1,2,1)
        plot(tgt,fitMinMSE,'o'); xlabel('Data'); ylabel('Fit')
        title(sprintf('Min-MSE model: lambda = %.3f, r^2 = %.2f',fitstats.LambdaMinMSE,r2))
        
        subplot(1,2,2)
        plot(fitMinMSE,tgt-fitMinMSE,'o')
        title(sprintf('Residuals of MinMSE model fit (r = %.2f)',r_res))
    end
    
    %*****using minimum constraints nonlinear solver****
    
elseif strcmp(fit.dependence,'fmincon')
    fprintf('\n Assuming independent explanatory variables: Using fmincon \n')
    param_in=all_param;
    pval_in=1;
    timings_in=all_timings;
    
    rptfit = 100; %repeat N times jittering the betas
    
    options = optimset('Display','off','MaxIter',100000,'TolFun',1e-10,'TolX',1e-10,...
        'DiffMaxChange',1e-2,'DiffMinChange',1e-4,'MaxFunEvals',10000,'LargeScale','off');
    warning off; %display,iter to see outputs
    
    for fill_bound= 1:proper_model_size
        lb(1,fill_bound) = -100;
        ub(1,fill_bound)= 100;
    end
    
    %     if exist('constant','var') && ~isnan(constant)
    %         inx = [inx mean(TargetBOLD)];
    %         lb = [lb 0];
    %         ub = [ub 1];
    %     end;
    
    result.modelName = 'fit_BOLD';
    result.inx = double(inx);
    result.lb  = lb;
    result.ub  = ub;
    b=inx;
    result.betas = b;
    [lse, BOLDpred, r2] = fit_BOLD(b, terms,tgt);
    
    for n=1:rptfit %jitter the betas
        [b,fval,exitflag,fitstats,lambda,grad,hessian] = fmincon(@fit_BOLD, double(b+0.1*randn(size(b))), [],[],[],[],lb,ub,[], options, terms,tgt);
        [lse2, BOLDpred, r2] = fit_BOLD(b, terms,tgt);
        if lse2<lse
            result.b = b;
            [lse, BOLDpred, r2,tstat] = fit_BOLD(b, terms,tgt);
        end;
    end;
    s=size(BOLDpred,1); df=s-2;
    p_value=tcdf(tstat,df,'upper');
    pval_in=p_value;
end

if ~isempty(strfind(fit.dependence,'fmincon')) || ~isempty(strfind(fit.dependence,'-el'))
    result.betas = b;
    result.BOLDpred = BOLDpred;
    result.rsqr = r2;
    result.bic = length(BOLDpred)*log(lse/length(BOLDpred)) + length(b)*log(length(BOLDpred));
    result.aic = length(BOLDpred)*log(lse/length(BOLDpred)) + 2*length(b);
    result.ROIs = param_in;
    result.pval=pval_in;
    result.timings=timings_in;
    result.model_in=0;
    %dummy
    if exist('constant','var') && ~isnan(constant)
        result.paramNames = [result.paramNames, {'constant'}];
    end
end
end

function b=lassocv(XTRAIN,ytrain,XTEST)
[lso,fitstats]=lasso(XTRAIN,ytrain,'CV',10);
b=XTEST*(lso(:,fitstats.IndexMinMSE));
end

function [lse, BOLDpred, BOLDr2,reg_tstat] =  fit_BOLD(x, terms,target)
for ind_fill= 1:size(terms,2)
    betas(ind_fill) = x(ind_fill);
end

BOLDpred=nansum(betas.*terms,2)/size(terms,2);

meanresult = mean(target,1);
lse = nansum((target-BOLDpred).^2); %sum least-squares error (SSR)
re=nansum((target-meanresult).^2);
BOLDr2 = 1-(lse/re);
s=size(BOLDpred,1); df=s-2;
reg_tstat=sqrt(BOLDr2)*sqrt(df/(1-BOLDr2));
end