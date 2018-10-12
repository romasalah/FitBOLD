function [result,tgt,shifts,fitstats] = BOLD_modelling(ROIs,timings,TargetBOLD,target_timing,~,...
    behav,dir,sub,shift,equal_regressors,whatcombine,howcombine,pooled,zscored,dependence,bimodal,MA,AR,disk_parallel,use_4D,social_ROIs,condition_before,model,FIR)
clear corct_BOLD
all_BOLD_combined=cell(0,0);
all_BOLD_singles=cell(0,0);
all_behav_combined=cell(0,0);
all_behav_singles=cell(0,0);

for look_duplicate=1:size(timings,2)
    if strcmp(timings{1,look_duplicate},target_timing{1,1})
        bimodal=0;
    end
end

if size(ROIs,2)==size(timings,2)
    combine_size=size(whatcombine,1);
    fprintf('\n ROIs and timings sizes are matching \n')
    fprintf('\n getting the corrected BOLD arrays (shifted and filtered) \n')
    fprintf('\n This script is using corr_BOLD function \n')
    
    %************loading target BOLD vector***********
    
    fprintf('\n Getting your target BOLD or target Behavioural vector \n')
    tgt=zeros(0,1);
    if behav==1
        [~,~,~,shifts,~,~,~,~,target_onsets,~,tgt_sess]=corr_BOLD(TargetBOLD,{target_timing{:,1}}',target_timing{1,1},dir,sub,equal_regressors,shift,0,0,disk_parallel,use_4D,MA,zscored,FIR);
        for z_tgt=1:size(tgt_sess,2)
                    tgt=vertcat(tgt,tgt_sess{1,z_tgt});
        end
        
        fprintf('\n loading target Behavioural vector was Successful \n')
    else
        [~,~,~,shifts,~,~,~,~,target_onsets,BOLD_tgt_sess]=corr_BOLD(TargetBOLD,{target_timing{:,1}}',target_timing{1,1},dir,sub,equal_regressors,shift,0,0,disk_parallel,use_4D,MA,zscored,FIR);
            for z_BOLD_tgt=1:size(BOLD_tgt_sess,2)
                tgt=vertcat(tgt,BOLD_tgt_sess{1,z_BOLD_tgt});
            end
        fprintf('\n loading target BOLD vector was Successful \n')
    end
     if FIR~=0
            filter= BOLD_filt(tgt);
      tgt=filtfilt(filter,tgt);
     end
    %%
    %*************getting and combining BOLD from different ROIs**************
    if exist('whatcombine','var') && whatcombine(1,1) ~=0
        ROIs_names=cell(1,size(whatcombine,1));
        timings_for_BOLD=cell(1,size(whatcombine,1));
        for ROI_names_i=size(whatcombine,2):size(whatcombine,2):numel(whatcombine)
            ROIs_names{1,ROI_names_i}=ROIs{ROI_names_i};
            timings_for_BOLD{1,ROI_names_i}=timings{ROI_names_i}(1,:);
        end
        for fill_rest_ROI_names= ROI_names_i+1:size(ROIs,2)
            ROIs_names{1,end+1}=ROIs{fill_rest_ROI_names}; %#ok<*AGROW>
            timings_for_BOLD{1,end+1}=timings{fill_rest_ROI_names};
        end
        non_empt_idx=find(~cellfun(@isempty,ROIs_names));
        ROIs_names={ROIs_names{1,[non_empt_idx]}};
        timings_for_BOLD={timings_for_BOLD{1,[non_empt_idx]}};
        %correct names to bilateral
        if ~isempty(ROIs_names)
            for corct_bi=1:size(ROIs_names,2)
                ROIs_names{1,corct_bi}=[ROIs_names{1,corct_bi}(1:end-10) 'Bi'];
            end
        end
        
        for combine_loop=1:size(whatcombine,1)
            combinebold1=['\n getting and combining BOLD from different ROIs \n'];
            combinebold2= ['Process number ' num2str(combine_loop) ' out of ' num2str(size(whatcombine,1)) ' combinations '];
            combinebold3=[' Doing ROIs number ' num2str(whatcombine(combine_loop,:))];
            fprintf(combinebold1);fprintf(combinebold2);fprintf(combinebold3);
            for roi_i=1:size(whatcombine,2)
                ROI_temp{1,roi_i}=ROIs{1,whatcombine(combine_loop,roi_i)};
            end
            if ~strcmp(timings{1,combine_loop},target_timing)
                equal_regressors=0;
            end
            [~,~,behav_onecombined,~,corct_BOLD_combined,~,~,~,onsets_filtered,BOLDpersession,behavpersessioncombined]=...
                corr_BOLD(ROI_temp,{timings{:,combine_loop}}',timings{1,combine_loop},dir,sub,equal_regressors,shift,0,howcombine,disk_parallel,use_4D,MA,zscored,FIR);
            all_BOLD_combined{1,combine_loop}=corct_BOLD_combined;
            all_behav_combined{1,combine_loop}=behav_onecombined;
            
            if ~strcmp(timings{1,combine_loop},target_timing)
                [BOLD_output,behav_output,~]=corct_timing(BOLDpersession,onsets_filtered,target_onsets,behavpersessioncombined);
                all_BOLD_combined{1,combine_loop}=BOLD_output;
                all_behav_combined{1,combine_loop}=behav_output;
                
                for ind_fill2= 1:size(BOLD_output,2)
%                     if zscored==1
%                         all_BOLD_combined{1,combined_ROIs_loop}{1,ind_fill2}=zscore(all_BOLD_combined{1,combined_ROIs_loop}{1,ind_fill2});
%                         all_behav_combined{1,combined_ROIs_loop}{1,ind_fill2}=zscore(all_behav_combined{1,combined_ROIs_loop}{1,ind_fill2});
%                     end
                end
                comb_temp=zeros(0,1);
                comb_temp_behav=zeros(0,1);
                for ind_fill2= 1:size(BOLD_output,2)
                    comb_temp=vertcat(comb_temp, all_BOLD_combined{1,combined_ROIs_loop}{1,ind_fill2});
                    comb_temp_behav=vertcat(comb_temp_behav, all_behav_combined{1,combined_ROIs_loop}{1,ind_fill2});
                end
                all_BOLD_combined{1,combined_ROIs_loop}=comb_temp;
                all_behav_combined{1,combined_ROIs_loop}=comb_temp_behav;
                clear comb_temp
                clear ROI_temp
                clear BOLDpersession
                clear onsets_filtered
                clear behavpersessioncombined
            else
                for ind_fill1= 1:size(BOLDpersession,2) 
                            all_BOLD_combined{1,combined_ROIs_loop}{1,ind_fill1}=BOLDpersession{1,ind_fill1};
                            all_behav_combined{1,combined_ROIs_loop}{1,ind_fill1}=behavpersessioncombined{1,ind_fill1};
                end
                comb_temp=zeros(0,1);
                comb_temp_behav=zeros(0,1);
                for ind_fill1= 1:size(BOLDpersession,2)
                    comb_temp=vertcat(comb_temp, all_BOLD_combined{1,combined_ROIs_loop}{1,ind_fill1});
                    comb_temp_behav=vertcat(comb_temp_behav, all_behav_combined{1,combined_ROIs_loop}{1,ind_fill1});
                end
                all_BOLD_combined{1,combined_ROIs_loop}=comb_temp;
                all_behav_combined{1,combined_ROIs_loop}=comb_temp_behav;
                clear comb_temp
                clear comb_temp_behav
                clear BOLDpersession
                clear onsets_filtered
                clear behavpersessioncombined
            end
            
        end
        
        idx=1:1:size(ROIs,2);
        whatcombinevector=(reshape(whatcombine,size(whatcombine,1)*size(whatcombine,2),1))';
        for find_rest=1:size(idx,2)
            if  ~isempty(find(whatcombinevector==idx(find_rest)))
                idx(find_rest)=0;
            end
        end
        rest_idx= nonzeros(idx);
        for single_roi_fill= 1:size(rest_idx,1)
            single_ROIs{1,single_roi_fill}=ROIs{1,rest_idx(single_roi_fill)};
        end
        
    else
        ROIs_names=ROIs;
        timings_for_BOLD=timings(1,:);
        single_ROIs=ROIs;
    end
    
    
    %*************Getting Sinlge ROIs BOLD vectors************
    if exist('single_ROIs','var')
        for single_ROIs_loop= 1:size(single_ROIs,2)
            gettingbold=[' Getting BOLD from ROI number ' num2str(single_ROIs_loop) ' out of ' num2str(size(ROIs,2)) ' ROIs '];
            fprintf(gettingbold);
            if ~strcmp(timings{1,single_ROIs_loop},target_timing)
                equal_regressors=0;
            end
            [~,~,~,~,~,~,~,~,onsets_filtered,BOLDpersession,behavpersessionsingles,session_num]=...
                corr_BOLD(single_ROIs{single_ROIs_loop},{timings{:,single_ROIs_loop}}',timings{1,single_ROIs_loop},dir,sub,equal_regressors,shift,0,0,disk_parallel,use_4D,MA,zscored,FIR);
            
            if ~strcmp(timings{1,single_ROIs_loop},target_timing)
                [BOLD_output,behav_output,~]=corct_timing(BOLDpersession,onsets_filtered,target_onsets,behavpersessionsingles);
                
                all_BOLD_singles{1,single_ROIs_loop}=BOLD_output;
                all_behav_singles{1,single_ROIs_loop}=behav_output;
                for ind_fill1= 1:size(BOLD_output,2)
%                     if zscored==1
%                         all_BOLD_singles{1,single_ROIs_loop}{1,ind_fill1}=zscore(all_BOLD_singles{1,single_ROIs_loop}{1,ind_fill1});
%                         all_behav_singles{1,single_ROIs_loop}{1,ind_fill1}=zscore(all_behav_singles{1,single_ROIs_loop}{1,ind_fill1});
%                     end
                end
                sing_temp=zeros(0,1);
                sing_temp_behav=zeros(0,1);
                for ind_fill1= 1:size(BOLD_output,2)
                    sing_temp=vertcat(sing_temp, all_BOLD_singles{1,single_ROIs_loop}{1,ind_fill1});
                    sing_temp_behav=vertcat(sing_temp_behav, all_behav_singles{1,single_ROIs_loop}{1,ind_fill1});
                end
                all_BOLD_singles{1,single_ROIs_loop}=sing_temp;
                BOLD_back=BOLDpersession;
                clear sing_temp
                clear sing_temp_behav
                clear BOLDpersession
                clear onsets_filtered
                clear behavpersessionsingles
            else
                
                for ind_fill1= 1:size(BOLDpersession,2)
                            all_BOLD_singles{1,single_ROIs_loop}{1,ind_fill1}=BOLDpersession{1,ind_fill1};
                            all_behav_singles{1,single_ROIs_loop}{1,ind_fill1}=behavpersessionsingles{1,ind_fill1};
                end
                sing_temp=zeros(0,1);
                sing_temp_behav=zeros(0,1);
                for ind_fill1= 1:size(BOLDpersession,2)
                    sing_temp=vertcat(sing_temp, all_BOLD_singles{1,single_ROIs_loop}{1,ind_fill1});
                    sing_temp_behav=vertcat(sing_temp_behav, all_behav_singles{1,single_ROIs_loop}{1,ind_fill1});
                end
                all_BOLD_singles{1,single_ROIs_loop}=sing_temp;
                all_behav_singles{1,single_ROIs_loop}=sing_temp_behav;
                BOLD_back=BOLDpersession;
                clear sing_temp
                clear sing_temp_behav
                clear BOLDpersession
                clear onsets_filtered
                clear behavpersessionsingles
            end
            
        end
        all_corct_BOLD=horzcat(all_BOLD_singles,all_BOLD_combined);
        all_corct_behav=horzcat(all_behav_combined,all_behav_singles);
    else
        all_corct_BOLD=all_BOLD_combined;
        all_corct_behav=all_behav_combined;
    end
    
    if bimodal==1
        all_corct_reg=horzcat(all_corct_BOLD,all_corct_behav);
    else
        all_corct_reg=all_corct_BOLD;
    end
    
else
    fprintf('\n ROIs and timings size mismatch \n make sure to put the first two parameters in curly braces \n')
    return
end

fprintf('\n Loading the BOLD arrays done \n')
proper_model_size=size(all_corct_reg,2);

if bimodal==1
    all_param=horzcat(ROIs_names,timings(1,:));
    all_timings=horzcat(timings_for_BOLD(1,:),timings(1,:));
else
    all_param=  ROIs_names;
    all_timings=timings_for_BOLD(1,:);
end
%%
%******Preparing the fit********
for ind_fill= 1:proper_model_size
    if FIR~=0
        filter= BOLD_filt([all_corct_reg{1,ind_fill}],2);
        terms(ind_fill,:)=filtfilt(filter,[all_corct_reg{1,ind_fill}]);
    else
    terms(ind_fill,:)=[all_corct_reg{1,ind_fill}];
    end
end
terms=terms';
nobs=size(terms,1);

clear CM
last_sub=0;
last_row=0;
if ~strcmp(timings{1,single_ROIs_loop},target_timing)
    BOLD_ref=BOLD_output;
else
    BOLD_ref=BOLD_back;
end
CV=zeros(0,1);
S_i=cell(size(session_num,1),1);
for const_m= 1:size(session_num,1)
    temp_ones=zeros(0,1);
    CM(1:size(terms,1),const_m)=zeros(size(terms,1),1);
    for sess_ones=last_sub+1:last_sub+session_num(const_m)
        temp_ones=vertcat(temp_ones,ones(size(BOLD_ref{1,sess_ones},1),1));
    end
    temp_ones_s=size(temp_ones,1);
    cv_sub=const_m-1+ones(size(temp_ones,1),1);
    CV=vertcat(CV,cv_sub);
    CM(last_row+1:(last_row+size(temp_ones,1)),const_m)= temp_ones;
    S_i{const_m,1}=last_row+1:(last_row+size(temp_ones,1));
    last_row=last_row+size(temp_ones,1);
    last_sub=last_sub+session_num(const_m);
    clear cv_sub
end
if zscored==0
    terms=horzcat(terms,CM);
end
nanterms=find(isnan(terms));
terms(nanterms)=0;
meanresult = mean(tgt,1);
re=nansum((tgt-meanresult).^2);
for fill_dum= 1:size(terms,2)
    inx(1,fill_dum) = 0.5;
end

%%
%Checking for Correlations between regressors
corrval_array=zeros(size(all_corct_reg,2),size(all_corct_reg,2));
pval_array=zeros(size(all_corct_reg,2),size(all_corct_reg,2));
for source_col=1:size(all_corct_reg,2)
    source_temp=all_corct_reg{1,source_col};
    rest_cols=find(1:size(all_corct_reg,2)~=source_col);
    for rest_cols_i= 1:size(rest_cols,2)
        rest_col_temp=all_corct_reg{1,rest_cols(rest_cols_i)};
        [corrval,pval]=corr(source_temp,rest_col_temp);
        corrval_array(source_col,rest_cols(rest_cols_i))=corrval;
        pval_array(source_col,rest_cols(rest_cols_i))=pval;
    end
end



if behav==1 tgt_type='behav'; else tgt_type='BOLD'; end
if howcombine==1 comb_type='Bi'; else comb_type='Uni';end
if pooled==1 pooled_type='pool'; else pooled_type='ind'; end
if strcmp(shift,'standard shift') shifttype='stdshft'; else shifttype='trushft'; end
if condition_before; which_con_before=condition_before; else which_con_before='NSCB';end
if social_ROIs==1; social_ROIs_name='SocROI'; else social_ROIs_name='NSocROI';end
if zscored; zscore_state='Zscr'; else zscore_state='intrcpt'; end
if AR~=0 AR_state= ['AR(n)']; else AR_state='';end
if isa(MA,'double') && MA~=0 MA_state= ['MA(' num2str(MA) ')']; elseif isa(MA,'double') && MA==0 MA_state=''; elseif ~isa(MA,'double') MA_state=MA.filterstructure; end
modelName=[model '_' tgt_type '_'  pooled_type '_' comb_type '_' ...
    shifttype '_' which_con_before '_' zscore_state '_' dependence '_' social_ROIs_name '_' AR_state '_' MA_state '_' num2str(size(sub,2)) 'Subjects'];
allvar_name=[ modelName '.mat'];
folder='/Users/NEMO_admin/Desktop/Omar/';
save(fullfile(folder,allvar_name))

%%
if strcmp(dependence,'stepwise')
    fprintf('\n Dependent explanatory variables: Using Stepwise regression \n')
    %%%*****Do Stepwise fit******
    [betas,se,pval,inmodel,fitstats,nxt]=stepwisefit(terms,tgt);
    inreg=find([inmodel==1]);
    if ~isempty(inreg)
        param_in={all_param{1,[inreg]}};
        timings_in={all_timings{1,[inreg]}};
        betas_in=betas(inreg,1);
        %          se_in=se(inreg,1);
        pval_in=pval(inreg,1);
        terms_in=terms(:,inreg);
    else
        inreg=find(pval==min(pval));
        inreg=inreg(1,1);
        terms_in=terms(:,inreg);param_in={all_param{1,[inreg]}}; timings_in={all_timings{1,[inreg]}};
        [betas_in2,~,~,~,fitstats]=regress(tgt,terms_in);
        betas_in=betas_in2(1,1);
        pval_in=fitstats(1,3);
        if isnan(pval_in)
            pval_in=1; %coverts the Nan P-value to P-value of 1 which will exclude the fit.
        end
    end
    %*****getting r-squared*****
    BOLDpred_stepwise=terms_in*betas_in;
    lse_stepwise = nansum((tgt-BOLDpred_stepwise).^2); %sum least-squares error (SSR)
    r2_stepwise = 1-(lse_stepwise/re);
    BOLDpred=BOLDpred_stepwise;
    r2=r2_stepwise;
    lse=lse_stepwise;
    b=betas_in;
    
    %%    %********Do Ridge regression*******
elseif strcmp(dependence,'ridge')
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
elseif (strfind(dependence,'fused'))
    fprintf('\n Dependent explanatory variables: \n Using Fused lasso followed by Autoregressive OLS regression \n')
    cd '/Applications/spm12/toolbox/spams-matlab-v2.6'
    start_spams
    addpath(genpath('/Applications/spm12/toolbox/jplv7'))
    fprintf('\n Dependent explanatory variables & autocorrelated observations \n Using Fused lasso regression \n')
    inx=inx';
    if zscored==0
        max_lambda=10000;
        nvar=size(terms,2);
        nobs=size(terms,1);
        lambda_1=nobs:max_lambda/29:max_lambda;
        lambda_2=nobs:max_lambda/29:max_lambda;
        lambda_3=nobs:max_lambda/29:max_lambda;
    else
        
        nvar=size(terms,2);
        nobs=size(terms,1);
        max_lambda=0.01*nobs;
        lambda_1=(max_lambda*(nvar/nobs))/9:(max_lambda*(nvar/nobs))/9:(max_lambda*(nvar/nobs));
        lambda_2=(max_lambda*(nvar/nobs))/9:(max_lambda*(nvar/nobs))/9:(max_lambda*(nvar/nobs));
        lambda_3=(max_lambda*(nvar/nobs))/9:(max_lambda*(nvar/nobs))/9:(max_lambda*(nvar/nobs));
        
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
    S=(1:size(sub,2))';
    CVS = cvpartition(S,'k',8);
    for i = 1:CVS.NumTestSets
        trsub = CVS.training(i);trIdx=[S_i{find(trsub==1),1}]';
        tesub = CVS.test(i);teIdx=[S_i{find(tesub==1),1}]';
        terms_train=terms(trIdx,:);tgt_train=tgt(trIdx,1);
        terms_test=terms(teIdx ,:);tgt_test=tgt(teIdx,1);
        j=size(lm,1);
        param=cell(1,j);out=cell(1,j);fitstats_lm=cell(1,j);lambda_BIC=cell(j,1);
        lso=cell(1,j);terms_train_const=cell(1,j); highest_likelihood_i=cell(1,j);
        terms_train_in =cell(1,j);out2=cell(1,j);sig_lag=cell(1,j);nlag=cell(1,j); nonempt3=cell(1,j);
        g=cell(1,j);
        parfor parlm=1:size(lm,1)
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
            if AR==0
            if zscored==0
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
                        if zscored==0
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
        empt_features=find(lmbc_features==0);
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
        terms_test_in = terms_test(:,lasso_train(:,1)~=0);
        param_in = {all_param{1,lasso_train(1:size(all_param,2),1)~=0}};
        timings_in={all_timings{1,lasso_train(1:size(all_param,2),1)~=0}};
        minlag=1;maxlag=10;min_dof=5;
        try
            if AR==0
                fitstats_train=fitglm(terms_train_in,tgt_train);
                betas=fitstats_train.Coefficients.Estimate;
%                 g2.lse=fitstats_train.SSE;
%                  g2.aic=fitstats_train.ModelCriterion.AIC;
%                 g2.bic=fitstats_train.ModelCriterion.BIC;
                %fitstats_train.lse = sum((fitstats_train.y-fitstats_train.yhat).^2);  
                
%                 fitstats_train.aic = length(fitstats_train.yhat)*log(fitstats_train.lse/length(fitstats_train.yhat)) + 2*fitstats_train.beta_size;
%                 fitstats_train.bic = length(fitstats_train.yhat)*log(fitstats_train.lse/length(fitstats_train.yhat)) + fitstats_train.beta_size*log(length(fitstats_train.yhat));
               
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
                    if zscored==0
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
        if AR~=0
        b=betas;
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
        
        test_pred=terms_test_in*betas(nlag_tr+1:size(terms_test_in,2)+nlag_tr,1);
        g2.test_yhat=test_pred;
        g2.test_y=tgt_test;
        
        
        if AR~=0
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
        
        if zscored==0
            [test_pred_r,test_pred_pval]=corr(zscore(test_pred),zscore(tgt_test));
            if AR~=0
                [test_pred_r_lagcorct,test_pred_pval_lagcorct]=corr(zscore(test_pred),zscore(tgt_test-lag_pred));
            end
        else
            [test_pred_r,test_pred_pval]=corr(test_pred,tgt_test);
            if AR~=0
                [test_pred_r_lagcorct,test_pred_pval_lagcorct]=corr(test_pred,tgt_test-lag_pred);
            end
        end
        g2.prediction_r2=(test_pred_r)^2;
        g2.predictionpval=test_pred_pval;
        if AR~=0
            fitstats_train.prediction_r2_lagcorct=(test_pred_r_lagcorct)^2;
            fitstats_train.prediction_pval_lagcorct=test_pred_pval_lagcorct;
            fitstats_train.test_lag_pred=lag_pred;
            result(i).prediction_r2_lagcorct=fitstats_train.prediction_r2_lagcorct;
            result(i).prediction_pval_lagcorct=test_pred_pval_lagcorct;
            result(i).lag_stats=out;
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
        result(i).test_yhat=test_pred;
        result(i).test_y=tgt_test;
        if AR~=0
        result(i).prediction_r2=(test_pred_r)^2;
        result(i).prediction_pval=(test_pred_pval);
        result(i).r2=fitstats_train.rsqr;
        result(i).pval=fitstats_train.pval;
        result(i).train_yhat=fitstats_train.yhat;
        result(i).train_y=fitstats_train.y;
        result(i).bic=fitstats_train.bic;
        result(i).aic=fitstats_train.aic;
        else
        result(i).prediction_r2= g2.prediction_r2;
        result(i).prediction_pval=g2.predictionpval;
        result(i).r2=fitstats_train.Rsquared.Ordinary;
        result(i).pval=fitstats_train.Coefficients.pValue;
        result(i).train_yhat=fitstats_train.Fitted.Response;
        result(i).train_y=fitstats_train.Variables.y;
        result(i).bic=fitstats_train.ModelCriterion.BIC;
        result(i).aic=fitstats_train.ModelCriterion.AIC;
        end
        clear fitstats_train
    end
    
    
    
    
    
    %%
    %%******Do the actual fitting with the best lambda********
    if AR~=0
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
    S=(1:size(sub,2))';
    CVB = cvpartition(S,'k',8);
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
        if AR==0
            fitstats_trainB=fitglm(terms_train_in,tgt_train);
% beta_size=size(all_corct_reg,2);
            betasB=fitstats_trainB.Coefficients.Estimate;
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
                    if zscored==0
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
        if AR~=0
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
        
        if AR~=0
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
        if zscored==0
            [test_predB_r,test_predB_pval]=corr(zscore(test_predB),zscore(tgt_testB));
            if AR~=0
                [test_predB_r_lagcorct,test_predB_pval_lagcorct]=corr(zscore(test_predB),zscore(tgt_testB-lag_predB));
            end
        else
            [test_predB_r,test_predB_pval]=corr(test_predB,tgt_testB);
            if AR~=0
                [test_predB_r_lagcorct,test_predB_pval_lagcorct]=corr(test_predB,tgt_testB-lag_predB);
            end
        end
        g3.prediction_r2=(test_predB_r)^2;
        g3.prediction_pval=test_predB_pval;
        if AR~=0
            fitstats_trainB.prediction_r2_lagcorct=(test_predB_r_lagcorct)^2;
            fitstats_trainB.prediction_pval_lagcorct=test_predB_pval_lagcorct;
            fitstats_trainB.lag_predB=lag_predB;
        end
        fitstats(iB).testB_yhat=g3.testB_yhat;
        fitstats(iB).testB_y=g3.testB_y;
        fitstats(iB).prediction_r2=g3.prediction_r2;
        fitstats(iB).prediction_pval=g3.prediction_pval;
        fitstats(iB).trainB=fitstats_trainB;
        if AR==0
        fitstats(iB).rsqr=fitstats_trainB.Rsquared.Ordinary;
        fitstats(iB).bic=fitstats_trainB.ModelCriterion.BIC;
        fitstats(iB).aic=fitstats_trainB.ModelCriterion.AIC;
        fitstats(iB).pval=fitstats_trainB.Coefficients.pValue;
        fitstats(iB).beta=fitstats_trainB.Coefficients.Estimate;
        fitstats(iB).y=fitstats_trainB.Variables.y;
        fitstats(iB).yhat=fitstats_trainB.Fitted.Response;
        end
        clear fitstats_trainB
        
        % % %
        % % %     %Pred,lse, R2, tstat, p_val of uncorrected Predictors
        % % %     lagged_BOLD_fit=ols(fitstats.y,fitstats.xmat(:,nlag+1:end));
        % % %     lagged_BOLDpred_fused=lagged_BOLD_fit.yhat;
        % % %     fitstats.lagged_yhat=lagged_BOLD_fit.yhat;
        % % %     meanresult = mean(lagged_BOLD_fit.y,1);
        % % %     re_lagged=sum((lagged_BOLD_fit.y-meanresult).^2);
        % % %     lagged_lse_fused = sum((lagged_BOLD_fit.y-lagged_BOLDpred_fused).^2); %sum least-squares error (SSR)
        % % %     fitstats.lagged_BOLDr2_fused = 1-(lagged_lse_fused/re_lagged);
        % % %     rt1=fitstats.lagged_BOLDr2_fused;
        % % %     dof1=size(lagged_BOLD_fit.y,1)-2;
        % % %     r_stat1= (rt1.*(sqrt(dof1)))./sqrt(1-rt1.^2);
        % % %     fitstats.lagged_BOLD_pval=tcdf(r_stat1,dof1,'upper');
        % % %     fitstats.lagged_BOLD_lse=lagged_lse_fused;
        % % %     fitstats.lagged_BOLD_tstat=r_stat1;
        % % %     fitstats.lagged_BOLD_bic = length(fitstats.lagged_yhat)*log(fitstats.lagged_BOLD_lse/length(fitstats.lagged_yhat)) + size(lagged_BOLD_fit.beta,1)*log(length(fitstats.lagged_yhat));
        % % %
        % % %     %Pred,lse, R2,tstat, p_val of lag correction terms
        % % %     lagterms_BOLD_fit=ols(fitstats.y,fitstats.xmat(:,1:nlag));
        % % %     %lag_betas=fitstats.beta(1:nlag,1);
        % % %     fitstats.lagterms_pred=lagterms_BOLD_fit.yhat;
        % % %     meanresult2 = mean(lagterms_BOLD_fit.y,1);
        % % %     re_lagterms=sum((lagterms_BOLD_fit.y-meanresult2).^2);
        % % %     lse_lagterms=sum((lagterms_BOLD_fit.y-lagterms_BOLD_fit.yhat).^2);
        % % %     lagterms_r2=1-(lse_lagterms/re_lagterms);
        % % %     fitstats.lagterms_r2=lagterms_r2;
        % % %     rt2=fitstats.lagterms_r2;
        % % %     if rt2>1 || rt2<-1
        % % %         rt2=0;
        % % %     end
        % % %     r_stat2= (rt2.*(sqrt(dof1)))./sqrt(1-rt2.^2);
        % % %     fitstats.lagterms_tstat=r_stat2;
        % % %     fitstats.lagterms_pval=tcdf(r_stat2,dof1,'upper');
        % % %     fitstats.lagterms_lse=lse_lagterms;
        % % %     fitstats.lagterms_bic = length(lagterms_BOLD_fit.yhat)*log(fitstats.lagterms_lse/length(lagterms_BOLD_fit.yhat)) + size(lagterms_BOLD_fit.beta,1)*log(length(lagterms_BOLD_fit.yhat));
        % % %
        % % %     %save results in the output stucture
        % % %     r2=fitstats.rsqr;
        % % %     BOLDpred=fitstats.yhat;
        % % %     fitstats.lse=sum((fitstats.y-fitstats.yhat).^2);
        % % %     fitstats.model_tstat= (r2.*(sqrt(dof1)))./sqrt(1-r2.^2);
        % % %     fitstats.model_pval=tcdf(fitstats.model_tstat,dof1,'upper');
        % % %
        % % %     lse=fitstats.lse;
        % % %     if isempty(param_in)
        % % %         pval_in=1;
        % % %     else
        % % %         for get_pval=1+nlag:size(param_in,2)+nlag
        % % %             p_value(get_pval-nlag,1)=fitstats.tprob(get_pval);
        % % %         end
        % % %         pval_in=p_value;
        % % %     end
    end
    %*******Do Lasso Regression*********
    cd '/Applications/spm12/toolbox/spams-matlab-v2.6'
    start_spams
    addpath(genpath('/Applications/spm12/toolbox/jplv7'))
    fprintf('\n Dependent explanatory variables & autocorrelated observations \n Using Fused lasso regression \n')
elseif ~isempty(strfind(dependence,'lasso')) && isempty(strfind(dependence,'fused'))
    if ~isempty(strfind(dependence,'ridge'))
        alpha=0.2;
    elseif ~isempty(strfind(dependence,'-el'))
        alpha=0.5;
    else
        alpha=1;
    end
    if pooled==1
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
            %         if zscored==0
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
    if zscored==0 %adds one constant term to the explanatory varaiables (irrespective of the number of subjects)
        %to calculate BIC and AIC properly
        b=betas(1:size(all_corct_reg,2)+1);
    end
    
    %Plot Lasso predictor weights*****
    if pooled==0
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
    
elseif strcmp(dependence,'fmincon')
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
    result.b = b;
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
if isempty(strfind(dependence,'fused'))
    result.b = b;
    result.BOLDpred = BOLDpred;
    result.r2 = r2;
    result.bic = length(BOLDpred)*log(lse/length(BOLDpred)) + length(b)*log(length(BOLDpred));
    result.aic = length(BOLDpred)*log(lse/length(BOLDpred)) + 2*length(b);
    result.paramNames = param_in;
    result.pval=pval_in;
    result.timings=timings_in;
    result.model_in=0; %dummy
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