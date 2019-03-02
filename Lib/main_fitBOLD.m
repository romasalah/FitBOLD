function [results,fitstats,modelstat]= main_fitBOLD(fit)
clear results
if ~isa(fit.preextracted,'char') && ~isa(fit.preextracted,'cell')
%% Defining the models***~
[fit]=spec_fitBOLD(fit);
%% Design the Model***
[fit]=design_fitBOLD(fit);
else
% load preextracted model
load(fit.preextracted);
fit.regress.threshold=0.05./size(fit.design,2);
clearvars -except fit
end
%% Model name and settings***
modelName= name_fitBOLD(fit);
%% Do the fitting and Draw results
if fit.pooled==1
    [results,~,fitstats]=fit_fitBOLD(fit);
    draw_fitBOLDp(results,fitstats,fit)
    
elseif fit.pooled==0
    sub_end=size(fit.sub,2);
    fit2=cell(sub_end,1);
    fitstats=cell(sub_end,1);
    hbar = parfor_progressbar(sub_end,'Doing Subjects in parallel...'); %create the progress bar
    parfor (parsub_id= 1:sub_end,fit.disk_parallel)
        fit2{parsub_id}=fit;
        fit2{parsub_id}.sub=disp(fit.sub(parsub_id));
        fprintf(['   \n Doing Subject: ' num2str(fit.sub(parsub_id)) '  ']);
        [results(parsub_id,1),~,~,fitstats{parsub_id,1}]=fit_fitBOLD(fit2{parsub_id});
        switch fit.dependence
            case 'fused-lasso'
                draw_fitBOLD(results(parsub_id,1),fitstats{parsub_id,1},fit,lambda_BIC);
            otherwise
                draw_fitBOLD(results(parsub_id,1),fitstats{parsub_id,1},fit,0);
                hbar.iterate(1);
        end
    end
       
    close(hbar);
    [rescell,sig_results,rescell_Bmodel]=orgres_fitBOLD(results,fitstats,fit,fit.sub,fit.pooled); %#ok<ASGLU>
%% Do group stats
    [modelstat,Best_modelstat]=groupstats_fitBOLD(fit,fitstats,results,rescell,rescell_Bmodel); %#ok<ASGLU>
end
%% Saving data
folder=fit.dest_dir;  % Desination folder
results_name=[modelName '_results.mat'];
MSg_name=[modelName '_model_groupstats.mat'];
BFg_name=[modelName '_Best_features_groupstats.mat'];
fitstats_name=[modelName '_fitstats.mat'];
BF_name=[modelName '_Best_features_stats.mat'];
FI_name=[modelName '_feature_iterations_stats.mat'];
fprintf('saving the results matrix and group statistics \n')
save(fullfile(folder,results_name),'results')

if  fit.pooled==1
    if ~isempty(strfind(fit.dependence,'fused'))
        save(fullfile(folder,MSg_name),'modelstat')
        save(fullfile(folder,BFg_name),'Best_modelstat')
    end
    try save(fullfile(folder,BF_name),'fitstats','-v7.3'); end
elseif fit.pooled==0
    save(fullfile(folder,fitstats_name),'rescell')
else
    try save(fullfile(folder,BF_name),'fitstats'); end
end

%% Iterative Regression after Clustering the subjects
if fit.pooled==0 && fit.clustered==1
    iterativeSubClusteredFit_fitBOLD(fit,results,n)
end
end