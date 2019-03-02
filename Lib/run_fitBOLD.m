%% User Specified
clear fit
%remember to remove the "/" at the end of the string.
fit.social={'active','passive','n-soc'}; % condition_before='active'; set it to zero to disable.
fit.models={'RL'}; %model='CS1'; %name of the model
%or scans from the whole experiment.
fit.designs='/Volumes/EVO 860 2TB/Designs';
fit.preextracted={...
    [fit.designs filesep 'FFX_CRminEV_IS_behav_behav_pool_Uni_trushft_active_Zscr_stepwise_bootstrap_NSocROI___40Subjects.mat'],...
    [fit.designs filesep 'FFX_CRminEV_IS_behav_behav_pool_Uni_trushft_passive_Zscr_stepwise_bootstrap_NSocROI___40Subjects.mat'],...
    [fit.designs filesep 'FFX_CRminEV_IS_behav_behav_pool_Uni_trushft_n-soc_Zscr_stepwise_bootstrap_NSocROI___40Subjects.mat'],...
    [fit.designs filesep 'FFX_CRminEV_RL_behav_all_pool_Uni_trushft_active_Zscr_stepwise_bootstrap_NSocROI___40Subjects.mat'],...
    [fit.designs filesep 'FFX_CRminEV_RL_behav_all_pool_Uni_trushft_passive_Zscr_stepwise_bootstrap_NSocROI___40Subjects.mat'],...
    [fit.designs filesep 'FFX_CRminEV_RL_behav_all_pool_Uni_trushft_n-soc_Zscr_stepwise_bootstrap_NSocROI___40Subjects.mat']};
fit.preextracted={...
    [fit.designs filesep 'FFX_CRminEV_RL_behav_all_pool_Uni_trushft_active_Zscr_stepwise_bootstrap_NSocROI___40Subjects.mat']};
fit.preextracted=0;
fit.disk_parallel=16; % the number of workers in multi-threading the reading process from disk. set to 0 to disable multi-threading.
fit.cpu_parallel=16;
fit.use_4D=1;
fit.social_ROIs=[0,1];
fit.use_subROIs= [1];
fit.sub=[ 44 62 63 65 69 70 71 74 77 79 81 83 84 90 92 93 94 95 96 98 101 102 103 104 106 110 120 124 125 126 129 136 155 166 171 173 174 175 176 177 ];
%% Clustering options
fit.cluster.mink=2;
fit.cluster.distance='Correlation';
fit.cluster.minvol=50;
fit.cluster.maxk=10;
%%
fit.resample.method='bootstrap';
fit.regress.validationfolds=5;
%% Almost Stable
fit.dest_dir='/Volumes/EVO 860 2TB/Res/fitBOLD41-testcluster';
fit.GLMs={'/Volumes/EVO 860 2TB/FFX_CRminEV_09Jan19'}; %SPM's GLM directory
fit.maskdir='/Volumes/EVO 860 2TB/fitBOLD/Masks';
fit.subROIsdir=[fit.dest_dir filesep 'subROIs'];
if ~exist(fit.subROIsdir,'dir'); mkdir(fit.subROIsdir); end
fit.dependence='stepwise'; %options 'fmincon', 'lasso', 'stepwise', 'lasso-elastic','lasso-ridge','lasso-fused, 'OLS''
fit.bimodal=0;
fit.AR=0;
fit.regress.incriterion='adjrsquared';
fit.regress.incriterion2='aic';
fit.regress.outcriterion='sse';
fit.regress.output='linear';
fit.regress.threshold=0.01;
fit.regress.useCImu=0;
fit.FIR=1; %remove low signal drift in BOLD using finite impulse response filtering. 
fit.clustered=0;
fit.pooled=[1]; %pooled=1; %pooling all the subjects into one vector and modelling them in one time. %set to 0 to do each individually or 1 to pool.
fit.zscored=[1];
fit.combine=[0];
%combine=0; %Combine BOLD from different ROIs and consider them as one ROI. e.g. Collapse_Bilaterality
%you have to specifiy the Combined ROIs in a vector under each model in the
fit.shifts={1}; %shift=0; %specifies of HRFs are shifted to sync with the behavioural task. set to "standard shift" as a string
%to do a shift of 2 scans for all timings. set to 0 to estimate the shift
%for each timing individually
fit.ME='/Users/NEMO_admin/Google Drive/Master(shared)/Data/SoDec-fMRI_fitBOLD.xlsx';
%% Never Change
fit.regress.failed=0;
fit.MA=0; %time interval for calculating the moving avergae. 1 time point is 10-15 seconds
fit.behav=1; %modelling happiness as a behavioural measure only. %behav=1; %if You want to model the behavioural measure, set to 1
%if you want to model BOLD signal in ROIs at the time of the behavioural task set to 0
fit.equal_regressors=30; %per session %set this number to number of trials you want to model. forces the corr_BOLD function to make all BOLD vectors of the same length as the target BOLD or bahevioural vector.
fit.use_marsbar=1;
%% Start
if ~isa(fit.preextracted,'char') && ~isa(fit.preextracted,'cell')
    routine=1;
else
    dbstop if error
    conds=fit.preextracted;
    for i=1:size(fit.preextracted,2)
    fit.preextracted=conds{1,i};
    res=main_fitBOLD(fit);
    close all
    end
    return
end
if routine==1
    dbstop if error
res=Batch_fitBOLD(fit);
end

