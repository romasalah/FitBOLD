dest_dir='/Volumes/SSD/fitBOLD12-(CV(IFLR-OLSAR(n)))x2-IS';
GLMs={'/Volumes/SSD/FFX_Rutledge_minEV_14aug2018_ISIfinal_Ch2/Subjects'}; %SPM's GLM directory
%remember to remove the "/" at the end of the string.
social={'n-soc','active','passive'}; % condition_before='active'; set it to zero to disable.
models={'IS'}; %model='CS1'; %name of the model
pooled=[1]; %pooled=1; %pooling all the subjects into one vector and modelling them in one time. %set to 0 to do each individually or 1 to pool.
zscored=[1];
combine=[0];
%combine=0; %Combine BOLD from different ROIs and consider them as one ROI. e.g. Collapse_Bilaterality
%you have to specifiy the Combined ROIs in a vector under each model in the
shifts={0}; %shift=0; %specifies of HRFs are shifted to sync with the behavioural task. set to "standard shift" as a string
%to do a shift of 2 scans for all timings. set to 0 to estimate the shift
%for each timing individually
scans={'behav'}; %which_scans='behav'; %Use Brain scans at the time of the modelled behavior only
%or scans from the whole experiment.
social_ROIs=[0];
MA=0; %time interval for calculating the moving avergae. 1 time point is 10-15 seconds
behav=1; %modelling happiness as a behavioural measure only. %behav=1; %if You want to model the behavioural measure, set to 1
%if you want to model BOLD signal in ROIs at the time of the behavioural task set to 0
equal_regressors=30; %per session %set this number to number of trials you want to model. forces the corr_BOLD function to make all BOLD vectors of the same length as the target BOLD or bahevioural vector.
dependence='lasso-fused'; %options 'fmincon', 'lasso', 'stepwise', 'lasso-elastic','lasso-ridge','lasso-fused'
bimodal=0;
disk_parallel=8; % the number of workers in multi-threading the reading process from disk. set to 0 to disable multi-threading.
use_4D=1;
AR=0; FIR=1;
% if pooled==1 dependence='lasso-elastic';end
%the current function gets the m
sub=[ 44 62 63 65 69 70 71 74 77 79 81 83 84 90 92 93 94 95 96 98 101 102 103 104 106 110 120 124 125 126 129 136 155 166 171 173 174 175 176 177 ];
% sub=sub(:,1:7);
%sub=[175]
for glm_i=1:size(GLMs,2)
    for pooled_i= 1:size(pooled,2)
        for models_i=1:size(models,2)
            for soc_stat= 2:2
                if soc_stat==1
                    condition_before=0;
                    for scans_i=1:size(scans,2)
                        for social_ROIs_i=1:size(social_ROIs,2)
                            for shifts_i=1:size(shifts,2)
                                for combine_i=1:size(combine,2)
                                    for zscore_i= 1:size(zscored,2)
                                        master_fit_BOLD(GLMs{1,glm_i},models{1,models_i},scans{1,scans_i},behav,...
                                            shifts{shifts_i},combine(1,combine_i),sub,pooled(1,pooled_i),...
                                            equal_regressors,dest_dir,condition_before,zscored(zscore_i),...
                                            social_ROIs(social_ROIs_i),dependence,bimodal,MA,AR,disk_parallel,use_4D,FIR)
                                        close all %keep this to avoid hundreds of figures on your screen
                                    end
                                end
                            end
                        end
                    end
                else
                    for soc_typ=1:size(social,2)
                        condition_before=social{1,soc_typ};
                        for scans_i=1:size(scans,2)
                            for social_ROIs_i=1:size(social_ROIs,2)
                                for shifts_i=1:size(shifts,2)
                                    for combine_i=1:size(combine,2)
                                        for zscore_i= 1:size(zscored,2)
                                            master_fit_BOLD(GLMs{1,glm_i},models{1,models_i},scans{1,scans_i},behav,...
                                                shifts{shifts_i},combine(1,combine_i),sub,pooled(1,pooled_i),...
                                                equal_regressors,dest_dir,condition_before,zscored(zscore_i),...
                                                social_ROIs(social_ROIs_i),dependence,bimodal,MA,AR,disk_parallel,use_4D,FIR)
                                            close all %keep this to avoid hundreds of figures on your screen
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end