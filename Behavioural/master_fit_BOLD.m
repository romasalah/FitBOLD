function [results,shifts,fitstats,modelstat]= master_fit_BOLD(glm_dir,model,which_scans,behav,shift,combine,...
    sub,pooled,equal_regressors,dest_dir,condition_before,zscored,social_ROIs,dependence,bimodal,MA,AR,disk_parallel,use_4D,FIR)
clear results
%****Defining the models***~
whatcombine=0;
howcombine=0;
switch model
    case 'RL'
        fprintf('\n Using Reinforcement Learning ROI \n ');
        if strcmp(which_scans,'behav')
            if behav==1
                if combine
                    whatcombine=[1 2; 3 4; 5 6; 7 8];%you need to specify this vector, each group of ROIs must be in the same row.
                    howcombine='mean'; %possible options for now are mean or median.
                end
                fprintf(' \n predicting a behavioural report from brain scans at time of the report');
                ROIs={'Amygdala Lt_roi.mat','Amygdala Rt_roi.mat',...
                    'N.Ac Lt_roi.mat','N.Ac Rt_roi.mat',...
                    'ACC Lt_roi.mat','ACC Rt_roi.mat',...
                    'Insula Lt_roi.mat','Insula Rt_roi.mat',...
                    'Subgenual vmPFC & mOFC Bi_roi.mat'};
                timings=cell(1,size(ROIs,2));
                for at_time_fill=1:size(ROIs,2)
                    timings{1,at_time_fill}='happiness';
                    if condition_before
                        timings{2,at_time_fill}= condition_before;
                    end
                end
            end
        elseif strcmp(which_scans,'behav')
            if behav==0
                fprintf(' \n Cannot use Brain scans at the time of the behvioural report to model ROIs at the same time \n the procedure is not valid statistically')
            end
        elseif strcmp(which_scans,'all')
            fprintf(' \n Using all the brain scans to model either another ROI or a behavioural report');
            if combine
                whatcombine=[1 2; 3 4; 5 6; 7 8;9 10];%you need to specify this vector, each group of ROIs must be in the same row.
                howcombine='mean'; %possible options for now are mean or median.
            end
            if social_ROIs==0
                ROIs={'Amygdala Lt_roi.mat','Amygdala Rt_roi.mat',...
                    'ACC Lt_roi.mat','ACC Rt_roi.mat',...
                    'N.Ac Lt_roi.mat','N.Ac Rt_roi.mat',...
                    'N.Ac Lt_roi.mat','N.Ac Rt_roi.mat',...
                    'ACC Lt_roi.mat','ACC Rt_roi.mat'};
                timings={'Min EV','Min EV','Min EV',...
                    'Min EV','Min EV','Min EV',...
                    'PE','PE','PE','PE'};
            elseif social_ROIs==1
                ROIs={'Amygdala Lt_roi.mat','Amygdala Rt_roi.mat',...
                    'ACC Lt_roi.mat','ACC Rt_roi.mat',...
                    'N.Ac Lt_roi.mat','N.Ac Rt_roi.mat',...
                    'N.Ac Lt_roi.mat','N.Ac Rt_roi.mat',...
                    'ACC Lt_roi.mat','ACC Rt_roi.mat',...
                    'BA5m Lt_roi.mat','BA5m Rt_roi.mat',...
                    'BA7m Lt_roi.mat','BA7m Rt_roi.mat',...
                    'BA23d_Lt_roi.mat','BA23d_Rt_roi.mat',...
                    'BA44v_Lt_roi.mat','BA44v_Rt_roi.mat',...
                    'BA31 Lt_roi.mat','BA31 Rt_roi.mat',...
                    'BA9m Lt_roi.mat','BA9m Rt_roi.mat',...
                    'BA10m Lt_roi.mat','BA10m Rt_roi.mat',...
                    'BA5m Lt_roi.mat','BA5m Rt_roi.mat',...
                    'BA7m Lt_roi.mat','BA7m Rt_roi.mat',...
                    'BA23d_Lt_roi.mat','BA23d_Rt_roi.mat',...
                    'BA44v_Lt_roi.mat','BA44v_Rt_roi.mat',...
                    'BA31 Lt_roi.mat','BA31 Rt_roi.mat',...
                    'BA9m Lt_roi.mat','BA9m Rt_roi.mat',...
                    'BA10m Lt_roi.mat','BA10m Rt_roi.mat'};
                
                timings={'Min EV','Min EV','Min EV',...
                    'Min EV','Min EV','Min EV',...
                    'PE','PE','PE','PE',...
                    'Min EV','Min EV','Min EV','Min EV','Min EV','Min EV','Min EV',...
                    'Min EV','Min EV','Min EV','Min EV','Min EV','Min EV','Min EV',...
                    'PE','PE','PE','PE','PE','PE','PE',...
                    'PE','PE','PE','PE','PE','PE','PE'};
            end
            for at_time_fill=1:size(ROIs,2)
                if condition_before
                    timings{2,at_time_fill}= condition_before;
                end
            end
        end
        TargetBOLD='Subgenual vmPFC & mOFC Bi_roi.mat';
        Target_timing{1,1}='happiness';
        if condition_before
            Target_timing{2,1}=condition_before;
        end
        
        
    case 'IS'
        
        fprintf(' \n Using Incentive Salience ROI \n ');
        if strcmp(which_scans,'behav')
            if behav==1
                if combine
                    whatcombine=[1 2; 3 4; 5 6;7 8;9 10];%you need to specify this vector, each group of ROIs must be in the same row.
                    howcombine='mean'; %possible options for now are mean or median.
                end
                fprintf('predicting a behavioural report from brain scans at time of the report');
                ROIs={'VPe Ant Lt_roi.mat','VPe Ant Rt_roi.mat',...
                    'VPe Post Lt_roi.mat','VPe Post Rt_roi.mat',...
                    'N.Ac Med Ant Lt_roi.mat','N.Ac Med Ant Rt_roi.mat',...
                    'BA32sg_Lt_roi.mat','BA32sg_Rt_roi.mat',...
                    'BA13_Lt_roi.mat','BA13_Rt_roi.mat'};
                timings=cell(1,size(ROIs,2));
                for at_time_fill=1:size(ROIs,2)
                    timings{1,at_time_fill}='happiness';
                    if condition_before
                        timings{2,at_time_fill}= condition_before;
                    end
                end
            end
        elseif strcmp(which_scans,'behav')
            if behav==0
                fprintf(' \n Cannot use Brain scans at the time of the behvioural report to model ROIs at the same time \n the procedure is not valid statistically')
            end
        elseif strcmp(which_scans,'all')
            fprintf(' \n Using all the brain scans to model either another ROI or a behavioural report');
            if combine
                whatcombine=[1 2; 3 4; 5 6;7 8;9 10];%you need to specify this vector, each group of ROIs must be in the same row.
                howcombine='mean'; %possible options for now are mean or median.
            end
            ROIs={ 'N.Ac Med Ant Lt_roi.mat','N.Ac Med Ant Rt_roi.mat',...
                'VPe Ant Lt_roi.mat','VPe Ant Rt_roi.mat',...
                'ACC Lt_roi.mat','ACC Rt_roi.mat',...
                'N.Ac Lt_roi.mat','N.Ac Rt_roi.mat',...
                'ACC Lt_roi.mat','ACC Rt_roi.mat'};
            timings={'Min EV','Min EV',...
                'Min EV','Min EV',...
                'PE','PE',...
                'PE','PE'};
            for at_time_fill=1:size(ROIs,2)
                if condition_before
                    timings{2,at_time_fill}= condition_before;
                end
            end
        end
        
        TargetBOLD='Subgenual vmPFC & mOFC Bi_roi.mat'; %respecify that
        Target_timing{1,1}='happiness';
        if condition_before
            Target_timing{2,1}=condition_before;
        end
        
    case 'CS'
        fprintf(' \n Using a Custom built combination of ROIs \n ')
        if strcmp(which_scans,'behav')
            if behav==1
                if combine
                    whatcombine=[1 2; 3 4; 5 6; 7 8; 9 10];%you need to specify this vector, each group of ROIs must be in the same row.
                    howcombine='mean'; %possible options for now are mean or median.
                end
                fprintf('\n predicting a behavioural report from brain scans at time of the report');
                ROIs={'N.Ac Med Ant Lt_roi.mat','N.Ac Med Ant Rt_roi.mat'...
                    'Insula Lt_roi.mat','Insula Rt_roi.mat',...
                    'BA32p_Lt_roi.mat','BA32p_Rt_roi.mat',...
                    'VPe Ant Lt_roi.mat','VPe Ant Rt_roi.mat',...
                    'VPe Post Lt_roi.mat','VPe Post Rt_roi.mat'};
                
                
                %                 'ACC Lt_roi.mat','ACC Rt_roi.mat',...
                %                 'Amygdala Lt_roi.mat','Amygdala Rt_roi.mat',...
                %                 'BA24rv_Rt_roi.mat','BA24rv_Lt_roi.mat',...
                %                 'BA32sg_Lt_roi.mat','BA32sg_Rt_roi.mat',...
                %                 'BA13_Lt_roi.mat','BA13_Rt_roi.mat'};
                timings=cell(1,size(ROIs,2));
                for at_time_fill=1:size(ROIs,2)
                    timings{1,at_time_fill}='happiness';
                    %                 if condition_before
                    %                     timings{2,at_time_fill}= condition_before;
                    %                 end
                end
            end
        elseif strcmp(which_scans,'behav')
            if behav==0
                fprintf('\n Cannot use Brain scans at the time of the behvioural report to model ROIs at the same time \n the procedure is not valid statistically')
            end
        elseif strcmp(which_scans,'all')
            fprintf('\n Using all the brain scans to model either another ROI or a behavioural report');
            ROIs={'N.Ac Med Ant Lt_roi.mat','N.Ac Med Ant Rt_roi.mat'...
                'Amygdala Lt_roi.mat','Amygdala Rt_roi.mat'...
                'ACC Lt_roi.mat','ACC Rt_roi.mat',...
                'Pos RPE -16_0_70_roi.mat','Pos RPE -20_2_-22_roi.mat'...
                'Pos RPE -23_2_0_roi.mat','Pos RPE 4_8_38_roi.mat',...
                'Pos RPE 25_3_0_roi.mat','Pos RPE 56_12_18_roi.mat',...
                'Neg RPE -24_-4_6_roi.mat','Neg RPE -30_14_64_roi.mat',...
                'Neg RPE -30_28_6_roi.mat','Neg RPE -42_-8_48_roi.mat',...
                'Neg RPE 4_2_60_roi.mat'};
            timings={'Min EV','Min EV',...
                'Min EV','Min EV',...
                'Min EV','Min EV',...
                'PE Pos','PE Pos','PE Pos','PE Pos','PE Pos','PE Pos'...
                'PE Neg','PE Neg','PE Neg','PE Neg','PE Neg'};
            for at_time_fill=1:size(ROIs,2)
                if condition_before
                    timings{2,at_time_fill}= condition_before;
                end
            end
        end
        TargetBOLD='Subgenual vmPFC & mOFC Bi_roi.mat'; %respecify that
        Target_timing{1,1}='happiness';
        if condition_before
            Target_timing{2,1}=condition_before;
        end
        
        
    case 'RLIS'
        
        fprintf(' \n Using Reinforcement learning and Incentive Salience ROI \n ');
        if strcmp(which_scans,'behav')
            if behav==0
                fprintf(' \n Cannot use Brain scans at the time of the behvioural report to model ROIs at the same time \n the procedure is not valid statistically')
            end
        elseif strcmp(which_scans,'all')
            fprintf(' \n Using all the brain scans to model either another ROI or a behavioural report');
            if combine
                whatcombine=[1 2; 3 4; 5 6;7 8;9 10];%you need to specify this vector, each group of ROIs must be in the same row.
                howcombine='mean'; %possible options for now are mean or median.
            end
            if social_ROIs==0
                ROIs={'Amygdala Lt_roi.mat','Amygdala Rt_roi.mat',...
                    'ACC Lt_roi.mat','ACC Rt_roi.mat',...
                    'N.Ac Med Ant Lt_roi.mat','N.Ac Med Ant Rt_roi.mat',...
                    'N.Ac Lt_roi.mat','N.Ac Rt_roi.mat',...
                    'ACC Lt_roi.mat','ACC Rt_roi.mat'...
                    'N.Ac Med Ant Lt_roi.mat','N.Ac Med Ant Rt_roi.mat'...
                    'Insula Lt_roi.mat','Insula Rt_roi.mat',...
                    'BA32sg_Lt_roi.mat','BA32sg_Rt_roi.mat',...
                    'VPe Ant Lt_roi.mat','VPe Ant Rt_roi.mat',...
                    'VPe Post Lt_roi.mat','VPe Post Rt_roi.mat'};
                timings={'Min EV','Min EV','Min EV',...
                    'Min EV','Min EV','Min EV',...
                    'PE','PE','PE','PE'...
                    'happiness','happiness',...
                    'happiness','happiness',...
                    'happiness','happiness',...
                    'happiness','happiness',...
                    'happiness','happiness'};
            elseif social_ROIs==1
                ROIs={'Amygdala Lt_roi.mat','Amygdala Rt_roi.mat',...
                    'ACC Lt_roi.mat','ACC Rt_roi.mat',...
                    'N.Ac Lt_roi.mat','N.Ac Rt_roi.mat',...
                    'N.Ac Lt_roi.mat','N.Ac Rt_roi.mat',...
                    'ACC Lt_roi.mat','ACC Rt_roi.mat',...
                    'BA5m Lt_roi.mat','BA5m Rt_roi.mat',...
                    'BA7m Lt_roi.mat','BA7m Rt_roi.mat',...
                    'BA23d_Lt_roi.mat','BA23d_Rt_roi.mat',...
                    'BA44v_Lt_roi.mat','BA44v_Rt_roi.mat',...
                    'BA31 Lt_roi.mat','BA31 Rt_roi.mat',...
                    'BA9m Lt_roi.mat','BA9m Rt_roi.mat',...
                    'BA10m Lt_roi.mat','BA10m Rt_roi.mat',...
                    'BA5m Lt_roi.mat','BA5m Rt_roi.mat',...
                    'BA7m Lt_roi.mat','BA7m Rt_roi.mat',...
                    'BA23d_Lt_roi.mat','BA23d_Rt_roi.mat',...
                    'BA44v_Lt_roi.mat','BA44v_Rt_roi.mat',...
                    'BA31 Lt_roi.mat','BA31 Rt_roi.mat',...
                    'BA9m Lt_roi.mat','BA9m Rt_roi.mat',...
                    'BA10m Lt_roi.mat','BA10m Rt_roi.mat'...
                    'N.Ac Med Ant Lt_roi.mat','N.Ac Med Ant Rt_roi.mat'...
                    'Insula Lt_roi.mat','Insula Rt_roi.mat',...
                    'BA32p_Lt_roi.mat','BA32p_Rt_roi.mat',...
                    'VPe Ant Lt_roi.mat','VPe Ant Rt_roi.mat',...
                    'VPe Post Lt_roi.mat','VPe Post Rt_roi.mat'};
                
                timings={'Min EV','Min EV','Min EV',...
                    'Min EV','Min EV','Min EV',...
                    'PE','PE','PE','PE',...
                    'Min EV','Min EV','Min EV','Min EV','Min EV','Min EV','Min EV',...
                    'Min EV','Min EV','Min EV','Min EV','Min EV','Min EV','Min EV',...
                    'PE','PE','PE','PE','PE','PE','PE',...
                    'PE','PE','PE','PE','PE','PE','PE'...
                    'happiness','happiness',...
                    'happiness','happiness',...
                    'happiness','happiness',...
                    'happiness','happiness',...
                    'happiness','happiness'};
            end
            for at_time_fill=1:size(ROIs,2)
                if condition_before
                    timings{2,at_time_fill}= condition_before;
                end
            end
        end
        
        TargetBOLD='Subgenual vmPFC & mOFC Bi_roi.mat'; %respecify that
        Target_timing{1,1}='happiness';
        if condition_before
            Target_timing{2,1}=condition_before;
        end
        
    case 'CS2' %time of options shown and time of rating of happiness are consdidered timings for activations of regions encoding pleasure directly
        %PE is not divided into positive and negative
        
        fprintf(' \n Using Reinforcement learning and Incentive Salience ROI \n ');
        if strcmp(which_scans,'behav')
            if behav==0
                fprintf(' \n Cannot use Brain scans at the time of the behvioural report to model ROIs at the same time \n the procedure is not valid statistically')
            end
        elseif strcmp(which_scans,'all')
            fprintf(' \n Using all the brain scans to model either another ROI or a behavioural report');
            if combine
                whatcombine=[1 2; 3 4; 5 6;7 8;9 10];%you need to specify this vector, each group of ROIs must be in the same row.
                howcombine='mean'; %possible options for now are mean or median.
            end
            if social_ROIs==0
                ROIs={'Amygdala Lt_roi.mat','Amygdala Rt_roi.mat',...%options shown
                    'ACC Lt_roi.mat','ACC Rt_roi.mat',...
                    'N.Ac Med Ant Lt_roi.mat','N.Ac Med Ant Rt_roi.mat'...
                    'VPe Ant Lt_roi.mat','VPe Ant Rt_roi.mat',...
                    'VPe Post Lt_roi.mat','VPe Post Rt_roi.mat'...
                    'Insula Lt_roi.mat','Insula Rt_roi.mat',...
                    'BA32sg_Lt_roi.mat','BA32sg_Rt_roi.mat',...
                    'N.Ac Lt_roi.mat','N.Ac Rt_roi.mat',...%Results shown
                    'ACC Lt_roi.mat','ACC Rt_roi.mat'...
                    'Amygdala Lt_roi.mat','Amygdala Rt_roi.mat',...%rating happiness
                    'N.Ac Med Ant Lt_roi.mat','N.Ac Med Ant Rt_roi.mat'...
                    'ACC Lt_roi.mat','ACC Rt_roi.mat',...
                    'Insula Lt_roi.mat','Insula Rt_roi.mat',...
                    'BA32sg_Lt_roi.mat','BA32sg_Rt_roi.mat',...
                    'VPe Ant Lt_roi.mat','VPe Ant Rt_roi.mat',...
                    'VPe Post Lt_roi.mat','VPe Post Rt_roi.mat'};
                timings={'Min EV','Min EV','Min EV',...
                    'Min EV','Min EV','Min EV',...
                    'Min EV','Min EV',...
                    'Min EV','Min EV',...
                    'Min EV','Min EV',...
                    'Min EV','Min EV',...
                    'PE','PE','PE','PE'...
                    'happiness','happiness',...
                    'happiness','happiness',...
                    'happiness','happiness',...
                    'happiness','happiness',...
                    'happiness','happiness',...
                    'happiness','happiness',...
                    'happiness','happiness'};
            elseif social_ROIs==1
                ROIs={'Amygdala Lt_roi.mat','Amygdala Rt_roi.mat',...%options shown
                    'ACC Lt_roi.mat','ACC Rt_roi.mat',...
                    'N.Ac Med Ant Lt_roi.mat','N.Ac Med Ant Rt_roi.mat'...
                    'VPe Ant Lt_roi.mat','VPe Ant Rt_roi.mat',...
                    'VPe Post Lt_roi.mat','VPe Post Rt_roi.mat'...
                    'Insula Lt_roi.mat','Insula Rt_roi.mat',...
                    'BA32sg_Lt_roi.mat','BA32sg_Rt_roi.mat',...
                    'BA5m Lt_roi.mat','BA5m Rt_roi.mat',... %Social ROIs/options shown
                    'BA7m Lt_roi.mat','BA7m Rt_roi.mat',...
                    'BA23d_Lt_roi.mat','BA23d_Rt_roi.mat',...
                    'BA44v_Lt_roi.mat','BA44v_Rt_roi.mat',...
                    'BA31 Lt_roi.mat','BA31 Rt_roi.mat',...
                    'BA9m Lt_roi.mat','BA9m Rt_roi.mat',...
                    'BA10m Lt_roi.mat','BA10m Rt_roi.mat',...
                    'N.Ac Lt_roi.mat','N.Ac Rt_roi.mat',...%Results shown
                    'ACC Lt_roi.mat','ACC Rt_roi.mat'...
                    'BA5m Lt_roi.mat','BA5m Rt_roi.mat',...%Social ROIs/results shown
                    'BA7m Lt_roi.mat','BA7m Rt_roi.mat',...
                    'BA23d_Lt_roi.mat','BA23d_Rt_roi.mat',...
                    'BA44v_Lt_roi.mat','BA44v_Rt_roi.mat',...
                    'BA31 Lt_roi.mat','BA31 Rt_roi.mat',...
                    'BA9m Lt_roi.mat','BA9m Rt_roi.mat',...
                    'BA10m Lt_roi.mat','BA10m Rt_roi.mat'...
                    'Amygdala Lt_roi.mat','Amygdala Rt_roi.mat',...%Happiness rating
                    'ACC Lt_roi.mat','ACC Rt_roi.mat',...
                    'N.Ac Med Ant Lt_roi.mat','N.Ac Med Ant Rt_roi.mat'...
                    'VPe Ant Lt_roi.mat','VPe Ant Rt_roi.mat',...
                    'VPe Post Lt_roi.mat','VPe Post Rt_roi.mat'...
                    'Insula Lt_roi.mat','Insula Rt_roi.mat',...
                    'BA32sg_Lt_roi.mat','BA32sg_Rt_roi.mat'};
                
                timings={'Min EV','Min EV','Min EV','Min EV','Min EV','Min EV','Min EV',...
                    'Min EV','Min EV','Min EV','Min EV','Min EV','Min EV','Min EV',...
                    'Min EV','Min EV','Min EV','Min EV','Min EV','Min EV','Min EV',...
                    'Min EV','Min EV','Min EV','Min EV','Min EV','Min EV','Min EV',...
                    'PE','PE','PE','PE',...
                    'PE','PE','PE','PE','PE','PE','PE',...
                    'PE','PE','PE','PE','PE','PE','PE'...
                    'happiness','happiness',...
                    'happiness','happiness',...
                    'happiness','happiness',...
                    'happiness','happiness',...
                    'happiness','happiness',...
                    'happiness','happiness',...
                    'happiness','happiness'};
            end
            for at_time_fill=1:size(ROIs,2)
                if condition_before
                    timings{2,at_time_fill}= condition_before;
                end
            end
        end
        
        TargetBOLD='Subgenual vmPFC & mOFC Bi_roi.mat'; %respecify that
        Target_timing{1,1}='happiness';
        if condition_before
            Target_timing{2,1}=condition_before;
        end
        
end

%*********Model name and settings********
glm_name=glm_dir(end-30:end-9);
if behav==1 tgt_type='behav'; else tgt_type='BOLD'; end
if combine==1 comb_type='Bi'; else comb_type='Uni';end
if pooled==1 pooled_type='pool'; else pooled_type='ind'; end
if strcmp(shift,'standard shift') shifttype='stdshft'; else shifttype='trushft'; end
if condition_before; which_con_before=condition_before; else which_con_before='NSCB';end
if social_ROIs==1; social_ROIs_name='SocROI'; else social_ROIs_name='NSocROI';end
if zscored; zscore_state='Zscr'; else zscore_state='intrcpt'; end
if AR~=0 AR_state= ['AR(n)']; else AR_state='';end
if isa(MA,'double') && MA~=0 MA_state= ['MA(' num2str(MA) ')']; elseif isa(MA,'double') && MA==0 MA_state=''; elseif ~isa(MA,'double') MA_state=MA.filterstructure; end
modelName=[glm_name '_' model '_' tgt_type '_' which_scans '_' pooled_type '_' comb_type '_' ...
    shifttype '_' which_con_before '_' zscore_state '_' dependence '_' social_ROIs_name '_' AR_state '_' MA_state];

%*****Start the fitting*****
dup_idx=find(strcmp(ROIs,TargetBOLD));
if ~isempty(dup_idx)
    dup_idx_timing=find(strcmp(timings{dup_idx},Target_timing));
    if ~isempty(dup_idx_timing)
        fprintf('\n I found the same ROI BOLD at the same time on both sides...\n of the modelling equation. This will give you a fake R squared of 1 \n removing from the left side of the equation.... \n')
        rest_idx=find(~strcmp(ROIs,TargetBOLD));
        ROIs={ROIs{1,[rest_idx]}};
        for timings_i= 1: size(timings,1)
            for rest_idx_i=1: size(rest_idx,2)
                timings_corct{timings_i,rest_idx_i}=timings{timings_i,rest_idx(rest_idx_i)};
            end
        end
    end
else timings_corct=timings;
end
model_size=size(ROIs,2);
what_stage_lbl='';
if pooled==1
    sub_end=1;
    [results,tgt,~,fitstats]=BOLD_modelling(ROIs,timings_corct,TargetBOLD,Target_timing,model_size,behav,glm_dir,...
        sub,shift,equal_regressors,whatcombine,howcombine,pooled,zscored,dependence,bimodal,MA,AR,disk_parallel,use_4D,social_ROIs,condition_before,model,FIR);
    if ~isempty(strfind(dependence,'fused'))
        if AR~=0
            model_b=find([results(:).prediction_r2_lagcorct]==max([results(:).prediction_r2_lagcorct]));
        else
            model_b=find([results(:).prediction_r2]==max([results(:).prediction_r2]));
        end
        rescell=cell(sub_end,20);
        rescell_Bmodel=cell(sub_end,20);
        for m = 1:size(results,2)
            figure('pos',[10 10 1200 900]);
            subplot(3,1,1)
            try  plot(results(m).train_y,'b');
                hold on; plot(results(m).train_yhat,'r'); end
            dataNames = {'Actual data','training data predictions'};
            legend(dataNames)
            if AR~=0
                rf=['Model No(' num2str(m) ') Training data OLSAR(' num2str(results(m).nlag) ') regression'];
            else
                rf=['Model No(' num2str(m) ') Training data OLSAR regression'];
            end
            title(rf, 'interpreter','none')
            ylabel(sprintf('Regression R2 = %.2f',results(m).r2));
            xlabel(['BIC=' num2str(results(m).bic) '  AIC=' num2str(results(m).aic)])
            subplot(3,1,2)
            plot(results(m).test_y,'b');
            hold on; plot(results(m).test_yhat,'r')
            dataNames = {'Actual data','test data predictions'};
            legend(dataNames)
            ylabel([sprintf('Test prediction R2 = %.2f',results(m).prediction_r2) sprintf('   p-value = %.2f',results(m).prediction_pval)]);
            rf=['Model No(' num2str(m) ') Test data predictions versus actual test data'];
            title(rf, 'interpreter','none')
            subplot(3,1,3)
            bar(results(m).betas(1:size(results(m).ROIs,2)))
            title([ 'Model No(' num2str(m) ') '  modelName ' - parameters'], 'interpreter','none')
            lbl=cell(1,size(results(m).ROIs,2));
            for x_lbl=1:size(results(m).ROIs,2)
                lbl{1,x_lbl}=[results(m).ROIs{1,x_lbl} ' at ' results(m).timings{1,x_lbl}];
               
                    if results(m).pval(x_lbl,1) <0.05
                        lbl{1,x_lbl}=[lbl{1,x_lbl} '*'];
                    end
 
            end
            set(gca,'XTick',1:size(results(m).ROIs,2));
            try xlabel_oblique(lbl); end
            
            %*****Plotting Best model cross validations
            figure('Name','Best model - Cross Validations','pos',[10 10 1200 900]);
            subplot(4,1,1)
            plot(fitstats(m).y,'b');
            hold on; plot(fitstats(m).yhat,'r')
            dataNames = {'Actual data','training data predictions'};
            legend(dataNames)
            if AR~=0
                rf=['Best model - Cross Validation No(' num2str(m)  ') Training data OLSAR(' num2str(fitstats(m).nlag) ') regression'];
            else
                rf=['Best model - Cross Validation No(' num2str(m)  ') Training data OLS regression'];
            end
            title(rf, 'interpreter','none')
            ylabel(sprintf('Regression R2 = %.2f',fitstats(m).rsqr));
            xlabel(['BIC=' num2str(fitstats(m).bic) '  AIC=' num2str(fitstats(m).aic)])
            
            subplot(4,1,2)
            plot(fitstats(m).testB_y,'b');
            hold on; plot(fitstats(m).testB_yhat,'r')
            dataNames = {'Actual data','test data predictions'};
            legend(dataNames)
            ylabel([sprintf('prediction R2 = %.2f',fitstats(m).prediction_r2) sprintf('   p-value = %.2f',fitstats(m).prediction_pval)]);
            rf=['Best model - Cross Validation No(' num2str(m) ') Test data predictions versus actual test data'];
            title(rf, 'interpreter','none')
            if AR~=0
                subplot(4,1,3)
                plot((fitstats(m).testB_y-fitstats(m).lag_predB),'b');
                hold on; plot(fitstats(m).testB_yhat,'r')
                dataNames = {'Actual data - lag corrected','test data predictions'};
                legend(dataNames)
                ylabel([sprintf('prediction R2 = %.2f',fitstats(m).prediction_r2_lagcorct) sprintf('   p-value = %.2f',fitstats(m).prediction_pval_lagcorct)]);
                rf=['Best model - Cross Validation No(' num2str(m) '): Test data predictions vs. lag(' num2str(fitstats(m).nlag) ') - corrected actual test data'];
                title(rf, 'interpreter','none')
                subplot(4,1,4)
            else
                subplot(4,1,3)
            end
            bar(fitstats(m).beta(1:size(results(model_b).ROIs,2)))
            title(['Best model - Cross Validation No(' num2str(m) modelName ' - parameters,'], 'interpreter','none')
            lbl=cell(1,size(results(model_b).ROIs,2));
            for x_lbl=1:size(results(model_b).ROIs,2)
                lbl{1,x_lbl}=[results(model_b).ROIs{1,x_lbl} ' at ' results(model_b).timings{1,x_lbl}];
                
                    if fitstats(m).pval(x_lbl,1) <0.05
                        lbl{1,x_lbl}=[lbl{1,x_lbl} '*'];
                    end;
            end
            set(gca,'XTick',1:size(results(model_b).ROIs,2));
            try xlabel_oblique(lbl); end
            
            
            %prepare group stats cells
            
            rescell{m,1}=m;               rescell_Bmodel{m,1}=m;
            rescell{m,2}=results(m).r2;   rescell_Bmodel{m,2}=fitstats(m).rsqr;
            rescell{m,3}=results(m).aic; rescell_Bmodel{m,3}=fitstats(m).aic ;
            rescell{m,4}=results(m).bic;  rescell_Bmodel{m,4}=fitstats(m).bic;
            rescell{m,5}=results(m).betas; rescell_Bmodel{m,5}=fitstats(m).beta;
            rescell{m,6}=results(m).ROIs;  rescell_Bmodel{m,6}=results(model_b).ROIs; %stable in the best model
            rescell{m,7}=results(m).pval'; rescell_Bmodel{m,7}=fitstats(m).pval'; %Transposing this pval vector makes a confusion
            rescell{m,8}=results(m).timings;rescell_Bmodel{m,8}=results(model_b).timings; %stable in the best model
            rescell{m,9}=results(m).prediction_r2;    rescell_Bmodel{m,9}=fitstats(m).prediction_r2;
            rescell{m,10}=results(m).prediction_pval; rescell_Bmodel{m,10}=fitstats(m).prediction_pval;
            
            rescell{m,12}=results(m).lambda1; %stable in the best model
            rescell{m,13}=results(m).lambda2;%stable in the best model
            rescell{m,14}=results(m).lambda3;%stable in the best model
            rescell{m,15}=results(m).lambda_BIC;%not done in the best model
            rescell{m,17}=results(m).model_in;%it's only one model
            if AR~=0
                rescell{m,11}=results(m).nlag;            rescell_Bmodel{m,11}=fitstats(m).nlag;
                rescell{m,21}=results(m).prediction_r2_lagcorct;   rescell_Bmodel{m,21}=fitstats(m).prediction_r2_lagcorct;
                rescell{m,22}=results(m).prediction_pval_lagcorct; rescell_Bmodel{m,22}=fitstats(m).prediction_pval_lagcorct;
            end
            if ~isempty({rescell{m,6}}) && ~isempty({rescell_Bmodel{m,6}})
                for comb_i=1:size(rescell{m,6},2)
                    rescell{m,19}{1,comb_i}=[rescell{m,6}{1,comb_i} ' at ' rescell{m,8}{1,comb_i}];
                end
                for comb_i2=1:size(rescell_Bmodel{m,6},2)
                    rescell_Bmodel{m,19}{1,comb_i2}=[rescell_Bmodel{m,6}{1,comb_i2} ' at ' rescell_Bmodel{m,8}{1,comb_i2}];
                end
            end
            
        end
        %%
    else
        reg_tstat=sqrt(results.r2)*sqrt(df/(1-results.r2));
        results.tstat=reg_tstat;
        try p_value=tcdf(reg_tstat,df,'upper'); results.pvalue=p_value; end
        
        for m = 1:size(results,2)
            figure('pos',[10 10 1200 900]);
            subplot(2,1,1)
            plot(fitstats(m).y,'b');
            hold on; plot(fitstats(m).yhat,'r')
            dataNames = {'Actual data','predicted data'};
            legend(dataNames)
            rf='Pooled data from all subjects';
            title(rf, 'interpreter','none')
            ylabel(sprintf('R2 = %.2f',results(m).r2));
            xlabel(['BIC=' num2str(results(m).bic) '  AIC=' num2str(results(m).aic)])
            subplot(2,1,2)
            bar(results(m).b(1:size(results(1,m).paramNames,2)))
            %         sj_s='';
            %     for ssj=1:size(sub,2)
            %         sub_lbl{1,ssj}=[sj_s num2str(sub(ssj))];
            %     end
            % %    lbl= horzcat(results(1,m).paramNames);
            %     set(gca,'XTick',(1:size(lbl,2)))
            title([modelName ' - parameters, Intercepts for each subject are not shown'], 'interpreter','none')
            lbl=cell(1,size(results.paramNames,2));
            for x_lbl=1:size(results.paramNames,2)
                lbl{1,x_lbl}=[results.paramNames{1,x_lbl} ' at ' results.timings{1,x_lbl}];
                try
                    if results.pval(x_lbl,1) <0.05
                        lbl{1,x_lbl}=[lbl{1,x_lbl} '*'];
                    end
                end
            end
            set(gca,'XTick',1:size(results.paramNames,2));
            try xlabel_oblique(lbl); end
            
            
            
        end
    end
    fprintf('Saving figures \n')
    h=findobj('type','figure'); % find the handles of the opened figures
    folder=dest_dir;  % Desination folder
    for k=1:numel(h)
        filename=sprintf('%d.jpg',k);
        filename2=sprintf('%d.fig',k);
        fileset=[glm_name '_' model '_' tgt_type '_' which_scans '_' pooled_type '_' comb_type '_' ...
            shifttype '_' which_con_before '_' zscore_state '_' dependence '_' social_ROIs_name '_' AR_state '_' MA_state];
        figure_name=[fileset '_' filename];
        figure_name2=[fileset '_' filename2];
        file=fullfile(folder,figure_name);
        file2=fullfile(folder,figure_name2);
        if ~isempty(h)
            saveas(h(k),file)
            saveas(h(k),file2)
        end
        close
    end
    %%
    %%
elseif pooled==0
    sub_end=size(sub,2);
    fitstats=cell(sub_end,1);
    rescell=cell(sub_end,20);
    for sub_id= 1:sub_end
        clear tgt
        fprintf(['   \n Doing Subject: ' num2str(sub(sub_id)) '  ']);
        evalc(' [results(sub_id,1),tgt,~,fitstats{sub_id,1}]=BOLD_modelling(ROIs,timings_corct,TargetBOLD,Target_timing,model_size,behav,glm_dir,sub(sub_id),shift,equal_regressors,whatcombine,howcombine,pooled,zscored,dependence,bimodal,MA,AR,disk_parallel,use_4D,social_ROIs,condition_before,model,FIR)');
        if ~isempty(strfind(dependence,'fused'))
            figure('pos',[10 10 700 1000]);subplot(4,1,1);
            plot(tgt,'b'); hold on; plot(fitstats{sub_id,1}.lagged_yhat,'r');
            title('Lagged predictions against actual ratings');
            dataNames_1 = {'Actual data','lagged predicted data'};
            legend(dataNames_1);
            ylabel([sprintf('R2 = %.2f',fitstats{sub_id,1}.lagged_BOLDr2_fused)...
                sprintf('   LSE = %.2f',fitstats{sub_id,1}.lagged_BOLD_lse)]);
            xlabel([sprintf('t-statistic = %.2f ',fitstats{sub_id,1}.lagged_BOLD_tstat)...
                sprintf('   p-value = %.2f',fitstats{sub_id,1}.lagged_BOLD_pval)]);
            subplot(4,1,2);
            plot(tgt,'b'); hold on; plot(fitstats{sub_id,1}.lagterms_pred,'r');
            title(['Autoregression(' num2str(fitstats{sub_id,1}.nlag) ') Lag terms predictions against actual ratings']);
            dataNames_1 = {'Actual data','predictions of lag terms'};
            legend(dataNames_1);
            ylabel([sprintf('R2 = %.2f',fitstats{sub_id,1}.lagterms_r2)...
                sprintf('   LSE = %.2f',fitstats{sub_id,1}.lagterms_lse)]);
            xlabel([sprintf('t-statistic = %.2f ',fitstats{sub_id,1}.lagterms_tstat)...
                sprintf('   p-value = %.2f',fitstats{sub_id,1}.lagterms_pval)]);
            
            subplot(4,1,[3 4])
            cVals = unique(fitstats{sub_id,1}.lambda_BIC(:,4));
            for i = 1:numel(cVals)                     % For every one of those unique values
                indices = find(fitstats{sub_id,1}.lambda_BIC(:,4) == cVals(i));         % Find the corresponding indices
                scatter3(fitstats{sub_id,1}.lambda_BIC(indices,1),fitstats{sub_id,1}.lambda_BIC(indices,2),fitstats{sub_id,1}.lambda_BIC(indices,3),100,'filled') % Plot
                hold on
            end
            lgd1=cell(1,0);
            for lgd=1:size(cVals,1)
                lgd1=horzcat(lgd1,{['BIC=  ' num2str(round(cVals(lgd,1)),2)]});
            end
            legend(lgd1)
            % plot([fitstats{sub_id,1}.lambda_BIC(:,1)], [fitstats{sub_id,1}.lambda_BIC(:,4)])
            xlabel('Fused lasso lambda');  ylabel('L1 norm lambda "Lasso"');zlabel('L2 norm lambda "Ridge"'); title('Flat Fused Lasso(FISTA optimization): BIC for each Regularization parameter Lambda');
            
        end
        figure('pos',[10 10 700 900]);
        subplot(2,1,1)
        plot(tgt,'b');
        hold on; plot(results(sub_id,1).BOLDpred,'r')
        dataNames = {'Actual data','predicted data'};
        if ~isempty(strfind(dependence,'fused'))
            dataNames{1,2}=['AR(' num2str(fitstats{sub_id,1}.nlag) ') corrected predictions'];
        end
        legend(dataNames)
        if fitstats{sub_id,1}.nlag >0
            txt=['  Sims Corrected p-value: ' num2str(round(fitstats{sub_id,1}.Sims_dof_correced_probability,3)) ' Granger marginal probability: ' num2str(round(fitstats{sub_id,1}.fprob,3))];
        else
            txt=[' Granger marginal probability ' num2str(round(fitstats{sub_id,1}.fprob,3))];
        end
        title(['Subject: ' num2str(sub(sub_id)) txt], 'interpreter','none')
        ylabel([sprintf('R2 = %.2f',results(sub_id,1).r2)...
            sprintf('   LSE = %.2f',fitstats{sub_id,1}.lse)]);
        xlabel([sprintf('t-statistic = %.2f ',fitstats{sub_id,1}.model_tstat)...
            sprintf('   p-value = %.2f',fitstats{sub_id,1}.model_pval)...
            sprintf('   F-statistic = %.2f ',fitstats{sub_id,1}.ftest)...
            sprintf('   F test p-val = %.2f',fitstats{sub_id,1}.fprob)]);
        subplot(2,1,2)
        bar(results(sub_id,1).b)
        title(['Parameters estimates'], 'interpreter','none');
        lbl=cell(1,size(results(sub_id,1).paramNames,2));
        for x_lbl=1:size(results(sub_id,1).paramNames,2)
            lbl{1,x_lbl}=[results(sub_id,1).paramNames{1,x_lbl} ' at ' results(sub_id,1).timings{1,x_lbl}];
            try
                if results(sub_id,1).pval(x_lbl,1) <0.05
                    lbl{1,x_lbl}=[lbl{1,x_lbl} '*'];
                end
            end
        end
        set(gca,'XTick',1:size(results(sub_id,1).paramNames,2));
        try xlabel_oblique(lbl); end
        drawnow
        fprintf('Saving figures \n')
        h=findobj('type','figure'); % find the handles of the opened figures
        folder=dest_dir;  % Desination folder
        
        for k=1:numel(h)
            filename=sprintf('%d.jpg',k);
            filename2=sprintf('%d.fig',k);
            fileset1=[glm_name '_' model '_' tgt_type '_' which_scans '_' pooled_type '_' comb_type '_' ...
                shifttype '_' which_con_before '_' zscore_state '_' num2str(sub(sub_id)) '_' dependence '_' social_ROIs_name '_' AR_state '_' MA_state ];
            figure_name=[fileset1 '_' filename];
            figure_name2=[fileset1 '_' filename2];
            file=fullfile(folder,figure_name);
            file2=fullfile(folder,figure_name2);
            saveas(h(k),file)
            saveas(h(k),file2)
            close
        end
        
        
        
        %  disp(reshape([results.r2],sub_id,[]));
        
        sig_results(sub_id)=results(sub_id);
        if isempty(find(results(sub_id).pval<0.05))
            sig_results(sub_id).aic=0;
            sig_results(sub_id).bic=0;
            sig_results(sub_id).BOLDpred=mean(tgt,1)+zeros(size(tgt,1),1);
            sig_results(sub_id).b=zeros(size(results(sub_id).b,1),1);
            sig_results(sub_id).paramNames='No significant terms';
            sig_results(sub_id).timings=' ';
        end
        if fitstats{sub_id,1}.model_pval>0.05
            sig_results(sub_id).model_in=-1;
        else
            sig_results(sub_id).model_in=1;
        end
        rescell{sub_id,1}=sub(sub_id);
        rescell{sub_id,2}=fitstats{sub_id,1}.rsqr;
        rescell{sub_id,3}=results(sub_id).aic ;
        rescell{sub_id,4}=results(sub_id).bic;
        rescell{sub_id,5}=results(sub_id).b;
        rescell{sub_id,6}=results(sub_id).paramNames;
        rescell{sub_id,7}=results(sub_id).pval;
        rescell{sub_id,8}=results(sub_id).timings;
        rescell{sub_id,9}=fitstats{sub_id,1};
        rescell{sub_id,10}=fitstats{sub_id,1}.rbar;
        rescell{sub_id,11}=fitstats{sub_id,1}.nlag;
        rescell{sub_id,12}=fitstats{sub_id,1}.lambda1;
        rescell{sub_id,13}=fitstats{sub_id,1}.lambda2;
        rescell{sub_id,14}=fitstats{sub_id,1}.lambda3;
        rescell{sub_id,15}=fitstats{sub_id,1}.lagged_BOLDr2_fused;
        rescell{sub_id,16}=fitstats{sub_id,1}.lagterms_r2;
        rescell{sub_id,17}=fitstats{sub_id,1}.lagged_BOLD_bic;
        rescell{sub_id,18}=fitstats{sub_id,1}.lagterms_bic;
        if ~isempty({rescell{sub_id,6}})
            for comb_i=1:size(rescell{sub_id,6},2)
                rescell{sub_id,19}{1,comb_i}=[rescell{sub_id,6}{1,comb_i} ' at ' rescell{sub_id,8}{1,comb_i}];
            end
        end
        
    end
end
%%
%%
%Do group stats

if (pooled==0 && ~isempty(strfind(dependence,'fused')))|| (pooled==1 && ~isempty(strfind(dependence,'fused')))
    if  (pooled==1 && ~isempty(strfind(dependence,'fused')))
        %     Best_models_in=find([results.model_in]~=-1);
        %     Best_modelstat.models_in=size(models_in,2);
        Best_modelstat.mean_r2= mean([fitstats.rsqr],2);
        Best_modelstat.median_r2= median([fitstats.rsqr],2);
        Best_modelstat.mean_prediction_r2= mean([fitstats.prediction_r2],2);
        Best_modelstat.median_prediction_r2= median([fitstats.prediction_r2],2);
        if AR~=0
            Best_modelstat.mean_prediction_r2_lagcorct= mean([fitstats.prediction_r2_lagcorct],2);
            Best_modelstat.median_prediction_r2_lagcorct= median([fitstats.prediction_r2_lagcorct],2);
        end
        Best_modelstat.BIC=sum([fitstats.bic],2);
        Best_modelstat.AIC=sum([fitstats.aic],2);
        Best_modelstat.zr2=zscore([fitstats.rsqr]); %which r2 to use?
    end
    models_in=find([results.model_in]~=-1);
    modelstat.models_in=size(models_in,2);
    modelstat.mean_r2= mean([results.r2],2);
    modelstat.median_r2= median([results.r2],2);
    modelstat.BIC=sum([results.bic],2);
    modelstat.AIC=sum([results.aic],2);
    modelstat.zr2=zscore([results.r2]);
    if pooled==0
        modelstat.mean_r2_adjusted=mean([rescell{:,10}],2);
        modelstat.median_r2_adjusted=median([rescell{:,10}],2);
        modelstat.BIC_lagged_BOLD=sum([rescell{:,17}],2);
        modelstat.BIC_lagterms=sum([rescell{:,18}],2);
        stat_ROIs=[sig_results.paramNames];
        modelstat.ROIs=unique(stat_ROIs);
        modelstat.zr2=zscore([results.r2]);
    end
end
%%
%%
if pooled==0 || (pooled==1 && ~isempty(strfind(dependence,'fused')))
    if ~isempty(strfind(dependence,'stepwise')) || ~isempty(strfind(dependence,'lasso'))
        all_fit_ROIs=cell(1,2);
        if ~isempty(strfind(dependence,'stepwise')) || ~isempty(strfind(dependence,'fused'))
            if pooled==1 && ~isempty(strfind(dependence,'fused'))
                what_stage_lbl={'Feature selection iterations: ', 'Best features Cross-validations: '};
                rescell_backup=rescell;
                for what_stage=1:2
                    
                    if what_stage==2
                        rescell=rescell_Bmodel;
                    end
                    
                    for n=1:size(rescell,1)
                        for ROIs_num=1:size(rescell{n,19},2)
                            %try
                                if rescell{n,7}(1,ROIs_num) <0.05 %only considers significant fittings
                                    all_fit_ROIs{what_stage}{1,end+1}=(rescell{n,19}{1,ROIs_num});
                                end;
                            %end
                        end
                    end
                    
                    unique_ROIs{what_stage}=unique(all_fit_ROIs{what_stage});
                    unique_ROI_array{what_stage}=cell(0,6);
                    for uri=1:size(unique_ROIs{what_stage},2)
                        if ~isempty(all_fit_ROIs{what_stage})
                            ROI_loc=strfind(all_fit_ROIs{what_stage},unique_ROIs{what_stage}{1,uri});
                            if ~isempty(ROI_loc)
                                non_empt_unique=find(~cellfun(@isempty,ROI_loc));
                                unique_ROI_array{what_stage}{end+1,1}=unique_ROIs{what_stage}{1,uri};
                                unique_ROI_array{what_stage}{end,2}=size(non_empt_unique,2);
                            end
                        end
                        for get_unique_r2=1:size(unique_ROI_array{what_stage},1)
                            for n=1:size(rescell,1)
                                if ~isempty({rescell{n,6}}) && ~isempty(rescell{n,5})
                                    r2_idx=strfind({rescell{n,19}{1,:}},unique_ROI_array{what_stage}{get_unique_r2,1});
                                    
                                    if ~isempty(r2_idx)
                                        nonemptidx=find(~cellfun(@isempty,r2_idx));
                                        %try
                                            if rescell{n,7}(1,nonemptidx) <0.05 %only considers significant fittings
                                                unique_ROI_array{what_stage}{get_unique_r2,3}(1,end+1)=(rescell{n,2}(1,1));
                                                unique_ROI_array{what_stage}{get_unique_r2,6}(1,end+1)=(rescell{n,5}(nonemptidx,1));
                                            end
                                        %end
                                    end
                                end
                            end
                            unique_ROI_array{what_stage}{get_unique_r2,2}=size(unique_ROI_array{what_stage}{get_unique_r2,3},2);
                        end
                        unique_ROI_array{what_stage}=sortrows(unique_ROI_array{what_stage},2);
                        
                        nonempt=find(~cellfun(@isempty,{unique_ROI_array{what_stage}{:,3}}));
                        unique_ROI_array{what_stage}=unique_ROI_array{what_stage}([nonempt],:);
                        if what_stage==1
                            modelstat.ROIs_frequency=unique_ROI_array{1};
                        elseif what_stage==2
                            Best_modelstat.ROIs_frequency=unique_ROI_array{2};
                        end
                        try
                            for n=1:size(unique_ROI_array{what_stage},1)
                                for ROIs_lbl_3d=1:unique_ROI_array{what_stage}{n,2}
                                    unique_ROI_array{what_stage}{n,4}(1,ROIs_lbl_3d)=unique_ROI_array{what_stage}{n,2};
                                end
                                unique_ROI_array{what_stage}{n,5}=1:unique_ROI_array{what_stage}{n,2};
                            end
                            %                                 figure;
                            %                                 stem3([unique_ROI_array{what_stage}{:,4}],[unique_ROI_array{what_stage}{:,5}],[unique_ROI_array{what_stage}{:,3}],'filled','LineWidth',10)
                            %                                 axis([0 10 0 size(rescell,1) 0 1])
                        end
                    end
                end
                rescell=rescell_backup;
            end
        end
        %%
        %%
        %*****R-squared descriptive statistics*****
        %DOF Adjusted R-squared
        
        if pooled==0
            figure('name',['R-squared group stats - ' modelName],'pos', [450 50 700 900]);
            subplot(3,1,1)
            bar([rescell{:,10}])
            title('Lag corrected fit - DOF Adjusted R-squared')
            sj=[rescell{:,1}];
            sj_s='';
            for ssj=1:size(sj,2)
                sub_lbl{1,ssj}=[sj_s num2str(sj(ssj))];
            end
            % xlabel('Subjects')
            set(gca,'XTick',(1:size(results,1)))
            try xlabel_oblique(sub_lbl); end
            subplot(3,1,2)
            hist([rescell{:,10}],size([rescell{:,10}],2)/5)
            title('DOF Adjusted R2 histogram/Subject ID')
            ylabel_meanr2=['mean R2= ' num2str(round(modelstat.mean_r2_adjusted,2)) ' & median R2= ' num2str(round(modelstat.median_r2_adjusted,2))];
            ylabel(ylabel_meanr2);
            if AR~=0
                subplot(3,1,3)
                histogram([rescell{:,11}])
                title('Auto regression lag Histogram')
            end
            %uncorrected BOLD and lag terms fit R-squared
            figure('name',['BOLD & lag terms R2 group stats - ' modelName],'pos', [450 50 700 900]);
            subplot(4,1,1)
            bar([rescell{:,15}])
            title('Uncorrected BOLD (OLS fit) - R2/Subject ID')
            sj=[rescell{:,1}];
            sj_s='';
            for ssj=1:size(sj,2)
                sub_lbl{1,ssj}=[sj_s num2str(sj(ssj))];
            end
            % xlabel('Subjects')
            set(gca,'XTick',(1:size(results,1)))
            try xlabel_oblique(sub_lbl); end
            subplot(4,1,2)
            hist([rescell{:,15}],size([rescell{:,15}],2)/5)
            title('Uncorrected BOLD (OLS fit) - R2 histogram')
            ylabel_meanr2=['mean R2= ' num2str(round(mean([rescell{:,15}]),2)) ' & median R2= ' num2str(round(median([rescell{:,15}]),2))];
            ylabel(ylabel_meanr2);
            xlabel(['BIC=' num2str(modelstat.BIC_lagged_BOLD)])
            subplot(4,1,3)
            bar([rescell{:,16}])
            title('Lag terms (OLS fit) - R2/Subject ID')
            sj=[rescell{:,1}];
            sj_s='';
            for ssj=1:size(sj,2)
                sub_lbl{1,ssj}=[sj_s num2str(sj(ssj))];
            end
            %  xlabel('Subjects')
            set(gca,'XTick',(1:size(results,1)))
            try xlabel_oblique(sub_lbl); end
            subplot(4,1,4)
            hist([rescell{:,16}],size([rescell{:,16}],2)/5)
            title('Lag terms (OLS fit) - R2 histogram')
            ylabel_meanr2=['mean R2= ' num2str(round(mean([rescell{:,16}]),2)) ' & median R2= ' num2str(round(median([rescell{:,16}]),2))];
            ylabel(ylabel_meanr2);
            xlabel(['BIC=' num2str(modelstat.BIC_lagterms)])
        end
        %%
        
        %*****Parameter/ROIs descriptive statistics*****
        if ~isempty(strfind(dependence,'stepwise')) || ~isempty(strfind(dependence,'lasso'))
            for what_stage=1:2
                if ~isempty([unique_ROI_array{what_stage}{:,2}])
                    C=zeros(1,0);
                    grp=zeros(1,0);
                    for get_box= 1:size(unique_ROI_array{what_stage},1)
                        C = [C unique_ROI_array{what_stage}{get_box,6} ];
                        grp = [grp ((get_box-1)+zeros(1,size(unique_ROI_array{what_stage}{get_box,6},2)))];
                    end
                    figure('name',['Parameters & ROIs for ' what_stage_lbl{what_stage} modelName],'pos',[10 10 900 1200]);
                    subplot(4,1,1)
                    if what_stage==1; bar([rescell{:,2}]); elseif what_stage==2 bar([fitstats.rsqr]); end
                    title([ what_stage_lbl{what_stage} 'Lag corrected fit - R2/Subject ID'])
                    sj=[rescell{:,1}];
                    sj_s='';
                    for ssj=1:size(sj,2)
                        sub_lbl{1,ssj}=[sj_s num2str(sj(ssj))];
                    end
                    % xlabel('Subjects')
                    set(gca,'XTick',(1:size(results,1)))
                    try xlabel_oblique(sub_lbl); end
                    
                    subplot(4,1,2)
                    title([what_stage_lbl{what_stage} 'Parameter estimates']);
                    boxplot(C,grp)
                    try xlabel_oblique({unique_ROI_array{what_stage}{:,1}}); end
                    subplot(4,1,4)
                    bar([unique_ROI_array{what_stage}{:,2}])
                    title([what_stage_lbl{what_stage} 'ROI frequency'])
                    set(gca,'XTick',(1:size([unique_ROI_array{what_stage}{:,2}],2)))
                    try xlabel_oblique({unique_ROI_array{what_stage}{:,1}}); end
                end
            end
            figure('name',['R2 group statistics for ' modelName],'pos',[10 10 900 1200]);
            subplot(2,1,1)
            histogram([rescell{:,2}])
            title('R2 histogram')
            ylabel_meanr2=['mean R2= ' num2str(round(modelstat.mean_r2,2)) ' & median R2= ' num2str(round(modelstat.median_r2,2))];
            ylabel(ylabel_meanr2);
            xlabel(['BIC=' num2str(modelstat.BIC)])
            if ~isempty(strfind(dependence,'fused'))
                title('Training set regression R2 histogram')
                subplot(2,1,2)
                histogram([results(:).prediction_r2])
                title('Test prediction R2 histogram')
                ylabel_meanr2=['mean R2= ' num2str(round(mean([results(:).prediction_r2],2),2)) ' & median R2= ' num2str(round(median([results(:).prediction_r2],2),2))];
                ylabel(ylabel_meanr2);
            end
            %Feature selection Lambda values
            figure('name',['Feature selection Lambda values ' modelName],'pos',[10 10 700 900]);
            subplot(5,1,1)
            histogram([rescell{:,12}])
            title('Autoregression Regularization "Fused" histogram')
            subplot(3,1,2)
            histogram([rescell{:,13}])
            title('L1 norm Regularization "Lasso" histogram')
            subplot(3,1,3)
            histogram([rescell{:,14}])
            title('L2 norm Regularization "Ridge" histogram')
        end
    end
end
fprintf('Saving figures \n')
h=findobj('type','figure'); % find the handles of the opened figures
folder=dest_dir;  % Desination folder

for k=1:numel(h)
    filename=sprintf('%d.jpg',k);
    filename2=sprintf('%d.fig',k);
    fileset=[glm_name '_' model '_' tgt_type '_' which_scans '_' pooled_type '_' comb_type '_' ...
        shifttype '_' which_con_before '_' zscore_state '_' dependence '_' social_ROIs_name '_' AR_state '_' MA_state];
    figure_name=[fileset '_' filename];
    figure_name2=[fileset '_' filename2];
    file=fullfile(folder,figure_name);
    file2=fullfile(folder,figure_name2);
    saveas(h(k),file)
    saveas(h(k),file2)
end

results_name=[modelName '_results.mat'];
MSg_name=[modelName '_model_groupstats.mat'];
BFg_name=[modelName '_Best_features_groupstats.mat'];
fitstats_name=[modelName '_fitstats.mat'];
BF_name=[modelName '_Best_features_stats.mat'];
FI_name=[modelName '_feature_iterations_stats.mat'];
fprintf('saving the results matrix and group statistics \n')

if ~isempty(strfind(dependence,'fused')) && pooled==1
    save(fullfile(folder,MSg_name),'modelstat')
    save(fullfile(folder,BFg_name),'Best_modelstat')
    save(fullfile(folder,BF_name),'fitstats')
    save(fullfile(folder,FI_name),'results')
else
    save(fullfile(folder,fitstats_name),'fitstats')
    save(fullfile(folder,results_name),'results')
end


end
