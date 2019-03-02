function [shifts,person_sess_onsets,person_sess_BOLD,corct_behav_data,session_num,sub_tgt_id,sub_coord]=corr_BOLD(fit,figr)
fprintf('\n Always set equal_regressors to 0 unless you are calling \n this function to model unequal number of trials \n')
clear CRB corct_BOLD behav_data x y
subjectsfolder=fit.GLMs; %WARNING: remeber to remove the "/" at the end
if exist('figr','var') && figr==1
    fig=figr;
else
    fig=0;
end
persons{1}=fit.sub;
person_sess_onsets=cell(1,0);
person_sess_BOLD=cell(1,0);
corct_behav_data=cell(1,0);
shifts=0;
session_num=0;
sub_coord=0;
if ~strcmp(fit.TargetBOLD,fit.ROIs)
    if exist('combine','var') && fit.combine(1,1) ~=0
        fprintf(' \n Combining all ROIs into one BOLD vector \n \n')
        [Subjects_BOLD,shifts,~,~,sub_coord]=shiftcombine_BOLD(fit.ROIs,fit.timings,1,fit.GLMs,fit.sub,fit.shifts ,fit.combine,fit);
        fit.ROIs={Subjects_BOLD{1,1}.ROI};
    else
        fprintf(' \n getting BOLD data and shifting the scans according HRF shape for the condition \n \n')
        evalc('[Subjects_BOLD,shifts,sub_coord]=extract_fitBOLD(fit)');
    end
end
% initialISIandSetupValuesFile = 'Rutledge_InitialISIs_final.txt';
%getting Onsets
fprintf('\n getting Condition onsets \n')
% Subjects=subdir(subjectsfolder);
fit.timings=cellstr(fit.timings);
fit.ROIs=cellstr(fit.ROIs);
fit.shifts =0;
Pooledcorrbehav=zeros(0,0);
pooledcorrBOLD=zeros(0,0);
sub_tgt_id=zeros(0,1);
for person_i = 1 : length( persons{1} )   
    corct_BOLD_subROIs=zeros(0,0);
    %corct_BOLD=zeros(0,1);
    cd(subjectsfolder)
    callperson = persons{1}(person_i); %calls a specific folder according to index
    personName = ['SODEC_FMRI_' num2str(callperson)]; %get the directory as a character array
    fprintf(' \n Checking new subject \n \n')
    shortName= ['\n' personName '\n \n'];
    fprintf(shortName);
    cd (personName)
    load('sessionData.mat');
    %loading all behavioural onsets and sorting them.
    ons_sess=cell(0,2);
    allons=zeros(0,0);
    for behav_sess= 1:size(sessionData,1)
        allons=horzcat(allons,sessionData{behav_sess,3}{1:end});
    end
    allons=allons';
    allons_sort=sort(allons);
    allons_sort=unique(allons_sort,'rows');
    for sess= 1:size(sessionData,1)
        ons_sess{1,sess}=(sessionData{sess,3})';
    end
    ons_behav=zeros(0,1);
    behav_data=zeros(0,1);
    corct_behav_data_V=zeros(0,1);
    
    for behav_sess= 1:size(sessionData,1)
        sess_ons=zeros(0,1);
        behav_data=zeros(0,1);
        for behav_look= 1: size(sessionData{behav_sess,2},2)
            if ~isempty(strfind(sessionData{behav_sess,2}{1,behav_look},fit.timings{1,1}))
                ons_behav=vertcat(ons_behav,transpose(sessionData{behav_sess,3}{1,behav_look}));
                behav_data=vertcat(behav_data,transpose(sessionData{behav_sess,6}{1,behav_look}));
                sess_ons=vertcat(sess_ons,transpose(sessionData{behav_sess,3}{1,behav_look}));
            end
        end
        person_sess_onsets{1,end+1}= sess_ons;
        corct_behav_data{1,end+1}=behav_data;
    end
    
    for i= size(corct_behav_data,2)-1 : size(corct_behav_data,2)
        if fit.zscored==1
            corct_behav_data{1,i}=zscore(corct_behav_data{1,i});
        end
    end
    
    %% Find a special Condition for each trial
    if ischar(fit.condition_before)
    fprintf('\n Found a special condition "must qualify a specific name for the previous trial"')
    for behav_sess3=1:size(sessionData,1)
        current_sess=size(person_sess_onsets,2)-size(sessionData,1)+behav_sess3;
        for get_behav_trial_sm=1:size(person_sess_onsets{1,current_sess},1)
            for get_behav_trial_lrg=1:size(allons_sort,1)
                if allons_sort(get_behav_trial_lrg,1)==person_sess_onsets{1,current_sess}(get_behav_trial_sm,1)
                    print_trialnm=['\n \n 1. found trial: ' num2str(get_behav_trial_sm) ' at  sorted index number ' num2str(get_behav_trial_lrg)];
                    fprintf(print_trialnm)
                    
                    %First trial in the experiment
                    if get_behav_trial_lrg==1
                        fprintf('\nthis is the very first trial, cannot know the state of a previous trial. Will exclude')
                        person_sess_onsets{1,current_sess}(get_behav_trial_sm,1)=9999;
                    else
                        
                        %all later trials in the experiment
                        prv_trial_ons=allons(get_behav_trial_lrg-1,1);
                        onsnum=['\n the onset index ' num2str(get_behav_trial_lrg-1) ' is at scan number ' num2str(prv_trial_ons)];
                        fprintf(onsnum)
                        for behav_sess= 1:size(sessionData,1)
                            for get_prv_trial_type_i=1:size(sessionData{1,3},2)
                                prv_trial_loc=zeros(0,1);
                                trial_name=0;
                                if ~isempty(find(sessionData{behav_sess,3}{1,get_prv_trial_type_i}==prv_trial_ons, 1))
                                    loc2=find(sessionData{behav_sess,3}{1,get_prv_trial_type_i}==prv_trial_ons);
                                    prv_trial_loc=vertcat(prv_trial_loc,[get_prv_trial_type_i,loc2(1,1)]);
                                    print_prvtrialnm=['\n 2. found trial onset ' num2str(prv_trial_ons) ' at session: ' num2str(behav_sess)...
                                        ' Con: ' num2str(prv_trial_loc(1,1)) ' scan: ' num2str(prv_trial_loc(1,2)) '\n'];
                                    fprintf(print_prvtrialnm)
                                    for trial_name_i=1:size(prv_trial_loc,1)
                                        trial_name=sessionData{behav_sess,2}{1,prv_trial_loc(trial_name_i,1)};
                                        if isempty(strfind(trial_name,fit.condition_before))
                                            person_sess_onsets{1,current_sess}(get_behav_trial_sm,1)=9999;
                                            corct_behav_data{1,current_sess}(get_behav_trial_sm,1)=9999;
                                            print_conname=['---Condition "' fit.condition_before '" is not there'];
                                            fprintf(print_conname)
                                        else
                                            fprintf('\n 3. found a matching previous ')
                                            trial_name
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
    
    
    corct_sess_ons=zeros(0,1);
    for behav_sess4=1:size(sessionData,1)
        current_sess=size(person_sess_onsets,2)-size(sessionData,1)+behav_sess4;
        non_nines_ons_behav=find(person_sess_onsets{1,current_sess}~=9999);
        person_sess_onsets{1,current_sess}=person_sess_onsets{1,current_sess}(non_nines_ons_behav);
        corct_sess_ons=vertcat(corct_sess_ons,person_sess_onsets{1,current_sess});
        non_nines_behav=find(corct_behav_data{1,current_sess}~=9999);
        corct_behav_data{1,current_sess}=corct_behav_data{1,current_sess}(non_nines_behav);
        corct_behav_data_V=vertcat(corct_behav_data_V,corct_behav_data{1,current_sess});
    end
    last_person_behav_data=corct_behav_data_V;
    ons_behav=corct_sess_ons;
    Pooledcorrbehav=vertcat(Pooledcorrbehav,corct_behav_data_V);
    end 
    scan_behav=round(ons_behav,-1)./2.5;
    %% Getting relevant BOLD and behavioural data vector
    if ~strcmp(fit.TargetBOLD,fit.ROIs)
        fprintf('\n Throwing unrelated scans \n \n ')
        nscans=size(Subjects_BOLD{1,person_i}.BOLD{1,1},1);
        nsubROIs=size(Subjects_BOLD{1,person_i}.BOLD,2);
        corct_BOLD=cell(1,nsubROIs);
        BOLD_behav=cell(1,nsubROIs);corct_BOLD_subROIs=cell(1,nsubROIs);
        if isempty(strfind(Subjects_BOLD{1,person_i}.ROI,fit.ROIs))==0
            for subROIs_loop=1:nsubROIs
                if ~isempty(strfind(Subjects_BOLD{1,person_i}.Conditions,fit.timings{1,1}))
                    if ~isempty(Subjects_BOLD{1,person_i}.BOLD{1,subROIs_loop})
                        subROI=Subjects_BOLD{1,person_i}.BOLD{1,subROIs_loop};
                        avgzeros_i= subROI==0;
                        subROI(avgzeros_i)=mean(mean(subROI,2),1);
                        BOLD_behav{1,subROIs_loop}=mean(subROI,2);
                    end
                end
                corct_BOLD_subROIs{1,subROIs_loop}=[corct_BOLD_subROIs{1,subROIs_loop},BOLD_behav{1,subROIs_loop}(scan_behav-shifts(person_i,subROIs_loop))];
                clear subROI
            end
           % corct_BOLD=vertcat(corct_BOLD,corct_BOLD_subROIs);
        end
        
        
        %% Session relevant BOLD
        for sess2_i= 1:size(sessionData,1)
            tmp=[];
            scansess_behav=round(person_sess_onsets{1,size(person_sess_BOLD,2)+1},-1)./2.5;
            for i=1:nsubROIs
                tmp(:,i)=BOLD_behav{1,i}(scansess_behav-shifts(person_i,i),:);
            end
            person_sess_BOLD{1,end+1}=tmp;
            % sub_tgt_id=vertcat(sub_tgt_id,zeros(size(person_sess_BOLD{1,end},1),1)+fit.sub(person_i));
        end
        session_num(person_i,1)=size(sessionData,1);
        % pooledcorrBOLD=vertcat(pooledcorrBOLD,corct_BOLD);
    end
    
    
end
end
