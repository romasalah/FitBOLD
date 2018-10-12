function [CRB,corct_BOLD,last_person_behav_data,shifts,corct_BOLD_equal,Pooledcorrbehav,...
    pooledcorrBOLD_equal,pooledcorrBOLD,person_sess_onsets,person_sess_BOLD,corct_behav_data,session_num]=corr_BOLD(ROIs,...
    Conditions,behav,dir,sub,equal_regressors,shift,figr,combine,disk_parallel,use_4D,MA,zscored,FIR)
fprintf('\n Always set equal_regressors to 0 unless you are calling \n this function to model unequal number of trials \n')
clear CRB corct_BOLD behav_data x y
subjectsfolder=dir; %WARNING: remeber to remove the "/" at the end
if ~exist('shift','var')
    shift=0;
end
if exist('figr','var') && figr==1
    fig=figr;
else
    fig=0;
end
if exist('sub','var')
    persons{1}=sub;
else
    fprintf(' \n No Subjects specified... doing all of them \n \n')
    persons{1} = [ 44 62 63 65 69 70 71 74 77 79 81 83 84 90 92 93 94 95 96 98 101 102 103 104 106 110 120 124 125 126 129 136 155 166 171 173 174 175 176 177 ];
    sub= persons{1};
end
person_sess_onsets=cell(1,0);
person_sess_behav=cell(1,0);
person_sess_BOLD=cell(1,0);
corct_behav_data=cell(1,0);
corct_behav_data_V=zeros(0,1);
if exist('combine','var') && combine(1,1) ~=0
    fprintf(' \n Combining all ROIs into one BOLD vector \n \n')
    [Subjects_BOLD,shifts,~]=shiftcombine_BOLD(ROIs,Conditions,1,dir,sub,shift,combine);
    ROIs={Subjects_BOLD{1,1}.ROI};
else
    fprintf(' \n getting BOLD data and shifting the scans according HRF shape for the condition \n \n')
    [Subjects_BOLD,shifts,~,Conditions_lookup]=BOLD(ROIs,Conditions,1,dir,sub,shift,disk_parallel,use_4D);
end
% initialISIandSetupValuesFile = 'Rutledge_InitialISIs_final.txt';
%getting Onsets
fprintf('\n getting Condition onsets \n')
Conditions_lookup=Conditions_lookup;
% Subjects=subdir(subjectsfolder);
Conditions=cellstr(Conditions);
ROIs=cellstr(ROIs);
shift=0;
CRB=zeros(0,1);

pooledcorrBOLD_equal=zeros(0,1);
Pooledcorrbehav=zeros(0,0);
pooledcorrBOLD=zeros(0,0);
pooled_onsets=zeros(0,0);

for person_i = 1 : length( persons{1} )
    corct_BOLD_equal=zeros(0,1);
    corct_BOLD=zeros(0,1);
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
            if ~isempty(strfind(sessionData{behav_sess,2}{1,behav_look},behav))
                ons_behav=vertcat(ons_behav,transpose(sessionData{behav_sess,3}{1,behav_look}));
                behav_data=vertcat(behav_data,transpose(sessionData{behav_sess,6}{1,behav_look}));
                sess_ons=vertcat(sess_ons,transpose(sessionData{behav_sess,3}{1,behav_look}));
            end
        end
        person_sess_onsets{1,end+1}= sess_ons;
        corct_behav_data{1,end+1}=behav_data;
    end
    
    for i= size(corct_behav_data,2)-1 : size(corct_behav_data,2)
%     if isa(MA,'double') && MA~=0   
if FIR~=0
    filter= BOLD_filt(corct_behav_data{1,i});
    corct_behav_data{1,i}=filtfilt(filter,corct_behav_data{1,i});
end
%         max=size(corct_behav_data{1,i},1);
%         corct_behav_data{1,i}=corct_behav_data{1,i}-movmean(corct_behav_data{1,i},round(MA*max));
%     end

    if zscored==1
            corct_behav_data{1,i}=zscore(corct_behav_data{1,i});
    end
    end
    
    
    for ROIs_loop= 1: size(ROIs,2)
        if isempty(strfind(Subjects_BOLD{1,person_i}(ROIs_loop).ROI,ROIs))==0
            for conds_loop= 1:size(Conditions,2)
                if ~isempty(strfind(Subjects_BOLD{1,person_i}(ROIs_loop).Conditions{conds_loop},behav))
                    if size(Conditions,1)>1
                        if ~isempty(Conditions{2,conds_loop})
                            fprintf('\n Found a special condition "must qualify a specific name for the previous trial"')
                            for behav_sess3=1:size(sessionData,1)
                                current_sess=size(person_sess_onsets,2)-size(sessionData,1)+behav_sess3;
                                for get_behav_trial_sm=1:size(person_sess_onsets{1,current_sess},1)
                                    for get_behav_trial_lrg=1:size(allons_sort,1)
%                                         if get_behav_trial_lrg>size(person_sess_onsets{1,current_sess},1)
%                                             get_behav_trial_lrg=get_behav_trial_lrg-size(person_sess_onsets{1,current_sess},1);
%                                         end
                                        if allons_sort(get_behav_trial_lrg,1)==person_sess_onsets{1,current_sess}(get_behav_trial_sm,1)
                                            print_trialnm=['\n \n 1. found trial: ' num2str(get_behav_trial_sm) ' at  sorted index number ' num2str(get_behav_trial_lrg)];
                                            fprintf(print_trialnm)
                                            if get_behav_trial_lrg==1
                                                fprintf('\nthis is the very first trial, cannot know the state of a previous trial. Will exclude')
                                                person_sess_onsets{1,current_sess}(get_behav_trial_sm,1)=9999;
                                            else
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
                                                                if isempty(strfind(sessionData{behav_sess,2}{1,prv_trial_loc(trial_name_i,1)},'Outcome'))
                                                                    fprintf('\n 2b. trial onset is at a 2nd level regressor "Outcome"\n')
                                                                end
                                                                if isempty(strfind(trial_name,Conditions{2,conds_loop}))
                                                                    person_sess_onsets{1,current_sess}(get_behav_trial_sm,1)=9999;
                                                                    corct_behav_data{1,current_sess}(get_behav_trial_sm,1)=9999;
                                                                    print_conname=['---Condition "' Conditions{2,conds_loop} '" is not there'];
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
    scan_behav=round(ons_behav,-1)./2.5;
    Pooledcorrbehav=vertcat(Pooledcorrbehav,corct_behav_data_V);
    
    
    
    %getting relevant BOLD and behavioural data vector
    
    fprintf('\n Throwing unrelated scans \n \n ')
    for ROIs_loop= 1: size(ROIs,2)
        if isempty(strfind(Subjects_BOLD{1,person_i}(ROIs_loop).ROI,ROIs))==0
            for conds_loop= 1:size(Conditions,2)
                if ~isempty(strfind(Subjects_BOLD{1,person_i}(ROIs_loop).Conditions{conds_loop},behav))
                  if ~isempty(Subjects_BOLD{1,person_i}(ROIs_loop).BOLD{conds_loop})
                            BOLD_behav=Subjects_BOLD{1,person_i}(ROIs_loop).BOLD{conds_loop};
                            if isa(MA,'double') && MA~=0
                                max=size( BOLD_behav,1);
                                BOLD_behav=BOLD_behav-movmean(BOLD_behav,round(MA*max));
                            
                            sd2=2*std(BOLD_behav);
                            n=size(BOLD_behav,1);
                            m=mean(BOLD_behav,1);
                            for corct_outliers=1:n
                                if BOLD_behav(corct_outliers)>sd2 || BOLD_behav(corct_outliers)< -1*sd2
                                BOLD_behav(corct_outliers)=m;
                                end
                            end
                            end
                            if zscored==1
                                BOLD_behav=zscore(BOLD_behav);
                            end
                                 if FIR~=0
                                  filter= BOLD_filt(BOLD_behav,'LSD');
                                  BOLD_behav=filtfilt(filter,BOLD_behav);
                                 end
                     corct_BOLD=vertcat(corct_BOLD,BOLD_behav(scan_behav,1));
                   end
                end
            end
        end
    end
    
    %session relevant BOLD
    for sess2_i= 1:size(sessionData,1)
        person_sess_BOLD{1,end+1}=BOLD_behav(round(person_sess_onsets{1,size(person_sess_BOLD,2)+1},-1)./2.5);
    end
%     if constant==1
%         for sess2_i= 1:size(sessionData,1)
%        ones=vertcat(ones,ones(size(person_sess_BOLD,1),1)); 
%         end
%     end
    session_num(person_i,1)=size(sessionData,1);
    
    pooledcorrBOLD=vertcat(pooledcorrBOLD,corct_BOLD);
    
    if equal_regressors~=0
        fprintf('\n applying a special condition of averaging all two trials together \n')
        fprintf('This is because each two values of each condition is correlated with one happiness rating \n')
        fprintf('please make sure to modify this condition for other modelling puposes \n')
        
        if size(corct_BOLD,1)==((equal_regressors*size(sessionData,1))*2)
            temp_BOLD_min=zeros(0,1);
            for avg_loop= 1:2:size(corct_BOLD,1)
                temp_BOLD_min(end+1,1)=(corct_BOLD(avg_loop,1)+corct_BOLD(avg_loop+1,1))./2;
            end
            corct_BOLD_equal=temp_BOLD_min;
            clear temp_BOLD_min
            pooledcorrBOLD_equal=vertcat(pooledcorrBOLD_equal,corct_BOLD_equal);
        elseif size(corct_BOLD,1)>(equal_regressors*size(sessionData,1)) && size(corct_BOLD,1)<((equal_regressors*size(sessionData,1))*2)
            corct_BOLD_equal=corct_BOLD(1:equal_regressors,1);
            pooledcorrBOLD_equal=vertcat(pooledcorrBOLD_equal,corct_BOLD_equal);
        else
            pooledcorrBOLD_equal=pooledcorrBOLD;
            corct_BOLD_equal=corct_BOLD;
        end
    else
        pooledcorrBOLD_equal=pooledcorrBOLD;
        corct_BOLD_equal=corct_BOLD;
    end
    if size(Conditions,1)>1
        if isempty(Conditions{2,conds_loop})
            fprintf('\n Correlating \n');
            if fig==1
                figure;
                scatter(corct_BOLD,behav_data)
                %     figure;
                %     histogram(shifts)
            end
            CRB=corr(corct_BOLD,behav_data);
        end
    end
    
end
end