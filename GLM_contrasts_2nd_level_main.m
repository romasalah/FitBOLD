clear, %close all

% ----------- Main settings -----------------------------------------------
runBatch = 1;

% determine task(s)
% tasks = {'FFXspecAndEst'};
durationsAll = [0];
FFXfolderNames = {'FFX_CREV_minEV_3sep2018'};
printCons = 0;
mainFolder = '/Volumes/MYBOOK1';
codeDir = '/Users/NEMO_admin/Google Drive/Master(shared)/Code/';

% insert file of Rutledge  ISIs when they are ready
initialISIandSetupValuesFile = 'Rutledge_InitialISIs_final.txt';
HRFderivatives = [1 1];
paraMod = 1;
NtrialTypes = 14; NparamMod = 13;

TR = 2.5;
imgType = 's6wuf';
processDuration = NaN;
usefig = 0; % set to 0 for -nojvm start

% All subjects:
% Subjects to do:
persons{1} = [ 44 55 62 63 65 69 70 71 74 77 79 81 83 84 90 92 93 94 95 96 98 101 102 103 104 106 110 120 124 125 126 129 136 155 166 171 173 174 175 176 177 ];
Ns = sum(cellfun(@length,persons)); % N subjects
% ------------ for each run type: ------------------
kinds = {''};

% order of the runs, for each subject:
runNames = {...
  'rutledge1','run1';...
  'rutledge2','run2';...
  %'SoDec1','run1';...
  %'SoDec2','run2';...
  };
Nses = size(runNames,1);
reportFile = [codeDir '.txt'];
crashedList = [];
doWhat = 'FFXspecAndEst';
if runBatch
    spm( 'defaults', 'fmri' );
    spm_jobman( 'initcfg' );
    if usefig;
        f_hdle = figure(99); f_hdle.Position = [130 877 560 420]; f_hdle.Color = [1 1 1]; g_hdle = gca; g_hdle.Visible='off'; drawnow; pause(.1)
    end
end
    counter = 0;
    for person_i = 1 : length( persons{1} )
        
        t1 = clock;
        % determine person name
        personNo = persons{1}(person_i);
        personName = ['SODEC_FMRI_' num2str( personNo ) ];
        fprintf('Doing subject: %s\n',personName)
        counter = counter + 1;
        
        %     try % so that if analysis for 1 subject crashes, the whole batch analysis does not stop
        % report
        if runBatch
            fid = fopen( reportFile, 'a' );
            if usefig;
                figure(99); cla; text(.01,.5,string,'fontsize',36,'color',[0 0 1],'interpreter','none'); g_hdle = gca; g_hdle.Visible='off'; drawnow
            end
        end
        
        % FFX folder
        FFXdirName = [mainFolder filesep FFXfolderNames{1} filesep personName filesep];
        if ~exist( FFXdirName, 'dir' ); mkdir(FFXdirName); fprintf(fid, 'Erstelle Ordner %s/r/n', FFXdirName ); end
        fprintf('FFX dir name: %s\n',FFXdirName)
        
        % ===================== what to do? =======================
        clear matlabbatch sessionData regs
        for ses = 1:Nses
            clear scanInput sessMultiFile sessMultiregFile
            
            % Get scans
            scanFolder = [mainFolder '/fMRIdataComplete' filesep personName filesep runNames{ses,1} filesep];
            fprintf('fMRI data folder: %s\n',scanFolder)
            scanListing = dir( [scanFolder imgType '*.nii'] );
            scanInput = cell( length( scanListing ), 1 );
            for scanFile = 1 : length( scanListing )
                scanInput{ scanFile } = [scanFolder scanListing( scanFile ).name ',1'];
            end
            scans = sort( scanInput );
            sessionData{ses,1} = scans;
            
            % Get behavioural data RUTLEDGE
            behavFile = dir([mainFolder '/fMRIdataComplete' filesep personName filesep '*Decision_Rutledge*' runNames{ses,2} '*.mat']);
            fname = [behavFile.folder filesep behavFile.name];
            fprintf('Behaviour file: %s\n',fname)
            initialISIandSetupValues = dlmread(initialISIandSetupValuesFile);
            idx = find(initialISIandSetupValues(:,1)==personNo);
            ISI = initialISIandSetupValues(idx,ses+1);
            if isempty(ISI); ISI = NaN; end
            ISIs(ses) = ISI;
            [onsets,durations,trial_type_names,param_mods,missingOnsets(person_i,ses,:)] = Rutledge_trialtypes_noCompReg_3sep2018(fname,ISI);
            disp('N onsets per condition:')
            fprintf('%d ',cellfun(@length,onsets))
            
            % if no clean initialISI, set missingOnsets for all trialtypes to 1
            if isnan(ISI)
                missingOnsets(person_i,ses,:) = 1;
            end
            
            sessionData{ses,2} = trial_type_names;
            sessionData{ses,3} = onsets;
            sessionData{ses,4} = durations;
            
            % Get realignment parameters
            runNameFull = [mainFolder '/fMRIdataComplete' filesep personName filesep runNames{ses,1}];
            multiregFile = dir( [runNameFull '/rp*.txt' ] );
            if length(multiregFile)>0
                sessMultiregFile = [runNameFull filesep multiregFile.name];
            end
            sessionData{ses,5} = sessMultiregFile;
            
            % Get parametric modulation variables
            sessionData{ses,6} = param_mods;
        end
        
        if length(sessionData(1,:)) ~= 6; warning('Input data missing!');
            sessionData
        end
        
        if ~runBatch
            sessionData
        else
            % clean the FFX directory
            owd = pwd; cd(FFXdirName); delete *; cd(owd);
        end
        
        % throw away session information for session with no ISI
        sessionData = sessionData(find(~isnan(ISIs)),:);
        
        % specify multi-session model
        matlabbatch{1}.spm.stats.fmri_spec.dir = { FFXdirName };
        matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
        matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 1;
        for ses = 1:size(sessionData,1) % add session-specific infos
            matlabbatch{1}.spm.stats.fmri_spec.sess(ses).scans = sessionData{ses,1};
            % load(sessionData{ses,3}) % loads onsets, durations and difficulties (-> parametric modulators) into workspace
            % matlabbatch{1}.spm.stats.fmri_spec.sess(ses).cond = struct('name', names, 'onset', onsets, 'duration', durations, 'tmod', {}, 'pmod', struct('name', {'difficulty'}, 'param', difficulties, 'poly', {1}));
            for c = 1:length(sessionData{ses,2})
                matlabbatch{1}.spm.stats.fmri_spec.sess(ses).cond(c).name = sessionData{ses,2}{c};
                matlabbatch{1}.spm.stats.fmri_spec.sess(ses).cond(c).onset = sessionData{ses,3}{c};
                %               matlabbatch{1}.spm.stats.fmri_spec.sess(ses).cond(c).duration = durationsAll(tt) * ones(length(sessionData{ses,3}{c}),1);% durations{c};
                matlabbatch{1}.spm.stats.fmri_spec.sess(ses).cond(c).duration = sessionData{ses,4}{c};
                matlabbatch{1}.spm.stats.fmri_spec.sess(ses).cond(c).tmod = 0;
                if paraMod & ~isempty(sessionData{ses,6}{c})
                    matlabbatch{1}.spm.stats.fmri_spec.sess(ses).cond(c).pmod = struct('name', {'Val'}, 'param', sessionData{ses,6}{c}, 'poly', {1});
                else
                    matlabbatch{1}.spm.stats.fmri_spec.sess(ses).cond(c).pmod = struct('name', {}, 'Val', {});
                end
            end
            matlabbatch{1}.spm.stats.fmri_spec.sess(ses).regress = struct('name', {}, 'val', {});
            %             for rrr = 1:size(sessionData{ses,6},2)
            %               NtimePoints = length(sessionData{ses,1}); % this is Nscans
            %               dif = size(sessionData{ses,6},1) - NtimePoints;
            %               if dif < 0; sessionData{ses,6} = [sessionData{ses,6}; zeros(ceil(-dif),5)]; end
            %               matlabbatch{1}.spm.stats.fmri_spec.sess(ses).regress(rrr).name = regressorNames{rrr};
            %               matlabbatch{1}.spm.stats.fmri_spec.sess(ses).regress(rrr).val = sessionData{ses,6}(1:NtimePoints,rrr);
            %             end
            matlabbatch{1}.spm.stats.fmri_spec.sess(ses).multi_reg = sessionData(ses,5);
            matlabbatch{1}.spm.stats.fmri_spec.sess(ses).hpf = 128;
        end
        matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
        matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = HRFderivatives;
        matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
        matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
        matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
        matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
        
        
                  % estimate model
                  matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = { [FFXdirName 'SPM.mat'] };
                  matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
                  save([FFXdirName '/sessionData'],'sessionData')
                  save([FFXdirName '/matlabbatch'],'matlabbatch')
                  if runBatch
                    spm_jobman( 'run', matlabbatch );
                  end
    end
    %   end
    
    clear all
    %%********session details**********
    %To be set by the experimenter
    maxnumofsessions=2;
    maxnumofuniquecovariates= 2;
    global numofmovregressors;
    numofmovregressors= 6;
    maxnumofconditions=40;
    subjectsfolder='/Volumes/SSD/FFX_Rutledge_minEV_2aug2018_ISIfinal/Subjects'; %WARNING: remeber to remove the "/" at the end
    maxmatfolder='/Volumes/SSD/FFX_Rutledge_minEV_2aug2018_ISIfinal/Subjects/SODEC_FMRI_70';
    secondlevelfolder='/Volumes/SSD/FFX_Rutledge_minEV_2aug2018_ISIfinal/Second Level';
    scriptsfolder='/Volumes/SSD/FFX_Rutledge_minEV_2aug2018_ISIfinal/Scripts';
    numofrepeatednamecharacters=29; %must be set for now
    Movement_regressor_id='R'; %used to identify movement regressors
    nummulticontrasts=3;
    Excluded_ANOVA='SODEC_FMRI_87';
    
    %*******make multiple condition contrasts by name*******To be set by the experimenter
    con_names= {'CR all', 'EV all', 'PE all'};
    find_con = {'CR', 'EV', 'PE'};
    exc_con = {'Val', 'Val', 'Val'};%will exclude regressors which have this string in their name
    signs= { 'plus' 'minus'};
    
    %usually stable
    maxnumofcontrasts=maxnumofconditions+maxnumofuniquecovariates;
    numsinglecontrasts=maxnumofcontrasts;
    required_contrasts=numsinglecontrasts+nummulticontrasts;
    global Constant_regressor_id;
    Constant_regressor_id='constant'; %used to identify constant regressors
    Movement_regressor_id_correct='Movement Regressor';
    numofconstants=1;
    maxnumofcontrasts=maxnumofconditions+maxnumofuniquecovariates;
    maxnumofregressorspersession=maxnumofcontrasts+numofconstants+numofmovregressors;
    numofunusedregressors=numofmovregressors+numofconstants;
    maxnumofregressors= numofconstants+((maxnumofconditions+numofmovregressors)*maxnumofsessions)+maxnumofuniquecovariates;
    non_unique_regressor_id= 'Ununique Regressor';
    missing_condition_id= 'Misssing Condition';
    remove_permutation_duplicate='minus';
    
    %dummies
    invalid_regressors=cell(1,1);
    scans_idx=cell(0,0);
    mov_reg_names=cell(1,numofmovregressors);
    global mov_reg_names2
    mov_reg_names2=cell(1,0);
    ANOVA_dir=[secondlevelfolder '/ANOVA/'];
    mkdir(ANOVA_dir);
    ttest_dir=[secondlevelfolder '/Paired-t-tests/'];
    mkdir(ttest_dir);
    %*********looks for conditions and exports them********
    cd(subjectsfolder)
    Subjects=subdir(subjectsfolder);
    condition_names=cell(106,length(Subjects));
    fprintf('\n Creating conditions array \n')
    for person_i = 1 : length( Subjects ) %loop through the folders
        callperson = Subjects(person_i); %calls a specific folder according to index
        personName = callperson{1,1}; %get the directory as a character array
        fprintf(' \n Checking new subject \n')
        shortname= erase(personName,subjectsfolder);
        shortName= ['\n' shortname '\n \n'];
        fprintf(shortName);
        % if exist( personName, 'dir' ) ~= 0 % sets a condition tostart the
        % script only if
        cd (personName)
        M_i=open('SPM.mat');
        S=M_i.SPM.Vbeta;
        nrows=size(S,2);
        f= cellstr(char(S(:).descrip));
        dataend=size(f,1);
        for u= 1: dataend
            f{u}=f{u}(29:end); %removes the first repeated characters in condition names
        end
        Remaining= 106- size(f,1);
        for R_i = 1: Remaining
            startzeros=dataend+R_i;
            ccd=cell(Remaining,1);
            % f(dataend:67,1)= ccd(:,1);
            f=vertcat(f(1:dataend,1),ccd(:,1));
        end
        condition_names(:,person_i)= f(:,1);
        d=char(cd);
        
        %*****Getting condition names from the best matrix you have*****
        if d(1,end-5:end) == maxmatfolder(1,end-5:end)
            fprintf('\n Getting condition names \n')
            conditions_reference=cell(1,maxnumofcontrasts);
            for i=1:maxnumofregressors
                conditions_reference{i}=condition_names{i};  % Remove the repeated characters in the contrast name.
            end
        end
    end
    cd ..
    %xlwrite('conditions.xls',C);
    
    fprintf('\n Conditions array created \n')
    
    %C=cell(maxnumofregressors,length(Subjects));
    % for t = 1:maxnumofconditions
    %  conditions_last_clean{=cell(conditions_last{);
    
    fprintf('\n\nPreparing the batch settings for \n \n');
    
    matching_array=cell(maxnumofregressors+1, length( Subjects )+1);
    for fill_matching_array = 1 : size(matching_array,1)
        for fill_matching_array_col = 2: length( Subjects )+1
            matching_array{fill_matching_array,fill_matching_array_col}=cell(1,0);
        end
    end
    matching_array{1,1}='Conditions/Names';
    CV=cell(maxnumofregressorspersession+nummulticontrasts,length( Subjects )+1);
    
    for person_i = 1 : length( Subjects ) %loop through the folders
        callperson = Subjects(person_i); %calls a specific folder according to index
        personName = callperson{1,1}; %get the directory as a character array
        cd (personName) % go to the directory
        matdir=[ personName '/SPM.mat']; %specifies the required matrix to open
        shortname= erase(personName,subjectsfolder);
        shortName= ['\n' shortname '\n \n'];
        fprintf(shortName);
        person_col=person_i+1;
        if exist( matdir ) ==0 %#ok<EXIST>
            fprintf('Cant find SPM.mat file \n \n');
            
        else
            % sets a condition to start the script only if the matrix exists
            s=load('SPM.mat');
            ss= s.SPM.Vbeta;
            
            nrows=size(ss,2);
            numofcontrasts=nrows-numofunusedregressors;
            %person_idx=find(Subjects{1,:},shortname);
            
            %identify non-unique regressors
            if isfield(s.SPM.xX,'xKXs') && ...
                    ~isempty(s.SPM.xX.xKXs) && isstruct(s.SPM.xX.xKXs)
                iX = 1;
                [nScan,nPar] = size(s.SPM.xX.xKXs.X);
            elseif isfield(s.SPM.xX,'X') && ~isempty(s.SPM.xX.X)
                iX = 0;
                [nScan,nPar] = size(s.SPM.xX.X);
            else
                error('Can''t find DesMtx in this structure!')
            end
            if iX,  est = spm_SpUtil('IsCon',s.SPM.xX.xKXs);
            else
                est = spm_SpUtil('IsCon',s.SPM.xX.X);
            end
            idx_ununique=zeros(1,0);
            for idx_ununique_i= 1:nPar
                if est(1,idx_ununique_i)==0
                    idx_ununique(1,end+1)= idx_ununique_i;
                end
            end
            %******Create the matching array*******
            matching_array{1,person_col}=cellstr(shortname);
            for condition_index= 1:size(conditions_reference,2) %sets a condition to start searching for from the reference array.
                for target_index = 1:nrows %sets a line in the target array to look in.
                    %repair subscripted condition name from the original matrix
                    l=char(ss(target_index).descrip);
                    l_last= l(numofrepeatednamecharacters:end);
                    % compare the two names from target array and the reference
                    % array
                    if strcmp(l_last,conditions_reference{condition_index}) == 1
                        t=['Condition ' conditions_reference{condition_index} ' found at row ' num2str(target_index) '\n'];
                        fprintf(t)
                        if condition_index==1
                            matching_array{condition_index+1,1} = conditions_reference{condition_index};
                            matching_array{condition_index+1,person_col}{1,end+1}= target_index;
                        elseif condition_index>1
                            repeat_res=zeros(1,size(matching_array,2));
                            for avoid_repeat_i= 2:size(matching_array,1)
                                repeat_res(1,avoid_repeat_i)=sum(strcmp(matching_array{avoid_repeat_i,1},conditions_reference{condition_index}));
                            end
                            if sum(repeat_res,2)==0
                                matching_array{condition_index+1,1} = conditions_reference{condition_index};
                                matching_array{condition_index+1,person_col}{1,end+1}= target_index;
                            else
                                first_matching_row=find(repeat_res,1);
                                matching_array{condition_index+1,person_col}{1,end+1}= target_index;
                            end
                        end
                    end
                end
            end
            
            %******Start dealing with non-interesting regressors******
            mov_reg_names=cell(1,numofmovregressors);
            for rr = 1: numofmovregressors
                mov_reg_names{1,rr}=[Movement_regressor_id num2str(rr)];
            end
            
            %start searching and correction in the reference array
            conditions_reference_mov_corct= conditions_reference;
            for reference_loop = 1:size(conditions_reference,2)
                for mov_reg_i = 1: size(mov_reg_names,2)
                    if strfind(conditions_reference{1,reference_loop},mov_reg_names{1,mov_reg_i})==1
                        conditions_reference_mov_corct{1,reference_loop}= [Movement_regressor_id_correct num2str(mov_reg_i)];
                    end
                end
            end
            
            %getting indices and corrects the name of Movement regressors in
            %the matching matrix
            
            mov_reg_i2= zeros(1,numofmovregressors);
            mov_reg_i3= zeros(1,numofmovregressors);
            for matching_loop = 1:size(matching_array,1)
                for mov_reg_i = 1: size(mov_reg_names,2)
                    if matching_loop < (maxnumofregressorspersession-numofconstants)+1 %the additional names row
                        if strfind(matching_array{matching_loop,1},mov_reg_names{1,mov_reg_i})==1
                            mov_reg_i2(1,mov_reg_i)= matching_loop;
                            matching_array{matching_loop,1}=[Movement_regressor_id_correct num2str(mov_reg_i)];
                            mov_reg_names2{1,end+1}=matching_array{matching_loop,1};
                        end
                    else
                        if strfind(matching_array{matching_loop,1},mov_reg_names{1,mov_reg_i})==1
                            fprintf('found')
                            mov_reg_i3(1,mov_reg_i)= matching_loop;
                            matching_array{matching_loop,1}=[Movement_regressor_id_correct num2str(mov_reg_i)];
                        end
                    end
                end
            end
            
            %corrects the name and contrast vector of ununique regressors
            matching_array_corct=matching_array;
            if isempty(idx_ununique)==0
                for matching_loop = 1:size(matching_array,1)
                    for ununique_repair = 1: size(idx_ununique,2)
                        for matching_col_i = 1:size(matching_array{matching_loop,person_col},2)
                            if isempty(matching_array{matching_loop,person_col}{1,matching_col_i})==0
                                ee=matching_array{matching_loop,person_col}{1,matching_col_i};
                                if ee == idx_ununique(1,ununique_repair)
                                    matching_array_corct{matching_loop,person_col}{1,matching_col_i}= [];
                                end
                            end
                        end
                    end
                end
            end
            
            % now cut out the empty conditions
            non_missing_idx=cell(1,0);
            for find_empty_i =1:size(matching_array_corct,1)
                if ~isempty(matching_array_corct{find_empty_i,1})==1
                    non_missing_idx{1,end+1}=find_empty_i;
                end
            end
            non_missing_idx=cell2mat(non_missing_idx);
            matching_array_corct=matching_array_corct(non_missing_idx,:);
            
            %*****Starting making the batch file******
            fprintf('\n\n Cooking the batch file \n')
            
            %specifies SPM.mat address
            all_con_batch{1}.spm.stats.con.spmmat = cellstr(matdir);
            
            %*******make contrast vectors******
            CV{1,1}=matching_array_corct{1,1};
            CV{1,person_col}=matching_array_corct{1,person_col};
            %invalid_regressors{:,1}=matching_array{:,1};
            fprintf('\n\n Arranging single condition contrasts \n')
            
            %Simple single vector contrasts
            scans_idx(1,:)=matching_array_corct(1,:); %includes person names in scans array
            for CV_condition= 2:size(matching_array_corct,1)
                CV{CV_condition,1}=matching_array_corct{CV_condition,1};
                reg_result=isreg(matching_array_corct{CV_condition,1});
                if reg_result==1
                    scans_idx{CV_condition,1}=matching_array_corct{CV_condition,1}; %removes the names of movment an contrast regressors from the scans array
                end
                CV{CV_condition,person_col}=zeros(1,nrows);
                
                %Arbitrary contrast vector for nonunique regressors
                if sum(cellfun(@isempty,matching_array_corct{CV_condition,person_col}))==size(matching_array_corct{CV_condition,person_col},2)
                    %matching_array_corct{CV_condition,1}=
                    %non_unique_regressor_id;
                    CV{CV_condition,person_col}(1,end)=1;
                    invalid_regressors{end+1,1}=shortname;
                    invalid_regressors{end,2}=CV_condition;
                    
                    %Arbitrary contrast vector for offset regressors
                elseif matching_array_corct{CV_condition,1}(end-7:end)== Constant_regressor_id(end-7:end)
                    for CV_index = 1:size(matching_array_corct{CV_condition,person_col},2)
                        if isempty(matching_array_corct{CV_condition,person_col}(1,CV_index))==0
                            location= matching_array_corct{CV_condition,person_col}(1,CV_index);
                            location=location{1,1};
                            CV{CV_condition,person_col}(1,location)= 1;
                        end
                    end
                    %Arbitrary contrast vector for movement regressors
                elseif matching_array_corct{CV_condition,1}(end-7:end)==Movement_regressor_id_correct(end-7:end)
                    for CV_index = 1:size(matching_array_corc{CV_condition,person_col},2)
                        if isempty(matching_array_corct{CV_condition,person_col}(1,CV_index))==0
                            location= matching_array_corct{CV_condition,person_col}(1,CV_index);
                            location=location{1,1};
                            CV{CV_condition,person_col}(1,1)= 1;
                        end
                    end
                    
                    %True contrast vectors
                elseif matching_array_corct{CV_condition,1}(end-7:end)==conditions_reference_mov_corct{1,CV_condition-1}(end-7:end) %the reference matrix doesn't have a row for person names, you must subtract 1.
                    for CV_index = 1:size(matching_array_corct{CV_condition,person_col},2)
                        if isempty(matching_array_corct{CV_condition,person_col}(1,CV_index))==0
                            location= matching_array_corct{CV_condition,person_col}(1,CV_index);
                            location=location{1,1};
                            CV{CV_condition,person_col}(1,location)= 1;
                            if CV_condition<=size(scans_idx,1)
                                reg_result=isreg(matching_array_corct{CV_condition,1});
                                if reg_result==1
                                    scans_idx{CV_condition,person_col}=CV_condition;
                                end
                            end
                        end
                    end
                end
            end
            %Trimming CV
            if person_i==1
                for CV_str_i = 1: size(CV,1)
                    if isempty(CV{CV_str_i,person_col})==0
                        CV_str{1,CV_str_i}=cellstr(CV{CV_str_i,1});
                    end
                end
                
                cut_CV_idx=find(~cellfun(@isempty,CV_str));
                CV=CV(cut_CV_idx,:); %cuts empty rows from CV
                CV_end_idx=cut_CV_idx(1,end);
            end
            % % %         %makes one-way multiple condition contrasts by name
            % % %         fprintf('\n\n Now doing one-way multiple condition contrasts \n')
            % % %         %             Starts finding conditions
            % % %         for con_names_i = 1:size(con_names,2)
            % % %             for regressor_sign_i = 1:2
            % % %                 regressor_sign= signs{1,regressor_sign_i};
            % % %                 regressor_func= str2func(char(signs(1,regressor_sign_i)));
            % % %                 %adding signs to the names
            % % %                 con_name_sign=con_names;
            % % %                 location2=cell(1,0);
            % % %                 for str_find_i = 1: size(matching_array_corct,1)
            % % %                     if isempty(strfind(matching_array_corct{str_find_i,1},find_con{con_names_i}))==0 && ...
            % % %                             isempty(strfind(matching_array_corct{str_find_i,1},exc_con{con_names_i}))==1
            % % %                         if isempty(matching_array_corct{str_find_i, person_col})==1
            % % %                             missing_condition = [ 'n\ Condition' matching_array_corct{str_find_i,1} 'is missing in subject' shortname 'n\'];
            % % %                             fprintf(missing_condition)
            % % %                         else
            % % %                             for CV_index = 1:size(matching_array{str_find_i,person_col},2)
            % % %                                 if isempty(matching_array_corct{str_find_i,person_col}(1,CV_index))==0
            % % %                                     location2(1,end+1)= matching_array_corct{str_find_i,person_col}(1,CV_index);
            % % %                                 end
            % % %                             end
            % % %                         end
            % % %                     end
            % % %                 end
            % % %                 %make the new names only with the first participant
            % % %                 if person_i==1
            % % %                     for sign_name_i = 1:size(con_names,2)
            % % %                         if  strcmp(regressor_sign,'minus')==1
            % % %                             con_name_sign{1,sign_name_i}= [regressor_sign ' ' con_names{1,sign_name_i} ];
            % % %                         end
            % % %                     end
            % % %                     CV{end+1,1}= con_name_sign{con_names_i};
            % % %                     CV{end,person_col}=zeros(1,nrows);
            % % %                     for fill_CV = 1: size(location2,2)
            % % %                         CV_location= location2(1,fill_CV);
            % % %                         CV_location= CV_location{1,1};
            % % %                         CV{end,person_col}(1,CV_location)= regressor_func(0,1);
            % % %                     end
            % % %                 else
            % % %                     for CV_sec_i= 1: con_names_i*regressor_sign_i
            % % %                         CV{CV_end_idx+CV_sec_i,person_col}=zeros(1,nrows);
            % % %                         for fill_CV = 1: size(location2,2)
            % % %                             CV_location= location2(1,fill_CV);
            % % %                             CV_location= CV_location{1,1};
            % % %                             CV{CV_end_idx+CV_sec_i,person_col}(1,CV_location)= regressor_func(0,1);
            % % %                             CV_end_multi=CV_end_idx+CV_sec_i;
            % % %                         end
            % % %                     end
            % % %                 end
            % % %             end
            % % %         end
            % % %         CV_end_multi=CV_end_idx+(con_names_i*regressor_sign_i);
            % % %         fprintf('\n\n making Paired-t-tests contrasts \n')
            % % %         %making Paired-t-tests contrasts
            % % %         if person_i==1
            % % %             for rename_minus= 1:size(signs,2)
            % % %                 con_name_sign{1,rename_minus}= ['minus ' con_names{1,rename_minus} ];
            % % %             end
            % % %
            % % %             idx_minus=cell(1,0);
            % % %             idx_pos=cell(1,0);
            % % %             fprintf('\n\n Calculating Permutations \n')
            % % %             %getting the condition to permute with
            % % %             for minus_i= 1:size(con_name_sign,2)
            % % %                 for get_con_i =1: size(CV,1)
            % % %                     if strfind(CV{get_con_i,1},con_name_sign{1,minus_i})==1
            % % %                         idx_minus{1,end+1}=get_con_i;
            % % %                     end
            % % %                 end
            % % %             end
            % % %             for pos_i=1:size(con_names,2)
            % % %                 for get_con_i2 =1: size(CV,1)
            % % %                     if strfind(CV{get_con_i2,1},con_names{1,pos_i})==1
            % % %                         idx_pos{1,end+1}=get_con_i2;
            % % %                     end
            % % %                 end
            % % %             end
            % % %             idx_minus=cell2mat(idx_minus);
            % % %             idx_pos=cell2mat(idx_pos);
            % % %             all_signs=[idx_minus idx_pos];
            % % %         end
            % % %         Permutations=nchoosek(all_signs,2);
            % % %         %remove empty
            % % %
            % % %         %Start making the contrast vectors
            % % %         for CV_paired_i= 1:size(Permutations,1)
            % % %             con_name_P=[CV{Permutations(CV_paired_i,1),1} '  ' CV{Permutations(CV_paired_i,2),1}];
            % % %             CV_target=cellstr(con_name_P);
            % % %             min_num=strfind(CV_target,remove_permutation_duplicate); %remove permutaion contrasts with double minuses to avoid duplicates
            % % %                             CV_P_o=Permutations(CV_paired_i,1);
            % % %                 CV_P_t=Permutations(CV_paired_i,2);
            % % %              CV_P_check= CV{CV_P_o,person_col}+...
            % % %                 CV{CV_P_t,person_col};
            % % %
            % % %             if person_i==1
            % % %                 if size(min_num{1,1},2)==2 || isempty(find(CV_P_check,1))==1
            % % %                     Permutations(CV_paired_i,:)=0;
            % % %                 end
            % % %             end
            % % %         end
            % % %          if person_i==1
            % % %             nn(:,1)=nonzeros(Permutations(:,1));
            % % %             nn(:,2)=nonzeros(Permutations(:,2));
            % % %             Permutations=nn;
            % % %         end
            % % %         for fill_CV_paired= 1: size(nn,1)
            % % %             if person_i==1
            % % %                 CV_P_one=Permutations(fill_CV_paired,1);
            % % %                 CV_P_two=Permutations(fill_CV_paired,2);
            % % %                 CV_P= CV{CV_P_one,person_col}+...
            % % %                 CV{CV_P_two,person_col};
            % % %                 CV{CV_end_multi+fill_CV_paired,person_col}=CV_P;
            % % %                 CV{CV_end_multi+fill_CV_paired,person_col}=CV_P;
            % % %                 con_name_P=[CV{Permutations(fill_CV_paired,1),1} '  ' CV{Permutations(fill_CV_paired,2),1}];
            % % %                 CV{CV_end_multi+fill_CV_paired,1}=cellstr(con_name_P);
            % % %             else
            % % %             CV_P_one=Permutations(fill_CV_paired,1);
            % % %             CV_P_two=Permutations(fill_CV_paired,2);
            % % %             CV_P= CV{CV_P_one,person_col}+...
            % % %                 CV{CV_P_two,person_col};
            % % %                 CV{CV_end_multi+fill_CV_paired,person_col}=CV_P;
            % % %                 con_name_P=[CV{Permutations(fill_CV_paired,1),1} '  ' CV{Permutations(fill_CV_paired,2),1}];
            % % %             end
            % % %         end
            %Re-trim CV
            if person_i==1
                for CV_str_i2 = 1: size(CV,1)
                    if isempty(CV{CV_str_i2,person_col})==0
                        CV_str{1,CV_str_i2}=cellstr(CV{CV_str_i2,1});
                    end
                end
                cut_CV_idx=find(~cellfun(@isempty,CV_str));
                CV=CV(cut_CV_idx,:); %cuts empty rows from CV
            end
            %CV_end_idx=CV_end_idx-size(CV_end_idx,2);
            fprintf('\n\n All Contrasts Prepared \n')
            all_con_batch{1}.spm.stats.con.delete = 1;
            %print names and contrast vectors to the batch
            for batch_loop = 2:size(CV,1)
                if isempty(CV{batch_loop,person_col})==0
                    all_con_batch{1}.spm.stats.con.consess{batch_loop-1}.tcon.name = char(CV{batch_loop,1});
                    all_con_batch{1}.spm.stats.con.consess{batch_loop-1}.tcon.convec = CV{batch_loop,person_col};
                end
            end
            
            %*******Saving the batch*****
            matlabbatch=all_con_batch;
            fid = fopen( 'matlabbatch_scripted.mat', 'w' );
            fprintf(' \n Saving the batch to a file \n\n ');
            save('matlabbatch_scripted.mat','matlabbatch'); %saves the batch file
            
            %********Start SPM job manager*********
            jobdir=[personName '/matlabbatch_scripted.mat'];
            nrun=1;
            jobfile = jobdir;
            jobs=repmat(jobfile, 1, nrun);
            inputs= cell(1, nrun);
            spm( 'defaults', 'FMRI' );
            spm_jobman('run',jobs,inputs{:});
            
        end
    end
    fprintf(' \n Contrasts were sucessfuly created \n')
    
    %****Start second level analysis****
    fprintf(' \n Starting the second level analysis \n')
    fprintf(' \n An ANOVA model will be created first  \n')
    clear matlabbatch %or change the name
    cd(ANOVA_dir);
    delete *.nii
    cd(subjectsfolder);
    %****get the directories of the valid scans ande exclude inavalid ones*****
    fprintf(' \n getting scans directories \n')
    Subjects=subdir;
    % Subject=Subjects{1:42}; %make sure no other folders are present except
    % subjects folders, otherwise use that.
    scans=cell(0,0);
    scans(:,1)=scans_idx(:,1);
    if size(scans_idx,2)~=length(Subjects)+1
        fprintf('Aborting: One or more subjects are missing')
    else
        for condition_scans_i= 2:size(scans_idx,1)
            for    person_scans_i= 2:size(scans_idx,2)
                scans{1,person_scans_i}=scans_idx{1,person_scans_i};
                if isempty(scans_idx{condition_scans_i,person_scans_i})==0
                    req_con=scans_idx{condition_scans_i,person_scans_i}-1;
                    condir=[Subjects{1,person_scans_i-1} '/con_*0' mat2str(req_con) '.nii'];
                    condir_s=dir(condir);
                    if isempty(condir_s)==1
                        scans{condition_scans_i,person_scans_i}=[];
                    else
                        condir_n=condir_s.name;
                        condir_f=condir_s.folder;
                        condir_all=[condir_f '/' condir_n ];
                        scans{condition_scans_i,person_scans_i}=condir_all;
                    end
                end
            end
        end
    end
    
    %*******make ANOVA and t-test factorial design batch****
    fprintf(' \n Cooking the batch file \n')
    matlabbatch{1}.spm.stats.factorial_design.dir={ANOVA_dir};
    for condition_cell_i=2:size(scans,1)
        matlabbatch{1}.spm.stats.factorial_design.des.anova.icell(condition_cell_i-1).scans=[];
        ttestbatch{condition_cell_i-1}.spm.stats.factorial_design.des.t1.scans=[];
        newdir=[ttest_dir scans{condition_cell_i,1}(1:end-6)];
        if exist(newdir,'dir')==0
            mknewdir=mkdir(newdir);
        end
        for fill_scans= 2:size(scans,2)
            if isempty(scans{condition_cell_i,fill_scans})==0 %&& strcmp(scans{1,fill_scans},Excluded_ANOVA)==0 used to exclude subjects
                matlabbatch{1}.spm.stats.factorial_design.des.anova.icell(condition_cell_i-1).scans{end+1,1}= scans{condition_cell_i,fill_scans};
                %results directory
                ttestbatch{condition_cell_i-1}.spm.stats.factorial_design.des.t1.scans{end+1,1}=scans{condition_cell_i,fill_scans};
                ttestbatch{condition_cell_i-1}.spm.stats.factorial_design.dir=cellstr(newdir);
            end
        end
    end
    
    %*******Saving ANOVA batch*****
    cd(ANOVA_dir);
    fid2 = fopen( 'matlabbatch_scripted_ANOVA.mat', 'w' );
    fprintf(' \n Saving the batch to a file \n\n ');
    save('matlabbatch_scripted_ANOVA.mat','matlabbatch');
    
    %*******Saving t-test design batch*****
    clear matlabbatch
    matlabbatch=ttestbatch;
    cd(ttest_dir);
    delete *.nii
    fid3= fopen( 'matlabbatch_scripted_ttest.mat', 'w' );
    fprintf(' \n Saving matlabbatch_scripted_ttest.mat \n\n ');
    save('matlabbatch_scripted_ttest.mat','matlabbatch');
    
    %******estimate ANOVA & P-t-tests model******
    ANOVA_matdir=[ANOVA_dir 'SPM.mat'];
    Pt_ANOVA_estbatch{1}.spm.stats.fmri_est.spmmat=cellstr(ANOVA_matdir);
    Pt_ANOVA_estbatch{1}.spm.stats.fmri_est.write_residuals = 0;
    Pt_ANOVA_estbatch{1}.spm.stats.fmri_est.method.Classical=1;
    
    for est_i=1:size(ttestbatch,2)
        batch_loc=est_i+1;
        est_dir=[ttestbatch{est_i}.spm.stats.factorial_design.dir{1,1}];
        cd(est_dir);
        delete *.nii
        if exist('SPM.mat')~=0
            delete SPM.mat
        end
        est_matdir=[est_dir '/SPM.mat'];
        Pt_ANOVA_estbatch{batch_loc}.spm.stats.fmri_est.spmmat{1,1}=est_matdir; %#ok<*SAGROW>
        Pt_ANOVA_estbatch{batch_loc}.spm.stats.fmri_est.write_residuals = 0;
        Pt_ANOVA_estbatch{batch_loc}.spm.stats.fmri_est.method.Classical=1;
        Pt_con_batch{est_i}.spm.stats.con.spmmat = cellstr(est_matdir);
        Pt_con_batch{est_i}.spm.stats.con.consess{1}.tcon.name = char(erase(est_dir,ttest_dir));
        Pt_con_batch{est_i}.spm.stats.con.consess{1}.tcon.convec = 1;
    end
    %*******Saving t-test design batch*****
    clear matlabbatch
    matlabbatch=Pt_ANOVA_estbatch;
    cd(ttest_dir);
    fid4= fopen( 'matlabbatch_scripted_Pt_ANOVA_est.mat', 'w' );
    fprintf(' \n Saving matlabbatch_scripted_Pt_ANOVA_est.mat \n\n ');
    save('matlabbatch_scripted_Pt_ANOVA_est.mat','matlabbatch');
    
    % %*******Saving Pt-test contrast batch*****
    clear matlabbatch
    matlabbatch=Pt_con_batch;
    cd(ttest_dir);
    fid5= fopen( 'matlabbatch_scripted_Pt_con.mat', 'w' );
    fprintf(' \n Saving matlabbatch_scripted_Pt_con.mat \n\n ');
    save('matlabbatch_scripted_Pt_con.mat','matlabbatch');
    
    %save and go back to root
    cd(ANOVA_dir)
    if exist('SPM.mat')~=0
        delete SPM.mat
    end
    cd(ttest_dir);
    delete *.nii
    if exist('SPM.mat')~=0
        delete SPM.mat
    end
    save('allvar.mat');
    cd(subjectsfolder);
    
    %******Run all jobs*****
    spm( 'defaults', 'FMRI' );
    nrun=1;
    jobdir=cell(0,0);
    jobdir{1}=[ ttest_dir 'matlabbatch_scripted_ttest.mat']; %run t-test design job
    jobdir{2}=[ ANOVA_dir 'matlabbatch_scripted_ANOVA.mat']; %run ANOVA design job
    jobdir{3}=[ ttest_dir 'matlabbatch_scripted_Pt_ANOVA_est.mat']; %run P-ttest and ANOVA estimation job
    jobdir{4}=[ ttest_dir 'matlabbatch_scripted_Pt_con.mat']; %run pt-test contrast job
    for job_i=1:size(jobdir,2)
        jobs = {jobdir{job_i}};
        spm_jobman('run', jobs);
    end
    
    %***Functions***
    function checkreg= isreg(x)
    global numofmovregressors
    global mov_reg_names2
    global Constant_regressor_id
    is_it_reg=zeros(1,numofmovregressors);
    for ref_reg= 1: numofmovregressors
        if strcmp(x, mov_reg_names2{1,ref_reg})==1
            is_it_reg(1,ref_reg)=1;
        end
    end
    if  sum(is_it_reg,2)==0 && strcmp(x, Constant_regressor_id)==0
        result=1;
    else
        result =0;
    end
    checkreg=result;
    end