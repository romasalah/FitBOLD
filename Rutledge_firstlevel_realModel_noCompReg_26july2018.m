
%% FFX batch code for SocialDecision experiments Summer 2017
% J Schultz based on code by Paul Jung

% copy of SoDec Code to translate to Rutledge

clear, %close all

% ----------- Main settings -----------------------------------------------
runBatch = 1;

% determine task(s)
% tasks = {'FFXspecAndEst','contrasts'};
% tasks = {'calculateModelFit'};
% tasks = {'checkRealignmentParameters'};
tasks = {'FFXspecAndEst'};
% tasks = {'contrasts'};

durationsAll = [0];
FFXfolderNames = {...
  'FFX_Rutledge_minEV_31july2018_ISIfinal',... % CR paramMod = CR; durations = 0; with decisions; 
%   'FFX_Rutledge_noCompRegs2',... % CR paramMod = CR-EV; no decisions; durations = RT
%   'FFX_Rutledge_noCompRegs3',... % CR paramMod = CR-EV; with decisions; durations = 0
  };
printCons = 0;
mainFolder = '/Volumes/MYBOOK1';
% mainFolder = '/Users/johannesschultz/Desktop/fakeData';
codeDir = '/Users/NEMO_admin/Google Drive/Master(shared)/Code/';

% insert file of Rutledge  ISIs when they are ready
initialISIandSetupValuesFile = 'Rutledge_InitialISIs_final.txt';
%initialISIandSetupValuesFile = 'RutledgeISIsNoDoubt.txt';
HRFderivatives = [0 0];
paraMod = 1;
NtrialTypes = 14; NparamMod = 13;

TR = 2.5;
imgType = 's6wuf';
processDuration = NaN;
usefig = 0; % set to 0 for -nojvm start

% All subjects:
% persons{1} = [ 44 62 63 65 69 70 71 74 77 78 79 81 83 84 87 90 92 93 94 95 96 98 101 102 103 104 106 110 120 121 124 125 126 129 136 155 166 171 173 174 175 176 177 ]; % all participants

% Subjects to do:
 persons{1} = [ 44 55 62 63 65 69 70 71 74 77 79 81 83 84 90 92 93 94 95 96 98 101 102 103 104 106 110 120 124 125 126 129 136 155 166 171 173 174 175 176 177 ]; % all participants
% persons{1} = [ 44 ]; % all participants
% persons{1} = [ 71 85 94 104 129 174 ]; % all participants with both runs with clean best ISIs
%  persons{1} = [ 55 63 65 70 71 74 83 85 87 94 98 101 102 104 110 124 129 136 171 174 ]; % all participants with clean best ISIs in at least 1 session and fMRIregressors
% persons{1} = [ 44 55 62 63 65 69 70 71 74 77 78 79 81 83 84 87 90 92 93 94 95 96 98 101 102 103 104 106 110 120 121 124 125 126 129 136 155 166 171 173 174 175 176 177 ]; % all participants with clean best ISIs in at least 1 session and fMRIregressors
% persons{1} = [ 174 ]; % all participants with clean best ISIs in at least 1 session and fMRIregressors
Ns = sum(cellfun(@length,persons)); % N subjects

%  ---------- contrasts ---------------------------------------------------
deleteCons = 1; % deletes all already calculated contrasts

% T-test over all conditions:
psfileNewNameRoot = 'SPMresults ';
conType = 'T';

conNo = 1;
Ncond = 7; Nbf = sum(HRFderivatives)+1; Nruns = 4; NregPerRun = Ncond*Nbf*(paraMod+2)+6;

% % contrasts (for 1 run; gets expanded for more runs below)
cons = {...
    'CR all',                               [zeros(1,0)  repmat([1 0],1,3) zeros(1,21+6)];...
    'EV all',                               [zeros(1,6)  repmat([1 0],1,3) zeros(1,15+6)];...
    'CR all - EV all',                      [zeros(1,0)  kron([1 1 1 -1 -1 -1],[1 0]) zeros(1,15+6)];...
    'PE all',                               [zeros(1,12) repmat([1 0],1,3) zeros(1,9+6)];...
    'Outcome self>other active',            [zeros(1,18) kron([1 0 0 0],[1 0]) zeros(1,1+6)];...
    'Outcome self>other passive',           [zeros(1,18) kron([0 1 0 0],[1 0]) zeros(1,1+6)];...
    'Outcome self<other active',            [zeros(1,18) kron([0 0 1 0],[1 0]) zeros(1,1+6)];...
    'Outcome self<other passive',           [zeros(1,18) kron([0 0 0 1],[1 0]) zeros(1,1+6)];...
    'Outcome self>other active-passive',    [zeros(1,18) kron([1 -1 0 0],[1 0]) zeros(1,1+6)];...
    'Outcome self<other active-passive',    [zeros(1,18) kron([0 0 1 -1],[1 0]) zeros(1,1+6)];...
    'PE active > passive',                  [zeros(1,12) 1 0 -1 0 0 0 zeros(1,9+6)];...
    'PE active > non-soc',                  [zeros(1,12) 1 0 0 0 -1 0 zeros(1,9+6)];...
    'Own > other Decision'                  [zeros(1,0)  kron([1 -2 1 1 -2 1],[1 0]) zeros(1,9+6)];...
    'Other > Own Decision'                  [zeros(1,0)  kron([-1 2 -1 -1 2 -1],[1 0]) zeros(1,9+6)];...
    };
conNamesANOVA = {...
    'CR active';...
    'CR active paraMod';...
    'CR passive';...
    'CR passive paraMod';...
    'CR n-soc';...
    'CR n-soc paraMod';...
    'EV active';...
    'EV active paraMod';...
    'EV passive';...
    'EV passive paraMod';...
    'EV n-soc';...
    'EV n-soc paraMod';...
    'PE active';...
    'PE active paraMod';...
    'PE passive';...
    'PE passive paraMod';...
    'PE n-soc';...
    'PE n-soc paraMod';...
    'Outcome self>other active';...
    'Outcome self>other active paraMod';...
    'Outcome self>other passive';...
    'Outcome self>other passive paraMod';...
    'Outcome self<other active';...
    'Outcome self<other active paraMod';...
    'Outcome self<other passive';...
    'Outcome self<other passive paraMod';...
    'happiness';...
    };
conValsANOVA = [eye(length(conNamesANOVA)) zeros(length(conNamesANOVA),6)];

% for c=1:length(conNamesANOVA)
%   i = size(cons,1);
%   cons{i+1,1} = conNamesANOVA{c};
%   cons{i+1,2} = conValsANOVA(c,:);
% end
%% Settings unlikely to change

% get matrix of optimal onset corrections:
% correctionFile = '/Users/johannesschultz/Documents/exp/KINO/Matlab_Skript/check_residuals/bestOnsetCorrectionMatrix_BA17and18';
% load(correctionFile) % gets resultsMatrix; subjNo x runNo, has stored optimal onsetTypes
kinds = {''};

% order of the runs, for each subject:
runNames = {...
  'rutledge1','run1';...
  'rutledge2','run2';...
  %'SoDec1','run1';...
  %'SoDec2','run2';...
  };
Nses = size(runNames,1);

%% Setup finished, get to work
for tt = 1:length(tasks)
  reportFile = [codeDir '.txt'];
  crashedList = [];
  doWhat = tasks{tt};
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
    
    try % so that if analysis for 1 subject crashes, the whole batch analysis does not stop
      % report
      if runBatch
        fid = fopen( reportFile, 'a' );
        string = sprintf('task %d/%d\ndoing %s\nsubj %d/%d\n%s\nLast model took %d seconds',...
          tt,length(tasks),doWhat,person_i,Ns,personName,processDuration);
        if usefig;
          figure(99); cla; text(.01,.5,string,'fontsize',36,'color',[0 0 1],'interpreter','none'); g_hdle = gca; g_hdle.Visible='off'; drawnow
        end
      end
      
      % FFX folder
      FFXdirName = [mainFolder filesep FFXfolderNames{1} filesep personName filesep];
      if ~exist( FFXdirName, 'dir' ); mkdir(FFXdirName); fprintf(fid, 'Erstelle Ordner %s/r/n', FFXdirName ); end
      fprintf('FFX dir name: %s\n',FFXdirName)
      
      % ===================== what to do? =======================
      switch doWhat
        case 'checkRealignmentParameters'
          figure;
          for ses = 1:Nses
            runNameFull = [mainFolder '/fMRIdataComplete' filesep personName filesep runNames{ses,1}];
            multiregFile = dir( [runNameFull '/rp*.txt' ] );
            if length(multiregFile)>0
              sessMultiregFile = [runNameFull filesep multiregFile.name];
            else
              disp('No realignment parameters file found!')
            end
            if Nses > 1
              subplot(round(sqrt(Nses)),round(sqrt(Nses)),ses)
            end
            rp = dlmread(sessMultiregFile);
            plot(rp)
            axis tight
            set(gca,'ylim',[-4 4])
            title([runNames{ses,1}])
          end
          if Nses>1
            suptitle([personName])
          else
            title([personName])
          end
          if isempty(find(abs(rp)>3))
            printfig(['realignmentParameters_' personName])
          else
            printfig(['realignmentParameters_' personName '_over3'])
          end
          close
          
        case 'calculateModelFit'
          error('deactivated')
          % multiply residual model mask and anatomical mask images
          % with each other and average all voxels
          clear matlabbatch
          matlabbatch{1}.spm.util.imcalc.input = {maskImg; ResMS; V1img};
          matlabbatch{1}.spm.util.imcalc.output = 'temp';
          matlabbatch{1}.spm.util.imcalc.outdir = {codeDir};
          matlabbatch{1}.spm.util.imcalc.expression = 'i1.*i2.*i3';
          matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
          matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
          matlabbatch{1}.spm.util.imcalc.options.mask = 0;
          matlabbatch{1}.spm.util.imcalc.options.interp = 1;
          matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
          spm_jobman( 'run', matlabbatch );
          img = spm_read_vols(spm_vol([codeDir 'temp.nii']));
          modelFit(ses,ISI) = mean(img(find(img)));
          
        case 'FFXspecAndEst'
          % ------------ for each run type: ------------------
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
            [onsets,durations,trial_type_names,param_mods,missingOnsets(person_i,ses,:)] = Rutledge_trialtypes_noCompReg_31july2018(fname,ISI);
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
          
          % ================== CONTRASTS ===========================
        case 'contrasts'
          % set the contrasts for each subject depending on the conditions available
          % first, determine for each contrast if the necessary conditions
          % have onsets
          for c = 1:size(cons,1)
            clear matlabbatch;
            matlabbatch{1}.spm.stats.con.spmmat = { [FFXdirName 'SPM.mat'] };
            matlabbatch{1}.spm.stats.con.delete = deleteCons;
            con1run = cons{c,2}; % get contrast values

            % find missing onsets
            initialISIandSetupValues = dlmread(initialISIandSetupValuesFile);
            idx = find(initialISIandSetupValues(:,1)==personNo);
            ISIs = initialISIandSetupValues(idx,2:end);
            clear missingOnsetsForCon
            for ses = 1:Nses
                behavFile = dir([mainFolder '/fMRIdataComplete' filesep personName filesep '*Decision_Rutledge*' runNames{ses,2} '*.mat']);
                fname = [behavFile.folder filesep behavFile.name];
                [onsets,durations,trial_type_names,values,missingOnsetsForCon(ses,1:2:(NtrialTypes+NparamMod))] = Rutledge_trialtypes_noCompReg_26july(fname,ISIs(ses));

                % if no clean initialISI, set missingOnsets for all trialtypes to 1
                if isnan(ISIs(ses))
                    missingOnsetsForCon(ses,:) = 1;
                end
            end
            
            % remove sessions without these trials
            runMultiplicator = ones(1,Nses);
            for ses = 1:Nses
                if sum(missingOnsetsForCon(ses,find(con1run(1:NtrialTypes+NparamMod))))
                    runMultiplicator(ses) = 0;
                end
            end
            runMultiplicator = runMultiplicator(find(~isnan(ISIs)));

            % build contrast vector
            conVec = [kron(runMultiplicator,con1run) zeros(1,length(runMultiplicator))];
            % remove missing conditions
            mm = [makerow([missingOnsetsForCon zeros(Nses,6)]') zeros(1,Nses)];
            conVec = conVec(~mm);
            if sum(abs(conVec)) % at least one non-zero value
                disp(['Running contrast: ' cons{c,1}])
                if c > 1;  matlabbatch{1}.spm.stats.con.delete = 0;  end
                switch conType
                    case 'F'
                        matlabbatch{1}.spm.stats.con.consess{1}.fcon.name = cons{c,1};
                        matlabbatch{1}.spm.stats.con.consess{1}.fcon.convec = conVec;
                        matlabbatch{1}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
                    case 'T'
                        matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = cons{c,1};
                        matlabbatch{1}.spm.stats.con.consess{1}.tcon.convec = conVec;
                        matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
                end
                if printCons
                    matlabbatch{2}.spm.stats.results.spmmat = { [FFXdirName 'SPM.mat'] };
                    matlabbatch{2}.spm.stats.results.conspec.titlestr = '';
                    matlabbatch{2}.spm.stats.results.conspec.contrasts = conNo(c);
                    matlabbatch{2}.spm.stats.results.conspec.threshdesc = 'none';
                    matlabbatch{2}.spm.stats.results.conspec.thresh = 0.001;
                    matlabbatch{2}.spm.stats.results.conspec.extent = 5;
                    matlabbatch{2}.spm.stats.results.conspec.mask = struct('contrasts', {}, 'thresh', {}, 'mtype', {});
                    matlabbatch{2}.spm.stats.results.units = 1;
                    matlabbatch{2}.spm.stats.results.print = true;
                end
                if runBatch
                    spm_jobman( 'run', matlabbatch );
                end
                if printCons
                    psfile = dir([FFXdirName filesep '*.' outputFileType]);
                    seps = findstr(FFXdirName,filesep);
                    newDir = FFXdirName;
                    newDir(seps([5 6])) = '_';
                    copyfile([FFXdirName filesep psfile(1).name],[newDir psfileNewNameRoot cons{c,1} '.' outputFileType])
                    delete(psfile(1).name)
                end
            else
                disp(['Contrast not estimated because of missing trials: ' cons{c,1}])
            end % not empty contrast
          end % each con
      end % doWhat
      % post batch cleanup
      switch doWhat
        case 'contrasts'
          % rename ps output file
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      t2 = clock;
      processDuration = floor(etime( t2,t1));
      if runBatch
        fprintf(fid, [datetime_js 'Dauer fuer Proband ' personName ' in Sekunden: %d\r\n'], processDuration );
        fprintf(fid, '\r\n' );
      end
      report{personNo} = [personName ' ran fine'];
      
    catch errorMsg
      crashedList = [crashedList personNo];
      reportedErrorMessages{personNo} = errorMsg;
      report{personNo} = [personName ' crashed: ' errorMsg.message];
      disp(errorMsg.stack(end))
      %                 keyboard
    end % try/catch
    % post batch cleanup
    switch doWhat
      case 'calculateModelFit'
        error('deactivated')
        figure('name',['ModelFit ' personName ', ' runNames{ses}]);
        plot(initialISI,modelFit(ses,initialISI),'o-')
        xlabel('ISIs [s]')
        ylabel('mean residuals in Area V1')
        [~,idx] = min(modelFit(ses,initialISI));
        bestFittingISI(personNo,ses) = initialISI(idx);
        owd = pwd; cd(codeDir)
        printfig
        close
        cd(owd)
    end
  end % fuer alle Probanden
  disp('....----====|||||====----....')
  if ~isempty(crashedList)
    disp(['Crashed: ' num2str(crashedList)])
  else
    disp('All jobs ran without crash')
  end
  switch doWhat
    case 'calculateModelFit'
      disp('Best fitting ISI in seconds for each subject and run:')
      dlmwrite('initialISIandSetupValues.txt',bestFittingISI)
  end
end % tasks

if runBatch
  fprintf(fid, 'Fertig!' );
  fclose( fid );
end
