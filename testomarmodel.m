% script to load and fit behavioural data from the modified Rutledge et al 2016 task

% clear all; close all
owd = pwd;
useTauForFMRIreg = 0;
showIndivmodelNameFitFigs = 1;

% definitions
% conditionNames = {'soc act','soc pass','non soc'};
% conditionNames = {'Own decision (social)','Other''s decision (social)','Own decision (non-social)'};
% binEdges = {-45:10:45,-100:20:100,-250:25:0,-45:10:45};


% get files
[files, pathname] = uigetfile('*Rutledge*.mat','Bitte Dateien auswaehlen:','MultiSelect', 'on');
if ~iscell(files); ff{1} = files; files = ff; clear ff; end

% sort or not?
if length(files) > 1
    %   sortOrNot = {'sort','don''t sort'};
    %   try selection = bttnChoiceDialog(sortOrNot, 'Sort according to subject number or not?', '','',[1 2]);
    %     sortOrNot = sortOrNot{selection};
    %   end
    sortOrNot = 'sort';
    switch sortOrNot
        case 'sort'
            files = sortrows(char(files));
        case 'don''t sort'
            files = char(files);
    end
else
    files = char(files);
end

% Zscore the happiness results separately for each run
zscoreSepEachRun = {'Yes','No'};
try selection = bttnChoiceDialog(zscoreSepEachRun,'',1,'Zscore happiness separately in each run?',[1 2]);
    zscoreSepEachRun = zscoreSepEachRun{selection};
end
zscoreSepEachRun = 'Yes';

% % Create fMRI regressors?
createfMRIreg = {'Yes','No'};
try selection = bttnChoiceDialog(createfMRIreg,'',1,'Create fMRI regressors?',[1 2]);
    createfMRIreg = createfMRIreg{selection};
end
if strcmpi(createfMRIreg,'Yes')
    toDoList = [1 2];
    spm12 % put spm functions in path
    [file, pathnameISI] = uigetfile('*ISI*.txt','Bitte Datei auswaehlen:');
    initialISIandSetupValues = dlmread('initialISIandSetupValues.txt');
    initialISIandSetupValues = dlmread([pathnameISI filesep file]);
else toDoList = 1;
end

% Show subject decisions?
nonComputAnalysis = {'Yes','No'};
try selection = bttnChoiceDialog(nonComputAnalysis,'',1,'Show subject decisions (non-computational)?',[1 2]);
    nonComputAnalysis = nonComputAnalysis{selection};
end

% =========================================================================

% load
for f = 1:size(files,1)
    Raw = load([pathname filesep deblank(files(f,:))]);
    Raw.results(1).subjName = deblank(files(f,:));
    ii = findstr(Raw.results(1).subjName,'_');
    temp = Raw.results(1).subjName(ii(end-1)+1:ii(end)-1);
    if strcmpi(temp(1),'M') || strcmpi(temp(1),'W') % then it's from Julia's experiment
        temp2 = str2num(temp(2:end));
        if strcmpi(temp(1),'M')
            temp2 = 100 + temp2;
        elseif strcmpi(temp(1),'W')
            temp2 = 200 + temp2;
        end
        temp = temp2;
    end
    if isempty(temp) % || isempty(str2num(temp))
        temp = Raw.results(1).subjName(ii(end-2)+1:ii(end-1)-1);
    end
    if ischar(temp); temp=str2num(temp); end
    %   if isempty(str2num(temp)) % it's not just a number
    %     temp2 = str2num(temp(2:end));
    %     if strcmpi(temp(1),'M')
    %       temp2 = 100 + temp2;
    %     elseif strcmpi(temp(1),'W')
    %       temp2 = 200 + temp2;
    %     end
    %     temp = temp2;
    %   else
    %     temp = str2num(temp);
    %   end
    for str_loop= 1:length(Raw.results)
        gl=[Raw.results.amountThisTrialSubject];
        Raw.results(str_loop).Agl_Abs=sum(gl(1:str_loop-1));
        riskchoices=Raw.results(str_loop).riskyOption;
        riskEV=(riskchoices*0.5);
        choicesEV=[riskEV Raw.results(str_loop).safeOption];
        minEV=min(choicesEV);
        Raw.results(str_loop).MinEV=minEV;
        if Raw.results(str_loop).play==1
            choiceEV=sum(riskchoices*0.5);
            Raw.results(str_loop).ChoiceEV=choiceEV;
            Raw.results(str_loop).Ch_diff=choiceEV-minEV;
            
        elseif Raw.results(str_loop).play==0
            choiceEV=Raw.results(str_loop).safeOption;
            Raw.results(str_loop).ChoiceEV=choiceEV;
            Raw.results(str_loop).Ch_diff=choiceEV-minEV;
        end
        Raw.results(str_loop).PE=(gl(str_loop)-choiceEV);
    end
    subjNumbers(f) = temp;
    R{f} = Raw.results;
    if strcmpi(zscoreSepEachRun,'Yes')
        % zscore happiness (separately for each run)
        happ = [R{f}.happiness];
        idx = find(~isnan(happ));
        happ = zscore(happ(idx));
        for i = 1:length(idx)
            R{f}(idx(i)).happiness = happ(i);
        end
    end
end

% combine results structures?
clear initialISIandSetup
if size(files,1) > 1
    combineOrNot = {'combine','don''t combine'};
    try selection = bttnChoiceDialog(combineOrNot, 'Combine data or not?', '','',[1 2]);
        combineOrNot = combineOrNot{selection};
    end
    combineOrNot = 'combine';
    switch combineOrNot
        case 'combine' % then decide how many files to group
            prompt={'How many files to combine'}; name='Combine';
            numlines=1; defaultanswer={'2'}; %change to 3 if comparing behavioral data
            N = inputdlg(prompt,name,numlines,defaultanswer); N=str2num(N{1});
            for g = 1:size(files,1) / N % number of groups
                idx = [1:N] + (g-1)*N;
                RR{g} = horzcat(R{idx});
                subjNumbTemp(g) = subjNumbers(idx(1));
                % if need to create fMRI regressors, load the initial ISI and setup time
                %                 if strcmpi(createfMRIreg,'Yes')
                %                     idx = find(initialISIandSetupValues(:,1,1) == subjNumbTemp(g));
                %                     if isempty(idx)
                %                         warning(sprintf('no ISI values found for subject %d',subjNumbTemp(g)));
                %                         initialISIandSetup(g,:) = NaN;
                %                     else
                %                         initialISIandSetup(g,:) = initialISIandSetupValues(idx,2:end); % go to the correct run and pick its number
                %                     end
                %                 end
            end
            R = RR;
            subjNumbers = subjNumbTemp;
        case 'don''t combine' % nothing to do
    end
end

% Use which modelName?
modelName = {'GuiltEnvy','ActPass','InequalSingle','NoInequal','Compare','MinEV+PE non-social'};
try selection = bttnChoiceDialog(modelName,'',1,'Use which modelName for computational analysis?',[2 2]);
    modelName = modelName{selection};
end

% % which trials to fit?
% whichTrials = {'all','social','non-social','social-active','social-passive','social-passive and non-social','social-active and non-social'};
% defSelection = whichTrials{1};
% try selection = bttnChoiceDialog(whichTrials,'',1,'Which trials?',[7 1]);
%   whichTrials = whichTrials{selection};
% end
whichTrials = 'social-active and non-social';

% % For the non-computational analysis, use only the trials directly before the happiness ratings?
% useAllTrialsForNonComp = {'No','Yes'};
% try selection = bttnChoiceDialog(useAllTrialsForNonComp,'',1,'Use all trials for non-computational analysis?',[1 2]);
%   useAllTrialsForNonComp = useAllTrialsForNonComp{selection};
% end
useAllTrialsForNonComp = 'no';

% Do the fitting and the behavioural data analysis
clear happiness* hap_idx* decisions
for i = 1:length(R)
    happinessTemp = [R{i}(:).happiness];
    if strcmpi(zscoreSepEachRun,'No')
        idx = find(~isnan(happinessTemp));
        happ = zscore(happinessTemp(idx));
        happinessTemp(idx) = happ;
    end
    
    % analysis of the behavioural data
    % get the happiness values for the different outcomes
    hap_typeNames_soc = {'win-win','win-lose','lose-lose','lose-win'}; % first subject, then partner
    hap_typeNames_nonsoc = {'win','lose'};
    hap_whoWon_soc = [1 1;1 0;0 0;0 1]; % first subject, then partner
    hap_whoWon_nonsoc = [1;0]; % only subject
    for ii = 1:length(hap_typeNames_soc)
        for cc = 1:2 % social active, social passive and non-social conditions
            hap_idx{ii,cc} = find([R{i}(:).subjectWon]==hap_whoWon_soc(ii,1) & [R{i}(:).partnerWon]==hap_whoWon_soc(ii,2) & [R{i}(:).condition]==cc);
            if strcmpi(useAllTrialsForNonComp,'yes')
                try happiness_soc(i,ii,cc) = nanmean([R{i}(hap_idx{ii,cc}).happiness R{i}(hap_idx{ii,cc}-1).happiness]); catch happiness_soc(i,ii,cc) = NaN; end
            else
                try happiness_soc(i,ii,cc) = nanmean([R{i}(hap_idx{ii,cc}).happiness]); catch happiness(i,ii,cc) = NaN; end
            end
        end
    end
    for ii = 1:length(hap_typeNames_nonsoc)
        hap_idx{ii,1} = find([R{i}(:).subjectWon]==hap_whoWon_nonsoc(ii) & [R{i}(:).condition]==3);
        if strcmpi(useAllTrialsForNonComp,'yes')
            try happiness_nonsoc(i,ii,1) = nanmean([R{i}(hap_idx{ii,1}).happiness R{i}(hap_idx{ii,1}-1).happiness]); catch happiness_nonsoc(i,ii,1) = NaN; end
        else
            try happiness_nonsoc(i,ii,1) = nanmean([R{i}(hap_idx{ii,1}).happiness]); catch happiness_nonsoc(i,ii,1) = NaN; end
        end
    end
    
    %     if strcmpi(nonComputAnalysis,'Yes')
    %         % look at the decisions that people made
    %         conditions = [R{i}(:).condition];
    %         %   Nbins = 7;
    %         colors = 'bgr';
    %         figure('name',sprintf('Subject %d',subjNumbers(i)));
    %         for c = 1:length(conditionNames)
    %             idx2 = find(conditions == c);
    %             riskyMean = mean(reshape([R{i}(idx2).riskyOption],2,[]));
    %             riskyDiff = diff(reshape([R{i}(idx2).riskyOption],2,[]));
    %             safe = [R{i}(idx2).safeOption];
    %             play = [R{i}(idx2).play];
    %             %     varNames = {'safe','riskyMean','riskyDiff','riskyMean-safe'};
    %             %     vars = {safe,riskyMean,riskyDiff,riskyMean-safe};
    %             varNames = {'EV_r_i_s_k_y-safe'};
    %             vars = {riskyMean-safe};
    %             for v = 1:length(vars)
    %                 [qq,ww,ee] = histcounts(vars{v},binEdges{v});
    %                 for iii = 1:length(binEdges{v})-1
    %                     decisions{v}(i,iii,c) = nanmean(play(find(ee==iii)));
    %                 end
    %                 ww = (ww(1:end-1)+ww(2:end))/2;
    %                 if length(vars)>1; subplot(2,2,v); end
    %                 hold on
    %                 plot(ww,squeeze(decisions{v}(i,:,c)),[colors(c) '.-'])
    %                 title(varNames{v})
    %                 %       ylabel('P(risky), \bar{x}\pm SEM','interpreter','tex')
    %                 ylabel('P(risky) [M+-SEM]')
    %             end
    %         end
    %         legend(conditionNames)
    %     end
    
    % now prepare the variables for the computational analysis
    choseSafe = [R{i}(:).play];
    %   choseGamble = 1-choseSafe;
    try, nonsocialTrial = ([R{i}(:).condition] == 3);
    catch; nonsocialTrial = ~isnan([R{i}(:).amountThisTrialPartner]); end
    try, socialActive = ([R{i}(:).condition] == 1);
    catch; socialActive = [R{i}.condition]; end
    try, socialPassive = ([R{i}(:).condition] == 2); end
    amountReceivedBySubject = [R{i}.amountThisTrialSubject];
%     CR = amountReceivedBySubject.*choseSafe;
%     EV = mean(reshape([R{i}(:).riskyOption],2,[])) .* (1-choseSafe);
%     RPE = EV-amountReceivedBySubject.*(1-choseSafe);
    
    matrix = zeros(length(R{i}),10);
    %     matrix(:,1) = [R{i}(:).condition];
         
   
   
    %     matrix(:,8) = [R{i}(:).amountThisTrialSubject];
    matrix(:,1) = [R{i}(:).Agl_Abs];
    matrix(:,2) = [R{i}(:).MinEV];
    matrix(:,3) = [R{i}(:).ChoiceEV];
    matrix(:,4) = [R{i}(:).Ch_diff];
    matrix(:,5) = [R{i}(:).PE];
    matrix(:,6) = happinessTemp;
    matrix(:,7) = [R{i}(:).safeOption];
    matrix(:,8:9) = reshape([R{i}(:).riskyOption],2,[])';
    matrix(:,10) = [R{i}(:).play];
    %     matrix(:,13) = isnan([R{i}(:).amountThisTrialPartner])+1;
    %     matrix(:,16) = [R{i}(:).amountThisTrialPartner];
    
    switch whichTrials
        case 'all'
            idx = 1:size(matrix,1);
        case 'social'
            idx = find(~nonsocialTrial);
        case 'non-social'
            idx = find(nonsocialTrial);
        case 'social-active'
            idx = find(socialActive);
        case 'social-passive'
            idx = find(socialPassive);
        case 'social-passive and non-social'
            idx = [find(socialPassive) find(nonsocialTrial)];
        case 'social-active and non-social'
            idx = [find(socialActive) find(nonsocialTrial)];
    end
    matrix = matrix(idx,:);
    
    alldata(i).socdata{1} = matrix;
end

% % show happiness ratings
% figure;
% happiness = happiness_soc;
% happiness(:,[5 6],3) = happiness_nonsoc;
% if i>1 % more than 1 subject
%   meansembar(happiness); drawnow
% else
%   bar(squeeze(happiness)'); drawnow
% end
% xlabel_oblique(conditionNames);
% legend([hap_typeNames_soc hap_typeNames_nonsoc])
% ylabel('Happiness (z-scored)')
%
% %non-computational
% if strcmpi(nonComputAnalysis,'Yes') && length(R)>1
%   % show decisions to play as function of diff variables and conditions
%   figure('name','Decision to play');
%   for v = 1:length(vars)
%     for c = 1:length(conditionNames)
%       ww = binEdges{v};
%       ww = (ww(1:end-1)+ww(2:end))/2;
%       if length(vars)>1; subplot(2,2,v); hold on; end
%       meansemplot2(ww,decisions{v})
%       set(gca,'ylim',[0 1])
%       ylabel('P(risky) [M&SEM]')
%       xlabel(varNames{v})
%       box on
%     end
%     legend(conditionNames,'location','SouthEast')
%   end
% end

clear result results
% alldata(:).socdata{1,1}=alldata(:).socdata{1,1}(:,1:6);
for whatToDo = toDoList % 1 = fit the data, 2 = create fMRI regressors.
    % first, fit behavioural data, then create fMRI regressors if desired
    
    for n=1:length(alldata), % n subjects
        fprintf(1,'%d,',n);
        for m=1:length(alldata(n).socdata), % m runs?
            % temp is the data, rows are trials, columns are variables:
            % 3: Agl_Abs
            % 4: MinEV
            % 5: ChoiceEV
            % 6: Ch_diff
            % 7: PE
            % 10: happiness score
            % 13: 1 if participant chose, 2 if other person chose
            % 16: other person reward
%             alldata(n).socdata{1,1}=alldata(n).socdata{1,1}(:,1:6);
            eval(sprintf('temp = alldata(%d).socdata{%d};',n,m)); 
%             temp(1,10:12)=nan; %collect matrices and toss first rating
            t2 = temp(:,6); rawhappy = t2(~isnan(t2)); happyind = find(~isnan(t2));
%             ; temp(isnan(temp(:,16)),16)=0; %other reward=0 when solo trial
            nrate=length(rawhappy); ntrial=size(temp,1); x=zeros(nrate,ntrial);
%             evmtx=x;  rpemtx=x;  othercertainmtx=x; otherrewardmtx=x; ipemtx=x; envymtx=x; guiltmtx=x; notenvymtx=x; notguiltmtx=x; youchosemtx=x; theychosemtx=x;
            Agl_Absmtx=x; MinEVmtx=x; ChoiceEVmtx=x; Ch_diffmtx=x; PEmtx=x; rewardmtx=x; certainmtx=x;
            
            switch whatToDo
                case 1 % fit the data
                    trialsToFit = happyind;
                case 2 % create fMRI regressors
                    trialsToFit = 1:size(matrix,1);
            end
            
            for q = 1:length(trialsToFit)
                temp2            = temp(1:trialsToFit(q),:); %clip out all trials up to rating
%                tempev           = mean(temp2(:,4:5),2) .* double(temp2(:,7)==1); %0 if no gamble or error, ev if gambled
                %                 temprpe          = temp2(:,8) .* double(temp2(:,7)==1) - tempev; %0 if no gamble or error, rpe if gambled
                tempreward       = temp2(:,3) .* double(temp2(:,10)==1); %gamble rewards (0 if chose certain)
                tempcertain      = temp2(:,7) .* double(temp2(:,10)==0); %0 if gambled, otherwise certain amount
                %                 tempothercertain = temp2(:,16) .* double(temp2(:,7)==0) .* double(temp2(:,13)==2); %other rewards
                %                 tempotherreward  = temp2(:,16) .* double(temp2(:,7)==1) .* double(temp2(:,13)==2);
                %                 tempipe          = (sign(temprpe))   .* double(temp2(:,1)<3) .* (tempotherreward-tempev); %sign or your rpe x sign of theirs
                %                 tempenvy         = (sign(temprpe)<0) .* double(temp2(:,1)<3) .* (tempotherreward-tempreward); %you lost and they got this much more than you
                %                 tempguilt        = (sign(temprpe)>0) .* double(temp2(:,1)<3) .* (tempreward-tempotherreward); %you won and you got this much more than them
                %                 tempnotenvy      = (sign(temprpe)<0) .* double(temp2(:,1)<3) .* double(~(tempotherreward-tempreward)); %you lost and they also lost dummy
                %                 tempnotguilt     = (sign(temprpe)>0) .* double(temp2(:,1)<3) .* double(~(tempreward-tempotherreward)); %you won and they also won dummy
                %                 tempyouchose     = double(temp2(:,13)==1);
                %                 temptheychose    = double(temp2(:,13)==2);
                %                 evmtx(q,1:length(tempev))                     = fliplr(transpose(tempev)) / 100; %convert to ?
                                certainmtx(q,1:length(tempcertain))           = fliplr(transpose(tempcertain)) / 100;
                %                 rpemtx(q,1:length(temprpe))                   = fliplr(transpose(temprpe)) / 100;
                                 rewardmtx(q,1:length(tempreward))             = fliplr(transpose(tempreward)) / 100;
                %                 ipemtx(q,1:length(tempipe))                   = fliplr(transpose(tempipe)) / 100;
                %                 othercertainmtx(q,1:length(tempothercertain)) = fliplr(transpose(tempothercertain)) / 100;
                %                 otherrewardmtx(q,1:length(tempotherreward))   = fliplr(transpose(tempotherreward)) / 100;
                %                 envymtx(q,1:length(tempenvy))                 = fliplr(transpose(tempenvy)) / 100;
                %                 guiltmtx(q,1:length(tempguilt))               = fliplr(transpose(tempguilt)) / 100;
                %                 notenvymtx(q,1:length(tempnotenvy))           = fliplr(transpose(tempnotenvy)); %dummy
                %                 notguiltmtx(q,1:length(tempnotguilt))         = fliplr(transpose(tempnotguilt)); %dummy
                %                 youchosemtx(q,1:length(tempyouchose))         = fliplr(transpose(tempyouchose)); %1 if you chose
                %                 theychosemtx(q,1:length(temptheychose))       = fliplr(transpose(temptheychose)); %1 if they chose
                Agl_Abs=double(temp2(:,1));
                MinEV=double(temp2(:,2));
                ChoiceEV=double(temp2(:,3));
                Ch_diff=double(temp2(:,4));
                PE=double(temp2(:,5));
                Agl_Absmtx(q,1:length(Agl_Abs))                     = fliplr(transpose(Agl_Abs)) / 100; %convert to ?
                MinEVmtx(q,1:length(MinEV))                    = fliplr(transpose(MinEV)) / 100; %convert to ?
                ChoiceEVmtx(q,1:length(ChoiceEV))                     = fliplr(transpose(ChoiceEV)) / 100; %convert to ?
                Ch_diffmtx(q,1:length(Ch_diff))                     = fliplr(transpose(Ch_diff)) / 100; %convert to ?
                PEmtx(q,1:length(PE))                     = fliplr(transpose(PE)) / 100; %convert to ?
            end;
%             clear mtx;  mtx.evmtx=single(evmtx); mtx.rpemtx=single(rpemtx); 
%             mtx.ipemtx=single(ipemtx); mtx.othercertainmtx=single(othercertainmtx); mtx.otherrewardmtx=single(otherrewardmtx);
%             mtx.envymtx=single(envymtx); mtx.guiltmtx=single(guiltmtx); mtx.youchosemtx=single(youchosemtx); mtx.theychosemtx=single(theychosemtx);
%             mtx.notenvymtx=single(notenvymtx); mtx.notguiltmtx=single(notguiltmtx);
            mtx.conds=temp(1:trialsToFit(end),1);
            mtx.Agl_Abs=single(Agl_Absmtx); mtx.MinEV=single(MinEVmtx); 
            mtx.ChoiceEV=single(ChoiceEVmtx); mtx.Ch_diff=single(Ch_diffmtx); 
            mtx.PE=single(PEmtx); mtx.rewardmtx=single(rewardmtx); mtx.certainmtx=single(certainmtx);
            eval(sprintf('alldata(%d).mtx{%d}=mtx;',n,m));
            eval(sprintf('alldata(%d).rawhappy{%d}=rawhappy;',n,m));
            eval(sprintf('temp = alldata(%d).socdata{%d};',n,m)); 
%             temp(1,10:12)=nan; %toss first rating
            t2 = temp(:,6); zhappy = t2(~isnan(t2));
            eval(sprintf('alldata(%d).zhappy{%d}=zhappy;',n,m));
            
            % do the fitting
            if whatToDo == 1 %2 & strcmpi(createfMRIreg,'Yes')
                if strcmpi(modelName,'ActPass')
                    results(n,1) = fit_happy_model_nogain_inequalityActPass(mtx,zhappy,NaN);
                elseif strcmpi(modelName,'GuiltEnvy')
                    results(n,1) = fit_happy_model_nogain_guiltenvy(mtx,zhappy,NaN);
                elseif strcmpi(modelName,'InequalSingle')
                    results(n,1) = fit_happy_model_nogain_inequalitySingle(mtx,zhappy,NaN);
                elseif strcmpi(modelName,'NoInequal')
                    results(n,1) = fit_happy_model_nogain_noInequal(mtx,zhappy,NaN);
                elseif strcmpi(modelName,'MinEV+PE non-social')
                    nonsoc_idx=find([Raw.results.condition]==3);
                    results(n,1) = fit_happy_model_Omar_1(mtx,zhappy,NaN);
                elseif strcmpi(modelName,'Compare')
                    results(n,1) = fit_happy_model_nogain_inequalityActPass(mtx,zhappy,NaN);
                    results(n,2) = fit_happy_model_nogain_guiltenvy(mtx,zhappy,NaN);
                    results(n,3) = fit_happy_model_nogain_inequalitySingle(mtx,zhappy,NaN);
                    results(n,4) = fit_happy_model_nogain_noInequal(mtx,zhappy,NaN);
                    
                end
                if showIndivmodelNameFitFigs
                    for m = 1:size(results,2)
                        figure;
                        subplot(2,1,1)
                        plot(zhappy,'b');
                        hold on; plot(results(n,m).happypred,'r')
                        dataNames = {'happiness','predicted happiness'};
                        legend(dataNames)
                        title([R{n}(1).subjName ' ' whichTrials], 'interpreter','none')
                        ylabel(sprintf('R2 = %.2f',results(n,m).r2));
                        subplot(2,1,2)
                        bar(results(n,m).b)
                        xlabel_oblique(results(n,m).paramNames)
                        title([results(n,m).modelName ' - parameters'], 'interpreter','none')
                        drawnow
                    end
                end
                %     regress_display(zhappy,results(n).happypred,'inputNames',dataNames)
            end
            
            % use tau to create fMRI regressors if desired
            if whatToDo == 2 & strcmpi(createfMRIreg,'Yes')
                idx = strmatch('decay',results(n,1).paramNames); % look for the decay parameter
                tau = results(n,1).b(idx); %decay constant
                if useTauForFMRIreg
                    [fMRIregressors,fMRIregressor_names] = create_fMRIregressors(tau,mtx,R{n},initialISIandSetup(n,:)); % use results of the fit to create fMRI regressors
                else
                    [fMRIregressors,fMRIregressor_names] = create_fMRIregressors(0,mtx,R{n},initialISIandSetup(n,:)); % use results of the fit to create fMRI regressors
                end
                q=R{n}(1).subjName; us = findstr(q,'_');
                subjNo = q(us(5)+1:us(6)-1);
                if double(subjNo(1)) > 57
                    subjNo = q(us(4)+1:us(5)-1);
                end
                for ses = 1:length(fMRIregressors)
                    dlmwrite(['fMRIregressors_' subjNo '_' num2str(ses) '.txt'],fMRIregressors{ses})
                    figure; imagesc(fMRIregressors{ses});
                    set(gca,'xtick',[1:size(fMRIregressors{ses},2)]+.4)
                    xlabel_oblique(fMRIregressor_names(1:size(fMRIregressors{ses},2)),15,12,1);
                    set(gca,'position',[0.1300    0.1800    0.7750    0.7450])
                    ylabel('Trial from start of run1')
                    title(['fMRIregressors for subject ' subjNo ', run ' num2str(ses)])
                end % ses
            end
        end
    end;
end;
% end

disp(reshape([results.r2],n,[]));

% combine the results in a matrix
clear resMtx
for m = 1:size(results,2)
    for n = 1:size(results,1)
        resMtx{m}(n,:) = [subjNumbers(n) results(n,m).r2 results(n,m).aic results(n,m).bic results(n,m).b]; % with aic and bic; but what do they represent?
        %   resMtx(n,:) = [subjNumbers(n) results(n).r2 results(n).b];
    end
    resMtxVarnames = [{'Subject number', 'R2', 'aic', 'bic'}, results(n,m).paramNames];
    % resMtxVarnames = [{'Subject number', 'R2'}, results(n).paramNames];
    
    % display
    if length(alldata) > 1
        % group results
        figure('name',['Parameters for modelName ' results(1,m).modelName]);
        showWhat = [2 5:size(resMtx{m},2)];
        plot([.5 length(showWhat)+.5],[0 0],'k')
        try notBoxPlot(resMtx{m}(:,showWhat))
        catch meansemplot(resMtx{m}(:,showWhat))
        end
        try; xlabel_oblique(resMtxVarnames(showWhat)); end
        
        % individual variation
        figure('name',['Individual variation of parameters for modelName ' results(1,m).modelName]);
        subplot(4,1,1)
        bar(resMtx{m}(:,2)); title('R2')
        subplot(4,1,2)
%         bar(nanmean(diff(happiness(:,2:3,1:2),[],3)')); title('Effect of own vs other''s decision on happiness in other loss trials')
%         subplot(4,1,3)
%         bar(resMtx{m}(:,end-1)); title(results(n,m).paramNames(end-1))
%         subplot(4,1,4)
%         bar(resMtx{m}(:,end)); title(results(n,m).paramNames(end))
        
        % write out data to analyse in SPSS or such
        dlmwrite(['RutledgeComputFitData_' results(n,m).modelName '.csv'],resMtx{m})
    end
    disp('Mean and median R2:')
    disp([mean(resMtx{m}(:,2)) median(resMtx{m}(:,2))])
end
