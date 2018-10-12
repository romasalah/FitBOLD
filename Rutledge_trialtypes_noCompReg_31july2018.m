function [ons,dur,tr_names,values,missing_onsets] = Rutledge_trialtypes_noCompReg_31july(fname,startISI)
% Analysis script for SoDec_fMRI
% rearranges the onsets of the important events of the fMRI experiment into
% a cell array according to their distinctive trial types
% Number of Trial types: 25

figFlag = 0;

if ~exist('fname','var')
  fname = '/Users/NEMO_admin/Google Drive/Master(shared)/Code/data/RutledgeBehav/socialDecision_Rutledge_behav_SoDec-behav_80_06-Jul-2017-10h16.mat';
end

load(fname);
for str_loop= 1:length(results)
        gl=[results.amountThisTrialSubject];
        results(str_loop).Agl_Abs=sum(gl(1:str_loop-1));
        riskchoices=results(str_loop).riskyOption;
        riskEV=(riskchoices*0.5);
        choicesEV=[riskEV results(str_loop).safeOption];
        minEV=min(choicesEV);
        results(str_loop).MinEV=minEV;
        if results(str_loop).play==1
            choiceEV=sum(riskchoices*0.5);
            results(str_loop).ChoiceEV=choiceEV;
            results(str_loop).Ch_diff=choiceEV-minEV;
            
        elseif results(str_loop).play==0
            choiceEV=results(str_loop).safeOption;
            results(str_loop).ChoiceEV=choiceEV;
            results(str_loop).Ch_diff=choiceEV-minEV;
        end
%         results(str_loop).PE=(gl(str_loop)-choiceEV);
end

% define basic variables to find trial types
cond = [results(:).condition]; % 1: social active; 2: social passive; 3: non-social
chRisky = [results(:).play];  chRisky(chRisky==0) = NaN;
chSafe = 1-[results(:).play];  chSafe(chSafe==0) = NaN;
CR = [results(:).safeOption] .* chSafe; % subject is certain to get this as the safe option was chosen
EV = mean(reshape([results(:).riskyOption],2,[])) .* chRisky; % average gain from risky options, which subject chose
PE = EV-[results.amountThisTrialSubject]; % prediction error from real outcome of risky option
outcomeDif = [results(:).amountThisTrialSubject]-[results(:).amountThisTrialPartner]; % the difference in the amount of money obtained by subject adn partner in a given trial
Agl_Abs=[results(:).Agl_Abs];
MinEV=[results(:).MinEV];
Ch_diff=[results(:).Ch_diff];

% set reference time if none available
offset = results(1).optionsShown;
if ~exist('startISI','var')
    startISI = 0;
end

% all kinds of durations and onsets
RT =     [results(:).RT]/1000; RT(isnan(RT)) = 3; % response time, for opponent this was always 3s
onsOpt = [results(:).optionsShown] - offset + startISI; % time when options were shown
onsDec = [[results(:).outcomeShown] + RT] - offset + startISI; % time when decision was taken or shown
onsOut = [results(:).outcomeShown] - offset + startISI; % time when the outcome was shown
onsHap = [results(:).ratingShown] - offset + startISI; % time when the happiness rating was given

% trial types: names, onsets, durations and parametric modulators (values)
% for certain rewards and expected value, onsets are the time at which the participant was presented the certain reward - we could use the decision time here, as only then was it certain that he would get this
% i = 1;  tr_names{i} = 'CR active';              idx = ~isnan(CR) & cond==1;         ons{i} = onsOpt(idx);   dur{i} = RT(idx);  values{i} = CR(idx);
% i = 2;  tr_names{i} = 'CR passive';             idx = ~isnan(CR) & cond==2;         ons{i} = onsOpt(idx);   dur{i} = RT(idx);  values{i} = CR(idx);
% i = 3;  tr_names{i} = 'CR n-soc';               idx = ~isnan(CR) & cond==3;         ons{i} = onsOpt(idx);   dur{i} = RT(idx);  values{i} = CR(idx);
% i = 4;  tr_names{i} = 'EV active';              idx = ~isnan(EV) & cond==1;         ons{i} = onsOpt(idx);   dur{i} = RT(idx);  values{i} = EV(idx);
% i = 5;  tr_names{i} = 'EV passive';             idx = ~isnan(EV) & cond==2;         ons{i} = onsOpt(idx);   dur{i} = RT(idx);  values{i} = EV(idx);
% i = 6;  tr_names{i} = 'EV n-soc';               idx = ~isnan(EV) & cond==3;         ons{i} = onsOpt(idx);   dur{i} = RT(idx);  values{i} = EV(idx);
% for prediction errors, onsets are the time at which the participant was presented with the outcome
i = 1;  tr_names{i} = 'PE active';              idx = ~isnan(PE) & cond==1;         ons{i} = onsOut(idx);   values{i} = PE(idx);
i = 2;  tr_names{i} = 'PE passive';             idx = ~isnan(PE) & cond==2;         ons{i} = onsOut(idx);   values{i} = PE(idx);
i = 3;  tr_names{i} = 'PE n-soc';               idx = ~isnan(PE) & cond==3;         ons{i} = onsOut(idx);   values{i} = PE(idx);
% for outcome differences, onsets are the time at which the participant was presented with the outcome
i = 4; tr_names{i} = 'Outcome self>other active';    idx = outcomeDif>0 & cond==1; ons{i} = onsOut(idx);   values{i} = outcomeDif(idx);
i = 5; tr_names{i} = 'Outcome self>other passive';   idx = outcomeDif>0 & cond==2; ons{i} = onsOut(idx);   values{i} = outcomeDif(idx);
i = 6; tr_names{i} = 'Outcome self<other active';    idx = outcomeDif<0 & cond==1; ons{i} = onsOut(idx);   values{i} = outcomeDif(idx);
i = 7; tr_names{i} = 'Outcome self<other passive';   idx = outcomeDif<0 & cond==2; ons{i} = onsOut(idx);   values{i} = outcomeDif(idx);
% % for decisions, onsets are the time at which the participant pressed the button indicating his decision
% i = 12; tr_names{i} = 'decision play active';   idx = ~isnan(chRisky) & cond==1;    ons{i} = onsOpt(idx);   values{i} = [];
% i = 13; tr_names{i} = 'decision play passive';  idx = ~isnan(chRisky) & cond==2;    ons{i} = onsOpt(idx);   values{i} = [];
% i = 14; tr_names{i} = 'decision play n-soc';    idx = ~isnan(chRisky) & cond==3;    ons{i} = onsOpt(idx);   values{i} = [];
% i = 15; tr_names{i} = 'decision safe active';   idx = ~isnan(chSafe) & cond==1;     ons{i} = onsOpt(idx);   values{i} = [];
% i = 16; tr_names{i} = 'decision safe passive';  idx = ~isnan(chSafe) & cond==2;     ons{i} = onsOpt(idx);   values{i} = [];
% i = 17; tr_names{i} = 'decision safe n-soc';    idx = ~isnan(chSafe) & cond==3;     ons{i} = onsOpt(idx);   values{i} = [];
% for happiness, onsets are the time at which the participant reported happiness - we don't have the onset of the happiness rating display
i = 8; tr_names{i} = 'happiness';                idx = ~isnan([results(:).happiness]);     ons{i} = onsHap(idx);                     values{i} = [results(idx).happiness];
% i = 9; tr_names{i} = 'Agl active';               idx = ~isnan(Agl_Abs) & cond==1;         ons{i} = onsOpt(idx);   dur{i} = RT(idx);  values{i} = Agl_Abs(idx);
% i = 10; tr_names{i} = 'Agl passive';              idx = ~isnan(Agl_Abs) & cond==2;         ons{i} = onsOpt(idx);   dur{i} = RT(idx);  values{i} = Agl_Abs(idx);
% i = 11; tr_names{i} = 'Agl non-social';           idx = ~isnan(Agl_Abs) & cond==3;         ons{i} = onsOpt(idx);   dur{i} = RT(idx);  values{i} = Agl_Abs(idx);
i = 9; tr_names{i} = 'Min EV active';            idx = ~isnan(MinEV) & cond==1;         ons{i} = onsOpt(idx);   dur{i} = RT(idx);  values{i} =MinEV(idx);
i = 10; tr_names{i} = 'Min EV passive';           idx = ~isnan(MinEV) & cond==2;         ons{i} = onsOpt(idx);   dur{i} = RT(idx);  values{i} = MinEV(idx);
i = 11; tr_names{i} = 'Min EV non-social';        idx = ~isnan(MinEV) & cond==3;         ons{i} = onsOpt(idx);   dur{i} = RT(idx);  values{i} = MinEV(idx);
i = 12; tr_names{i} = 'Ch_diff active';           idx = ~isnan(Ch_diff) & cond==1;         ons{i} = onsOpt(idx);   dur{i} = RT(idx);  values{i} = Ch_diff(idx);
i = 13; tr_names{i} = 'Ch_diff passive';          idx = ~isnan(Ch_diff) & cond==2;         ons{i} = onsOpt(idx);   dur{i} = RT(idx);  values{i} =Ch_diff(idx);
i = 14; tr_names{i} = 'Ch_diff non-social';       idx = ~isnan(Ch_diff) & cond==3;         ons{i} = onsOpt(idx);   dur{i} = RT(idx);  values{i} = Ch_diff(idx);
% durations are all 0 for now.
for c = 1:length(tr_names)
  dur{c} = zeros(1,length(ons{c}));
end

% remove NaNs
for t = 1:length(ons)
    ons{t} = ons{t}(~isnan(ons{t}));
    dur{t} = dur{t}(~isnan(dur{t}));
end

% in case no onset available, replace with dummy very large onset
missing_onsets = zeros(1,length(ons));
for t = 1:length(ons)
    if isempty(ons{t})
        missing_onsets(t) = 1;
        ons{t} = max([ons{:}])+1;
        dur{t} = 0;
    end
end
disp(cellfun(@length,ons))

% to check: onsets should increase linearly, and as they are in seconds, if
% divided by 60 the max value should be around 20 as this was roughly the
% duration of the scan
if figFlag; figure; plot(sort([ons{:}]/60),'.'); ylabel('Time in min.'); xlabel('Onsets of all types'); 
  figure
  for c=1:length(ons)
    subplot(length(ons),1,c)
    stem(ons{c},ones(1,length(ons{c})))
  end
end
