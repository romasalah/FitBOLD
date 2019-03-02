function [fit]=design_fitBOLD(fit)
all_BOLD_combined=cell(0,0);
all_BOLD_singles=cell(0,0);
all_behav_combined=cell(0,0);
all_behav_singles=cell(0,0);

for look_duplicate=1:size(fit.timings,2)
    if strcmp(fit.timings{1,look_duplicate},fit.Target_timing{1,1})
        fit.bimodal=0;
    end
end

if size(fit.ROIs,2)==size(fit.timings,2)
    combine_size=size(fit.whatcombine,1);
    
    %************loading target BOLD vector***********
    
    fprintf('\n Getting your target BOLD or target Behavioural vector \n')
    tgt=zeros(0,1);
    fitTarget=fit; fitTarget.timings=fit.Target_timing; fitTarget.ROIs=fit.TargetBOLD;
    if fit.behav==1
       evalc('[~,target_onsets,~,tgt_sess,session_num,sub_tgt_id,~]=corr_BOLD(fitTarget,0)');
        for z_tgt=1:size(tgt_sess,2)
            tgt=vertcat(tgt,tgt_sess{1,z_tgt});
        end
        
        fprintf('\n loading target Behavioural vector was Successful \n')
    else
        evalc('[shifts,target_onsets,BOLD_tgt_sess]=corr_BOLD(fitTarget,0)')
        for z_BOLD_tgt=1:size(BOLD_tgt_sess,2)
            tgt=vertcat(tgt,BOLD_tgt_sess{1,z_BOLD_tgt});
        end
        fprintf('\n loading target BOLD vector was Successful \n')
    end
    %*************getting and combining BOLD from different ROIs**************
    if exist('fit.whatcombine','var') && fit.whatcombine(1,1) ~=0
        ROIs_names=cell(1,size(fit.whatcombine,1));
        timings_for_BOLD=cell(1,size(fit.whatcombine,1));
        for ROI_names_i=size(fit.whatcombine,2):size(fit.whatcombine,2):numel(fit.whatcombine)
            ROIs_names{1,ROI_names_i}=fit.ROIs{ROI_names_i};
            timings_for_BOLD{1,ROI_names_i}=fit.timings{ROI_names_i}(1,:);
        end
        for fill_rest_ROI_names= ROI_names_i+1:size(fit.ROIs,2)
            ROIs_names{1,end+1}=fit.ROIs{fill_rest_ROI_names}; %#ok<*AGROW>
            timings_for_BOLD{1,end+1}=fit.timings{fill_rest_ROI_names};
        end
        non_empt_idx=find(~cellfun(@isempty,ROIs_names));
        ROIs_names={ROIs_names{1,[non_empt_idx]}};
        timings_for_BOLD={timings_for_BOLD{1,[non_empt_idx]}};
        %correct names to bilateral
        if ~isempty(ROIs_names)
            for corct_bi=1:size(ROIs_names,2)
                ROIs_names{1,corct_bi}=[ROIs_names{1,corct_bi}(1:end-10) 'Bi'];
            end
        end
        
        for combine_loop=1:size(fit.whatcombine,1)
            combinebold1=['\n getting and combining BOLD from different ROIs \n'];
            combinebold2= ['Process number ' num2str(combine_loop) ' out of ' num2str(size(fit.whatcombine,1)) ' combinations '];
            combinebold3=[' Doing ROIs number ' num2str(fit.whatcombine(combine_loop,:))];
            fprintf(combinebold1);fprintf(combinebold2);fprintf(combinebold3);
            for roi_i=1:size(fit.whatcombine,2)
                ROI_temp{1,roi_i}=fit.ROIs{1,fit.whatcombine(combine_loop,roi_i)};
            end
            if ~strcmp(fit.timings{1,combine_loop},fit.Target_timing)
                equal_regressors=0;
            end
            fitcomb=fit;
            fitcomb.timings={fit.timings{:,combine_loop}}'; fitcomb.ROIs=ROI_temp;
            evalc('[~,onsets_filtered,BOLDpersession,behavpersessioncombined,session_num,sub_tgt_id,sub_coord{:,combine_loop}]=corr_BOLD(fitcomb,0)');
            all_BOLD_combined{1,combine_loop}=corct_BOLD_combined;
            all_behav_combined{1,combine_loop}=behav_onecombined;
            
            if ~strcmp(fit.timings{1,combine_loop},fit.Target_timing)
                [BOLD_output,behav_output,~]=corct_timing(BOLDpersession,onsets_filtered,target_onsets,behavpersessioncombined);
                all_BOLD_combined{1,combine_loop}=BOLD_output;
                all_behav_combined{1,combine_loop}=behav_output;
                comb_temp=zeros(0,1);
                comb_temp_behav=zeros(0,1);
                for ind_fill2= 1:size(BOLD_output,2)
                    comb_temp=vertcat(comb_temp, all_BOLD_combined{1,combined_ROIs_loop}{1,ind_fill2});
                    comb_temp_behav=vertcat(comb_temp_behav, all_behav_combined{1,combined_ROIs_loop}{1,ind_fill2});
                end
                all_BOLD_combined{1,combined_ROIs_loop}=comb_temp;
                all_behav_combined{1,combined_ROIs_loop}=comb_temp_behav;
                clear comb_temp
                clear ROI_temp
                clear BOLDpersession
                clear onsets_filtered
                clear behavpersessioncombined
            else
                for ind_fill1= 1:size(BOLDpersession,2)
                    all_BOLD_combined{1,combined_ROIs_loop}{1,ind_fill1}=BOLDpersession{1,ind_fill1};
                    all_behav_combined{1,combined_ROIs_loop}{1,ind_fill1}=behavpersessioncombined{1,ind_fill1};
                end
                comb_temp=zeros(0,1);
                comb_temp_behav=zeros(0,1);
                for ind_fill1= 1:size(BOLDpersession,2)
                    comb_temp=vertcat(comb_temp, all_BOLD_combined{1,combined_ROIs_loop}{1,ind_fill1});
                    comb_temp_behav=vertcat(comb_temp_behav, all_behav_combined{1,combined_ROIs_loop}{1,ind_fill1});
                end
                all_BOLD_combined{1,combined_ROIs_loop}=comb_temp;
                all_behav_combined{1,combined_ROIs_loop}=comb_temp_behav;
                clear comb_temp
                clear comb_temp_behav
                clear BOLDpersession
                clear onsets_filtered
                clear behavpersessioncombined
            end
            
        end
        
        idx=1:1:size(fit.ROIs,2);
        whatcombinevector=(reshape(fit.whatcombine,size(fit.whatcombine,1)*size(fit.whatcombine,2),1))';
        for find_rest=1:size(idx,2)
            if  ~isempty(find(whatcombinevector==idx(find_rest)))
                idx(find_rest)=0;
            end
        end
        rest_idx= nonzeros(idx);
        for single_roi_fill= 1:size(rest_idx,1)
            single_ROIs{1,single_roi_fill}=fit.ROIs{1,rest_idx(single_roi_fill)};
        end
        
    else
        ROIs_names=fit.ROIs;
        timings_for_BOLD=fit.timings(1,:);
        single_ROIs=fit.ROIs;
    end
    
    
    %*************Getting Sinlge ROIs BOLD vectors************
    if exist('single_ROIs','var')
        for single_ROIs_loop= 1:size(single_ROIs,2)
            gettingbold=[' Getting BOLD from ROI number ' num2str(single_ROIs_loop) ' out of ' num2str(size(fit.ROIs,2)) ' ROIs \n '];
            fprintf(gettingbold);
            if ~strcmp(fit.timings{1,single_ROIs_loop},fit.Target_timing)
                equal_regressors=0;
            end
            fitsing=fit;
            fitsing.timings={fit.timings{:,single_ROIs_loop}}';
            fitsing.ROIs=single_ROIs{single_ROIs_loop};
            evalc('[~,onsets_filtered,BOLDpersession,behavpersessionsingles,session_num,~,sub_coord{single_ROIs_loop}]=corr_BOLD(fitsing,0)');

          %  if ~strcmp(fit.timings{1,single_ROIs_loop},fit.Target_timing)
                [BOLD_output,behav_output,~]=corct_timing(BOLDpersession,onsets_filtered,target_onsets,behavpersessionsingles);
                
                all_BOLD_singles{1,single_ROIs_loop}=BOLD_output;
                all_behav_singles{1,single_ROIs_loop}=behav_output;
                for ind_fill1= 1:size(BOLD_output,2)
                    %                     if fit.zscored==1
                    %                         all_BOLD_singles{1,single_ROIs_loop}{1,ind_fill1}=zscore(all_BOLD_singles{1,single_ROIs_loop}{1,ind_fill1});
                    %                         all_behav_singles{1,single_ROIs_loop}{1,ind_fill1}=zscore(all_behav_singles{1,single_ROIs_loop}{1,ind_fill1});
                    %                     end
                end
                sing_temp=zeros(0,1);
                sing_temp_behav=zeros(0,1);
                for ind_fill1= 1:size(BOLD_output,2)
                    sing_temp=vertcat(sing_temp, all_BOLD_singles{1,single_ROIs_loop}{1,ind_fill1});
                    sing_temp_behav=vertcat(sing_temp_behav, all_behav_singles{1,single_ROIs_loop}{1,ind_fill1});
                end
                all_BOLD_singles{1,single_ROIs_loop}=sing_temp;
                BOLD_back=BOLDpersession;
                clear sing_temp
                clear sing_temp_behav
                clear BOLDpersession
                clear onsets_filtered
                clear behavpersessionsingles
%             else
%                 for ind_fill1= 1:size(BOLDpersession,2)
%                     all_BOLD_singles{1,single_ROIs_loop}{1,ind_fill1}=BOLDpersession{1,ind_fill1};
%                     all_behav_singles{1,single_ROIs_loop}{1,ind_fill1}=behavpersessionsingles{1,ind_fill1};
%                 end
%                 sing_temp=zeros(0,1);
%                 sing_temp_behav=zeros(0,1);
%                 for ind_fill1= 1:size(BOLDpersession,2)
%                     sing_temp=vertcat(sing_temp, all_BOLD_singles{1,single_ROIs_loop}{1,ind_fill1});
%                     sing_temp_behav=vertcat(sing_temp_behav, all_behav_singles{1,single_ROIs_loop}{1,ind_fill1});
%                 end
%                 all_BOLD_singles{1,single_ROIs_loop}=sing_temp;
%                 all_behav_singles{1,single_ROIs_loop}=sing_temp_behav;
%                 BOLD_back=BOLDpersession;
%                 clear sing_temp
%                 clear sing_temp_behav
%                 clear BOLDpersession
%                 clear onsets_filtered
%                 clear behavpersessionsingles
%             end
            
        end
        all_corct_BOLD=horzcat(all_BOLD_singles,all_BOLD_combined);
        all_corct_behav=horzcat(all_behav_combined,all_behav_singles);
    else
        all_corct_BOLD=all_BOLD_combined;
        all_corct_behav=all_behav_combined;
    end
    
    if fit.bimodal==1
        all_corct_reg=horzcat(all_corct_BOLD,all_corct_behav);
    else
        all_corct_reg=all_corct_BOLD;
    end
    
else
    fprintf('\n ROIs and timings size mismatch \n make sure to put the first two parameters in curly braces \n')
    return
end

fprintf('\n Loading the BOLD arrays done \n')
proper_model_size=size(all_corct_reg,2);

if fit.bimodal==1
    all_param=horzcat(ROIs_names,fit.timings(1,:));
    all_timings=horzcat(timings_for_BOLD(1,:),fit.timings(1,:));
else
    all_param=  ROIs_names;
    all_timings=timings_for_BOLD(1,:);
end
%%
%******Preparing the fit********
terms=zeros(0,0);
for ind_fill= 1:proper_model_size
    terms=[terms all_corct_reg{1,ind_fill}];
end
%terms=terms';
nobs=size(terms,1);

clear CM
last_sub=0;
last_row=0;
if ~strcmp(fit.timings{1,single_ROIs_loop},fit.Target_timing)
    BOLD_ref=BOLD_output;
else
    BOLD_ref=BOLD_back;
end
CV=zeros(0,1);
S_i=cell(size(session_num,1),1);
for const_m= 1:size(session_num,1)
    temp_ones=zeros(0,1);
    CM(1:size(terms,1),const_m)=zeros(size(terms,1),1);
    for sess_ones=last_sub+1:last_sub+session_num(const_m)
        temp_ones=vertcat(temp_ones,ones(size(BOLD_ref{1,sess_ones},1),1));
    end
    temp_ones_s=size(temp_ones,1);
    cv_sub=const_m-1+ones(size(temp_ones,1),1);
    CV=vertcat(CV,cv_sub);
    CM(last_row+1:(last_row+size(temp_ones,1)),const_m)= temp_ones;
    S_i{const_m,1}=last_row+1:(last_row+size(temp_ones,1));
    last_row=last_row+size(temp_ones,1);
    last_sub=last_sub+session_num(const_m);
    clear cv_sub
end
fit.S_i=S_i;
if fit.zscored==0
    terms=horzcat(terms,CM);
end
nanterms=find(isnan(terms));
terms(nanterms)=0;
meanresult = mean(tgt,1);
re=nansum((tgt-meanresult).^2);
for fill_dum= 1:size(terms,2)
    inx(1,fill_dum) = 0.5;
end
MEdat_rlv=cell(size(fit.sub,2),0);
if isa('fit.ME','char') && fit.pooled==1
    
    MEdat=readtable(fit.ME,'FileType','spreadsheet');
    for iii= 2:size(MEdat,1)
        ME_ids(iii,1)=str2num(MEdat{iii,1}{1,1});
    end
    for iii2=1:size(fit.sub,2)
        gs=find(ME_ids==fit.sub(1,iii2));
        if ~isempty(gs)
            gs=nonzeros(gs);
            idsub(1,iii2)=gs(1,1);
        end
    end
    idsub=nonzeros(idsub);
    MEdat_rlv=MEdat{idsub,1:end};
    gender=zeros(size(MEdat_rlv,1),1);
    for gender_i=1:size(MEdat_rlv,1)
        if MEdat_rlv{gender_i,2}(1,1)=='W'
            gender(gender_i,1)=0.5;
        elseif MEdat_rlv{gender_i,2}(1,1)=='M'
            gender(gender_i,1)=-0.5;
        else
            gender(gender_i,1)=0;
        end
    end
    terms(:,end+1)=zeros(size(terms,1),1);
    for get_sub=1:size(MEdat_rlv,1)
        thissub=str2num(MEdat_rlv{get_sub,1});
        subdata= find(sub_tgt_id==thissub);
        terms(subdata,end)=gender(get_sub,1);
    end
    fit.ROIs{1,size(fit.ROIs,2)+1}='Gender';
    fit.timings{1,size(fit.timings,2)+1}='All timings';
end

            ROIs2=fit.ROIs; 
            %             if isa('ME','char') && fit.pooled==1
            %                 catvars=[size(terms,2)];
            %             else
            %                 catvars=[];
            %             end
            lbl={};
            for ROIs_i=1:size(all_corct_reg,2)
                for subROIs_lbl_i=1:size(all_corct_reg{1,ROIs_i},2)
                    lbl=horzcat(lbl,[fit.ROIs{ROIs_i}(1:end-8) '_sub_' num2str(subROIs_lbl_i) ' at ' fit.timings{1,ROIs_i}]);
                end
            end
            lbl=horzcat(lbl,fit.ROIs{1,size(all_corct_reg,2)+1:end});
            lbl{1,size(lbl,2)+1}='Happiness_Ratings';
            lbl_renamed=cellfun(@(x) replace(x,'.',''),lbl,'un',0);
            lbl_renamed=cellfun(@(x) replace(x,' ','_'),lbl_renamed,'un',0);
            lbl_renamed=cellfun(@deblank,lbl_renamed,'un',0);
            design=[ array2table(terms) array2table(tgt)];
            design.Properties.VariableNames=lbl_renamed;
            fit.terms=terms; fit.design=design; fit.tgt=tgt; fit.sub_coord=sub_coord;
            allvar_name=[ name_fitBOLD(fit) num2str(fit.cluster.mink) '-' num2str(fit.cluster.maxk) 'k-clusters_' num2str(size(fit.sub,2)) 'Subjects.mat'];
            fit.quadsize=(size(fit.design,2)-2)*2+1;
folder='/Volumes/EVO 860 2TB/Designs';
try save(fullfile(folder,allvar_name)); end

end