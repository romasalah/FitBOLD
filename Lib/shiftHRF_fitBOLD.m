%% shiftHRF_fitBOLD
function [allshifts]=shiftHRF_fitBOLD(fit,SODECallspmloadeddata,mask)
Subjects=[fit.sub];
conds_id= cellstr(fit.timings{1,:});
prefix='SODEC_FMRI_';
%% *****Dummies****
j=length( Subjects );
callperson=cell(1,j);personName=cell(1,j);
marsbar_loc_mean=cell(1,j);
marsbar_sess_mean=cell(1,j);
conds_derv=cell(1,j);conds_betas=cell(1,j);mean_HRF=cell(1,j);
scan_shift=cell(1,j);cond_loc=cell(1,j);
cond=cell(1,j);marsbar_beta=cell(1,j);correct_x=cell(1,j);
derv_look=cell(1,j);Peak_idx=cell(1,j);Peak_scan=cell(1,j);
for subROIs_i=1:size(mask,2)
    allsubROIs{subROIs_i}=load(mask{subROIs_i});
end
parfor (person_i = 1 : length( Subjects ),fit.disk_parallel)
   % for person_i = 1 : length( Subjects )
    if size(SODECallspmloadeddata{person_i}.SPM.xBF.bf,2)==3
        conds_derv{person_i}={'*bf(1)','*bf(2)','*bf(3)'};
    elseif size(SODECallspmloadeddata{person_i}.SPM.xBF.bf,2)==1
        conds_derv{person_i}={'*bf(1)'};
    else
        fprintf('Cannot know the derivatives situation')
    end
    conds_betas{person_i}=cell(size(conds_id,2),size(conds_derv{person_i},2)+1);
    mean_HRF{person_i}=cell(1,size(conds_id,2));
    cd(fit.GLMs)
    callperson{person_i} = [prefix num2str(Subjects(person_i))]; %calls a specific folder according to index
    personName{person_i} = callperson{person_i}; %get the directory as a character array
    cd (personName{person_i})
    for subROIs_i=1:size(mask,2)
        cursubROI=allsubROIs{subROIs_i};
        cond_loc{person_i}=zeros(1,0);
        % conds_betas.condition=[conds_id{conds_loop}];
        for derv_loop= 1:size(conds_derv{person_i},2)
            cond{person_i}=cellstr(conds_id);
            derv_look{person_i}=cellstr([conds_derv{person_i}{derv_loop}]);
            %             rep=['\nCondition: ' cond{person_i}{1,1} '  HRF Derivative: ' derv_look{person_i}{1,1} '\n'];
            %             fprintf(rep)
            %             BOLD{person_i}(masks_loop).Conditions{1,conds_loop}=cond{person_i};
            %             fprintf('\n \n  New HRF shape/scan shift \n')
            
            for mat_loop= 1:size(SODECallspmloadeddata{person_i}.SPM.Vbeta,2)
                if ~isempty(strfind(SODECallspmloadeddata{person_i}.SPM.Vbeta(mat_loop).descrip,cond{person_i}{1,1}))...
                        && ~isempty(strfind(SODECallspmloadeddata{person_i}.SPM.Vbeta(mat_loop).descrip,derv_look{person_i}{1,1}))...
                        && isempty(strfind(SODECallspmloadeddata{person_i}.SPM.Vbeta(mat_loop).descrip,'xVal'))
                    cond_loc{person_i}(1,end+1)=mat_loop;
                    %                     reg=['found at regressor number: ' num2str(mat_loop) '\n'];
                    %                     fprintf(reg)
                    marsbar_beta{person_i}=cell(1,size(cond_loc{person_i},2));
                    for loc_loop= 1:size(cond_loc{person_i},2)
                        marsbar_beta{person_i}{loc_loop}=getdata(cursubROI.roi,SODECallspmloadeddata{person_i}.SPM.Vbeta(cond_loc{person_i}(loc_loop)).fname);
                        marsbar_loc_mean{person_i}(loc_loop)=nanmean(marsbar_beta{person_i}{loc_loop},2);
                    end
                    marsbar_sess_mean{person_i}=nanmean(marsbar_loc_mean{person_i},2);
                end
            end
            conds_betas{person_i}{subROIs_i}{1,derv_loop}=marsbar_sess_mean{person_i};
        end
        %         if marsbar_sess_mean{person_i}==0
        %             conditions_lookup{person_i}=Conditions_lookup{person_i};
        %             fprintf('\n HRF beta estimates are Zeros, \n which mostly means that I cannot find the condition timing your HRF on')
        %             fprintf('\n I have collected a list of conditions in this GLM in Variable called Conditions_lookup \n')
        %             fprintf(' please make sure to have a look at it \n')
        %             %                         conditions_lookup{person_i};
        %             %                         return
        %         end
        %Get HRF shape and Peak of each condition
        mean_HRF{person_i}{subROIs_i}=SODECallspmloadeddata{person_i}.SPM.xBF.bf*([conds_betas{person_i}{subROIs_i}{1,1:end}]');
        Peak_idx{person_i}{subROIs_i}=find(mean_HRF{person_i}{subROIs_i}==max(abs(mean_HRF{person_i}{subROIs_i})));
        correct_x{person_i}{subROIs_i}=[0:1/(SODECallspmloadeddata{person_i}.SPM.xBF.T):(size(SODECallspmloadeddata{person_i}.SPM.xBF.bf,1)/(SODECallspmloadeddata{person_i}.SPM.xBF.T))];
        Peak_scan{person_i}{subROIs_i}=correct_x{person_i}{subROIs_i}(Peak_idx{person_i}{subROIs_i});
        scan_shift{person_i}{subROIs_i}=round(Peak_scan{person_i}{subROIs_i});
        if isempty(scan_shift{person_i}{subROIs_i}) || scan_shift{person_i}{subROIs_i}<1 || scan_shift{person_i}{subROIs_i}>4
            fprintf('\n Scan shift is far from what we expect biologically.\n It will be reverted back to a standard shift of 2 scans \n')
            scan_shift{person_i}{subROIs_i}=2;
        end
    end
end

for subROIs_i=1:size(mask,2)
    dirsep=strfind(mask{subROIs_i},filesep);
    charsep=strfind(mask{subROIs_i},'_');
    subROIlbl=mask{subROIs_i}(dirsep(end)+1:charsep(end)-1);
    subROIlbl=strrep(subROIlbl,'_',' ');
    figure;
   
for i=1:size(mean_HRF,2)
    plot(mean_HRF{i}{subROIs_i}); hold on;
    allshifts(i,subROIs_i)=[scan_shift{i}{subROIs_i}];
end
 title([' Mean HRF shape over all social conditions in ' subROIlbl ' during ' fit.timings{1,1} ' for ' num2str(size(fit.sub,2)) ' Subjects'])
end

xlabel('microtime resolution = 1/16 seconds.')

end
