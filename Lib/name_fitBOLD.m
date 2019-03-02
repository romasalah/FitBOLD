function modelName=name_fitBOLD(fit)
di=strfind(fit.GLMs,'/');
di=di(end,end);
glm_name=fit.GLMs(di+1:end-8);
if fit.behav==1 tgt_type='behav'; else tgt_type='BOLD'; end
if fit.combine==1 comb_type='Bi'; else comb_type='Uni';end
if fit.pooled==1 pooled_type='pool'; else pooled_type='ind'; end
if strcmp(fit.shifts,'standard shift') shifttype='stdshft'; else shifttype='trushft'; end
if fit.condition_before; which_con_before=fit.condition_before; else which_con_before='NSCB';end
if fit.social_ROIs==1; social_ROIs_name='SocROI'; elseif fit.social_ROIs==2; social_ROIs_name='ToM'; else social_ROIs_name='NSocROI';end
if fit.zscored; zscore_state='Zscr'; else zscore_state='intrcpt'; end
if fit.AR~=0 AR_state= ['AR(n)']; else AR_state='';end
if isa(fit.MA,'double') && fit.MA~=0 MA_state= ['MA(' num2str(fit.MA) ')']; elseif isa(fit.MA,'double') && fit.MA==0 MA_state=''; elseif ~isa(fit.MA,'double') MA_state=fit.MA.filterstructure; end
modelName=[glm_name '_' fit.models '_' tgt_type '_' fit.scans '_' pooled_type '_' comb_type '_' ...
    shifttype '_' which_con_before '_' zscore_state '_' fit.dependence '_' fit.resample.method '_' social_ROIs_name '_' AR_state '_' MA_state];
end