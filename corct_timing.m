function [BOLD_output,behav_output,scans_output]=corct_timing(BOLD,timings,target_timings,behavioural)
scans_output_1=cell(1,size(target_timings,2));
BOLD_output_1=zeros(size(target_timings,1),1);
behav_output_1=zeros(size(target_timings,1),1);
BOLD_output=zeros(1,0);
behav_output=zeros(1,0);
scans_output=cell(1,0);

for sess_i=1:size(timings,2)
    for tgt_i= 1:size(target_timings{1,sess_i},1)
        tgt_time_sort=sort(target_timings{1,sess_i});
        sort_ons=find(tgt_time_sort==target_timings{1,sess_i}(tgt_i,1));
        if sort_ons==1
           temp_ons=zeros(0,1);
            BOLD_ons=zeros(0,1);
            behav_ons=zeros(0,1);
            for tim_i= 1:size(timings{1,sess_i},1)
                if timings{1,sess_i}(tim_i,1)<=target_timings{1,sess_i}(tgt_i,1)
                    temp_ons(end+1,1)=timings{1,sess_i}(tim_i,1);
                    BOLD_ons(end+1,1)=BOLD{1,sess_i}(tim_i,1);
                   try behav_ons(end+1,1)=behavioural{1,sess_i}(tim_i,1);end
                else
                    BOLD_ons(end+1,1)=mean([BOLD{1,sess_i}],1);
                    try behav_ons(end+1,1)=mean([behavioural{1,sess_i}],1); end
                end
            end
            BOLD_m=nanmean(BOLD_ons,1);
            behav_m=nanmean(behav_ons,1);
            scans_output_1{1,tgt_i}=temp_ons;
            BOLD_output_1(tgt_i,1)=BOLD_m;
            behav_output_1(tgt_i,1)=behav_m;
            clear temp_ons
            clear BOLD_ons
            clear behav_ons

        else
            temp_ons=zeros(0,1);
            BOLD_ons=zeros(0,1);
            behav_ons=zeros(0,1);
            for tim_i= 1:size(timings{1,sess_i},1)
                if timings{1,sess_i}(tim_i,1)<=target_timings{1,sess_i}(tgt_i,1) && timings{1,sess_i}(tim_i,1)>tgt_time_sort(sort_ons-1,1)
                    temp_ons(end+1,1)=timings{1,sess_i}(tim_i,1);
                    BOLD_ons(end+1,1)=BOLD{1,sess_i}(tim_i,1);
                    try behav_ons(end+1,1)=behavioural{1,sess_i}(tim_i,1); end
               else
                    BOLD_ons(end+1,1)=mean([BOLD{1,sess_i}],1);
                   try behav_ons(end+1,1)=mean([behavioural{1,sess_i}],1);end
                end
            end
            BOLD_m=nanmean(BOLD_ons,1);
            behav_m=nanmean(behav_ons,1);
            scans_output_1{1,tgt_i}=temp_ons;
            BOLD_output_1(tgt_i,1)=BOLD_m;
            behav_output_1(tgt_i,1)=behav_m;
            clear temp_ons
            clear BOLD_ons
            clear behav_ons
        end
    end
    BOLD_output{1,end+1}=BOLD_output_1;
    behav_output{1,end+1}=behav_output_1;
    scans_output{1,end+1}=scans_output_1;
   clear scans_output_1
   clear BOLD_output_1
   clear behav_output_1
end
end