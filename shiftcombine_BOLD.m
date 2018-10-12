function [Subjects_BOLD,shifts,pooledBOLD]=shiftcombine_BOLD(ROIs,timings,marsbar,dir,sub,shift,combine)

for ROI_loop=1:size(ROIs,2)
ROI_shifts{1,ROI_loop}=BOLD(ROIs{ROI_loop},timings,marsbar,dir,sub,shift);
end

for person_i =1:size(sub,2)
 BOLD_mat=zeros(0,0);
BOLD_vec=[ROI_shifts{1,ROI_loop}{1,person_i}.BOLD{1,1}];
BOLD_mat=horzcat(BOLD_mat,BOLD_vec);
clear BOLD_vec
if strcmp(combine,'mean')
    finalBOLD=mean(BOLD_mat,2);
elseif strcmp(combine,'median')
    finalBOLD=median(BOLD_mat,2);
end
clear BOLD_mat
tempstr=ROI_shifts{1,1}{1,person_i};
for ROI_names_shifts=1:size(ROIs,2)
    if ROI_names_shifts >1
tempstr_name=[tempstr.ROI(1:end-4) ' and ' ROI_shifts{1,ROI_names_shifts}{1,person_i}.ROI(1:end-4)];
tempstr.ROI=tempstr_name;
    end
tempstr.shifts(ROI_names_shifts,1)=ROI_shifts{1,ROI_names_shifts}{1,person_i}.shifts;
end
tempstr.BOLD{1,1}=finalBOLD;
combined_ROIs{1,person_i}=tempstr;
clear tempstr
end
Subjects_BOLD=combined_ROIs;
Subjects_shifts=zeros(0,1);
pool_sub.BOLD=zeros(0,0);
for i= 1:size(sub,2)
        for i3= 1:size(ROIs,2)
            if ~isempty(Subjects_BOLD{1,i}.shifts{1,1})
                Subjects_shifts(end+1,1)=Subjects_BOLD{1,i}.shifts{i3,1};
                pool_sub.BOLD=vertcat(pool_sub.BOLD,Subjects_BOLD{1,i}.BOLD{1,1});
            end
    end  
end
 shifts=Subjects_shifts;
 pooledBOLD=pool_sub;
end