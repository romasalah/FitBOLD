%% LSD_fitBOLD
function [unclustered_BOLD_filt]=LSD_fitBOLD(unclustered_BOLD)
n=size(unclustered_BOLD,1);
for voxels_i=1:size(unclustered_BOLD,2)
    sd2=2*std(unclustered_BOLD(:,voxels_i));
    m=mean(unclustered_BOLD(:,voxels_i));
    for corct_outliers=1:n
        if unclustered_BOLD(corct_outliers,voxels_i)>m+sd2 || unclustered_BOLD(corct_outliers,voxels_i)< m-1*sd2
            unclustered_BOLD(corct_outliers,voxels_i)=m;
        end
    end
end
filter= BOLD_filt(unclustered_BOLD(:,1),'LSD');
for t1=1:size(unclustered_BOLD,2)
    unclustered_BOLD_filt(:,t1)=filtfilt(filter,unclustered_BOLD(:,t1));
end
end
