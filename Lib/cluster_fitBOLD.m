   %% cluster_fitBOLD
function [clustered_BOLD,sub_coord,lbls]=cluster_fitBOLD(unclustered_BOLD_filt,coord,mm3,fit)
krange=fit.cluster.mink:fit.cluster.maxk;
preclustered=dir(fit.subROIsdir);
Subjects=fit.sub;
lbls={};
basename=[fit.ROIs(1:end-8) '_subROI_All' num2str(size(Subjects,2))];
savename_roi=[ basename 'Subjects.mat'];
findpreclustered=strfind({preclustered(:).name},[fit.ROIs(1:end-8) '_subROI_All']);
preclusteredloc=find(~cellfun(@isempty,findpreclustered));

if ~isempty(preclusteredloc)
    load([preclustered(preclusteredloc).folder filesep  preclustered(preclusteredloc).name])
    for subROIs_i= 1:size(sub_coord,2)
    lbls{subROIs_i} = [fit.subROIsdir filesep fit.ROIs(1:end-8) '_subROI_' num2str(subROIs_i) '_' num2str(size(Subjects,2)) 'Subjects.mat'];
    end
 else
    
    if fit.cpu_parallel>0 || fit.disk_parallel>0; cluspar=1 ; else cluspar=0; end
    for k_i=1:size(krange,2)
        if size(unclustered_BOLD_filt,2)<fit.cluster.minvol; k=0 ; else; k=krange(k_i); end
    [allclustered_BOLD{k_i},~,allsub_coord{k_i},allsilval{k_i}]=Cluster_ROIs(unclustered_BOLD_filt,k,coord{1,1},fit.cluster.distance,1,0,cluspar,fit);
    silperk(k_i,1)=krange(k_i);
    silperk(k_i,2)=mean(allsilval{k_i});
    end
    bestk=silperk(find(silperk(:,2)==max(silperk(:,2))),1);
    bestk_i=find(krange==bestk);
     clustered_BOLD=allclustered_BOLD{bestk_i};
     sub_coord=allsub_coord{bestk_i};
     figure; h=barh(allsilval{bestk_i});
     CA=gca; CA.XLim=[-1,1];
%     for gg= 1:size(clustered_BOLD,2)
%         gg2(:,gg)=mean(clustered_BOLD{1,gg},2);
%     end
%     figure; plot(gg2);
    savefig_fitBOLD(fit.ROIs,fit.subROIsdir)
    save([fit.subROIsdir filesep savename_roi],'sub_coord','clustered_BOLD');
    save([fit.dest_dir filesep 'mm3.mat'],'mm3')
    
    for subROIs_i= 1:size(sub_coord,2)
        img2 = zeros(91,109,91);
        lbl_roi = [fit.subROIsdir filesep fit.ROIs(1:end-8) '_subROI_' num2str(subROIs_i) '_' num2str(size(Subjects,2)) 'Subjects.nii'];
        lbl_roi_mat=[fit.subROIsdir filesep fit.ROIs(1:end-8) '_subROI_' num2str(subROIs_i) '_' num2str(size(Subjects,2)) 'Subjects.mat'];
        for fill_3d=1:size(sub_coord{1,subROIs_i},2)
            img2( sub_coord{1,subROIs_i}(1,fill_3d), sub_coord{1,subROIs_i}(2,fill_3d), sub_coord{1,subROIs_i}(3,fill_3d))=1;
        end
        img_roi = maroi_pointlist(struct('XYZ',sub_coord{1,subROIs_i},'mat',mm3,'descrip','', 'label', lbl_roi),'vox');   
        saveroi(img_roi,lbl_roi_mat)
        save_as_image(img_roi,lbl_roi);
        lbls{subROIs_i}=lbl_roi_mat;
        clear lbl_roi
    end
 end
end