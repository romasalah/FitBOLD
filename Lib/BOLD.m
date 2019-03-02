%This Script needs MarsBar tool for functionality
function [Subjects_BOLD,scan_shift,sub_coord]=extract_fitBOLD(fit)
clear BOLD Subject_BOLD mean_HRF SPM
maskdir=fit.maskdir;
use_subROIs=fit.use_subROIs;
datadir=fit.GLMs;
cd(datadir)
if fit.use_marsbar==1
    use_marsbar=1;
    marsbar_dir='/Applications/spm12/toolbox/marsbar';
end
Subjects=[fit.sub];
marsbar_ROIs=fit.ROIs;
conds_id= fit.timings{1,:};
%% *****Dummies****
zBOLD_array_vox=cell(1,size(Subjects,2));
Subjects_BOLD=cell(1,size(Subjects,2));
nscans=[];
addpath(marsbar_dir);
conds_id=cellstr(conds_id);
prefix='SODEC_FMRI_';
j=length( Subjects );
BOLD=cell(1,j);
callperson=cell(1,j);personName=cell(1,j);Doing_subject=cell(1,j);matdir=cell(1,j);
BOLD_array_means=cell(1,j);Conditions_lookup=cell(1,j);
conds_derv=cell(1,j);conds_betas=cell(1,j);mean_HRF=cell(1,j);marsbar_BOLD=cell(1,j); all_BOLD_means=cell(1,j); shifts_all=cell(1,j);
BOLD_array_means_shifted=cell(1,j);
shifts_conds=cell(1,j);
unclustered_BOLD=[];
sub_coord=0;
unclustered_BOLD_all=[];
%% load all Subjects data and masks
if ~exist('SODECallspmloadeddata','var')
    global SODECallspmloadeddata
    %loop through the folders
    for add_files= 1 : length( Subjects )
        callperson{add_files} = [prefix num2str(Subjects(add_files))]; %calls a specific folder according to index
        personName{add_files} = callperson{add_files}; %get the directory as a character array
        matdir{add_files}=[ datadir '/' personName{add_files} '/SPM.mat']; %specifies the required matrix to open
        %     poolobj = gcp;
        %     addAttachedFiles(poolobj,matdir{add_files})
        SODECallspmloadeddata{add_files}=load(matdir{add_files});
    end
end

cd(maskdir)
mask=load(marsbar_ROIs);



%% Start extracting Unclustered BOLD
parfor (person_i = 1 : length( Subjects ),fit.disk_parallel)
    %for person_i = 1 : length( Subjects )
    BOLD{person_i}=struct('ROI',marsbar_ROIs,'Conditions','','BOLD','','shifts','');
    callperson{person_i} = [prefix num2str(Subjects(person_i))]; %calls a specific folder according to index
    personName{person_i} = callperson{person_i}; %get the directory as a character array
    all_BOLD_means{person_i}=cell(1,size(marsbar_ROIs,2));
%% Get raw BOLD.
if fit.use_4D==1
    sp=strfind(SODECallspmloadeddata{person_i}.SPM.xY.VY(1).fname,filesep);
    dest=SODECallspmloadeddata{person_i}.SPM.xY.VY(1).fname(1:sp(1,end));
    img=[dest '4D.nii'];
    [marsbar_BOLD{person_i},~,coord{person_i},mm3map{person_i}]=getdata(mask.roi,img);
else
    marsbar_BOLD{person_i}=getdata(mask.roi,SODECallspmloadeddata{person_i}.SPM.xY.VY);
end
%% filter BOLD
[ marsbar_BOLD_filt{person_i}]=LSD_fitBOLD(marsbar_BOLD{person_i});
%% get mean activation or zscore BOLD for clustering
if use_subROIs==1
    for voxels2_i=1:size(marsbar_BOLD_filt{person_i},2)
        zBOLD_array_vox{person_i}{1,1}(:,voxels2_i)=zscore(marsbar_BOLD_filt{person_i}(:,voxels2_i));
    end
else
    for voxels_loop=1:size(marsbar_BOLD_filt{person_i},1)
        if fit.zscored==1
            BOLD_array_means{person_i}{1,1}(voxels_loop,1)=zscore(nanmean(marsbar_BOLD_filt{person_i}(voxels_loop,:),2));
        else
            BOLD_array_means{person_i}{1,1}(voxels_loop,1)=nanmean(marsbar_BOLD_filt{person_i}(voxels_loop,:),2);
        end
    end
end
end
global mm3
mm3=mm3map{1,1};
%% Cluster BOLD

if use_subROIs && ~strcmp(fit.TargetBOLD,fit.ROIs)
    for cat=1:size( zBOLD_array_vox,2)
      unclustered_BOLD_all=vertcat(unclustered_BOLD_all, zBOLD_array_vox{cat}{1,1});
    end
    [clustered_BOLD,sub_coord,lbls]=cluster_fitBOLD(unclustered_BOLD_all,coord,mm3,fit);
    for person_i = 1 : length( Subjects )
        nscans(person_i)=size(zBOLD_array_vox{person_i}{1,1},1);
        for subset_ROI=1:size(clustered_BOLD,2)
            BOLD_array_means{person_i}{1,subset_ROI}=clustered_BOLD{1,subset_ROI}(1:nscans(person_i),:);
        end
    end
end

%% Get shifts
if fit.shifts==1
scan_shift=shiftHRF_fitBOLD(fit,SODECallspmloadeddata,lbls);
else
   scan_shift= zeros(1,size(Subjects,2))+2;
end

%% do the shifts

%parfor (person_i = 1 : length( Subjects ),fit.cpu_parallel)
for person_i = 1 : length( Subjects )
BOLD_array_means_shifted{person_i}=cell(0,0);
for subset_ROI=1:size(BOLD_array_means{person_i},2)
BOLD_array_means_shifted{person_i}{1,subset_ROI}=vertcat(BOLD_array_means{person_i}{1,subset_ROI}(1+scan_shift(person_i,subset_ROI):end,:),zeros(scan_shift(person_i),size(BOLD_array_means{person_i}{1,subset_ROI},2)));
end
BOLD{person_i}.shifts=scan_shift(person_i);
BOLD{person_i}.Conditions=fit.timings{1,1};
BOLD{person_i}.BOLD=BOLD_array_means_shifted{person_i};
Subjects_BOLD{1,person_i}=BOLD{person_i};
end
end
