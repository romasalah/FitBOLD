%This Script needs MarsBar tool for functionality
function [Subjects_BOLD,shifts,pooledBOLD,conditions_lookup]=BOLD(ROIs,Conditions,marsbar,dir,subjects,shift,disk_parallel,use_4D)
clear BOLD Subject_BOLD mean_HRF SPM
maskdir='/Users/NEMO_admin/Desktop/Omar/Masks/';

if exist('dir','var')
    datadir=dir;
else
    fprintf('No directory specified')
    return
end
cd(datadir)

if marsbar==1
    use_marsbar=1;
    marsbar_dir='/Applications/spm12/toolbox/marsbar';
end

if ~exist('subjects','var')
    Subjects = [44 62 63 65 69 70 71 74 77 79 81 83 84 90 92 93 94 95 96 98 101 102 103 104 106 110 120 124 125 126 129 136 155 166 171 173 174 175 176 177];
else
    Subjects=[subjects];
end
% if ~isempty(strfind(ROIs,'all_roi'))
%     %*****Directories*****
%     marsbar_ROIs={'Frontal Medial Orbital (minimum)_0_35_-20_roi.mat','Frontal_medial_Orbital_1_48_16_roi.mat','Frontal_medial_Orbital_3_16_-30_roi.mat'...
%         ,'Negative RPE ALE Mask_-24_-4_6_roi.mat','Negative RPE ALE Mask_-30_14_64_roi.mat','Negative RPE ALE Mask_-30_28_6_roi.mat'...
%         ,'Negative RPE ALE Mask_-42_-8_48_roi.mat','Negative RPE ALE Mask_4_2_60_roi.mat','Orbitofrontal_bilateral_AAl_2_42_-13_roi.mat'...
%         ,'Positive RPE ALE mask_-16_0_70_roi.mat','Positive RPE ALE mask_-20_2_-22_roi.mat','Positive RPE ALE mask_-23_2_0_roi.mat'...
%         ,'Positive RPE ALE mask_4_8_38_roi.mat','Positive RPE ALE mask_25_3_0_roi.mat','Positive RPE ALE mask_56_12_18_roi.mat'...
%         ,'VENTRAL_PALLIDUM_-10_2_-8_roi.mat','VENTRAL_PALLIDUM_11_3_-8_roi.mat'};
% else
%     marsbar_ROIs=ROIs;
% end
marsbar_ROIs=ROIs;
%****Names****
% conds_prefix='spm_spm:beta (0041) - Sn(1) ';
% if ~isempty(strfind(Conditions,'all'))
%     conds_id={'PE Pos active','Min EV zero active','Min EV Neg active','Ch_diff active','PE Neg active','PE Pos active','happiness'};
% else
    conds_id= Conditions{1,:};
% end
conds_types={'','xVal^'};

%*****Dummies****

marsbar_beta=0;
marsbar_loc_mean=0;
BOLD_array_means=0;
BOLD=0;
Subjects_BOLD=cell(1,size(Subjects,2));
marsbar_sess_mean=0;
Conditions_lookup=0;

%************************ Start*******************
addpath(marsbar_dir);
marsbar_ROIs=cellstr(marsbar_ROIs);
conds_id=cellstr(conds_id);
prefix='SODEC_FMRI_';
j=length( Subjects );
BOLD=cell(1,j);
callperson=cell(1,j);personName=cell(1,j);Doing_subject=cell(1,j);matdir=cell(1,j);
BOLD_array_means=cell(1,j);Conditions_lookup=cell(1,j);marsbar_loc_mean=cell(1,j);
marsbar_sess_mean=cell(1,j);
conds_derv=cell(1,j);conds_betas=cell(1,j);mean_HRF=cell(1,j);marsbar_BOLD=cell(1,j);
scan_shift=cell(1,j);cond_loc=cell(1,j); all_BOLD_means=cell(1,j); shifts_all=cell(1,j);
cond=cell(1,j);marsbar_beta=cell(1,j);correct_x=cell(1,j);BOLD_array_means_conds=cell(1,j);
derv_look=cell(1,j);Peak_idx=cell(1,j);Peak_scan=cell(1,j);shifts_conds=cell(1,j);

%loop through the folders
for add_files= 1 : length( Subjects ) 
    callperson{add_files} = [prefix num2str(Subjects(add_files))]; %calls a specific folder according to index
    personName{add_files} = callperson{add_files}; %get the directory as a character array
    matdir{add_files}=[ datadir '/' personName{add_files} '/SPM.mat']; %specifies the required matrix to open
%     poolobj = gcp;
%     addAttachedFiles(poolobj,matdir{add_files})
    data{add_files}=load(matdir{add_files});
end
for masks_loop1=1:size(marsbar_ROIs,2)
    cd(maskdir)
marsbar_masks{masks_loop1}=load(marsbar_ROIs{masks_loop1});
end
        %load(matdir{person_i});
parfor (person_i = 1 : length( Subjects ),disk_parallel)
        BOLD{person_i}=struct('ROI',marsbar_ROIs,'Conditions','','BOLD','','shifts','');
        callperson{person_i} = [prefix num2str(Subjects(person_i))]; %calls a specific folder according to index
        personName{person_i} = callperson{person_i}; %get the directory as a character array
        Doing_subject{person_i}=[' \n \n Doing Subject: ' personName{person_i} '\n \n'];
        fprintf(Doing_subject{person_i})
        Conditions_lookup{person_i}= char(data{person_i}.SPM.Vbeta.descrip);
        
        if size(data{person_i}.SPM.xBF.bf,2)==3
            conds_derv{person_i}={'*bf(1)','*bf(2)','*bf(3)'};
        elseif size(data{person_i}.SPM.xBF.bf,2)==1
            conds_derv{person_i}={'*bf(1)'};
        else
            fprintf('Cannot know the derivatives situation')
        end
        conds_betas{person_i}=cell(size(conds_id,2),size(conds_derv{person_i},2)+1);
        mean_HRF{person_i}=cell(1,size(conds_id,2));
        %Using Marsbar Library
        if use_marsbar
            
            all_BOLD_means{person_i}=cell(1,size(marsbar_ROIs,2));
            shifts_all{person_i}=zeros(0,1);
            for masks_loop=1:size(marsbar_ROIs,2)
                cd(maskdir)
                
                fprintf(marsbar_ROIs{masks_loop});
                if use_4D==1
                sp=strfind(data{person_i}.SPM.xY.VY(1).fname,filesep);
                dest=data{person_i}.SPM.xY.VY(1).fname(1:sp(1,end));
                img=[dest '4D.nii'];
                marsbar_BOLD{person_i}=getdata(marsbar_masks{masks_loop}.roi,img);
                else
                marsbar_BOLD{person_i}=getdata(marsbar_masks{masks_loop}.roi,data{person_i}.SPM.xY.VY);
                %Get raw BOLD.
                end
                for voxels_loop=1:size(marsbar_BOLD{person_i},1)
                    BOLD_array_means{person_i}(voxels_loop,1)=nanmean(marsbar_BOLD{person_i}(voxels_loop,:),2);
                end
                %return to the data
                cd(datadir)
                cd (personName{person_i})
                %Get betas of HRF+derivatives.
                scan_shift{person_i}=0;
                
                for conds_loop= 1:size(conds_id,2)
                    cond_loc{person_i}=zeros(1,0);
                    conds_betas{person_i}{conds_loop,1}=[conds_id{conds_loop}];
                    for derv_loop= 1:size(conds_derv{person_i},2)
                        cond{person_i}=cellstr([conds_id{conds_loop}]);
                        derv_look{person_i}=cellstr([conds_derv{person_i}{derv_loop}]);
                        rep=['\nCondition: ' cond{person_i}{1,1} '  HRF Derivative: ' derv_look{person_i}{1,1} '\n'];
                        fprintf(rep)
                        BOLD{person_i}(masks_loop).Conditions{1,conds_loop}=cond{person_i};
                        fprintf('\n \n  New HRF shape/scan shift \n')
                        
                        for mat_loop= 1:size(data{person_i}.SPM.Vbeta,2)
                            if ~isempty(strfind(data{person_i}.SPM.Vbeta(mat_loop).descrip,cond{person_i}{1,1}))...
                                    && ~isempty(strfind(data{person_i}.SPM.Vbeta(mat_loop).descrip,derv_look{person_i}{1,1}))...
                                    && isempty(strfind(data{person_i}.SPM.Vbeta(mat_loop).descrip,'xVal'))
                                    cond_loc{person_i}(1,end+1)=mat_loop;
                                    reg=['found at regressor number: ' num2str(mat_loop)];
                                    fprintf(reg)
                                marsbar_beta{person_i}=cell(1,size(cond_loc{person_i},2));
                                for loc_loop= 1:size(cond_loc{person_i},2)
                                    marsbar_beta{person_i}{loc_loop}=getdata(marsbar_masks{masks_loop}.roi,data{person_i}.SPM.Vbeta(cond_loc{person_i}(loc_loop)).fname);
                                    marsbar_loc_mean{person_i}(loc_loop)=nanmean(marsbar_beta{person_i}{loc_loop},2);
                                end
                                fprintf('\n')
                                marsbar_sess_mean{person_i}=nanmean(marsbar_loc_mean{person_i},2);
                                % %                             else
                                % %                                 cannotfind=['\n Cannot find the Condition ' cond{1,1} ' Please recheck the name \n'];
                                % %                                 fprintf(cannotfind)
                            end
                        end
                        conds_betas{person_i}{conds_loop,derv_loop+1}=marsbar_sess_mean{person_i};
                    end
                    if marsbar_sess_mean{person_i}==0
                       conditions_lookup{person_i}=Conditions_lookup{person_i};
                        fprintf('\n HRF beta estimates are Zeros, \n which mostly means that I cannot find the condition timing your HRF on')
                        fprintf('\n I have collected a list of conditions in this GLM in Variable called Conditions_lookup \n')
                        fprintf(' please make sure to have a look at it \n')
%                         conditions_lookup{person_i};
%                         return 
                    end
                    %Get HRF shape and Peak of each condition
                    mean_HRF{person_i}{1,masks_loop}{1,conds_loop}=data{person_i}.SPM.xBF.bf*([conds_betas{person_i}{conds_loop,2:end}]');
                    Peak_idx{person_i}=find(mean_HRF{person_i}{1,masks_loop}{1,conds_loop}==max(abs(mean_HRF{person_i}{1,masks_loop}{1,conds_loop})));
                    correct_x{person_i}=[0:1/(data{person_i}.SPM.xBF.T):(size(data{person_i}.SPM.xBF.bf,1)/(data{person_i}.SPM.xBF.T))];
                    Peak_scan{person_i}=correct_x{person_i}(Peak_idx{person_i});
                    scan_shift{person_i}=round(Peak_scan{person_i});
                    if isempty(scan_shift{person_i}) || scan_shift{person_i}<1 || scan_shift{person_i}>4
                        fprintf('\n Scan shift is far from what we expect biologically.\n It will be reverted back to a standard shift of 2 scans \n')
                        scan_shift{person_i}=2;
                        %                         round(6/SPM.xY.RT)
                    end
                    if shift=='standard shift'
                        fprintf('\n \n  Using a standard shift of 2 scans \n')
                        scan_shift{person_i}=2;
                    else
                    end
                    BOLD_array_means_conds{person_i}=zeros(0,0);
                    mean_HRF{person_i}{1,masks_loop}{2,conds_loop}=scan_shift{person_i};
                    BOLD_array_means_conds{person_i}(:,conds_loop)=vertcat(BOLD_array_means{person_i}(1+scan_shift{person_i}:end),zeros(scan_shift{person_i},1));
                    shifts_conds{person_i}=[mean_HRF{person_i}{1,masks_loop}{2,:}]';
                    BOLD{person_i}(masks_loop).shifts{1,conds_loop}=shifts_conds{person_i};
                    BOLD{person_i}(masks_loop).BOLD{1,conds_loop}=BOLD_array_means_conds{person_i}(:,conds_loop);
                end
            end
            all_BOLD_means{person_i}{1,masks_loop}=BOLD_array_means_conds{person_i};
            shifts_all{person_i}=vertcat(shifts_all{person_i}(:,1),shifts_conds{person_i}(:,1));
           
            %alternative method to marsbar
        else
            % % %             sessions=size(SPM.Sess,1);
            % % %             nscans=SPM.nscan;
            % % %             for masks_loop=1:length(ROIs)
            % % %                 spm_get_data(SPM.xY.VY,XYZ)
            % % %                 target_mask=[maskdir ROIs{masks_loop}];
            % % %                 current_mask=load_nii(target_mask);
            % % %                 mask_voxels=current_mask.img;
            % % %                 mask_ROI=find(mask_voxels~=0);
            % % %                 BOLD_array=cell(1,sessions+1);
            % % %                 BOLD_array_means=cell(1,sessions+1);
            % % %
            % % %                 for sess_loop=1:sessions
            % % %                     BOLD_array{sess_loop+1}=zeros(nscans(1,sess_loop),size(mask_ROI,1));
            % % %                     BOLD_array_means{sess_loop+1}=zeros(nscans(1,sess_loop),1);
            % % %                     for scan_loop=1:nscans(1,sess_loop)
            % % %                         ss=SPM.xY.VY(scan_loop).private.dat.fname;
            % % %                         ss2=load_nii(ss);
            % % %                         img=ss2.img;
            % % %                         for voxel_loop= 1:size(mask_ROI,1)
            % % %                             BOLD_array{1}(1,end+1)=mask_ROI(voxel_loop);
            % % %                             BOLD_array_means{1}=mask_ROI(voxel_loop);
            % % %                             BOLD_array{sess_loop+1}(scan_loop,voxel_loop)=img(mask_ROI(voxel_loop));
            % % %                         end
            % % %                         BOLD_array_means{sess_loop+1}(scan_loop,1)=mean(BOLD_array{sess_loop+1}(scan_loop,:));
            % % %                     end
            % % %                 end
            % % %                 mask_BOLD_means= vertcat(BOLD_array_means{1,1:sess_loop+1});
            % % %             end
            
            
        end
        Subjects_BOLD{1,person_i}=BOLD{person_i};
%     end
end
%Histogram all the BOLD shifts
Subjects_shifts=zeros(0,1);
pool_sub.BOLD=zeros(0,0);
for i= 1:size(subjects,2)
    for i2= 1:size(marsbar_ROIs,2)
        for i3= 1:size(conds_id,2)
            if ~isempty(Subjects_BOLD{1,i}(i2).shifts{size(conds_id,2)})
                Subjects_shifts(end+1)=Subjects_BOLD{1,i}(i2).shifts{size(conds_id,2)}(i3,1);
                pool_sub.BOLD=vertcat(pool_sub.BOLD,Subjects_BOLD{1,i}(i2).BOLD{1,1});
            end
        end
    end
end
shifts=Subjects_shifts';
pooledBOLD=pool_sub;
conditions_lookup=Conditions_lookup;

end
