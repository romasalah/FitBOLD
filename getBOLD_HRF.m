%get a matrix of required voxels' indices.
%add the option to extract these indices from a contrast image.
%This Script needs MarsBar tool for functionality
clear all
use_marsbar=1;
marsbar_dir='/Applications/spm12/toolbox/marsbar';
marsbar_ROIs={'Frontal Medial Orbital (minimum)_0_35_-20_roi.mat','Frontal_medial_Orbital_1_48_16_roi.mat','Frontal_medial_Orbital_3_16_-30_roi.mat'...
    ,'Negative RPE ALE Mask_-24_-4_6_roi.mat','Negative RPE ALE Mask_-30_14_64_roi.mat','Negative RPE ALE Mask_-30_28_6_roi.mat'...
    ,'Negative RPE ALE Mask_-42_-8_48_roi.mat','Negative RPE ALE Mask_4_2_60_roi.mat','Orbitofrontal_bilateral_AAl_2_42_-13_roi.mat'...
    ,'Positive RPE ALE mask_-16_0_70_roi.mat','Positive RPE ALE mask_-20_2_-22_roi.mat','Positive RPE ALE mask_-23_2_0_roi.mat'...
    ,'Positive RPE ALE mask_4_8_38_roi.mat','Positive RPE ALE mask_25_3_0_roi.mat','Positive RPE ALE mask_56_12_18_roi.mat'...
    ,'VENTRAL_PALLIDUM_-10_2_-8_roi.mat','VENTRAL_PALLIDUM_11_3_-8_roi.mat'};
maskdir='/Users/NEMO_admin/Desktop/Omar/Masks/';
% ROIs={'Ventral Pallidum Bilateral.nii','Pos RPE 160 70.nii'};
% onsetsfile= 'Rutledge_InitialISIs_final.txt';
% onsets=dlmread(onsetsfile);
conds_prefix='spm_spm:beta (0041) - Sn(1) ';
conds_id={'PE Pos active*bf','Min EV zero active*bf','Min EV Neg active*bf','Ch_diff active*bf','PE Neg active*bf','PE Pos active*bf','happiness*bf'};
conds_types={'','xVal^'};
conds_derv={'(1)','(2)','(3)'};
conds_betas=cell(size(conds_id,2),size(conds_derv,2)+1);
marsbar_beta=0;
marsbar_mean=0;

%************************ Start*******************
datadir='/Volumes/MYBOOK1/second level completed/FFX_Rutledge_minEV_31july2018_ISIfinal/Subjects';
cd(datadir)
prefix='SODEC_FMRI_';
% these conditions are added to get the rellevant scans for each condition
Subjects = [44 55 62 63 65 69 70 71 74 77 79 81 83 84 90 92 93 94 95 96 98 101 102 103 104 106 110 120 124 125 126 129 136 155 166 171 173 174 175 176 177];
for person_i = 1 : length( Subjects ) %loop through the folders
    cd(datadir)
    callperson = [prefix num2str(Subjects(person_i))]; %calls a specific folder according to index
    personName = callperson; %get the directory as a character array
    cd (personName) % go to the directory
    matdir=[ personName '/SPM.mat']; %specifies the required matrix to open
    person_col=person_i+1;
    if exist( matdir ) ==0 %#ok<EXIST>
        fprintf('Cant find SPM.mat file \n \n');
    else
        load(matdir);
        mean_HRF=cell(size(SPM.xBF.bf,1),size(conds_id,2));
        %Using Marsbar Library
        if use_marsbar
            addpath(marsbar_dir);
            all_BOLD_means=zeros(sum(SPM.nscan,2),size(marsbar_ROIs,2));
            %         perosn_onsets = find(initialISIandSetupValues(:,1)==Subjects(person_i));
            %         ISI = initialISIandSetupValues(idx,ses+1);
            % sets a condition to start the script only if the matrix exists
            for masks_loop=1:length(marsbar_ROIs)
                cd(maskdir)
                load(marsbar_ROIs{masks_loop})
                marsbar_BOLD=getdata(roi,SPM.xY.VY);
                %Get raw BOLD.
                for voxels_loop=1:size(marsbar_BOLD,1)
                    BOLD_array_means(voxels_loop,1)=mean(marsbar_BOLD(voxels_loop,:),2);
                end
                all_BOLD_means(1:voxels_loop,masks_loop)=BOLD_array_means(1:voxels_loop,1);
                %Get betas of HRF+derivatives.
                for conds_loop= 1:size(conds_id,2)
                    conds_betas{conds_loop,1}=[conds_prefix conds_id{conds_loop}];
                    for derv_loop= 1:size(conds_derv,2)
                        cond=[conds_prefix conds_id{conds_loop} conds_derv{derv_loop}];
                        for mat_loop= 1:size(SPM.Vbeta,2)
                            if strcmp(SPM.Vbeta(mat_loop).descrip,cond)==1
                                marsbar_beta=getdata(roi,SPM.Vbeta(mat_loop).private);
                                marsbar_mean=mean(marsbar_beta,2);
                            end
                        end
                        conds_betas{conds_loop,derv_loop}=marsbar_mean;
                        %Get HRF shape and Peak of each condition
                        mean_HRF{conds_loop}=SPM.xBF.bf*(conds_betas{conds_loop}');
                        find(mean_HRF==max(HRF)) % assumes that the peak is the best time to get the relevant BOLD from
                        correct_x=[0:1/16:12.8]
                    end
                end
            end 
            
            
            
            
            
        else
            sessions=size(SPM.Sess,1);
            nscans=SPM.nscan;
            for masks_loop=1:length(ROIs)
                spm_get_data(SPM.xY.VY,XYZ)
                target_mask=[maskdir ROIs{masks_loop}];
                current_mask=load_nii(target_mask);
                mask_voxels=current_mask.img;
                mask_ROI=find(mask_voxels~=0);
                BOLD_array=cell(1,sessions+1);
                BOLD_array_means=cell(1,sessions+1);
                
                for sess_loop=1:sessions
                    BOLD_array{sess_loop+1}=zeros(nscans(1,sess_loop),size(mask_ROI,1));
                    BOLD_array_means{sess_loop+1}=zeros(nscans(1,sess_loop),1);
                    for scan_loop=1:nscans(1,sess_loop)
                        ss=SPM.xY.VY(scan_loop).private.dat.fname;
                        ss2=load_nii(ss);
                        img=ss2.img;
                        for voxel_loop= 1:size(mask_ROI,1)
                            BOLD_array{1}(1,end+1)=mask_ROI(voxel_loop);
                            BOLD_array_means{1}=mask_ROI(voxel_loop);
                            BOLD_array{sess_loop+1}(scan_loop,voxel_loop)=img(mask_ROI(voxel_loop));
                        end
                        BOLD_array_means{sess_loop+1}(scan_loop,1)=mean(BOLD_array{sess_loop+1}(scan_loop,:));
                    end
                end
                mask_BOLD_means= vertcat(BOLD_array_means{1,1:sess_loop+1});
            end
        end
    end
end
% [fit_results]= fit_BOLD_Omar_1(RegionsBOLD,TargetBOLD,constant);
%assumes that each row of the cell array is a different subject