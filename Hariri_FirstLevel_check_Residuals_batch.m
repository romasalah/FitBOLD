%-----------------------------------------------------------------------
% Job saved on 29-Mar-2018 15:45:30 by cfg_util (rev $Rev: 6460 $)
% spm SPM - SPM12 (6470)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------

mainFolder = '/Volumes/MYBOOK1/FFX_Rutledge_minEV_30july2018_ISIfinal';
%E:\LOXY_MRI\FirstLevel_Hariri_mitsmoothing_4x4x4_mitResponseButton';
%E:\LOXY_MRI\FirstLevel_Hariri_mitsmoothing_4x4x4
%E:\LOXY_MRI\FirstLevel_Hariri_mitsmoothing_2x2x2_changed_implizitMask'
%E:\LOXY_MRI\FirstLevel_mitsmoothing_2x2x2_mitResponseButton;
%'E:\LOXY_MRI\FirstLevel_Hariri_mitsmoothing_2x2x2';
%
codeDir = '/Volumes/MYBOOK1/FFX_Rutledge_minEV_30july2018_ISIfinal/Model fits/'; %outputFile
anatMaskDir = '/Volumes/SSD/FFX_Rutledge_noCompRegs_19dec_ISIfinal/masks/';
anatMask{1} = [anatMaskDir 'Whole brain mask (grey).nii'];

persons{1} = [ 44 62 63 65 69 70 71 74 77 79 81 83 84 90 92 93 94 95 96 98 101 102 103 104 106 110 120 124 125 126 129 136 155 166 171 173 174 175 176 177 ];
    %1 2 4 5 6 7 8 9 10 12 14 15 16 17 18 22 23 24 26 27 28 29 30 32 33 34 35 36 38 39 41 42 45 46 47 48 51 52 54 57 58 59 60 61 65 67 68 70 71 72 73 74 75 76 ];
%101 102 103 105 106 107 108 110 112 114 115 116 117 118 120 121 122 123 126 127 129 130 131 134 135 136 137 139 140 141 142 143 144 145 146 147 150 151 152 153 154 155 156 157 158 159 161 162 164 165 166 167 168 171 172 173 174 175 176 177 178 179 180


for person_i = 1 : length( persons{1} )

personNo = persons{1}(person_i);
personName = ['SODEC_FMRI_' num2str( personNo ) ];
fprintf('Doing subject: %s\n',personName)



FFXdirName = [mainFolder filesep personName filesep];   
ResMS = [FFXdirName 'ResMS.nii'];
maskImg = [FFXdirName 'mask.nii'];
clear matlabbatch
matlabbatch{1}.spm.util.imcalc.input = {maskImg; ResMS; anatMask{1}};
matlabbatch{1}.spm.util.imcalc.output = 'temp';
matlabbatch{1}.spm.util.imcalc.outdir = {FFXdirName};
matlabbatch{1}.spm.util.imcalc.expression = 'i1.*i2.*i3';
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
spm_jobman( 'run', matlabbatch );
img = spm_read_vols(spm_vol([mainFolder filesep personName filesep 'temp.nii']));
modelFit(personNo) = mean(img(find(img)));
end
%Warning: The images do not all have the same dimensions. - using 1st image.
%Warning: The images do not all have same orientation and/or voxel sizes. - using 1st image.