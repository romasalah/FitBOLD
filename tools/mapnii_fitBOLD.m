function allsubrois=mapnii_fitBOLD(drawFlag)
addpath('/Users/NEMO_admin/Downloads/xjview96/xjview')
clear fnames
c=cd;
files=dir([c filesep '*.nii']);
fnames={files(1:end).name};
masksfiles=dir([c filesep '*All*.mat']);
masks={masksfiles(1:end).name};
load('mm3.mat');
allsubrois={};
for i=1:size(masks,2)
    niitmp=load(masks{1,i});
    for subroi_i=1:size(niitmp.sub_coord,2)
        subname= strrep(masks{1,i},'All',[num2str(subroi_i) '_']);
        allsubrois{1,end+1}= subname;
        coordinmm= round(mean(niitmp.sub_coord{1,subroi_i},2));
        allsubrois{2,end} = vox2coordSPM(coordinmm,mm3);
    end
    clear niitmp coordinmm
end
if drawFlag
input='';
for i=1:size(fnames,2)
    input=[input ''' , ''' fnames{i}];
end
input=[input(5:end) ''''];
input=['xjview(' input ')'];
eval(input)
end
end

