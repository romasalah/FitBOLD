function img=fourD(threeD,name)
cont=1;
matlabbatch{1}.spm.util.cat.vols = threeD;
if nargin<2
    name='4D.nii';
end
sp=strfind(threeD{1,end},filesep);
dest=threeD{1,end}(1:sp(1,end));
img=[dest name];
if exist(img,'file')==2
    prompt='There is a 4D volume for the this Subject, Overwrite it? 1=Yes 0=No \n';
    resp=input(prompt,'s');
    if isempty(resp) || resp==1
        cont=1;
    else
    cont=0;
    end
end
if cont==1
matlabbatch{1}.spm.util.cat.name = name;
matlabbatch{1}.spm.util.cat.dtype = 0;
dir=cd;
jobfile = [ dir '/4D_job.mat'];
save(jobfile,'matlabbatch')
jobs = repmat(jobfile, 1, 1);
inputs = cell(0, 1);
spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:});
else
    return
end
delete(jobfile)
end