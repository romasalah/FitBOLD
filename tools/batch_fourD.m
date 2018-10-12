subjects = [44 62 63 65 69 70 71 74 77 79 81 83 84 90 92 93 94 95 96 98 101 102 103 104 106 110 120 124 125 126 129 136 155 166 171 173 174 175 176 177];
prefix='SODEC_FMRI_';
datadir='/Volumes/SSD/FFX_Rutledge_minEV_14aug2018_ISIfinal_Ch2/Subjects';
for add_files= 1 : length( Subjects ) 
    callperson{add_files} = [prefix num2str(Subjects(add_files))]; %calls a specific folder according to index
    personName{add_files} = callperson{add_files}; %get the directory as a character array
    matdir{add_files}=[ datadir '/' personName{add_files} '/SPM.mat']; %specifies the required matrix to open
    data{add_files}=load(matdir{add_files});
end
for person=1:length(subjects)
    sp=strfind(data{1,person}.SPM.xY.VY(end).fname,filesep);
    dest=data{1,person}.SPM.xY.VY(1).fname(1:sp(1,end));
    dest_4D{person,1}=fourD({data{1,person}.SPM.xY.VY(:).fname}');
end
