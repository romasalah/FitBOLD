%set your ROIs Directory
close all
ROI_dir='/Users/NEMO_admin/Desktop/Omar/Masks';
cd(ROI_dir);
clear all_corr
sub=[ 44 62 63 65 69 70 71 74 77 79 81 83 84 90 92 93 94 95 96 98 101 102 103 104 106 110 120 124 125 126 129 136 155 166 171 173 174 175 176 177 ];
files=dir(ROI_dir);
all_corr=cell(0,6);
data_dir={'/Volumes/MYBOOK1/FFX_Rutledge_minEV_14aug2018_ISIfinal_Ch2/Subjects',...
    '/Volumes/MYBOOK1/second level completed/FFX_Rutledge_minEV_31july2018_ISIfinal/Subjects',...
    '/Volumes/MYBOOK1/second level completed/FFX_Rutledge_minEV_31july2018_ISIfinal/Subjects',...
    '/Volumes/MYBOOK1/second level completed/FFX_Rutledge_minEV_27july2018_ISIfinal/Subjects',...
    '/Volumes/MYBOOK1/second level completed/FFX_Rutledge_noCompRegs_20june2018_ISIfinal/Subjects',...
    '/Volumes/SSD/FFX_Rutledge_minEV_6aug2018_ISIfinal_Ch2/Subjects'};
data_dir={'/Volumes/MYBOOK1/FFX_Rutledge_minEV_14aug2018_ISIfinal_Ch2/Subjects'};
condition={'happiness',...
    'PE Neg active','PE Neg passive','PE Neg non-soc'...
    'PE Pos active','PE Pos passive','PE Pos non-soc',...
    'Outcome self>other','Outcome self<other',...
    'Agl active','Agl passive','Agl non-soc'...
    'Min EV Neg active','Min EV Neg passive','Min EV Neg non-soc'...
    'Ch_diff active','Ch_diff passive','Ch_diff non-soc'...
    'CR active','CR passive','CR non-soc',...
    'EV active','EV passive','EV non-soc'};
condition={'happiness','Min EV'};
for dir_loop= 1:size(data_dir,2)
    d_dir=data_dir{dir_loop};
    newglm=['\n \n Testing a new GLM loacted at \n' d_dir];
    fprintf(newglm)
    for con_loop= 1: size(condition,2)
        con=condition{con_loop};
        newcon=['\n \n Starting a new condition \n \n' con];
        fprintf(newcon)
        for ROI_loop= 1: size(files,1)
            if ~isempty(strfind(files(ROI_loop).name,'.mat'))
                try
                    [x,CRB,y]=corr_BOLD(files(ROI_loop).name,con,con,d_dir,sub,'standard shift',1);
                    all_corr{end+1,1}=CRB;all_corr{end,2}=y; all_corr{end,3}=x; all_corr{end,4}=files(ROI_loop).name; all_corr{end,5}=d_dir; all_corr{end,6}=con;
                end
            end
        end
    end
end
cd(ROI_dir);
save('Big_correlation_results.mat','all_corr')
