function namessign=orthviews_fitBOLD(folder,type,imagespec,sign,drawF)
if folder==0
    folder=pwd;
end
cd(folder);
timingsref={'Min EV','PE','happiness'};
%cd(fit.dest_dir)
query= ['*' type '*results.mat'];
load('mm3.mat');
addpath('/Applications/spm12/toolbox/NIfTI')
files= dir(query);
for load_i=1: size(files,1)
    try all_results(load_i)=load(files(load_i).name); catch failed_loads; end
end
conds={'active','passive','n-soc'};
names=erase_fitBOLD(files,type,'results.mat',conds);
numres=size(all_results,2)+1;
%spm_orthviews('reset')
for template_i=1:numres
    template(template_i)=spm_vol('/Volumes/EVO 860 2TB/fitBOLD/templates/single_subj_T1w_250um.nii');
end
switch imagespec
    case 'linear'
        imagespec=1;
        powerlbl='linearily';
    case 'quad'
        powerlbl='quadratically';
        imagespec=2;
end
mn = numres;
n  = round(mn^0.4);
m  = ceil(mn/n);
w  = 1/n;
h  = 1/m;
ds = (w+h)*0.02;
clear global st
for results_i=1:numres
    
    %% load and position a new template
    
    i = 1-h*(floor((results_i-1)/n)+1);
    j = w*rem(results_i-1,n);
    if drawF
        spm_orthviews('image',template(results_i),[(j+ds/2)/1.5 i+ds/2 w-ds h-ds]);
        spm_orthviews('mapping', 'histeq');
        hold on;
    end
    if results_i<numres
        current=(all_results(results_i).results);
        coords= current.sub_coord;
        allcoords=[];subcoords=[];
        ints=[];
        %% get initial ROIs
        for unpackv_i=1:size(coords,2)
            subsizes(unpackv_i)=size(coords{1,unpackv_i},2);
            subcoords=horzcat(subcoords,coords{1,unpackv_i});
        end
        subcoordsquad=horzcat(subcoords,subcoords);
        subcoordsquad_in={subcoordsquad{1,(sum(current.betas,1)~=0)}};
        maxsubsizes=max(subsizes);
        %% get ROIs to a specific power
        subcoordspowered={subcoordsquad_in{1,find(current.linearityidx(imagespec,:)==1)}};
        weightspowered=  current.coefci{1,find(current.linearityidx(imagespec,:)==1)};
        namespowered={current.predictorROIs{1,find(current.linearityidx(imagespec,:)==1)}};
        for get_timing=1:size(namespowered,2);
            timings{get_timing}=namespowered{1,get_timing}(strfind(namespowered{1,get_timing},'at')+3:end);
            timings{get_timing}=strrep(timings{get_timing},'_',' ');
        end
        for t_idx=1:size(timingsref,2)
            timingsidx{t_idx}=strcmp(timings,timingsref{t_idx});
        end
        for unpacksubrois=1:size(subcoordspowered,2)
            allcoords=horzcat(allcoords,subcoordspowered{1,unpacksubrois});
            ints=vertcat(ints,zeros(size(subcoordspowered{1,unpacksubrois},2),1)+weightspowered(unpacksubrois));
        end
        %% get signed intensities
        namessign(results_i).model=namespowered;
        switch sign
            case 'pos'
                signlbl='Positive';
                signev=1;
                intsign=ints(ints>0);
                coordsign=allcoords(:,find(ints>0));
                namessign(results_i).ROIs={namespowered{1,find(weightspowered>0)}};
                
            case 'neg'
                signlbl='Negative';
                signev=-1;
                intsign=ints(ints<0);
                coordsign=allcoords(:,find(ints<0));
                namessign(results_i).ROIs={namespowered{1,find(weightspowered<0)}};
            otherwise
                intsign=ints;
                coordsign=allcoords;
                signev=0;
                signlbl='Positive and Negative';
        end
        %% Divide by timings
        %% draw blobs
        %maxroi=coordsign(:,abs(ints)==max(abs(intsign)));
        %maxvox=maxroi(:,1);
        if drawF
            blb_lbl=char(names(results_i).name);
            spm_orthviews('addBlobs',results_i,coordsign,intsign,mm3,blb_lbl);
            title(char(names(results_i).name))
            spm_orthviews('Reposition',[0 0 0])
%             if results_i>1
%             spm_orthviews('Xhairs','off')
%             end
            spm_orthviews('hld','Sinc')
        end
    else
        if drawF
            x=-10:0.01:10;
            y= signev*x.^imagespec;
            subplot('Position',[(j+20*ds/2)/1.5 i+3*ds/2 w-ds*12 h-ds*5])
            plot(x,y);
            title([powerlbl ' behaving ROIs, ' signlbl ' betas'], 'interpreter','none')
        end
    end
end
if drawF
lastidx=strfind(files(load_i).name,'shft_');
figname=char(files(load_i).name(1:lastidx+3));
set(gcf,'Name',figname)
results_i=results_i+1;
set(0,'units','pixels')
screenres = get(0,'screensize');
pos=get(gcf,'position');
if pos(1,4)<1000
    ndim=2.5*pos(1,[3,4]);
    set(gcf,'position',[pos(1,[1,2]),ndim(1,[1,2])])
end
savefig(figname)
end
end