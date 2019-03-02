function names=erase_fitBOLD(files,type,mtt,conds)
names=files;
for i=1:size(names,1)
    names(i).name=erase(names(i).name,mtt);
    names(i).name=erase(names(i).name,'Zscr_');
    names(i).name=erase(names(i).name,'stepwise_');
    names(i).name=erase(names(i).name,'__');
    names(i).name=erase(names(i).name,type);
    for i2= 1:size(conds,2)
    initloc=strfind(names(i).name,conds{1,i2});
    if ~isempty(initloc)
    names(i).name=names(i).name(initloc(1,1):end);
    end
    end
end
end