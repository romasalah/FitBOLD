function startpool_fitBOLD(workers)
profiles={'fullpool1','fullpool2'};

if exist('poolobj_fitBOLD','var')
    for profile_i=1:size(profiles,2)
        if ~strcmp(poolobj_fitBOLD,profiles{profile_i})
            global poolobj_fitBOLD
           poolobj_fitBOLD=parpool(profiles{profile_i},workers);
        end
    end
end
end