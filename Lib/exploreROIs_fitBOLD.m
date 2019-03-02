function [resultROIs,constr,restROIs]=exploreROIs_fitBOLD(query)
warning off
drawF=0;
ROIs=get_ROIs(drawF);
conds={'active','passive','n-soc'};
models={'_NSoc','_ToM','_Soc'};

for id=1:size(models,2)
    models_i{id}=~cellfun(@isempty,strfind(ROIs.Properties.RowNames,models{id}));
    models_i_lbl=models;
end
for con_i=1:size(conds,2)
    conds_i{con_i}=~cellfun(@isempty,strfind(ROIs.Properties.RowNames,conds{con_i}));
    conds_i_lbl=conds;
end
signs_i{1}=~cellfun(@isempty,strfind(ROIs.Properties.VariableNames,'Pos'));
signs_i{2}=~cellfun(@isempty,strfind(ROIs.Properties.VariableNames,'Neg'));
lin_i{1}=~cellfun(@isempty,strfind(ROIs.Properties.VariableNames,'Linear'));
lin_i{2}=~cellfun(@isempty,strfind(ROIs.Properties.VariableNames,'Quadratic'));
lin_i_lbl={'Linear','Quadratic'};
signs_i_lbl={'positive','negative'};
%% Set all possible constraints
cons=1:4; %Constraints
consarray{1}=signs_i;
consarray{2}=lin_i;
consarray{3}=conds_i;
consarray{4}=models_i;
sizes={size(signs_i,2), size(lin_i,2), size(conds,2), size(models,2)};
pos={'col','col','row','row'};
lbl={'signs_i','lin_i','conds_i','models_i'};
consarray(2,cons)=sizes;
consarray(3,cons)=pos;
consarray(4,cons)=lbl;

consperms=perms(cons);
for perm_i=1:size(consperms,1)
    consarrayperm{perm_i}=consarray(1:end,consperms(perm_i,:));
    elements{perm_i} = {1:consarrayperm{perm_i}{2,1}, 1:consarrayperm{perm_i}{2,2}, 1:consarrayperm{perm_i}{2,3}, 1:consarrayperm{perm_i}{2,4}}; %cell array with N vectors to combine
    for const_i=1:size(elements{perm_i},2)
        combinations = cell(1, const_i); %set up the varargout result
        [combinations{:}] = ndgrid(elements{perm_i}{1:const_i});
        combinations = cellfun(@(x) x(:), combinations,'uniformoutput',false); %there may be a better way to do this
        constraints_i{perm_i}{1,const_i} = [combinations{:}]; % NumberOfCombinations by N matrix. Each row is unique.
        constraints_i{perm_i}{3,const_i}=consarrayperm{perm_i}(3,1:const_i);
        constraints_i{perm_i}{4,const_i}=consarrayperm{perm_i}(4,1:const_i);
    end
end
cases={}; allconst_i={};allpos={};allconst2={};
for unpack_const=1:size(constraints_i,2)
    allconst_i(unpack_const,cons)=constraints_i{1,unpack_const}(1,cons);
    for const_size=cons
        allpos{unpack_const,const_size}=constraints_i{1,unpack_const}{3,const_size};
        alllbl{unpack_const,const_size}=constraints_i{1,unpack_const}{4,const_size};
        for get_const2=1:const_size
            cases{unpack_const,const_size}{get_const2}=constraints_i{unpack_const}{1,const_size}(get_const2,:);
######YOua re not getting all the cases            
whichcase=cases{unpack_const,const_size}{get_const2};
            tmp={consarrayperm{unpack_const}{1,1:const_size}};
            for whichcase_i=1:const_size
                allconst2{unpack_const,const_size}{get_const2}{whichcase_i}=tmp{whichcase_i}{whichcase(whichcase_i)};
            end
        end
    end
end
idx=allconst2;
%% Find ROIs which satisfy the constraints
for const_perm=1:size(idx,1)
    for const_num=1:size(idx,2)
        for const_set_i=1:const_num
            critlbl=cell(1,const_num);
            joinedcritlbl='';
            rows=[];cols=[];
            for const_i =1:const_set_i
                
                curconst=idx{const_perm,const_num}{const_set_i}{const_i};
                curpos=allpos{const_perm,const_num}{const_i};
                curlbl=alllbl{const_perm,const_num}{const_i};
                critlbl{const_i}=eval([curlbl '_lbl{cases{unpack_const,const_size}{get_const2}(const_i)}']);
                if strcmp(curpos,'row')
                    rows=[rows; curconst']; cols=[cols; ones(1,size(ROIs,2))];
                else
                    cols=[cols; curconst]; rows=[rows; ones(1,size(ROIs,1))];
                end
                joinedcritlbl=[joinedcritlbl ' & ' char(critlbl{1,const_i})];      
            end
            rowsloc=find(strcmp({allpos{const_perm,const_num}{1:const_set_i}},'row')==1);
            colsloc=find(strcmp({allpos{const_perm,const_num}{1:const_set_i}},'col')==1);
            if length(rowsloc)==1
                rows_set=find(rows(rowsloc(1),:)==1);
            elseif isempty(rowsloc)
                rows_set=1:size(ROIs,1);
            else
            rows_set=find(rows(rowsloc(1),:)==(rows(rowsloc(2),:)) & rows(rowsloc(1),:)~=0 );
            end
            if length(colsloc)==1
                cols_set=find(cols(colsloc(1),:)==1);
            elseif isempty(colsloc)
                 cols_set=1:size(ROIs,2);
            else
            cols_set=find(cols(colsloc(1),:)==(cols(colsloc(2),:))& cols(colsloc(1),:)~=0);
            end
            restROIs{const_perm,const_num}{const_set_i}=ROIs{rows_set,cols_set};
            constr{const_perm,const_num}{const_set_i}=joinedcritlbl;
        end
    end
end
%% get queired constraint
nr=size(constr,1);
nc=size(constr,2);

for q_const=1:size(query,2)
            if iscell(query{q_const}) && size(query{q_const},2)>1
                expr='1';
                for i=1:size(query{q_const},2)
                    expr= strcat(expr,[' & contains(constr{lbl1,lbl2},query{q_const}{' num2str(i) '})']);
                end
            else
                expr='contains(constr{lbl1,lbl2},query{q_const}{1})';
            end

    ROIs_allperm{q_const}={};
    for lbl1=1:nr
        for lbl2=1:nc
            ROIs_oneperm={};
fprintf(num2str(sum(eval(expr))>0))
             if sum(eval(expr))>0
                 tgtsize=find(eval(expr)==1);
                 roispack=restROIs{lbl1,lbl2}{1,tgtsize};
                 for unp_i1=1:size(roispack,1)
                     for unp_i2=1:size(roispack,2)
                     ROIs_oneperm=horzcat(ROIs_oneperm,roispack{unp_i1,unp_i2});
                     end
                 end
             end
        end
        ROIs_allperm{q_const}= horzcat(ROIs_allperm{q_const},ROIs_oneperm);
        resultROIs{q_const}=unique(ROIs_allperm{q_const});
    end
end
end
