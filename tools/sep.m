clear lowdata 
clear highdata
xx=zeros(0,0);
for t4=1:size(alldata,2)
for t2=1:120
xx(end+1,:)=alldata(t4).socdata{1,1}(t2,:);
end
end
lowdata=zeros(0,0);
highdata=zeros(0,0);
for t5= 1: size(alldata,2)
    gg=find(lowvec==t5);
    if ~isempty(gg)
        for t6=1:120
            lowdata(end+1,:)=alldata(t5).socdata{1,1}(t6,:);
        end
    elseif isempty(gg)
        for t7=1:120
            highdata(end+1,:)=alldata(t5).socdata{1,1}(t7,:);
        end
    end
end
for t8=1:size(lowdata,2)
figure
histfit(lowdata(:,t8),100)
figure
qqplot(lowdata(:,t8))
end
for t9=1:size(highdata,2)
figure
histfit(highdata(:,t9),100)
figure
qqplot(highdata(:,t9))
end
