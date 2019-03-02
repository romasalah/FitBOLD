function [subROI,subROIz,sub_coord,silval]=Cluster_ROIs(voxels,k,coord,distance,silh,p,parallelFlag,fit)
if k==0 %uses n as cluster size in voxels. useful for clustering small ROIs
    n=round(size(voxels,2)./10);
    if n>fit.cluster.maxk; n=fit.cluster.maxk; end
else n=k;
end
subROI=cell(1,n);
subROIz=cell(1,n);
sub_coord=cell(1,n);
if n>1
for u=1:size(voxels,2)
    zvoxels(:,u)=zscore(voxels(:,u));
end
v=zvoxels';
if parallelFlag==1
    stream = RandStream('mlfg6331_64');  % Random number stream
    options = statset('UseParallel',1,'UseSubstreams',1,...
        'Streams',stream);
else options='default';
end
cidx=kmeans(v,n,'Options',options,'Replicate',8,'distance',distance);
if silh==1
   silval=silhouette(v,cidx,distance);
end

for i=1:n
subROI{1,i}=voxels(:,cidx==i);
subROIz{1,i}=zvoxels(:,cidx==i);
sub_coord{1,i}=coord(:,cidx==i);
if p; figure;plot(subROIz{1,i}); end
end
else 
    subROI{1,1}=mean(voxels,2);
    subROIz=zscore(subROI{1,1});
    sub_coord{1,1}=coord;
end
end