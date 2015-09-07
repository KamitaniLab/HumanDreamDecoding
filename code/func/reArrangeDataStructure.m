function [D_new]=reArrangeDataStructure(D,params,fieldnames)

% select voxels within the specified ROI
D_new.data=D.data(:,D.metaData.(params.useRois)==1);

% get label and label definition
labelIdx=find(~cellfun(@isempty,strfind(fieldnames,'Synset_')));
for j=1:length(labelIdx)
    D_new.labels(:,j)=D.data(:,D.metaData.(fieldnames{labelIdx(j)})==1);
end

% get xyz coordinate
D_new.labels_type=fieldnames(labelIdx)';
xyz=[D.metaData.X;D.metaData.Y;D.metaData.Z];
D_new.xyz=xyz(:,D.metaData.(params.useRois)==1);

end

