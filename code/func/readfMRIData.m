function [iDs,iDp,props]=readfMRIData(p)
%
% integrate fMRI data for the specified subject
%

%% set parameters
fMRIdataGroupIdx=1;
labelGroupIdx4SleepData=[2:7];
labelGroupIdx4PerceptionData=[2:4];

switch p.subjectID
    case 1
        nSleepDataSession=26;
        nPerceptionDataSession=32;
    case 2
        nSleepDataSession=14;
        nPerceptionDataSession=40;
    case 3
        nSleepDataSession=15;
        nPerceptionDataSession=40;
end
%% read data
fprintf('Read files...\n')
% sleep data
Ds=cell(nSleepDataSession,1);
for i=1:nSleepDataSession
    fprintf('Sleep data loading(%d/%d)\n',i,nSleepDataSession)
    Ds{i}=readHDF5AsStruct([p.datdir,'SleepDataSubject',num2str(p.subjectID),'Session',num2str(i)]);
end
% perception data
Dp=cell(nPerceptionDataSession,1);
for i=1:nPerceptionDataSession
    fprintf('Perception data loading(%d/%d)\n',i,nPerceptionDataSession)
    Dp{i}=readHDF5AsStruct([p.datdir,'PerceptionDataSubject',num2str(p.subjectID),'Session',num2str(i)]);
end    


%% integrate fMRI data
fprintf('Integrate fMRI data...\n')
fMRIgroupName=['group',num2str(fMRIdataGroupIdx)];
thresVals=100;

% integrate data
data_tmp=shiftdim(Ds{1}.(fMRIgroupName).data,1);
tmp_mask=data_tmp(:,:,:,1)>thresVals;
tmp_mask=tmp_mask(:);
tmp_nVox=sum(tmp_mask);
clear data_tmp

iDs=integrateData(Ds,nSleepDataSession,fMRIgroupName,tmp_mask);
iDp=integrateData(Dp,nPerceptionDataSession,fMRIgroupName,tmp_mask);

% extract ROI mask
% the ROI masks are common across files [Ds(i) and Dp(i)]
roiNames=fieldnames(Ds{1}.(fMRIgroupName).props.roi);
nRois=length(roiNames);
roiMask=zeros(nRois,tmp_nVox);
for i=1:nRois
    tmp=Ds{1}.(fMRIgroupName).props.roi.(roiNames{i})(:);
    roiMask(i,:)=tmp(tmp_mask);
end
thresVox=any([iDs.data;iDp.data]>thresVals);
props.roiMask=roiMask(:,thresVox);
props.roiNames=roiNames;

% reserve non-zero voxel data
iDs.data=iDs.data(:,thresVox);
iDp.data=iDp.data(:,thresVox);

% reserve xyz coordinate
iDs.xyz=Ds{1}.(fMRIgroupName).props.location(:,tmp_mask);
iDp.xyz=Dp{1}.(fMRIgroupName).props.location(:,tmp_mask);
iDs.xyz=iDs.xyz(:,thresVox);
iDp.xyz=iDp.xyz(:,thresVox);

%% integrate labels
fprintf('Integrate labels...\n')
% sleep data
nlabelSleepData=length(labelGroupIdx4SleepData);
nlabelPerceptionData=length(labelGroupIdx4PerceptionData);

% sleep data
cnt=0;
iDs.labels=cell(1,nlabelSleepData);
for j=1:nlabelSleepData
    label_tmp=cell(1,nSleepDataSession);
    labelGroupName=['group',num2str(labelGroupIdx4SleepData(j))];
    for i=1:nSleepDataSession
        label_tmp{i}=Ds{i}.(labelGroupName).data';
    end
    iDs.labels{j}=[label_tmp{:}]';
    if strcmp(Ds{1}.(labelGroupName).props.title,'Base synset labels')
        for k=1:length(Ds{1}.(labelGroupName).props.labelNames)
            cnt=cnt+1;
            iDs.labels_type{cnt}=['Synset_',Ds{1}.(labelGroupName).props.labelNames{k}];
            iDs.labels_def{cnt}=Ds{1}.(labelGroupName).props.description;
        end
    else
        cnt=cnt+1;
        iDs.labels_type{cnt}=Ds{1}.(labelGroupName).props.title;
        iDs.labels_def{cnt}=Ds{1}.(labelGroupName).props.description;
    end
end
iDs.labels=[iDs.labels{:}];

% perception data
cnt=0;
iDp.labels=cell(1,nlabelPerceptionData);
for j=1:nlabelPerceptionData
    label_tmp=cell(1,nlabelPerceptionData);
    labelGroupName=['group',num2str(labelGroupIdx4PerceptionData(j))];
    for i=1:nPerceptionDataSession
        label_tmp{i}=Dp{i}.(labelGroupName).data';
    end
    iDp.labels{j}=[label_tmp{:}]';
    if strcmp(Dp{1}.(labelGroupName).props.title,'Base synset labels')
        for k=1:length(Dp{1}.(labelGroupName).props.labelNames)
            cnt=cnt+1;
            iDp.labels_type{cnt}=['Synset_',Dp{1}.(labelGroupName).props.labelNames{k}];
            iDp.labels_def{cnt}=Dp{1}.(labelGroupName).props.description;
        end
    else
        cnt=cnt+1;
        iDp.labels_type{cnt}=Dp{1}.(labelGroupName).props.title;
        iDp.labels_def{cnt}=Dp{1}.(labelGroupName).props.description;
    end
end
iDp.labels=[iDp.labels{:}];

end

%% functions
function iD=integrateData(D,nDataSession,fMRIgroupName,tmp_mask)
iD.data=cell(nDataSession,1);
for i=1:nDataSession
    data_tmp=shiftdim(D{i}.(fMRIgroupName).data,1);
    n=size(data_tmp,4);
    nVox=sum(tmp_mask);
    iD.data{i}=zeros(nVox,n);
    for j=1:n
        tmp=data_tmp(:,:,:,j);
        iD.data{i}(:,j)=tmp(tmp_mask);
    end
end
iD.data=[iD.data{:}]';
end











