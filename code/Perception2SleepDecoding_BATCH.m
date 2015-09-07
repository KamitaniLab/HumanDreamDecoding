% function Perception2SleepDecoding_BATCH
clear all
% location of the downloaded data
p.datdir = '../data/preproc/';
% add path to classifier and functions
addpath('./libsvm-mat-3.0.1/');
addpath('./func/');

%% parameter settings
% file names
p.sbjNames={...
    'Subject1'
    'Subject2'
    'Subject3'
    };
p.sleepfiles={...
    'PreprocessedSleepDataSubject1.h5'
    'PreprocessedSleepDataSubject2.h5'
    'PreprocessedSleepDataSubject3.h5'
    };
p.perceptionfiles={...
    'PreprocessedPerceptionDataSubject1.h5'
    'PreprocessedPerceptionDataSubject2.h5'
    'PreprocessedPerceptionDataSubject3.h5'
    };
p.propertyfiles={...
    'propsSubject1.h5'
    'propsSubject2.h5'
    'propsSubject3.h5'
    };
% voxel selection
params.useRois = 'HVC'; % select ROIs from ['V1','V2','V3','LOC','FFA','PPA','LVC','HVC'];
params.num_comp=1000; % # of voxels
p.nSubjects=length(p.sleepfiles);

%% load fmri data
fprintf('Load fMRI data======\n'),tic
pDs=cell(p.nSubjects,1);
pDp=cell(p.nSubjects,1);
props=cell(p.nSubjects,1);
fprintf('Sleep data...\n')
for i=1:p.nSubjects
    pDs{i}=readHDF5AsStruct([p.datdir,p.sleepfiles{i}]);
end
fprintf('Perception data...\n')
for i=1:p.nSubjects
    pDp{i}=readHDF5AsStruct([p.datdir,p.perceptionfiles{i}]);
end
fprintf('Property data...\n')
for i=1:p.nSubjects
    props{i}=readHDF5AsStruct([p.datdir,p.propertyfiles{i}]);
    props{i}.fieldnamesSleep=fieldnames(pDs{i}.metaData);
    props{i}.fieldnamesPerception=fieldnames(pDp{i}.metaData);
end;toc

%% arrange fMRI data, labels, properties
% data structure is rearranged for classification
fprintf('Arrange data structure======\n'),tic
for i=1:p.nSubjects
    pDs{i}=reArrangeDataStructure(pDs{i},params,props{i}.fieldnamesSleep);
    pDp{i}=reArrangeDataStructure(pDp{i},params,props{i}.fieldnamesPerception);
end;toc

%% perform decoding analysis
fprintf('Perform pairwise decoding======\n')
accCell=cell(p.nSubjects,1);
for i = 1:p.nSubjects;tic
    accCell{i}=zeros(1,length(props{i}.synsetPairs));
    for j = 1:length(props{i}.synsetPairs); % perform decoding for all pairs
        % specify pairs
        params.class1 = find(ismember(pDs{i}.labels_type,props{i}.synsetPairs{j,1}));
        params.class2 = find(ismember(pDs{i}.labels_type,props{i}.synsetPairs{j,2}));

        % perform decoding
        results = performPerception2SleepDecoding(pDs{i},pDp{i},params);
        
        % reserve accuracy in matrix
        accCell{i}(j) = results.corrRate;
    end;toc
end

%% check results 
fprintf('Visualize results======\n')
% get tested pairs' accuracy
accuracies=[accCell{:}];

% visualize results
fprintf('Mean accuracy(%s)======\n',params.useRois)
for i = 1:p.nSubjects
    subplot(2,2,i);hist(accCell{i},0:5:100), axis square, xlim([0,100]), xlabel('Accuracy'),ylabel('Frequency'),title(sprintf('Subject%d: %s: %d pairs',i,params.useRois,length(accCell{i})))
    fprintf('%s = %.1f%%\n',p.sbjNames{i},mean(accCell{i}))
end
subplot(224);hist(accuracies,0:5:100), axis square, xlim([0,100]), xlabel('Accuracy'),ylabel('Frequency'),title(sprintf('Pooled: %s: %d pairs',params.useRois,length(accuracies)))
fprintf('Pooled   = %.1f%%\n',mean(accuracies))

%%
