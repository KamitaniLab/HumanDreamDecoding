% function convertRawDatatoPreprocessedData
%
% [note]
%   Run this script first to create and save preprocessed data.
%   it takes about 30 minutes to run this script.
%
% [main variables]
%   Ds                           -  raw data structure of sleep experiment (all voxels with non-zero values)
%   Ds,data                      -  2D matrix of fMRI data ([time(sample) x space(voxel)] format)
%   Ds.xyz                       -  xyz coordinate of each voxels the individual brain [3 x voxel]
%   Ds.labels                    -  condition labels of each sample ([time x 1] format)
%   Ds.labels_type               -  label names. each corresponds to each column of .labels field
%   Ds.labels_def                -  definition of label value in .labels field
%
%   Dp                           -  raw data structure of perception experiment (all voxels with non-zero values)
%   Dp,data                      -  2D matrix of fMRI data ([time(sample) x space(voxel)] format)
%   Dp.xyz                       -  xyz coordinate of each voxels the individual brain [3 x voxel]
%   Dp.labels                    -  condition labels of each sample ([time x 1] format)
%   Dp.labels_type               -  label names. each corresponds to each column of .labels field
%   Dp.labels_def                -  definition of label value in .labels field
%
%   props                        -  property structure
%   props.synsetPairs            -  synsetPairs list for pairwise decoding analysis
%
%   pDs                          -  Preprocessed data structure of sleep experiment (only voxels within pre-defined ROIs)
%   pDs.data                     -  2D matrix of [fMRIdata, labels, design]
%   pDs.metaData                 -  label index which is used to extract corresponding
%                                   data from the .data field
%   pDs.metaData.voxelData       -  index for voxelData-columns in .data field 
%   pDs.metaData.(X,Y,Z)         -  xyz coordinates of corresponding voxeldata
%   pDs.metaData.(ROI:e.g., V1)  -  index for voxels within each ROI
%   pDs.metaData.label           -  index for label-columns in .data field 
%   pDs.metaData.EEG_sleep_score -  index for eeg score label in .data field
%   pDs.metaData.Synset_XXX      -  index for each object(synset) category in .data field
%   pDs.metaData.design          -  index for design-columns in .data field 
%   pDs.metaData.session_number  -  index for session number of fMRI data in .data field
%   pDs.metaData.sample_number   -  index for sample number of fMRI data in .data field
%   pDs.metaDefinition           -  definition of label and design values
% 
%   pDp                          -  Preprocessed data structure of perception experiment (only voxels within pre-defined ROIs)
%   pDp.data                     -  2D matrix of [fMRIdata, labels, design]
%   pDp.metaData                 -  label index which is used to extract corresponding
%                                   data from the .data field
%   pDp.metaData.voxelData       -  index for voxelData-columns in .data field 
%   pDp.metaData.(ROI:e.g., V1)  -  index for voxels within each ROI
%   pDp.metaData.(X,Y,Z)         -  xyz coordinates of corresponding voxeldata
%   pDp.metaData.label           -  index for label-columns in .data field 
%   pDp.metaData.Synset_XXX      -  index for each object(synset) category in .data field
%   pDp.metaData.design          -  index for design-columns in .data field 
%   pDp.metaData.session_number  -  index for session number of fMRI data in .data field
%   pDp.metaData.sample_number   -  index for sample number of fMRI data in .data field
%   pDp.metaDefinition           -  definition of label and design values
%
%   [Example]
%    % To get fMRIdata, labels, or design in pDp structure
%      voxelData = pDp.data(:,pDp.metaData.voxelData==1);
%      labels    = pDp.data(:,pDp.metaData.label==1);
%      design    = pDp.data(:,pDp.metaData.design==1);
%
%
%
% [Note]
%  In the preprocessing stage the following procedures were applied
%   -voxel selection within pre-defined ROIs
%   -shift perception data for hemodynamic delay
%   -reduce outliers
%   -linear detrending
%   -feature normalization
%   -averaging to create samples
%   -select synset pairs for pairwise decoding analysis
%   -exclude awake samples (EEG sleep score = 0) from sleep data
%   -exclude irrevant labels
%   -match label order of sleep and perception data
%
% 
%   Written by Tomoyasu Horikawa horikawa-t@atr.jp
%
%clear all
%% initial settings
% CHANGE PATH
cd('/home/mu/project/open_data/science2013/code/')
% add path to functions
addpath('./func/');
% location of the donwloaded data
p.datdir = '../data/raw/';
% location of the save directory
p.savdir = setdir('../data/preproc/');
% # of subject
nSubjects = 3;

%% Read & preprocess fmri data: create decoding sample (it takes about an hour)
for i = 1:nSubjects
    p.subjectID=i;

    fprintf('Read data for Subject%d======\n',i);tic
    [Ds,Dp,props_org] = readfMRIData(p);toc
    % if you want whole brain data, save Ds and Dp
    
    fprintf('Preprocess data======\n');tic
    [pDs,pDp,props] = preprocfMRIData(Ds,Dp,props_org);toc
    
    fprintf('Save preprocessed data======\n')
    writeHDF5FromStruct([p.savdir,'PreprocessedSleepDataSubject',num2str(i),'.h5'],pDs);
    writeHDF5FromStruct([p.savdir,'PreprocessedPerceptionDataSubject',num2str(i),'.h5'],pDp);
    writeHDF5FromStruct([p.savdir,'propsSubject',num2str(i),'.h5'],props);
end

%%
