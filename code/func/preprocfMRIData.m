function [pDs,pDp,props_new]=preprocfMRIData_test(Ds,Dp,props)
%
% preprocess fMRI data to create features for decoding analysis
%

%% voxel selection within pre-defined ROIs
fprintf('Voxel selection...\n')
useVox=any(props.roiMask,1);
props.roiMask=props.roiMask(:,useVox);
Ds.data=Ds.data(:,useVox);
Dp.data=Dp.data(:,useVox);
Ds.xyz=Ds.xyz(:,useVox);
Dp.xyz=Dp.xyz(:,useVox);
props.xyz=Ds.xyz;

%% shift data for hemodynamic delay
fprintf('Shift data...\n')
params.nshift= 1;% [volume] default = 1 = 3 sec

% sleep data is not shifted
% perception data
[Dp,params]=shiftData(Dp,params);

%% reduce outliers
fprintf('Reduce outliers...\n')
params.std_thres=3;
params.num_its = 10;
params.remove  = 1;
params.method  = 2;
params.max_val = inf;
params.min_val = 100;

% sleep data
Ds=reduceOutliers(Ds,params);
% perception data
Dp=reduceOutliers(Dp,params);

%% feature matching
% this is necessary because reduceOutliers may change the # of features
[nouse,idx_s,idx_p]=intersect(Ds.xyz',Dp.xyz','rows');

% sleep data
Ds.data=Ds.data(:,idx_s);
Ds.xyz=Ds.xyz(:,idx_s);
% perception data
Dp.data=Dp.data(:,idx_p);
Dp.xyz=Dp.xyz(:,idx_p);

% roi mask data
[nouse,idx_r]=intersect(props.xyz',Ds.xyz','rows');
props.xyz=props.xyz(:,idx_r);
props.roiMask=props.roiMask(:,idx_r);

%% linear detrending
fprintf('Detrending...\n')
params.sub_mean  = 0;
params.method    = 'linear';

% sleep data
Ds=detrenddata(Ds,params);
% perception data
Dp=detrenddata(Dp,params);

%% feature normalization
fprintf('Normalize data...\n')
% sleep data
params.analysisPeriod=40; % [volume] default = 40 = 2 min before awakening for each period
params.normalizeBaselineDuration=10; %[volume] default = 10 = 30 sec duration for normalization
params.normalizationOnsetfromAwake=20; % [vlume] defualt = 20 = 1 min before awakening
params.zero_thres     = 1;

[Ds,params]=normalizeSleepData(Ds,params);

% perception data
params.base_labels    = 'all';
params.norm_mode      = 0;
params.zero_thres     = 1;

Dp=normalizeData(Dp,params);

%% average to create samples
fprintf('Average data...\n')
params.twShift=1;
params.winDur=3; % volume default = 3 = 9s before awakening

% sleep data
Ds=averageSleepData(Ds,params);

% perception data
Dp=averagePerceptionData(Dp);

%% select calculate pairs
nSampleThreshold=10;
synsetLabel=find(~cellfun(@isempty,strfind(Ds.labels_type,'Synset_')));
sum(Ds.labels(:,synsetLabel));
nSynset=length(synsetLabel);
nminSamp=zeros(nSynset);
for i=1:(nSynset-1)
    for j=(i+1):nSynset
        lab1=Ds.labels(:,synsetLabel(i))==1;
        lab2=Ds.labels(:,synsetLabel(j))==1;
        overlap=(lab1+lab2)==2;
        lab1(overlap)=0;
        lab2(overlap)=0;
        nminSamp(j,i)=min(sum([lab1,lab2]));
    end
end
[tri,i,j]=getTri(nminSamp);
usePairsIdx=(tri>=nSampleThreshold);

% rename synset names
synsetNames=Ds.labels_type(synsetLabel);
for ix=1:nSynset
    synsetNames{ix}=strrep(synsetNames{ix},'[','_');
    synsetNames{ix}=strrep(synsetNames{ix},':','_');
    synsetNames{ix}=synsetNames{ix}(1:end-1);
end

props.synsetPairs=[synsetNames(i(usePairsIdx));synsetNames(j(usePairsIdx))]';
props.synsetNames=synsetNames;

%% exclude awake samples from sleep data
awakeTrials=Ds.labels(:,ismember(Ds.labels_type,'EEG sleep score'))==0;
Ds.data(awakeTrials,:)=[];
Ds.labels(awakeTrials,:)=[];

%% create data matrix
% sleep data
synsetLabelIdx=find(~cellfun(@isempty,strfind(Ds.labels_type,'Synset_')));
eegScoreIdx=find(~cellfun(@isempty,strfind(Ds.labels_type,'EEG sleep score')));
sessionNumberIdx=find(~cellfun(@isempty,strfind(Ds.labels_type,'Session number')));
sleep_fmri=[...
    Ds.data,... % fMRI data
    ];
sleep_labels=[...
    Ds.labels(:,synsetLabelIdx),... % base synset labels (data)
    Ds.labels(:,eegScoreIdx),... % EEG sleep socre (labels)
    ];
sleep_design=[...
    Ds.labels(:,sessionNumberIdx),... % session number (design)
    (1:size(Ds.data,1))',... % sample number (design)
    ];
sleep_data=[sleep_fmri,sleep_labels,sleep_design];

% perception data
synsetLabelPerception=find(~cellfun(@isempty,strfind(Dp.labels_type,'Synset_')));
synsetNamesPerception=Dp.labels_type(synsetLabelPerception);
[ismem,synserMatchOrd]=ismember(Ds.labels_type(synsetLabel),synsetNamesPerception);
sessionNumberIdx=find(~cellfun(@isempty,strfind(Dp.labels_type,'Session number')));
perception_fmri=[...
    Dp.data,... % fMRI data
    ];
perception_labels=[...
    Dp.labels(:,synsetLabelPerception(synserMatchOrd)),... % base synset labels (data)
    ];
perception_design=[...
    Dp.labels(:,sessionNumberIdx),... % session number (design)
    (1:size(Dp.data,1))',... % sample number (design)
    ];
perception_data=[perception_fmri,perception_labels,perception_design];

%% create feature meta data
% sleep data
voxelData=[...
    ones(1,size(sleep_fmri,2));... % voxelData
    Ds.xyz;... % xyz coordinate
    props.roiMask... % roi mask
    ];
voxelData_definition={'voxelData','X','Y','Z',props.roiNames{:}};

labelData=[...
    ones(1,size(sleep_labels,2));... % label
    diag(ones(1,size(sleep_labels,2)));... % index for each label
    ];
label_definition={'label',props.synsetNames{:},'EEG_sleep_score'};

designData=[...
    ones(1,size(sleep_design,2));... % design
    diag(ones(1,size(sleep_design,2)));... % index for each design
    ];
design_definition={'design','session_number','sample_number'};

sleepMetaData=blkdiag(voxelData,labelData,designData);
sleepDefinition={voxelData_definition{:},label_definition{:},design_definition{:}};

% perception data
voxelData=[...
    ones(1,size(perception_fmri,2));... % voxelData
    Dp.xyz;... % xyz coordinate
    props.roiMask... % roi mask
    ];
voxelData_definition={'voxelData','X','Y','Z',props.roiNames{:}};

labelData=[...
    ones(1,size(perception_labels,2));... % label
    diag(ones(1,size(perception_labels,2)));... % index for each label
    ];
label_definition={'label',props.synsetNames{:}};

designData=[...
    ones(1,size(perception_design,2));... % design
    diag(ones(1,size(perception_design,2)));... % index for each design
    ];
design_definition={'design','session_number','sample_number'};

perceptionMetaData=blkdiag(voxelData,labelData,designData);
perceptionDefinition={voxelData_definition{:},label_definition{:},design_definition{:}};

%% create feature definition data
sleepMetaDefinition=cell(length(sleepDefinition),1);
for i=1:length(sleepDefinition)
    switch sleepDefinition{i}
        case {'voxelData','label','design'}
            sleepMetaDefinition{i}=['0 = not ',sleepDefinition{i},', 1 = ',sleepDefinition{i}];
        case {'X','Y','Z'}
            sleepMetaDefinition{i}=['Value = ',sleepDefinition{i},' coordinate'];
        case props.roiNames
            sleepMetaDefinition{i}=['0 = not ',sleepDefinition{i},' voxel, 1 = ',sleepDefinition{i},' voxel'];
        case 'session_number'
            sleepMetaDefinition{i}='Number = Sesssion number';
        case 'sample_number'
            sleepMetaDefinition{i}='Number = Sample number';
        otherwise
            sleepMetaDefinition{i}=['0 = absent, 1 = present'];
    end
end
% sleepMetaDefinition

perceptionMetaDefinition=cell(length(perceptionDefinition),1);
for i=1:length(perceptionDefinition)
    switch perceptionDefinition{i}
        case {'voxelData','label','design'}
            perceptionMetaDefinition{i}=['0 = not ',perceptionDefinition{i},', 1 = ',perceptionDefinition{i}];
        case {'X','Y','Z'}
            perceptionMetaDefinition{i}=['Value = ',perceptionDefinition{i},' coordinate'];
        case props.roiNames
            perceptionMetaDefinition{i}=['0 = not ',perceptionDefinition{i},' voxel, 1 = ',perceptionDefinition{i},' voxel'];
        case 'session_number'
            perceptionMetaDefinition{i}='Number = Sesssion number';
        case 'sample_number'
            perceptionMetaDefinition{i}='Number = Sample number';
        otherwise
            perceptionMetaDefinition{i}=['0 = absent, 1 = present'];
    end
end
% perceptionMetaDefinition

%% renew data structure
clear pD*

pDs.data=sleep_data;
for i=1:size(sleepMetaData,1)
pDs.metaData.(sleepDefinition{i})=sleepMetaData(i,:);
end
pDs.metaDefinition=sleepMetaDefinition;

pDp.data=perception_data;
for i=1:size(perceptionMetaData,1)
pDp.metaData.(perceptionDefinition{i})=perceptionMetaData(i,:);
end
pDp.metaDefinition=perceptionMetaDefinition;

props_new=props;

%% end
fprintf('Preprocessing end\n')

end




%% functions
% shift data ----------------------------------------
function [D,params]=shiftData(D,params)
nshift   = getFieldDef(params,'nshift',1);

sessionIdx=D.labels(:,ismember(D.labels_type,'Session number'));
unisession=unique(sessionIdx);
num_breaks = length(unisession);

inds_del_data = zeros(num_breaks*nshift,1);
inds_del_labels = zeros(num_breaks*nshift,1);

for itr=1:num_breaks
    idx=find(sessionIdx==unisession(itr));
    bi=idx(1);
    ei=idx(end);
    
    inds_del_data(nshift*(itr-1)+1:nshift*itr) = bi:bi+nshift-1;
    inds_del_labels(nshift*(itr-1)+1:nshift*itr) = ei-nshift+1:ei;
end


D.data(inds_del_data,:) = [];
D.labels(inds_del_labels,:) = [];
end

% reduce outliers ----------------------------------------
function [D,params]=reduceOutliers(D,params)
app_dim   = getFieldDef(params,'app_dim',1);
std_thres = getFieldDef(params,'std_thres',3);
num_its   = getFieldDef(params,'num_its',10);
max_val   = getFieldDef(params,'max_val',inf);
min_val   = getFieldDef(params,'min_val',-inf);

if min_val>-inf || max_val<inf
    method = getFieldDef(params,'method',3);
else
    method = getFieldDef(params,'method',1);
end

if method>=2
    remove = getFieldDef(params,'remove',1);
else
    remove = getFieldDef(params,'remove',0);
end

sessionIdx=D.labels(:,ismember(D.labels_type,'Session number'));
unisession=unique(sessionIdx);
num_breaks = length(unisession);

ind_all = sparse(size(D.data,1),size(D.data,2));     % keep num of outliers
for itb=1:num_breaks
    idx=sessionIdx==unisession(itb);
    data_temp = D.data(idx,:);
    data_size = size(data_temp);
    ind_m1    = zeros(data_size);   % indexes of outliers found in method1 (max std deviation)
    ind_m2    = zeros(data_size);   % indexes of outliers found in method2 (constant min_val, max_val)
    if method==1 || method==3
        for its=1:num_its       % do num_its iterations
            mu = mean(data_temp,app_dim);       % mean
            sd = std(data_temp,0,app_dim);      % standard deviation
            % Find and clip values OVER threshold:
            if app_dim==1,      thres_mat = repmat(mu+sd*std_thres,data_size(1),1);     % make threshold matrix
            else                thres_mat = repmat(mu+sd*std_thres,1,data_size(2));     end
            ind_out            = data_temp>thres_mat;       % find values over thres_mat
            data_temp(ind_out) = thres_mat(ind_out);        % set those to thres_mat
            ind_m1(ind_out)    = ind_m1(ind_out)+1;
            % Find and clip values UNDER threshold:
            if app_dim==1,      thres_mat = repmat(mu-sd*std_thres,data_size(1),1);     % make threshold matrix
            else                thres_mat = repmat(mu-sd*std_thres,1,data_size(2));     end
            ind_out            = data_temp<thres_mat;       % find values under thres_mat
            data_temp(ind_out) = thres_mat(ind_out);        % set those to thres_mat
            ind_m1(ind_out)    = ind_m1(ind_out)+1;
        end
    end
    if method>=2
        if max_val<inf
            % Over max_val:
            ind_out            = data_temp>max_val;
            data_temp(ind_out) = max_val;
            ind_m2(ind_out)    = ind_m2(ind_out)+1;
        end
        if min_val>-inf
            % Below min_val:
            ind_out            = data_temp<min_val;
            data_temp(ind_out) = min_val;
            ind_m2(ind_out)    = ind_m2(ind_out)+1;
        end
    end
    
    D.data(idx,:)  = data_temp;
    ind_all(idx,:) = sparse(ind_m1+ind_m2);
end

if remove
    ind_remove_chans           = find(sum(ind_all,1));
    D.data(:,ind_remove_chans) = [];
    D.xyz(:,ind_remove_chans) = [];
    params.ind_remove_chans=ind_remove_chans;
end
end

% detrend data ----------------------------------------
function D=detrenddata(D,params)
sub_mean  = getFieldDef(params,'sub_mean',0);
method    = getFieldDef(params,'method','linear');

sessionIdx=D.labels(:,ismember(D.labels_type,'Session number'));
unisession=unique(sessionIdx);
num_breaks = length(unisession);
for itb=1:num_breaks
    idx=sessionIdx==unisession(itb);
    
    data_temp1 = D.data(idx,:);
    data_temp2 = detrend(data_temp1,method);
    
    if sub_mean==0
        data_temp2 = data_temp2 + repmat(mean(data_temp1,1),size(data_temp1,1),1);
    end
    
    D.data(idx,:) = data_temp2;
end
end

% normalization for sleep data--------------------------------------
% sleep data
function [D,params]=normalizeSleepData(D,params)
zero_thres = getFieldDef(params,'zero_thres',1);
analysisPeriod = getFieldDef(params,'analysisPeriod',40);
normalizeBaselineDuration = getFieldDef(params,'normalizeBaselineDuration',10);
normalizationOnsetfromAwake = getFieldDef(params,'normalizationOnsetfromAwake',20);

% get analyzed sample index
SleepAwakePeriods = cpDetect(D.labels(:,ismember(D.labels_type,'Experimental instructions')));
AwakePeriods = SleepAwakePeriods(:,2:2:end);
AnalyzeWindow = zeros(2,size(AwakePeriods,2));
NormalizeWindow = zeros(2,size(AwakePeriods,2));
visualDreamIdx = D.labels(AwakePeriods(1,:),ismember(D.labels_type,'Visual dream report'))>0;
ExcludePeriodIdx = D.labels(AwakePeriods(1,:),ismember(D.labels_type,'Excluded samples'))>0;
cnt=0;
for i=1:size(AwakePeriods,2)
    if visualDreamIdx(i)>0 && ExcludePeriodIdx(i)==0
        cnt=cnt+1;
        AnalyzeWindow(:,cnt) = [AwakePeriods(1,i)-analysisPeriod;AwakePeriods(1,i)-1];
        NormalizeWindow(:,cnt) = [AwakePeriods(1,i)-normalizationOnsetfromAwake-normalizeBaselineDuration;AwakePeriods(1,i)-normalizationOnsetfromAwake-1];
    end
end
AnalyzeWindow(:,cnt+1:end) = [];
NormalizeWindow(:,cnt+1:end) = [];

% Sleep data extraction and normalization
data_tmp = cell(1,size(AnalyzeWindow,2));
labels_tmp = cell(1,size(AnalyzeWindow,2));
sleeplabels_tmp = cell(1,size(AnalyzeWindow,2));
zero_ind = cell(1,size(AnalyzeWindow,2));
sleepLabelIdx = ismember(D.labels_type,'EEG sleep score');
triallabels = cell(1,size(AnalyzeWindow,2));
for i = 1:size(AnalyzeWindow,2)
    baseline = mean(D.data(NormalizeWindow(1,i):NormalizeWindow(2,i),:),1);
    baseMat = repmat(baseline,analysisPeriod,1);
    zero_ind{i}  =  (abs(baseline) <=  zero_thres)';
    data_tmp{i} = 100*(D.data(AnalyzeWindow(1,i):AnalyzeWindow(2,i),:)-baseMat)'./baseMat';
    labels_tmp{i} = repmat(D.labels(AnalyzeWindow(2,i)+1,:),analysisPeriod,1)';
    sleeplabels_tmp{i} = repmat(D.labels(AnalyzeWindow(2,i),sleepLabelIdx),analysisPeriod,1)';
    triallabels{i} = ones(1,analysisPeriod)*i;
end
data_tmp2 = [data_tmp{:}]';
labels_tmp2 = [labels_tmp{:}]';
sleeplabels_tmp2 = [sleeplabels_tmp{:}]';
params.triallabels = [triallabels{:}]';
zero_ind2 = [zero_ind{:}]';
if any(sum(zero_ind2)), data_tmp2(:,sum(zero_ind2)>0)  =  0;     end
D.data = data_tmp2;
D.labels = labels_tmp2;
D.labels(:,ismember(D.labels_type,'EEG sleep score')) = sleeplabels_tmp2;

end


% normalization for perception data----------------------------------------
function D=normalizeData(D,params)
base_labels = getFieldDef(params,'base_labels',1);
zero_thres = getFieldDef(params,'zero_thres',1);
norm_mode       = getFieldDef(params,'norm_mode',0);


sessionIdx=D.labels(:,ismember(D.labels_type,'Session number'));
unisession=unique(sessionIdx);
num_breaks = length(unisession);

for itb=1:num_breaks
    idx=find(sessionIdx==unisession(itb));
    
    % Pull section (run) out:
    data_temp = D.data(idx,:);
    
    % Find indexes of base condition:
    if isstr(base_labels) && strcmp(base_labels,'all')
        ind_use = 1:numel(idx);
    elseif isnumeric(base_labels)
        ind_use = find(ismember(D.labels(idx),base_labels));
    end
    
    % Calc baseline:
    baseline = mean(data_temp(ind_use,:),1);
    sd = std(data_temp(ind_use,:),[],1);
    
    % Find baseline indexes ~= 0 (to avoid dividing by them):
    if norm_mode == 0 || norm_mode == 1
        zero_ind  = find(abs(baseline) <= zero_thres);
    elseif norm_mode == 3
        zero_ind  = find(abs(sd) <= zero_thres);
    else
        zero_ind = [];
    end
    
    num_zeros = numel(zero_ind);
    if num_zeros>1      % If there are some zero values in the baseline
        baseline(zero_ind) = zero_thres;    % set them to zero_thres to avoid dividing by them
        
    end
    
    % mean mat
    baseline_mat = repmat(baseline,size(data_temp,1),1);
    % sd mat
    sd_mat = repmat(sd,size(data_temp,1),1);
    
    switch norm_mode
        case 0
            % percent-signal change
            data_temp    = 100 * (data_temp - baseline_mat) ./ baseline_mat;
        case 1
            % division by mean only
            data_temp    = 100 * data_temp ./ baseline_mat;
        case 2
            % subtraction of mean only
            data_temp    = data_temp - baseline_mat;
        case 3
            % z-score
            data_temp    = (data_temp - baseline_mat) ./ sd_mat;
    end
    
    % Set zero_ind data to zero:
    if num_zeros>1,      data_temp(zero_ind) = 0;     end
    % Put normalized section (run) back:
    D.data(idx,:) = data_temp;
end
end

% averaging sleep data ----------------------------------------
function D=averageSleepData(D,params)
twShift = getFieldDef(params,'twShift',1);
winDur = getFieldDef(params,'winDur',3);
if ~isfield(params,'triallabels')
    error('missing trial labels.')
else
    triallabels=params.triallabels;
end

% extract analysis window and average them
inds_blocks=cpDetect(triallabels);
data_tmp=zeros(size(inds_blocks,2),size(D.data,2));
labels_tmp=zeros(size(inds_blocks,2),size(D.labels,2));
for i=1:size(inds_blocks,2)
    idx=inds_blocks(2,i)-(twShift+winDur-1:-1:twShift)+1;
    data_tmp(i,:)=mean(D.data(idx,:));
    labels_tmp(i,:)=D.labels(idx(end),:);
end
D.data=data_tmp;
D.labels=labels_tmp;
end



% averaging perception data ----------------------------------------
function D=averagePerceptionData(D)
analyzeLabels=[4]; % 4 is labels for stimuls periods
% averaged data
conditionlabelIdx=ismember(D.labels_type,'Experimental conditions');
inds=ismember(D.labels(:,conditionlabelIdx),analyzeLabels);
cp=cpDetect(inds);
use_periods=cp(:,2:2:end);
data_tmp=zeros(size(use_periods,2),size(D.data,2));
labels_tmp=zeros(size(use_periods,2),size(D.labels,2));
for itr=1:size(use_periods,2)
    ind=use_periods(1,itr):use_periods(2,itr);
    data_tmp(itr,:)=mean(D.data(ind,:),1);
    labels_tmp(itr,:)=D.labels(ind(1),:);
end
D.data=data_tmp;
D.labels=labels_tmp;
end



