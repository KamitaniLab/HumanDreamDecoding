function [results, pars] = libsvm_h(D, pars)
%libsvm_h - performs multi SVM using libsvm-mat, either train or test
%[results, pars] = libsvm_h(D, pars)
%
% Inputs:
%	D.data    - 2D matrix of data
%	D.labels  - labels matching samples of 'data'; only [] for test mode
%
% Optional:
%	pars.model- SVM 'model', including weights; optional for training
%	pars.mode   - train=1 (make weights) or test=2 mode (use weights)
%   pars.conds - class label, which is necessary to relabel correctly
%
%
% LibSVM pars fields:
%	.kernel - 0=linear, 1=poly, 2=rbf, 3=sigmoid; default=0
%	.cost   - C of C-SVC, epsilon-SVR, and nu-SVR; default=1
%	.gamma  - set gamma in kernel function; default 1/k
%	.coef   - set coef0 in kernel function; default=0
%	.degree - set degree in kernel function; default=3
%	.prob   - output probabilities as decVals? 0-no, 1-yes; default=1
%
% Outputs:
%	results - struct contain ANY result as a field, typically:
%       .mode   - name of this function
%       .pars   - parameters used (minus weights)
%		.decVals- decision values (raw classifier output)
%		.preds	- labels predicted by the models
%	pars    - modified pars, new weights will be added here
%
% Note: modelSwitch will make the remaining fields of results.
%
% Example:
%	>> % dataTrain, dataTest - [nSigs x nSamples]
%	>> pars.verbose = 0;	% to avoid printing (new)
%	>> % 1 - train mode:
%	>> [weights, results, pars] = libsvm_h(dataTrain, pars, labelsTrain, 1);
%	>> % 2 - test  mode:
%	>> [weights, results, pars] = libsvm_h(dataTest, pars, labelsTest, 2);
%
% Calls: svmtrain, svmpredict (help svmtrain for more info)
% Requires: libsvm-mat-3.01, (c) 2000-2005 Chih-Chung Chang & Chih-Jen Lin
% Info: http://www.csie.ntu.edu.tw/~cjlin/libsvm
% Status: basic testing
%


%% Check and get pars:
if exist('D','var')==0    || isempty(D),        error('Wrong args');        end
if exist('pars','var')==0 || isempty(pars),     pars = [];                  end

pars    = getFieldDef(pars,mfilename,pars);    % unnest, if need
model   = getFieldDef(pars,'model',[]);
mode    = getFieldDef(pars,'mode',1);
conds   = getFieldDef(pars,'conds',[]);
normMeanMode  = getFieldDef(pars,'normMeanMode','feature');
normScaleMode = getFieldDef(pars,'normScaleMode','feature');
normMean  = getFieldDef(pars,'normMean',0);
normScale = getFieldDef(pars,'normScale',1);
normMode = getFieldDef(pars,'normMode','test');


if     mode==1 && isempty(D.labels),    error('must have ''labels'' for train');
elseif mode==2 && isempty(model),       error('must have ''model'' for test');      end

if ~isempty(conds) && max(conds)~=length(unique(conds)) || isempty(find(conds==0,1))==0
  fprintf('\nWarning: unrecommended format of ''labels''');
  fprintf('\n rename ''labels'' in ''libsvm''\n');
    
  [D.labels, labels_new, labels_old] = reIndex(D.labels,1:length(conds),conds);
end

%% SVM pars:
kernel = getFieldDef(pars,'kernel',0);      % linear
gamma  = getFieldDef(pars,'gamma',0);       % NOTE: gamma=0 defaults to 1/k
prob   = getFieldDef(pars,'prob',1);
cost   = getFieldDef(pars,'cost',1);
coef   = getFieldDef(pars,'coef',0);
degree = getFieldDef(pars,'degree',3);

ops1 = sprintf('-t %d -c %g -r %g -d %g -q ', kernel, cost, coef, degree);
if prob,        ops1 = [ops1 ' -b 1'];   ops2 = '-b 1';     end
if gamma,       ops1 = [ops1 ' -g ' num2str(gamma)];        end


%% Test mode:
if mode==2
    if strcmp(normMode,'test')
        if size(D.data,1) == 1 && strcmp(normMeanMode,'feature')
            fprintf('\nWARNINIG: data sample size is 1. this normalization convnert all features into 0.\n');
        end
        data = normFeature(D.data,normMeanMode,normScaleMode);
    elseif strcmp(normMode,'training')
        data = normFeature(D.data,normMeanMode,normScaleMode,normMean,normScale);
    else
        error('normalization mode error');
    end  
    
    [preds, nouse, dec_vals] = svmpredict(D.labels,data,model,ops2);
    
    
%% Train mode:
else
    
    [data normMean normScale] = normFeature(D.data,normMeanMode,normScaleMode);
    
    model                    = svmtrain(D.labels,data,ops1);
    [preds, nouse, dec_vals] = svmpredict(D.labels,data,model,ops2);
    
    pars.model = model;
end

%% Retrun results:
if exist('labels_old','var') && isempty(labels_old)==0  
    D.labels = reIndex(D.labels,labels_old,labels_new);
    preds    = reIndex(preds,labels_old,labels_new);
end


%% Return results:
results.model    = mfilename;
results.preds    = preds;
results.labels   = D.labels;
results.dec_vals = dec_vals;
results.weights  = model.SVs' * model.sv_coef;
results.bias     = -model.rho;

