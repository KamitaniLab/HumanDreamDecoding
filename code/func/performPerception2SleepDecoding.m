function results=performPerception2SleepDecoding(Ds,Dp,params)
%
% perform dream decoding analysis using perception-trained decoder
%

%% sample selection
% sleep data
if params.class1==params.class2 % 1 vs others
    label1_ind=Ds.labels(:,params.class1)==1;
    label2_ind=Ds.labels(:,params.class2)==0;
else % 1 vs 1
    label1_ind=Ds.labels(:,params.class1)==1;
    label2_ind=Ds.labels(:,params.class2)==1;
end
% remove overlaps (= samples labeled with both class)
overlap=label1_ind==label2_ind;
label1_ind(overlap)=0;
label2_ind(overlap)=0;
Ds.labels=sum([label1_ind,label2_ind*2],2);
results.sample=[sum(label1_ind),sum(label2_ind)];

% perception data
if params.class1==params.class2 % 1 vs others
    label1_ind=Dp.labels(:,params.class1)==1;
    label2_ind=Dp.labels(:,params.class2)==0;
else % 1 vs 1
    label1_ind=Dp.labels(:,params.class1)==1;
    label2_ind=Dp.labels(:,params.class2)==1;
end
Dp.labels=sum([label1_ind,label2_ind*2],2);

% use samples labeled with either the pairs
Dp.data=Dp.data(Dp.labels>0,:);
Dp.labels=Dp.labels(Dp.labels>0);

%% voxel selection using paired t-test
params.conds=[1,2]; % used labels
[Dp,params]=selectTopTvals(Dp,params);

%% z-normalization for each sample
Dp.data=zscore(Dp.data')';
Ds.data=zscore(Ds.data')';

%% classification
P.mode=1;
[r]=libsvm_h(Dp,P);

[nouse,ord]=ismember(Dp.xyz',Ds.xyz','rows');
w=zeros(size(Ds.data,2),1);
w_tmp=uninorm(r.weights);
w(ord)=w_tmp;
bias=r.bias/sqrt(sum(r.weights.^2));
w=[w;bias];
dcv=addBias(Ds.data)*w;
if Dp.labels(1)==1
    r.preds=(dcv<0)+1;
else
    r.preds=(dcv>0)+1;
end
r.labels=Ds.labels;

results.preds=r.preds;
results.labels=r.labels;
results.decvals=dcv;
results.weights=w;
results.corrRate=mean(r.preds(r.labels>0)==r.labels(r.labels>0))*100;
%% end
end
