function [D,params]=selectTopTvals(D,params)
conds       = getFieldDef(params,'conds',[1,2]);
num_comp   = getFieldDef(params,'num_comp',size(D.data,2));
ind1=find(D.labels==conds(1));
ind2=find(D.labels==conds(2));

% class1
n1=length(ind1);
mu1=mean(D.data(ind1,:),1);
var1=var(D.data(ind1,:),[],1);
% class 2
n2=length(ind2);
mu2=mean(D.data(ind2,:),1);
var2=var(D.data(ind2,:),[],1);

% calc tvals
s=sqrt((var1*(n1-1)+var2*(n2-1))/(n1+n2-2));
tvals=(mu1-mu2)./(s*sqrt(1/n1+1/n2));
tvals =abs(tvals)';

% Select top N within range
[tvals, inds_tvals] = sort(tvals,'descend');

num_comp            = min(abs(num_comp),length(inds_tvals));
inds_tvals          = inds_tvals(1:abs(num_comp));
D.data              = D.data(:,inds_tvals);
D.xyz              = D.xyz(:,inds_tvals);

params.tvals      = tvals(1:abs(num_comp));
params.inds_tvals = inds_tvals;

end