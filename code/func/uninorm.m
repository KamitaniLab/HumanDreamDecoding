function n=uninorm(X,dir)
% uninorm -- caclculate unit norm vector of input X 
% 
% [Inputs]
%   -X: data
%   -dir:direction [col=1(default), row=2]
% 
% 
% [Outputs] 
%   -n: unit norm vector matrix
% 
% 
% 
% Written by Tomoyasu Horikawa horikawa-t@atr.jp 2011/10/07
% 
% X=randn(10,5);
% n=uninorm(X)
% sum(n.*n,1)
% n=uninorm(X,2)
% sum(n.*n,2)
% 
% 
%%
if nargin < 2
    dir=1;
end
if ~isreal(X) 
    warning('Computing distance table using imaginary inputs. Results may be off.');
end

if dir==1
    X=X';
end
% n=bsxfun(@rdivide,X,sqrt(diag(X*X')));
% unit norm
n=X;
for itr=1:size(n,1)
    n(itr,:)=X(itr,:)/sqrt(X(itr,:)*X(itr,:)');
end

if dir==1
    n=n';
end

%%

