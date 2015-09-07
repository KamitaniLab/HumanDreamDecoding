function X_1=addBias(X)
% addBias - add the bias term
%
% [Input]
%   -X: training data [N x M]
%
%
% [Output]
%	-X_1: training data with bias term [N x M+1] 
%
% [Example]
%   X=randn(100,1);
%   X_1=addBias(X)
% 
% 
% 
% 
%
% [Related function]
%  
%
% Created  By: Tomoyasu Horikawa horikawa-t@atr.jp 2009/10/11
%
%

%% add the bias term
X_1=[X,ones(size(X,1),1)];
