function [data normMean normScale] = normFeature(data,meanMode,scaleMode,normMean,normScale)
%
% normFeature - normalize feature vector by specified mode
%
% [input]
%  - data: [samples x features]
%  - meanMode: mean normalization mode [all, feature, sample, or none]
%  - scaleMode: scale normalization mode [all, feature, sample, or none]
%  - normMean: if supplied, normalized by this value(s). should match
%  feature dim
%  - normScale: if supplied, normalized by this value(s). should match
%  feature dim
% 
% [output]
%  - data: normalized data [samples x features]
%  - normMean: mean used for normalization [scalar, or vector]
%  - normScale: scale used for normalization [scalar, or vector] 
%
% [note]
%  - mostly identical what normalize_feature.m in SLR toolbox does.
%
% 09/09/10 written by Yoichi Miyawaki (yoichi_m@atr.jp)

if nargin < 3 || nargin > 5

  help normFeature;
  error;
  return;

end


if nargin < 5
  
  switch scaleMode
    case 'all'
      %normScale = max(abs(data(:)));
      normScale = std(data(:));
    case 'feature'
      %normScale = max(abs(data),[],1);
      normScale = std(data,0,1);      
    case 'sample'
      %normScale = max(abs(data),[],1);
      normScale = std(data,0,2);      
    case 'none'
      normScale = 1;
  end

end


if nargin < 4

  switch meanMode
    case 'all'
      normMean = mean(data(:));
    case 'feature'
      normMean = mean(data,1);
    case 'sample'
      normMean = mean(data,2);
    case 'none'
      normMean = 0;
  end
  
end


data = repadd(data,-normMean);
data = repmultiply(data,1./normScale);  


% end