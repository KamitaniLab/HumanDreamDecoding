function [cpm cpcomp cplen]=cpDetect(input)
% cpDetect -- detect the changepoint of components in input vector
% this function can be used for converting design matrix to ***_inds data/
% 
% [Input]
%   -input: vector
% 
% [Output]
%   -cpm    : change point index matrix
%   -cpcomp : top components of each cp interval
%   -cplen  : lengh of each cp interval
% 
% e.g.
%  input=[1,1,1,2,2,3,3,3,3,3,1,1,1,1,2]
% [cpm cpcomp]=cpDetect(input)
% cpm =
% 
%      1     4     6    11    15
%      3     5    10    14    15
% cpcomp =
%      1     2     3     1     2
% 
% input=['r','r','r','b','b','b','r','r','r']
% cpm =
%      1     4     7
%      3     6     9
% cpcomp =
%      rbr
% 
% input=[{'ra'},{'ra'},{'rrr'},{'br'},{'br'},{'rb'},{'rb'}]
% cpm =
%      1     3     4     6
%      2     3     5     7
% cpcomp = 
%     'ra'    'rrr'    'br'    'rb'
% 
% created  by HORIKAWA tomoyasu 09/09/03
% modified by HORIKAWA tomoyasu 09/10/22 add cpcomp
a=input(1);
cpm=zeros(2,1);
count=1;
if isnumeric(a)||islogical(a)
    for i=1:length(input)
        if input(i)~=a
            cpm(2,count)=i-1;
            cpm(1,count+1)=i;
            count=count+1;
            a=input(i);
        end
    end
elseif ischar(a)||iscell(a)
    for i=1:length(input)
        if ~strcmp(input(i),a)
            cpm(2,count)=i-1;
            cpm(1,count+1)=i;
            count=count+1;
            a=input(i);
        end
    end
 
end
cpm(1,1)=1;
cpm(end)=length(input);
cpcomp=input(cpm(1,:));
cplen=cpm(2,:)-cpm(1,:)+1;
