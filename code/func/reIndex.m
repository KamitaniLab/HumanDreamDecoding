function [list_new, new_indx, old_indx] = reIndex(list, new_indx, old_indx)
%reIndex - re-indexes a list, replacing old_indx items with new_indx items
%[list_new, new_indx, old_indx] = reIndex(list, new_indx, old_indx)
%
% Replaces items of 'list' (any matrix) specified in old_indx with corresponding items in
% new_indx.  If old_indx is absent or [], items are found with with unique(list), which
% sorts them.  If new_indx is absent or [], then it uses [1 2 3 ...] as the new index. 
%
% Exp:
%	>> list = [-1 0 153 0; 0 153 153 -1]
%    -1     0   153     0
%     0   153   153    -1
%	>> reIndex(list,[0 3 2])
%     0     3     2     3
%     3     2     2     0
%	>> reIndex(list)
%     1     2     3     2
%     2     3     3     1
%
% Created  By: Alex Harner (1),     alexh@atr.jp      06/10/27
% Modified By: Alex Harner (1),     alexh@atr.jp      06/10/30
% Modified By: Satoshi MURATA (1),  satoshi-m@atr.jp  08/10/09
% (1) ATR Intl. Computational Neuroscience Labs, Decoding Group


%% Check and get pars:
if exist('list','var')==0 || isempty(list)
    list_new = [];
    new_indx = [];
    old_indx = [];
    return;
end

if exist('old_indx','var')==0 || isempty(old_indx)
    old_indx = unique(list)';
end
if exist('new_indx','var')==0 || isempty(new_indx)
    new_indx = 1:length(old_indx);
end

if length(old_indx)~=length(new_indx)
    error('''new_indx'' and ''old_indx'' must be same length');
end


%% Re-indexes:
list_new = list;
for iti=1:length(old_indx)
    list_new(list==old_indx(iti)) = new_indx(iti);
end
